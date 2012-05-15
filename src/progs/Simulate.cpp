/**
 * \file  Simulate.cpp
 * \brief TODO
 *
 * \authors Philip J. Uren, Andrew D. Smith
 *
 * \section copyright Copyright Details
 * Copyright (C) 2012
 * University of Southern California,
 * Philip J. Uren, Andrew D. Smith
 *
 * \section license License Details
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \section bugs Known Bugs
 *
 * \section history Revision History
 **/

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cassert>

// From piranha common
#include "NegativeBinomialRegression.hpp"
#include "TruncatedPoissonRegression.hpp"
#include "PoissonRegression.hpp"
#include "ZTNBRegression.hpp"
#include "DistributionMixtureModel.hpp"
#include "RegressionMixtureModel.hpp"
#include "FittingMethod.hpp"
#include "ZTNB.hpp"

// From smithlab_cpp
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "RNG.hpp"

// From BAMTools
#include "api/BamWriter.h"
#include "api/SamHeader.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"

using std::vector;
using std::string;
using std::endl;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::stringstream;
using std::accumulate;
using std::tr1::unordered_map;

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamWriter;

template<typename T, size_t N>
T * end(T (&ra)[N]) {
    return ra + N;
}

/******************************************************************************
 * Static functions for sampling from the Poisson, NB, ZTP, and ZTNB distros
 */

/**
 * \brief TODO
 */
static double
calcMu (const vector<double>& covariates, const vector<double>& coefficients) {
  //assert(len(covariates) == len(regressionParams))

  double sum = 0;
  for (size_t i=0; i<covariates.size(); i++) {
    sum += covariates[i] * coefficients[i];
  }
  return exp(sum);
}

/**
 * \brief TODO
 */
static size_t
simulateReadCountPoisson(const double lambda) {
  gsl_rng *r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, rand());
  size_t res = gsl_ran_poisson(r, lambda);
  gsl_rng_free(r);
  return res;
}

/**
 * \brief TODO
 */
static size_t
simulateReadCountNB(const double mu,
                    const double alpha) {
  gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(rng, rand());
  // to convert from alpha and mu to p and r...
  const double r = 1/alpha;
  const double p = r / (r + mu);
  size_t res = gsl_ran_negative_binomial(rng, p, r);
  gsl_rng_free(rng);
  return res;
}

/**
 * \brief             TODO
 * \param truncatedAt TODO
 * \param dropped     A return parameter. Will be increased by the number of
 *                    sampled counts that were dropped because they didn't meet
 *                    the truncation point
 */
static size_t
simulateReadCountTruncatedPoisson(const double lambda, const int truncateAt,
                                  int& dropped) {
  int res = -1;
  dropped -= 1;
  while ((res == -1) || (res <= truncateAt)) {
    dropped += 1;
    res = simulateReadCountPoisson(lambda);
  }
  return res;
}

/**
 * \brief             TODO
 * \param truncatedAt TODO
 * \param dropped     A return parameter. Will be increased by the number of
 *                    sampled counts that were dropped because they didn't meet
 *                    the truncation point
 */
static size_t
simulateReadCountTruncatedNB(const double mu, const double alpha,
                                       const int truncateAt, int& dropped) {
  int res = -1;
  dropped -= 1;
  while ((res == -1) || (res <= truncateAt)) {
    dropped += 1;
    res = simulateReadCountNB(mu, alpha);
  }
  return res;
}


/******************************************************************************
 * Static functions for sampling from the Poisson, NB, ZTP, and ZTNB
 * regression distributions
 */

/**
 * \brief TODO
 */
static size_t
simulateReadCountPoissonRegression(const vector<double>& covariates,
                                   const vector<double>& coefficients) {
  const double mu = calcMu(covariates, coefficients);
  return simulateReadCountPoisson(mu);
}

/**
 * \brief TODO
 */
static size_t
simulateReadCountNBRegression(const vector<double>& covariates,
                                          const vector<double>& coefficients,
                                          const double alpha) {
  const double mu = calcMu(covariates, coefficients);
  return simulateReadCountNB(mu, alpha);
}

/**
 * \brief             TODO
 * \param truncatedAt TODO
 * \param dropped     is a return parameter for the number of sampled counts
 *                    that were dropped because they didn't meet the truncation
 *                    point
 */
static size_t
simulateReadCountTruncatedPoissonRegression(const vector<double>& covariates,
                                   const vector<double>& coefficients,
                                   const int truncateAt,
                                   int& dropped) {
  const double mu = calcMu(covariates, coefficients);
  return simulateReadCountTruncatedPoisson(mu, truncateAt, dropped);
}

/**
 * \brief             TODO
 * \param truncatedAt TODO
 * \param dropped     is a return parameter for the number of sampled counts
 *                    that were dropped because they didn't meet the truncation
 *                    point
 */
static size_t
simulateReadCountTruncatedNBRegression(const vector<double>& covariates,
                                   const vector<double>& coefficients,
                                   const double alpha,
                                   const int truncateAt,
                                   int& dropped) {
  const double mu = calcMu(covariates, coefficients);
  return simulateReadCountTruncatedNB(mu, alpha, truncateAt, dropped);
}



/******************************************************************************
 * Static functions for the simulation of read counts from either simple
 * distributions, or regression -based distributions. Supported distros
 * are currently: Poisson, negative binomial, zero-truncated Poisson,
 * zero-truncated negative binomial (all simple distributions), and the
 * corresponding regression based models for the regression read counts
 */

/**
 * \brief simulate scores for the count regions using a simple distribution
 *        type
 */
void
simulateScoresSimple(vector<GenomicRegion>& countRegions,
                     const string& distType,
                     const double maxMu,
                     const double maxAlpha) {
  // just allow zero-truncation
  const int truncateAt = 0;

  // how many counts we dropped because they don't meet the truncation point
  int dropped = 0;

  // we might not need alpha, but we'll just go ahead and pick one anyway
  double rmax = (double) RAND_MAX;
  const double mu = (rand() / rmax) * maxMu;
  const double alpha = (rand() / rmax) * maxAlpha;

  for (size_t i=0; i<countRegions.size(); i++) {
    size_t score;
    if (distType == "Poisson") {
      score = simulateReadCountPoisson(mu);
    }
    else if (distType == "NegativeBinomial") {
      score = simulateReadCountNB(mu, alpha);
    }
    else if (distType == "ZeroTruncatedPoisson") {
      score = simulateReadCountTruncatedPoisson(mu, truncateAt, dropped);
    }
    else if (distType == "ZeroTruncatedNegativeBinomial") {
      score = simulateReadCountTruncatedNB(mu, alpha, truncateAt, dropped);
    }
    else {
      throw SMITHLABException("Unknown simple distribution type");
    }
    countRegions[i].set_score(float(score));
  }

  cerr << "generated " << countRegions.size() << " regions with counts "
       << "drawn from " << distType << " with mu = " << mu;
  if ((distType == "NegativeBinomial") ||
     (distType == "ZeroTruncatedNegativeBinomial")) {
    cerr << " and alpha = " << alpha;
  }
  cerr << endl;
  if (dropped != 0) {
    cerr << "truncated " << dropped << " elements" << endl;
  }
}

/**
 * \brief simulate ead counts drawn from a distribution that is a regression
 *        model, either the Poisson, negative binomial, zero-truncated
 *        Poisson or zero-truncated negative binomial distributions
 */
void
simulateScoresRegression(vector<GenomicRegion>& countRegions,
                             vector< vector<double> >& covariates,
                             vector<double> coefficients,
                             const string& distType,
                             const double maxMu,
                             const double maxAlpha) {
  // just allow zero-truncation
  const int truncateAt = 0;

  // how many counts we dropped because they don't meet the truncation point
  int dropped = 0;

  // we might not need alpha, but we'll just go ahead and pick one anyway
  double rmax = (double) RAND_MAX;
  const double alpha = (rand() / rmax) * maxAlpha;
  if ((distType == "NegativeBinomialRegression") ||
      (distType == "ZeroTruncatedNegativeBinomialRegression"))
    cerr << "selected alpha: " << alpha << endl;

  for (size_t i=0; i<countRegions.size(); i++) {
    size_t score;
    if (distType == "PoissonRegression") {
      score = simulateReadCountPoissonRegression(covariates[i],
                                                 coefficients);
    }
    else if (distType == "NegativeBinomialRegression") {
      score = simulateReadCountNBRegression(covariates[i],
                                            coefficients, alpha);
    }
    else if (distType == "ZeroTruncatedPoissonRegression") {
      score = simulateReadCountTruncatedPoissonRegression(covariates[i],
                                                          coefficients,
                                                          truncateAt,
                                                          dropped);
    }
    else if (distType == "ZeroTruncatedNegativeBinomialRegression") {
      score = simulateReadCountTruncatedNBRegression(covariates[i],
                                                     coefficients, alpha,
                                                     truncateAt,
                                                     dropped);
    }
    else {
      throw SMITHLABException("Unknown regression distribution type");
    }
    countRegions[i].set_score(float(score));
  }
}

/**
 * \brief Simulate regions from the given set of chroms. Regions have a
 *        fixed bin size, and will always be offset from start of the chrom
 *        by a multiple of the bin size.
 */
static vector<GenomicRegion>
simulateRegions(const size_t num, size_t binSize,
                const vector<string>& chroms) {
  const size_t MAX_BIN_START = 1000000;
  const size_t MAX_POOL = MAX_BIN_START / binSize;
  double rmax = (double) RAND_MAX;

  vector<size_t> numPerChrom;
  for (size_t i=0; i<num; i++) numPerChrom.push_back(num/chroms.size());
  numPerChrom[numPerChrom.size()-1] += num % chroms.size();

  vector<GenomicRegion> res;
  for (size_t i=0; i<chroms.size(); i++) {
    for (size_t j=0; j<numPerChrom[i]; j++) {
      size_t s = size_t ((rand() / rmax) * (MAX_POOL)) * binSize;
      size_t e = s + binSize;
      res.push_back(GenomicRegion(chroms[i], s, e, "X", 0, '+'));
    }
  }

  sort(res.begin(), res.end());
  return res;
}


/**
 * \brief TODO
 */
static void
GenomicRegionToBamAlignment(const unordered_map<string, size_t> &chrom_lookup,
                            const GenomicRegion &r, BamAlignment &ba) {

  const unordered_map<string, size_t>::const_iterator
    the_chrom(chrom_lookup.find(r.get_chrom()));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no id for chrom: " + r.get_chrom());

  /* Fields corresponding directly to GenomicRegion format fields : */
  ba.RefID = the_chrom->second; // ID number for reference sequence
  ba.Position = r.get_start();   // position (0-based) where alignment starts
  ba.Length = r.get_width();     // length of query sequence
  ba.Name = r.get_name();        // read name
  ba.MapQuality = r.get_score(); // mapping quality score
  ba.QueryBases = string(r.get_width(), 'A');
  ba.Qualities = string(r.get_width(), 'I');
  ba.SetIsReverseStrand(r.neg_strand());     // sets value of "alignment
                                             // mapped to reverse strand" flag
}


/**
 * \brief TODO
 *
 */
static void
MakeBamHeader(const vector<GenomicRegion> regions, RefVector &refs,
              unordered_map<string, size_t> &chrom_lookup,
              string &header) {
  // figure out the 'size' of the chromosomes
  chrom_lookup.clear();
  typedef unordered_map<string, size_t>::const_iterator miterator;
  for (size_t i=0; i<regions.size(); i++) {
    const miterator the_chrom(chrom_lookup.find(regions[i].get_chrom()));
    if ((the_chrom == chrom_lookup.end()) ||
        (the_chrom->second < regions[i].get_end())) {
      chrom_lookup[regions[i].get_chrom()] = regions[i].get_end();
    }
  }

  vector<string> chroms;
  vector<size_t> chromlen;
  for (miterator it = chrom_lookup.begin(); it != chrom_lookup.end(); ++it) {
    chroms.push_back(it->first);
  }
  sort(chroms.begin(), chroms.end());
  for (size_t i=0; i<chroms.size(); i++) {
    chromlen.push_back(chrom_lookup[chroms[i]]);
  }

  // Mandatory first line
  header += "@HD\tVN:1.0\tSO:unsorted\n";

  // create a BAM header (@SQ) entry for each chrom in the BEDTools genome file.
  for (size_t i = 0; i < chroms.size(); ++i) {
    chrom_lookup[chroms[i]] = i;

    // add to the header text
    header += ("@SQ\tSN:" + chroms[i] + "\tLN:" + toa(chromlen[i]) + "\n");

    // create a chrom entry and add it to
    // the RefVector
    BamTools::RefData chrom;
    chrom.RefName = chroms[i];
    chrom.RefLength = chromlen[i];
    refs.push_back(chrom);
  }
}

/**
 * \brief TODO
 */
static void
writeBAM(const vector<GenomicRegion>& regions, const string& outfn) {
  RefVector rv;
  string bamHeader;
  unordered_map<string, size_t> chrom_lookup;
  MakeBamHeader(regions, rv, chrom_lookup, bamHeader);

  SamHeader sh;
  sh.SetHeaderText(bamHeader);
  assert(sh.IsValid(true));

  BamWriter bw;
  BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
  bw.SetCompressionMode(compressionMode);
  bw.Open(outfn, sh, rv);

  for (size_t i=0; i<regions.size(); i++) {
    BamAlignment ba;
    GenomicRegionToBamAlignment(chrom_lookup, regions[i], ba);
    bw.SaveAlignment(ba);
  }
  bw.Close();
}

static void
writeBED(const vector<GenomicRegion>& regions, const string& outfn) {
  ofstream outfh(outfn.c_str());
  for (size_t i=0; i<regions.size(); i++) outfh << regions[i] << endl;
}

/**
 * \brief TODO
 */
static void
generateReads(GenomicRegion &r, size_t readLength,
              size_t &readCounter, Runif &rng,
              vector<GenomicRegion> &reads) {
  for (size_t j=0; j<r.get_score(); j++) {
    // generate a read in this bin
    size_t start = rng.runif(r.get_start(), r.get_end());
    size_t end = start + readLength;
    stringstream name;
    name << "read_" << readCounter;
    ++readCounter;
    GenomicRegion read(r.get_chrom(), start, end,
                    name.str(), 0, '+');
    reads.push_back(read);
  }
}

/**
 * \brief TODO
 */
static void
writeResponseFile(const string& outfn, vector<GenomicRegion> regions,
                  string formatString) {
  Runif rng;
  size_t readLength = 36;
  size_t readCounter = 0;
  vector<GenomicRegion> outputRegions;
  if ((formatString == "BED_UNBINNED") || (formatString == "BAM_UNBINNED")) {
    cerr << "creating individual reads" << endl;
    for (size_t i=0; i<regions.size(); i++) {
      generateReads(regions[i], readLength, readCounter, rng, outputRegions);
    }
    cerr << "writing output file" << endl;
    if (formatString == "BED_UNBINNED") writeBED(outputRegions, outfn);
    else writeBAM(outputRegions, outfn);
  } else if (formatString == "BED_BINNED") {
    writeBED(regions, outfn);
  } else {
    stringstream ss;
    ss << "Unknown output format for response file: " << formatString;
    throw SMITHLABException(ss.str());
  }
}

/**
 * \brief Todo
 */
int
main(int argc, const char* argv[]) {
  try {
    // command line args are put into these
    string distributionType = "ZeroTruncatedNegativeBinomial";
    size_t numPoints;
    //size_t numComponents = 1;
    bool VERBOSE = false;
    string covariateFilenames, countFilename;
    string formatString = "BED_BINNED";
    size_t randSeed = size_t(time(NULL));
    size_t binSize = 1;

    // constants used in this function
    const double MAX_MU = 300;
    const double MAX_ALPHA = 10;
    const double MAX_COEFF = 1.5;
    const double MAX_COVAR_VAL = 2;
    const char* DEFAULT_CHROMS[] = {"chr1", "chr2", "chr3", "chr4", "chr5"};


    /****************** COMMAND LINE OPTIONS ********************/
    string about = "Simulate Version 1.0 -- A program to simulate input data "
                   "for the Piranha peak caller.  Written by Philip J. Uren "
                   "and Andrew D. Smith.";
    OptionParser opt_parse(strip_path(argv[0]), about, "");
    opt_parse.add_opt("distribution", 'd', "what kind of distribution to use "
                      "for generating read counts? Currently supports "
                      "Poisson, NegativeBinomial, ZeroTruncatedPoisson, "
                      "ZeroTruncatedNegativeBinomial, PoissonRegression, "
                      "NegativeBinomialRegression, "
                      "ZeroTruncatedPoissonRegression, "
                      "ZeroTruncatedNegativeBinomialRegression",
                      false, distributionType);
    opt_parse.add_opt("numPoints", 'n', "number of locations with at least a "
                      "single read mapping.",  true, numPoints);
    opt_parse.add_opt("VERBOSE", 'v', "output additional messages about run "
                      "to stderr if set", false, VERBOSE);
    opt_parse.add_opt("seed", 's', "seed (positive integer) for random "
                      "numbers (defaults to system time)", false, randSeed);
    // not implemented yet
    /*opt_parse.add_opt("components", 'k', "number of mixture components", false,
                      numComponents); */
    opt_parse.add_opt("countFilename", 'r', "filename to write counts to",
                      true, countFilename);
    opt_parse.add_opt("covFilename", 'c', "filename(s) to write covariates "
                      "to (if you have more than 1, wrap in quotes)", false,
                      covariateFilenames);
    opt_parse.add_opt("format", 'f', "format for the counts file. Options are "
                      "BED_BINNED (default), BED_UNBINNED, BAM_UNBINNED",
                      false, formatString);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);

    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }

    /****************** END COMMAND LINE OPTIONS *****************/
    // so we don't get the same stuff every time...
    srand(randSeed);

    // first we need to split the covariates filenames
    vector<string> covariateFilenamesV;
    smithlab::split_whitespace(covariateFilenames, covariateFilenamesV);

    // now make the regions
    vector<GenomicRegion> regions;
    regions = simulateRegions(numPoints, binSize,
        vector<string>(DEFAULT_CHROMS, end(DEFAULT_CHROMS)));

    // now simulate the read counts for simple distributions
    if ((distributionType == "Poisson") ||
        (distributionType == "NegativeBinomial") ||
        (distributionType == "ZeroTruncatedPoisson") ||
        (distributionType == "ZeroTruncatedNegativeBinomial")) {
      if (VERBOSE) cerr << "generating data" << endl;
      simulateScoresSimple(regions, distributionType, MAX_MU, MAX_ALPHA);
      if (VERBOSE) cerr << "writing data" << endl;
      writeResponseFile(countFilename, regions, formatString);
    // or for regression dists
    } else if ((distributionType == "PoissonRegression") ||
             (distributionType == "NegativeBinomialRegression") ||
             (distributionType == "ZeroTruncatedPoissonRegression") ||
             (distributionType == "ZeroTruncatedNegativeBinomialRegression")) {
      // randomly generate the coefficients
      vector<double> coefficients;
      cerr << "selected coefficients: ";
      for (size_t i=0; i<covariateFilenamesV.size(); i++) {
        double rmax = (double) RAND_MAX;
        double v = MAX_COEFF * rand() / rmax;
        cerr << v;
        if (i != covariateFilenamesV.size() - 1) cerr << ", ";
        else cerr << endl;
        coefficients.push_back(v);
      }

      // randomly generate the covariates
      vector< vector<double> > covariateVals;
      for (size_t i=0; i<regions.size(); i++) {
        covariateVals.push_back(vector<double>());
        for (size_t j=0; j<covariateFilenamesV.size(); j++) {
          double rmax = (double) RAND_MAX;
          covariateVals[i].push_back(MAX_COVAR_VAL * rand() / rmax);
        }
      }

      simulateScoresRegression(regions, covariateVals, coefficients,
                               distributionType, MAX_MU, MAX_ALPHA);
      writeResponseFile(countFilename, regions, formatString);

      for (size_t i=0; i<covariateFilenamesV.size(); i++) {
        ofstream outfh(covariateFilenamesV[i].c_str());
        for (size_t j=0; j<regions.size(); j++) {
          regions[j].set_score(covariateVals[j][i]);
          outfh << regions[j] << endl;
        }
      }
    }
  } catch (const SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
