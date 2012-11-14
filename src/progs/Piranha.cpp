/**
 * \file  Piranha.cpp
 * \brief Piranha is a program to identify peaks in read count profiles
 *        produced from CLIP/RIP-seq experiments.
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

// From piranha common
#include "config.hpp"
#include "NegativeBinomialRegression.hpp"
#include "TruncatedPoissonRegression.hpp"
#include "PoissonRegression.hpp"
#include "ZTNBRegression.hpp"
//#include "CovariateList.hpp"
#include "DistributionMixtureModel.hpp"
#include "RegressionMixtureModel.hpp"
#include "FittingMethod.hpp"
#include "ZTNB.hpp"
#include "RegressionBuilder.hpp"
#include "ReadBinner.hpp"

// From smithlab_cpp
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"

// From BAMTools --  if we have BAM support
#ifdef BAM_SUPPORT
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#endif

using std::vector;
using std::string;
using std::endl;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::stringstream;
using std::tr1::unordered_map;
using std::transform;

// if we have BAM support
#ifdef BAM_SUPPORT
using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;
#endif

#ifdef BAM_SUPPORT
enum inputType{BED, BAM};
#endif
#ifdef NO_BAM_SUPPORT
enum inputType{BED};
#endif

const size_t ALREADY_BINNED = 0;

/**
 * \brief functor for sorting genomic regions
 */
struct regionComp {
  bool operator()(const GenomicRegion& a,
          const GenomicRegion& b) const {
    return a < b;
  }
};

/**
 * \brief returns true if the two GenomicRegions cover exactly the same genomic
 *        location, cares only about chrom, start and end.
 */
static bool
sameRegion(const GenomicRegion &r1, const GenomicRegion &r2) {
  if ((r1.get_chrom() != r2.get_chrom()) ||
      (r1.get_start() != r2.get_start()) ||
      (r1.get_end() != r2.get_end())) {
        return false;
  }
  return true;
}

#ifdef BAM_SUPPORT
/**
 * \brief convert a BamAlignment object to a GenomicRegion object
 * \param chrom_lookup[in]   map of reference id to chrom name from BAM file
 * \param ba[in]             the BamAlignment
 * \param r[out]             this GenomicRegion will be populated to match ba
 */
static void
BamAlignmentToGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
             const BamAlignment &ba, GenomicRegion &r) {

  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));

  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  const string name(ba.Name);
  const float score = ba.MapQuality;
  const char strand = (ba.IsReverseStrand() ? '-' : '+');
  const string seq = ba.QueryBases;
  const string scr = ba.Qualities;
  r = GenomicRegion(chrom, start, end, name, score, strand);
}
#endif

/**
 * \brief guess whether this is a BAM file or a BED file -- doesn't check
 *        whether the file exists or not. Very basic check of file extension.
 * \param filename of the file to check
 */
static inputType
guessFileType(const string &filename, const bool verbose = false) {
  string tmp = smithlab::strip(filename);
  string ext = tmp.substr(tmp.size()-3);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

  if ((verbose) && (ext == "BAM"))
      cerr << "guessed that " << filename << " is BAM format" << endl;
  else if (verbose)
      cerr << "guessed that " << filename << " is BED format" << endl;

// only guess BAM format if we have compiled with BAM support
#ifdef BAM_SUPPORT
  if (ext == "BAM") return BAM;
#endif

  return BED;
}



/**
 * \brief load the response counts from the given filename
 * \param filename[in]   the filename that contains the input in either BAM
 *                       or BED format
 * \param response[out]  we'll put the response counts into this vector
 * \param regions[out]   we'll put the bins (genomicRegions) into this vector
 * \param binSize[in]    if we need to bin the reads, what is the bin size?
 *                       Can be set to ALREADY_BINNED if no binning should be
 *                       done (default).
 * \param SORT[in]       The input needs to be sorted; if this is true, we
 *                       sort the input ourselves (faster to leave this false
 *                       and just have the regions pre-sorted though)
 * \param verbose[in]    if true, print status messages to stderr
 */
static void
loadResponses(const string &filename, vector<double> &response,
              vector<GenomicRegion> &regions,
              const size_t binSize = ALREADY_BINNED,
              const bool SORT = false,
              const bool verbose = false) {
  regions.clear();
  response.clear();
  inputType type = guessFileType(filename, verbose);

#ifdef BAM_SUPPORT
  // Cry on BAM input with already binned reads -- can't fit the read counts
  // into the map quality attribute in the SAM format, so we can't do this.
  if ((type == BAM) && (binSize == ALREADY_BINNED)) {
    stringstream ss;
    ss << "Responses cannot be provided in BAM format if they are not raw "
       << "reads. Did you mean to use the -b option?";
    throw SMITHLABException(ss.str());
  }
#endif

  // first we load all the data into the vector of regions...
  if (type == BED) {
    ReadBEDFile(filename, regions);
  }
#ifdef BAM_SUPPORT
  else if (type == BAM) {
    size_t skippedInvalidChromID = 0;
    BamReader reader;
    // TODO -- need to check that filename is valid; BAMTools seems not to
    reader.Open(filename);

    // Get header and reference
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    unordered_map<size_t, string> chrom_lookup;
    for (size_t i = 0; i < refs.size(); ++i)
      chrom_lookup[i] = refs[i].RefName;

    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
      if (bam.RefID == -1) {
        skippedInvalidChromID += 1;
      } else {
        GenomicRegion r;
        BamAlignmentToGenomicRegion(chrom_lookup, bam, r);
        regions.push_back(r);
      }
    }

    if (skippedInvalidChromID > 0)
      cerr << "WARNING: Skipped " << skippedInvalidChromID << " in BAM file "
           << filename << " due to unknown reference IDs" << endl;

  }
#endif

  // do we need to sort things?
  if (SORT) {
    sort(regions.begin(), regions.end(), regionComp());
  }
  if (!check_sorted(regions)) {
    stringstream ss;
    ss << filename << " is not sorted. Required sort order is chrom, end, "
          "start, strand. You must either sort them yourself, or use the "
          "-s option";
    throw SMITHLABException(ss.str());
  }

  // now, do we need to bin the reads?
  if (binSize != ALREADY_BINNED) {
    vector<GenomicRegion> binned;
    ReadBinner b(binSize);
    b.binReads(regions, binned);
    swap(binned, regions);
  }

  // put the scores into the vector and clear them from the regions
  for (size_t i=0; i<regions.size(); i++) {
    response.push_back(regions[i].get_score());
    regions[i].set_score(0);
  }
}



/**
 * \brief load covariates for Piranha from the given filenames.
 * \param filenames[in]   filenames of the input files that contain the
 *                        covariate values and bins -- can be BED format
 * \param sites[in]       the response sites -- these will be checked to make
 *                        sure they match up with the bins defined in the
 *                        covariate files.
 * \param covariates[out] A matrix M which will be filled in such that M[i][j]
 *                        is the value of the ith covariate for the jth bin.
 * \param NORMALISE[in]   if true (the default), we normalise all covariate
 *                        values so they are between zero and one (exclusive).
 * \param SORT[in]        The input needs to be sorted; if this is true, we
 *                        sort the input ourselves (faster to leave this false
 *                        and just have the regions pre-sorted though)
 * \param VERBOSE[in]     if true, print status messages to stderr
 * \throws SMITHLABException if there is any response bin missing a covariate
 *                           bin.
 */
static void
loadCovariates(const vector<string> &filenames,
               const vector<GenomicRegion> &sites,
               vector<vector<double> > &covariates,
               const bool SORT = false,
               const bool NORMALISE_COVARS=true,
               const bool VERBOSE=false) {
  if (filenames.size() <= 0) {
    stringstream ss;
    ss << "failed to load covariates. Reason: didn't get any input filenames";
    throw SMITHLABException(ss.str());
  }

  // load covariates
  covariates.resize(filenames.size());
  for (size_t i = 0; i < covariates.size(); i++) {
    // first load the covariate file
    vector<GenomicRegion> covTmp;
    ReadBEDFile(filenames[i], covTmp);

    // do we need to sort things?
    if (SORT) {
      sort(covTmp.begin(), covTmp.end(), regionComp());
    }
    if (!check_sorted(covTmp)) {
      stringstream ss;
      ss << filenames[i] << " is not sorted. Required sort order is chrom, "
         << "end, start, strand. You must either sort them yourself, or use "
         << "the -s option";
      throw SMITHLABException(ss.str());
    }

    // tell the user how many values we found...
    if (VERBOSE)
      cerr << "loaded " << covTmp.size() << " elements from "
           << filenames[i] << endl;

    // cry if it doesn't match what we expect
    if (covTmp.size() < sites.size()) {
      stringstream ss;
      ss << "failed loading covariate from " << filenames[i] << ". Reason: "
         << "expected " << sites.size() << " regions, but found "
         << covTmp.size() << endl;
      throw SMITHLABException(ss.str());
    }

    // pull out the values and match up the regions
    size_t covIndex = 0;
    for (size_t j=0; j<sites.size(); j++) {
      while ((covIndex < covTmp.size()) && (covTmp[covIndex] < sites[j])) {
        covIndex += 1;
      }
      if ((covIndex >= covTmp.size()) || (sites[j] < covTmp[covIndex])) {
        stringstream ss;
        ss << "failed loading covariate from " << filenames[i] << ". Reason: "
           << "could not find a value to associate with the response region "
           << sites[j];
        throw SMITHLABException(ss.str());
      }
      else if (sameRegion(covTmp[covIndex], sites[j])) {
        covariates[i].push_back(covTmp[covIndex].get_score());
      }
    }
  }

  // normalise covariates?
  if (NORMALISE_COVARS) {
    for (size_t i = 0; i < covariates.size(); i++) {
      double maxVal = *std::max_element(covariates[i].begin(),
                                        covariates[i].end());
      double minVal = *std::min_element(covariates[i].begin(),
                                        covariates[i].end());
      for (size_t j = 0; j < covariates[i].size(); j++) {
        covariates[i][j] = (covariates[i][j] - minVal) / (maxVal - minVal);
        // avoid having anything at exactly 0 or 1
        if (covariates[i][j] == 0) covariates[i][j] = 0.01;
        if (covariates[i][j] == 1) covariates[i][j] = 0.99;
      }
    }
  }
}

/**
 * \brief Use a single regression model to find peaks and output each
 *        input region with a p-value
 */
static void
FindPeaksSingleComponentRegression(const bool VERBOSE, const bool FITONLY,
                                   const vector<GenomicRegion> &sites,
                                   const vector<double> &responses,
                                   const vector< vector<double> > &covariates,
                                   const FittingMethod ftmthd,
                                   const string &distType,
                                   const string& modelfn, ostream& ostrm) {
  if (covariates.size() < 1) {
    stringstream ss;
    ss << "Peak finding using regression method " << distType << "failed. "
       << "Reason: No covariates found";
    throw SMITHLABException(ss.str());
  }

  // a little sanity check on what was passed to the function, make sure
  // all the dimensions of vectors make sense.
  const size_t n_sites = responses.size();
  if (n_sites != sites.size()) {
    stringstream ss;
    ss << "Failed peak finding using " << distType << ". Reason: number "
       << "of responses doesn't match number of sites.";
    throw SMITHLABException(ss.str());
  }
  for(size_t i = 0; i < covariates.size(); i++) {
    if (covariates[i].size() != n_sites) {
      stringstream ss;
      ss << "Failed peak finding using " << distType << ". Reason: number "
         << "of values for covariate " << i << " does not match number of "
         << "sites";
      throw SMITHLABException(ss.str());
    }
  }

  // now we're going to flip the covariates matrix
  Matrix m(covariates);
  const vector< vector<double> > covariates_t =
    m.transpose().asVectorOfVector();

  // Load or fit the model
  RegressionModel *distro =
    RegressionBuilder::buildRegression(distType, covariates_t.size(), ftmthd);
  if (modelfn.empty())
    distro->estimateParams(responses, covariates_t, VERBOSE);
  else distro->load(modelfn);

  // Give our output as the scored sites, or just the model?
  if (FITONLY) {
    distro->save(ostrm);
  } else {
    for (size_t i=0; i < n_sites; i++) {
      GenomicRegion tmp(sites[i]);
      tmp.set_score(responses[i]);
      ostrm << tmp << '\t' << distro->pvalue(responses[i], covariates_t[i])
            << endl;
    }
  }

  // cleanup and we're done
  delete distro;
}


/**
 * \brief Use a single simple distribution to find peaks and output each
 *        input region with a p-value
 * \param VERBOSE   output additional run details to stderr if true
 * \param FITONLY   If true, just output the model, otherwise we
 *                  use the model to score the sites and output those
 * \param distType  the distribution type to use
 * \param modelfn   if not empty, load the model from this file rather than
 *                  fitting
 * \param sites     genomic regions representing the bins being considered
 * \param responses read counts for each of the sites
 * \param outfn     write regions to this output stream
 */
static void 
FindPeaksSingleComponentSimple(const bool VERBOSE, const bool FITONLY,
                               const string &distType, const string &modelfn,
                               const vector<GenomicRegion> &sites,
                               const vector<double> &responses,
                               ostream& ostrm) {

  if (VERBOSE)
    cerr << "Simple 1-component fitting using " << distType << endl;

  const size_t n_sites = responses.size();
  if (n_sites != sites.size()) {
    stringstream ss;
    ss << "Failed peak finding using " << distType << ". Reason: number "
       << "of responses doesn't match number of sites.";
    throw SMITHLABException(ss.str());
  }
  
  // Load or fit the model
  Distribution *distro = Distribution::create(distType);
  if (modelfn.empty())
    distro->estimateParams(responses, VERBOSE);
  else distro->load(modelfn);
  
  // Give our output as the scored sites, or just the model?
  if (FITONLY) {
    distro->save(ostrm);
  } else {
    for (size_t i=0; i < n_sites; i++) {
      GenomicRegion tmp(sites[i]);
      tmp.set_score(responses[i]);
      ostrm << tmp << '\t' << distro->pvalue(responses[i]) << endl;
    }
  }

  // cleanup
  delete distro;
}


/**
 * \brief TODO
 */
int 
main(int argc, const char* argv[]) {
  try {
    string modelfn = "", outfn;
    bool VERBOSE = false;
    bool FITONLY = false;
    bool NO_NORMALISE_COVARS = false;
    bool SORT = false;

    size_t numComponents = 1;
    string distType = "DEFAULT";
    string fittingMethodStr = defaultFittingMethod.toString();
    size_t binSize = ALREADY_BINNED;
    
    /****************** COMMAND LINE OPTIONS ********************/
    string about = "Piranha Version 1.1.3 -- A program for finding peaks in "
                   "high throughput RNA-protein interaction data (e.g. RIP- "
                   "and CLIP-Seq). Written by Philip J. Uren and "
                   "Andrew D. Smith. See README.TXT for further details";
    OptionParser opt_parse(strip_path(argv[0]), about, "[*.bed] *.bed");
    opt_parse.add_opt("output", 'o', "Name of output file, STDOUT if omitted", 
                      false, outfn);
    opt_parse.add_opt("sort", 's', "indicates that input is unsorted and "
                      "Piranha should sort it for you", false, SORT);
    opt_parse.add_opt("bin_size", 'b', "indicates that input is raw reads and "
                      "should be binned into bins of this size",
                      false, binSize);
    opt_parse.add_opt("fit", 'f', "Fit only, output model to file", 
                      false, FITONLY);
    opt_parse.add_opt("dist", 'd', "Distribution type. Currently supports "
                                   "Poisson, NegativeBinomial, "
                                   "ZeroTruncatedPoisson, "
                                   "ZeroTruncatedNegativeBinomial "
                                     "(default with no covariates), "
                                   "PoissonRegression, "
                                   "NegativeBinomialRegression, "
                                   "ZeroTruncatedPoissonRegression, "
                                   "ZeroTruncatedNegativeBinomialRegression "
                                     "(default with covariates)",
                      false, distType);
    opt_parse.add_opt("fitMethod", 't', "component fitting method",
                      false, fittingMethodStr);
    opt_parse.add_opt("model", 'm', "Use the specified model file instead of "
                      "fitting to input data", false, modelfn);
    opt_parse.add_opt("VERBOSE", 'v', "output additional messages about run "
                      "to stderr if set", false, VERBOSE);
    opt_parse.add_opt("no_normalisation", 'n', "don't normalise covariates",
                      false, NO_NORMALISE_COVARS);

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
    if (leftover_args.empty()) {
      // TODO handle input via stdin?
      cerr << "must provide input files" << endl;
      return EXIT_FAILURE;
    }
    /****************** END COMMAND LINE OPTIONS *****************/
    
    // warn the user that this is a little silly, but let them do it anyway
    if ((modelfn != "") && (FITONLY)) {
      cerr << "warning: loading model from file, and not scoring input -- "
	   << "will just copy model from input to output streams" << endl;
    }

    // decide where our output is going, stdout or file?
    std::ofstream of;
    if (!outfn.empty()) of.open(outfn.c_str());
    ostream ostrm(outfn.empty() ? cout.rdbuf() : of.rdbuf());

    // go ahead and load the input files
    vector<GenomicRegion> sites;
    vector<double> responses;
    vector< vector<double> > covariates;
    loadResponses(leftover_args.front(), responses, sites, binSize, SORT);
    if (leftover_args.size() > 1) {
      vector<string> cfnames =
          vector<string>(leftover_args.begin() + 1, leftover_args.end());
      loadCovariates(cfnames, sites, covariates, SORT, !NO_NORMALISE_COVARS);
    }

    // tell the user what options were set and what files we found
    if (VERBOSE) {
      cerr << "loaded " << responses.size() <<  " elements from "
           << leftover_args.front() << endl;
      if (leftover_args.size() > 1) {
        vector<string> cfnames =
            vector<string>(leftover_args.begin() + 1, leftover_args.end());
        cerr << "loaded " << covariates.size() << " covariates from ";
        for (size_t i=0; i<cfnames.size(); i++) cerr << cfnames[i] << "\t";
        cerr << endl;
      } else cerr << "no covariates found" << endl;
      if (!outfn.empty()) cerr << "writing output to " << outfn << endl;
      else cerr << "writing output to stdout" << endl;
      if (binSize == ALREADY_BINNED) cerr << "input is already binned" << endl;
      else cerr << "binning input into bins of size " << binSize << endl;
      if (FITONLY) cerr << "Not scoring input regions" << endl;
      else cerr << "scoring input regions" << endl;
      if (modelfn != "") cerr << "loading model from " << modelfn << endl;
      cerr << "model type? " << distType << endl;
      if (NO_NORMALISE_COVARS) cerr << "normalise covariates? no" << endl;
      else cerr << "normalise covariates? yes" << endl;
      if (SORT) cerr << "sort input files? yes" << endl;
      else cerr << "sort input files? no" << endl;
    }

    // if we're fitting only a single component, rather than a mixture..
    if (numComponents == 1) {
      if (leftover_args.size() == 1) {
        // just one file given, must be simple distribution
        if (distType == "DEFAULT")
          distType = "ZeroTruncatedNegativeBinomial";
        FindPeaksSingleComponentSimple(VERBOSE, FITONLY, distType, modelfn,
                                       sites, responses, ostrm);
      } else {
        // more than one input file given, must be regression
        if (distType == "DEFAULT")
          distType = "ZeroTruncatedNegativeBinomialRegression";
        FindPeaksSingleComponentRegression(VERBOSE, FITONLY,
                                           sites, responses, covariates,
                                           FittingMethod(fittingMethodStr),
                                           distType, modelfn, ostrm);
      }
    } else {
      // do the full mixture thing
      cerr << "mixtures not supported in this version" << endl;
      /*FindPeaksMixture(leftover_args, ostrm, modelfn, distType, numComponents,
		       FITONLY, FittingMethod(fittingMethodStr),
		       ignoreCovariates, !NO_NORMALISE_COVARS, VERBOSE);*/
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
