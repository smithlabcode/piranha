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
#include <map>

// From piranha common
#include "config.hpp"
#include "NegativeBinomialRegression.hpp"
#include "TruncatedPoissonRegression.hpp"
#include "PoissonRegression.hpp"
#include "ZTNBRegression.hpp"
#include "Utility.hpp"
#include "DistributionMixtureModel.hpp"
#include "RegressionMixtureModel.hpp"
#include "FittingMethod.hpp"
#include "ZTNB.hpp"
#include "RegressionBuilder.hpp"
#include "GenomicRegionAggregator.hpp"
#include "ReadBinner.hpp"
#include "FDR.hpp"

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
using std::map;
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
enum aggregationType {NO_AGGREGATION, SUMMITS_AGGREGATION,
                      CLUSTERS_AGGREGATION};

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
 * \brief TODO
 * \param allResponses[in/out]
 * \param foregroundResponses[out]
 */
static void
splitResponses(vector<double> &allResponses,
                    vector<GenomicRegion> &allSites,
                    vector<double> &foregroundResponses,
                    vector<GenomicRegion> &foregroundSites,
                    const double thresh,
                    const bool VERBOSE) {
  // do nothing if the split puts everything into the background
  if (thresh != 1) {
    const size_t N = allResponses.size();
    if (N != allSites.size()) {
      stringstream ss;
      ss << "Failed to split responses, sites and responses have different "
            "dimension";
      throw SMITHLABException(ss.str());
    }

    // here we calculate what is the smallest response value such that <thresh>
    // proportion of the responses are less than it (i.e. background)
    double absThresh = thresh * N, soFar = 0, responseThresh=0;
    map<double, size_t> hist = countItemsOrdered<double>(allResponses);
    map<double, size_t>::iterator prev = hist.end();
    for (map<double, size_t>::iterator v = hist.begin(); v != hist.end(); ++v) {
      soFar += v->second;
      if (soFar >= absThresh) {
        if (prev == hist.end()) {
          stringstream ss;
          ss << "Failed to split responses, smallest response accounts for more "
             << "than " << (thresh * 100) << "% of all responses. Try "
             <<  "increasing the threshold";
          throw SMITHLABException(ss.str());
        }
        responseThresh = prev->first;
        break;
      }
      prev = v;
    }
    if (VERBOSE) cerr << "Selected foreground count-threshold as "
                      << responseThresh << endl;

    // now build a vector of the indices that go into the foreground
    std::tr1::unordered_set<size_t> fg_indices;
    for (size_t i=0; i<allResponses.size(); i++) {
      if (allResponses[i] > responseThresh) fg_indices.insert(i);
    }

    // new we split the responses and the sites based on these indices
    copyRemoveByIndex<double>(allResponses, foregroundResponses,
                              fg_indices);
    copyRemoveByIndex<GenomicRegion>(allSites, foregroundSites,
                                     fg_indices);
  }

  if (VERBOSE) cerr << "fg/bg split was " << foregroundResponses.size() << "/"
                    << allResponses.size() << endl;
}


/**
 * \brief TODO
 * \param ...
 */
static void
splitResponsesAndCovariates(vector<double> &allResponses,
                            vector< vector<double> > &allCovariates,
                            vector<GenomicRegion> &allSites,
                            vector<double> &foregroundResponses,
                            vector< vector<double> > &foregroundCovariates,
                            vector<GenomicRegion> &foregroundSites,
                            const double thresh,
                            const bool VERBOSE) {
  // do nothing if the split puts everything into the background
  if (thresh != 1) {
    const size_t N = allResponses.size();
    if (N != allSites.size()) {
      stringstream ss;
      ss << "Failed to split responses and covariates, sites and responses "
            "have different dimension";
      throw SMITHLABException(ss.str());
    }

    // here we calculate what is the smallest response value such that <thresh>
    // proportion of the responses are less than it (i.e. background)
    double absThresh = thresh * N, soFar = 0, responseThresh=0;
    map<double, size_t> hist = countItemsOrdered<double>(allResponses);
    map<double, size_t>::iterator prev = hist.end();
    for (map<double, size_t>::iterator v = hist.begin(); v != hist.end(); ++v) {
      soFar += v->second;
      if (soFar >= absThresh) {
        if (prev == hist.end()) {
          stringstream ss;
          ss << "Failed to split responses, smallest response accounts for more "
             << "than " << (thresh * 100) << "% of all responses. Try "
             <<  "increasing the threshold";
          throw SMITHLABException(ss.str());
        }
        responseThresh = prev->first;
        break;
      }
      prev = v;
    }
    if (VERBOSE) cerr << "Selected foreground count-threshold as "
                      << responseThresh << endl;

    // now build a vector of the indices that go into the foreground
    std::tr1::unordered_set<size_t> fg_indices;
    for (size_t i=0; i<allResponses.size(); i++) {
      if (allResponses[i] > responseThresh) fg_indices.insert(i);
    }

    // new we split the responses and the sites based on these indices
    copyRemoveByIndex<double>(allResponses, foregroundResponses,
                              fg_indices);
    copyRemoveByIndex<GenomicRegion>(allSites, foregroundSites,
                                     fg_indices);
    for (size_t i = 0; i < allCovariates.size(); i++) {
      foregroundCovariates.push_back(vector<double>());
      copyRemoveByIndex<double>(allCovariates[i], foregroundCovariates[i],
                                fg_indices);
    }
  }

  if (VERBOSE) cerr << "fg/bg split was " << foregroundResponses.size() << "/"
                    << allResponses.size() << endl;
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
 * \param binSize[in]     if we need to bin the covariate, what is the bin size?
 *                        This will cause the covariates to be treated as a
 *                        set of reads basically. Can be set to ALREADY_BINNED
 *                        if no binning should be done (default).
 * \param SORT[in]        The input needs to be sorted; if this is true, we
 *                        sort the input ourselves (faster to leave this false
 *                        and just have the regions pre-sorted though)
 * \param LOG[in]         transform covariate values into log (must not
 *                        contain any zeros)
 * \param NORMALISE[in]   if true (the default), we normalise all covariate
 *                        values so they are between zero and one (exclusive).
 * \param VERBOSE[in]     if true, print status messages to stderr
 *
 * \throws SMITHLABException if there is any response bin missing a covariate
 *                           bin.
 */
static void
loadCovariates(const vector<string> &filenames,
               const vector<GenomicRegion> &sites,
               vector<vector<double> > &covariates,
               const size_t binSize = ALREADY_BINNED,
               const bool SORT=false,
               const bool LOG=false,
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

    // now, do we need to bin the covariate?
    // if we're binning, we'll add missing bins as well
    if (binSize != ALREADY_BINNED) {
      cerr << "binning" << endl;
      vector<GenomicRegion> binned;
      ReadBinner b(binSize);
      b.binReads(covTmp, binned, sites, 1);
      swap(binned, covTmp);
    }

    // cry if the sizes don't match what we expect -- we assume from here
    // on in that there are at least as many covariate values as sites.
    // TODO: I'm not sure this is needed, given the checks in the block of
    //       code below.
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

  if (LOG) {
    for (size_t i = 0; i < covariates.size(); i++) {
      for (size_t j = 0; j < covariates[i].size(); j++) {
        if (covariates[i][j] == 0) {
          stringstream ss;
          ss << "can't take log of covariates, some values at zero";
          throw SMITHLABException(ss.str());
        }
        // TODO fix magic number
        else if (covariates[i][j] == 1) covariates[i][j] = 0.00001;
        else {
          covariates[i][j] = log(covariates[i][j]);
          // avoid having anything at exactly 0 or 1
          if (covariates[i][j] == 0) covariates[i][j] = 0.01;
          if (covariates[i][j] == 1) covariates[i][j] = 0.99;
        }
      }
    }
  }
}


/**
 * \brief TODO
 */
static void
outputSites(const vector<double> &fgResponses,
                    const vector<double> &bgResponses,
                    const vector<GenomicRegion> &fgSites,
                    const vector<GenomicRegion> &bgSites,
                    const vector<double> &fgPvals,
                    const vector<double> &bgPvals,
                    const double pThresh,
                    aggregationType aggMethod,
                    ostream &ostrm) {
  // we first build a vector of the signfiicant sites
  // because we want to post-process to merge nearby sites or identify
  // cluster summits and this makes it easier.
  vector<GenomicRegion> sigSites;
  vector<double> sig_ps;
  size_t fg_index = 0, bg_index = 0;
  while ((fg_index < fgSites.size()) || (bg_index < bgSites.size())) {
    if ((fg_index >= fgSites.size()) ||
        (bgSites[bg_index] < fgSites[fg_index])) {
      // background site comes first; either because we ran out of
      // foreground sites or because of ordering. Score and output it.
      const double p = bgPvals[bg_index];
      if (p <= pThresh) {
        GenomicRegion tmp(bgSites[bg_index]);
        tmp.set_score(bgResponses[bg_index]);
        ostrm << tmp << '\t' << p << endl;
      }
      bg_index += 1;
    } else {
      // foreground site comes first
      const double p = fgPvals[fg_index];
      if (p <= pThresh) {
        GenomicRegion tmp(fgSites[fg_index]);
        tmp.set_score(fgResponses[fg_index]);
        sigSites.push_back(tmp);
        sig_ps.push_back(p);
      }
      fg_index += 1;
    }
  }

  // swap out the scores for the pvals, do the aggregation, then put them back.


  // TODO pass agg method and distance, cry if agg method is none and
  // distance isn't zero
  // TODO visitor pattern?
  // TODO would be __much__ better if we didn't have to copy these first...
  if (aggMethod == NO_AGGREGATION) {
    for (size_t i=0; i<sigSites.size(); i++)
      ostrm << sigSites[i] << '\t' << sig_psp << endl;
  } else if (aggMethod == SUMMITS_AGGREGATION) {
    GenomicRegionAggregator(10).aggregate(sigSites.begin(), sigSites.end(),
         ClusterSummitPrinter(ostrm, ClusterSummitPrinter::CLUSTER_MIN_SCORE));
  } else if (aggMethod == CLUSTERS_AGGREGATION) {
    GenomicRegionAggregator(10).aggregate(sigSites.begin(), sigSites.end(),
         ClusterLimitsPrinter(ostrm));
  }


}




/**
 * \brief Use a single regression model to find peaks and output each
 *        input region with a p-value
 * \param VERBOSE       TODO
 * \param FITONLY       TODO
 * \param sites         TODO
 * \param responses     TODO
 * \param covariates    TODO
 * \param ftmthd        TODO
 * \param distType      TODO
 * \param modelfn       TODO
 * \param pThresh       output sites that meet this threshold on significance
 * \param ostrm         TODO
 */
static void
FindPeaksSingleComponentRegression(const bool VERBOSE, const bool FITONLY,
                                   const bool NO_PVAL_CORRECT,
                                   const vector<GenomicRegion> &sites,
                                   const vector<GenomicRegion> &fgSites,
                                   const vector<double> &responses,
                                   const vector<double> &fgResponses,
                                   const vector< vector<double> > &covariates,
                                   const vector< vector<double> > &fgCovariates,
                                   const FittingMethod ftmthd,
                                   const string &distType,
                                   const string& modelfn, const double pThresh,
                                   ostream& ostrm) {
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

  // now we're going to flip the covariates matrix for the bg sites (we'll
  // use these to fit our model)
  Matrix m(covariates);
  const vector< vector<double> > covariates_t =
    m.transpose().asVectorOfVector();
  // Load or fit the model
  RegressionModel *distro =
    RegressionBuilder::buildRegression(distType, covariates_t.size(), ftmthd);
  if (modelfn.empty()) distro->estimateParams(responses, covariates_t, VERBOSE);
  else distro->load(modelfn);

  // Give our output as the scored sites, or just the model?
  if (FITONLY) {
    distro->save(ostrm);
  } else {
    // we first compute all the p-values, because we probably want to correct
    // them and to do that we need them all.
    vector<double> fg_pvals, pvals;
    for (size_t i=0; i < n_sites; i++)
      pvals.push_back(distro->pvalue(responses[i], covariates_t[i]));
    if (fgCovariates.size() > 0) {
      Matrix m_fg(fgCovariates);
      const vector< vector<double> > covariates_t_fg =
        m_fg.transpose().asVectorOfVector();
      for (size_t i=0; i < fgResponses.size(); i++)
        fg_pvals.push_back(distro->pvalue(fgResponses[i], covariates_t_fg[i]));
    }
    if (!NO_PVAL_CORRECT) {
      FDR::correctP(fg_pvals);
      FDR::correctP(pvals);
    }
    outputSites(fgResponses, responses, fgSites, sites, fg_pvals, pvals,
                pThresh, ostrm);
  }

  // cleanup and we're done
  delete distro;
}


/**
 * \brief Use a single simple distribution to find peaks and output each
 *        input region with a p-value
 * \param VERBOSE          output additional run details to stderr if true
 * \param FITONLY          If true, just output the model, otherwise we
 *                         use the model to score the sites and output those
 * \param NO_PVAL_CORRECT  don't correct p-values for multiple hypothesis
 *                         testing. By default we correct using B&H
 * \param distType         the distribution type to use
 * \param modelfn          if not empty, load the model from this file rather
 *                         than fitting
 * \param sites            genomic regions representing the bins being considered
 * \param fgSites          TODO
 * \param responses        read counts for each of the sites
 * \param fgResponses      TODO
 * \param pThresh     output sites that meet this threshold on significance
 * \param ostrm       write regions to this output stream
 */
static void 
FindPeaksSingleComponentSimple(const bool VERBOSE, const bool FITONLY,
                               const bool NO_PVAL_CORRECT,
                               const string &distType, const string &modelfn,
                               const vector<GenomicRegion> &sites,
                               const vector<GenomicRegion> &fgSites,
                               const vector<double> &responses,
                               const vector<double> &fgResponses,
                               const double pThresh,
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
    // we have to pre-compute all of the p-values so we can adjust them
    // for multiple hypothesis testing
    vector<double> fg_pvals, bg_pvals;
    for (size_t i=0; i<responses.size(); i++)
      bg_pvals.push_back(distro->pvalue(responses[i]));
    for (size_t i=0; i<fgResponses.size(); i++)
      fg_pvals.push_back(distro->pvalue(fgResponses[i]));
    if (!NO_PVAL_CORRECT) {
      FDR::correctP(fg_pvals);
      FDR::correctP(bg_pvals);
    }

    // output the scored sites
    outputSites(fgResponses, responses, fgSites, sites, fg_pvals, bg_pvals,
                pThresh, ostrm);
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
    const double DEFAULT_PTHRESH = 0.01;
    const double DEFAULT_BGTHRESH = 0.99;

    string modelfn = "", outfn;
    bool VERBOSE = false;
    bool FITONLY = false;
    bool NO_NORMALISE_COVARS = false;
    bool SORT = false;
    bool NO_PVAL_CORRECT = false;
    bool LOG_COVARS = false;

    size_t numComponents = 1;
    string distType = "DEFAULT";
    double pthresh = DEFAULT_PTHRESH;
    double bgThresh = DEFAULT_BGTHRESH;
    string fittingMethodStr = defaultFittingMethod.toString();
    size_t binSize_response = ALREADY_BINNED;
    size_t binSize_covars = ALREADY_BINNED;
    size_t binSize_both = ALREADY_BINNED;
    
    /************************* COMMAND LINE OPTIONS **************************/
    string about = "Piranha Version 1.2.0 -- A program for finding peaks in "
                   "high throughput RNA-protein interaction data (e.g. RIP- "
                   "and CLIP-Seq). Written by Philip J. Uren and "
                   "Andrew D. Smith. See README.TXT for further details";
    OptionParser opt_parse(strip_path(argv[0]), about, "[*.bed] *.bed");
    opt_parse.add_opt("output", 'o', "Name of output file, STDOUT if omitted", 
                      false, outfn);
    opt_parse.add_opt("sort", 's', "indicates that input is unsorted and "
                      "Piranha should sort it for you", false, SORT);
    opt_parse.add_opt("p_threshold", 'p', "significance threshold for sites",
                      false, pthresh);
    opt_parse.add_opt("no_pval_correct", 'c', "don't correct p-values for "
                                              "multiple hypothesis testing. "
                                              "We correct by default using B&H.",
                      false, NO_PVAL_CORRECT);
    opt_parse.add_opt("background_thresh", 'a', "indicates that this "
                                                "proportion of the lowest "
                                                "scores should be considered "
                                                "the background. Default is "
                                                "0.99",
                      false, bgThresh);
    opt_parse.add_opt("bin_size_reponse", 'b', "indicates that the response "
                                       "(first input file) is raw reads and "
                                       "should be binned into bins of this size",
                      false, binSize_response);
    opt_parse.add_opt("bin_size_covars", 'i', "indicates that the covariates "
                                              "(all except first file) are "
                                              "raw reads and should be binned "
                                              "into bins of this size",
                      false, binSize_covars);
    opt_parse.add_opt("bin_size_both", 'z', "synonymous with -b x -i x for any x",
                      false, binSize_both);
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
    opt_parse.add_opt("log_covars", 'l', "convert covariates to log scale",
                      false, LOG_COVARS);

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
    
    // Check that the options the user picked all make sense
    if ((modelfn != "") && (FITONLY)) {
      cerr << "WARNING: loading model from file, and not scoring input -- "
	   << "will just copy model from input to output streams" << endl;
    }
    if ((bgThresh > 1) || (bgThresh < 0)) {
      stringstream ss;
      ss << "background threshold is out-of-range. Must be between 0 "
         << "and 1 (inclusive). Found: " << bgThresh;
      throw SMITHLABException(ss.str());
    }
    if ((pthresh > 1) || (pthresh < 0)) {
      stringstream ss;
      ss << "p-value threshold is out-of-range. Must be between 0 "
         << "and 1 (inclusive). Found: " << pthresh;
      throw SMITHLABException(ss.str());
    }
    if (binSize_both != ALREADY_BINNED) {
      if ((binSize_response != ALREADY_BINNED) ||
          (binSize_covars != ALREADY_BINNED)) {
        stringstream ss;
        ss << "the -z option is mutually exclusive with -b and -i options";
        throw SMITHLABException(ss.str());
      }
      binSize_response = binSize_both;
      binSize_covars = binSize_both;
    }
    if (LOG_COVARS && !NO_NORMALISE_COVARS) {
      cerr << "WARNING: normalisation and log-conversion of covariates are "
           << "mutually exclusive. Ignoring normalisation flag" << endl;
      NO_NORMALISE_COVARS = true;
    }

    // decide where our output is going, stdout or file?
    std::ofstream of;
    if (!outfn.empty()) of.open(outfn.c_str());
    ostream ostrm(outfn.empty() ? cout.rdbuf() : of.rdbuf());

    // go ahead and load the input files
    vector<GenomicRegion> sites;
    vector<double> responses;
    vector< vector<double> > covariates;
    loadResponses(leftover_args.front(), responses, sites,
                  binSize_response, SORT);
    if (leftover_args.size() > 1) {
      vector<string> cfnames =
          vector<string>(leftover_args.begin() + 1, leftover_args.end());
      loadCovariates(cfnames, sites, covariates, binSize_covars, SORT,
                     LOG_COVARS, !NO_NORMALISE_COVARS, VERBOSE);
    }

    // tell the user what options were set and what files we found
    if (VERBOSE) {
      cerr << "Piranha was run with the following options: " << endl;
      cerr << "Significance threshold: " << pthresh << endl;
      if (!outfn.empty()) cerr << "writing output to " << outfn << endl;
      else cerr << "writing output to stdout" << endl;

      // binning our data
      if (binSize_response == ALREADY_BINNED)
        cerr << "responses are already binned" << endl;
      else
        cerr << "binning responses into bins of size "
             << binSize_response << endl;
      if (binSize_covars == ALREADY_BINNED)
        cerr << "covariates are already binned" << endl;
      else
        cerr << "binning covariates into bins of size "
             << binSize_response << endl;

      if (FITONLY) cerr << "Not scoring input regions" << endl;
      else cerr << "scoring input regions" << endl;
      if (modelfn != "") cerr << "loading model from " << modelfn << endl;
      cerr << "model type? " << distType << endl;
      if (NO_NORMALISE_COVARS) cerr << "normalise covariates? no" << endl;
      else cerr << "normalise covariates? yes" << endl;
      if (SORT) cerr << "sort input files? yes" << endl;
      else cerr << "sort input files? no" << endl;
      cerr << "loaded " << responses.size() <<  " elements from "
           << leftover_args.front() << endl;
      if (leftover_args.size() > 1) {
        vector<string> cfnames =
            vector<string>(leftover_args.begin() + 1, leftover_args.end());
        cerr << "loaded " << covariates.size() << " covariates from ";
        for (size_t i=0; i<cfnames.size(); i++) cerr << cfnames[i] << "\t";
        cerr << endl;
      } else cerr << "no covariates found" << endl;
    }

    // if we're fitting only a single component, rather than a mixture..
    if (numComponents == 1) {
      if (leftover_args.size() == 1) {
        // just one file given, must be simple distribution
        if (distType == "DEFAULT")
          distType = "ZeroTruncatedNegativeBinomial";
        vector<GenomicRegion> fgSites;
        vector<double> fgResponses;
        splitResponses(responses, sites, fgResponses, fgSites, bgThresh,
                       VERBOSE);
        FindPeaksSingleComponentSimple(VERBOSE, FITONLY, NO_PVAL_CORRECT,
                                       distType, modelfn, sites, fgSites,
                                       responses, fgResponses, pthresh, ostrm);
      } else {
        // more than one input file given, must be regression
        if (distType == "DEFAULT")
          distType = "ZeroTruncatedNegativeBinomialRegression";
        vector<GenomicRegion> fgSites;
        vector<double> fgResponses;
        vector< vector<double> > fgCovariates;
        splitResponsesAndCovariates(responses, covariates, sites, fgResponses,
                                    fgCovariates, fgSites, bgThresh, VERBOSE);
        FindPeaksSingleComponentRegression(VERBOSE, FITONLY, NO_PVAL_CORRECT,
                                           sites, fgSites, responses,
                                           fgResponses, covariates, fgCovariates,
                                           FittingMethod(fittingMethodStr),
                                           distType, modelfn, pthresh, ostrm);
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
