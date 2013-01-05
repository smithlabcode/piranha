/**
  \file  GenomicRegionAggregator.hpp
  \brief This header defines functions for combining genomic regions

  \authors Philip J. Uren, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Andrew D. Smith

  \section license License Details
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  \section bugs Known Bugs

  \section history Revision History
**/


#ifndef GREGAGG_HPP
#define GREGAGG_HPP

#include <vector>
#include <string>
#include <iterator>
#include <numeric>
#include <list>

#include "GenomicRegion.hpp"

using namespace std;

class GenomicRegionAggregator {
public :
  /**
   * \brief Construct a cluster aggregator with the given distance
   */
  GenomicRegionAggregator(size_t clusterDist) : clusterDistance(clusterDist) {};

  /**
   *  \brief Aggregate a set of genomic regions where only some of them are
   *         considered significant. Only significant regions are allowed to
   *         extend clusters (if they fall within the required distance), but
   *         all bins within the cluster are available for the function f
   *         (which will be called on each cluster)
   *  \param f function to call on each cluster. Must take start and end
   *           iterators for the genomic regions, pvalues and covariates
   *           covered by the cluster
   */
  template<typename Func>
  void aggregate(
      const std::vector<GenomicRegion>::const_iterator r_start,
      const std::vector<GenomicRegion>::const_iterator r_end,
      const std::vector<double>::const_iterator p_start,
      const std::vector<double>::const_iterator p_end,
      const std::vector< vector<double>::const_iterator > c_starts,
      const std::vector< vector<double>::const_iterator > c_ends,
      const double alpha,
      const Func f) const {
    // just to cut down on typing..
    typedef std::vector<GenomicRegion>::const_iterator gIter;
    typedef std::vector<double>::const_iterator dIter;

    // check length and number of covariate iterators
    if (c_starts.size() != c_ends.size()) {
      std::stringstream ss;
      ss << "aggregating genomic regions failed -- didn't get the same number "
         << "of start iterators for covariates as we got end iterators";
      throw SMITHLABException(ss.str());
    }
    for (size_t i=1; i<c_starts.size(); i++) {
      if (std::distance(c_starts[0], c_ends[0]) !=
          std::distance(c_starts[i], c_ends[i])) {
        std::stringstream ss;
        ss << "aggregating genomic regions failed -- didn't get the same "
           << "number of covariates values for all covariates";
        throw SMITHLABException(ss.str());
      }
    }
    // Check that we have an equal number of regions, pvals and covariate
    // vectors
    if ((std::distance(r_start, r_end) != std::distance(p_start, p_end)) ||
        ((c_starts.size() > 0) &&
         (std::distance(r_start, r_end) != std::distance(c_starts[0], c_ends[0])))) {
      std::stringstream ss;
      ss << "aggregating genomic regions failed -- didn't get the same number "
         << "of regions as pvalues or covariate matrices";
      throw SMITHLABException(ss.str());
    }


    std::string clusterChrom = "";
    size_t clusterStartCoord = 0, clusterEndCoord = 0;

    gIter clusterStart_region = r_start, clusterEnd_region = r_start;
    gIter prev_region;
    dIter clusterStart_pval = p_start, clusterEnd_pval = p_end;
    std::vector<dIter> clusterStart_covars, clusterEnd_covars;
    for (size_t i=0; i<c_starts.size(); ++i) {
      clusterStart_covars.push_back(c_starts[i]);
      clusterEnd_covars.push_back(c_ends[i]);
    }

    gIter it_region = r_start;
    dIter it_pval = p_start;
    vector<dIter> it_covars;
    for (size_t i=0; i<c_starts.size(); ++i) it_covars.push_back(c_starts[i]);

    bool first = true;
    for (gIter it_region = r_start; it_region != r_end; ++it_region) {
      // make sure the regions are sorted
      if (!first) {
        if (it_region < prev_region) {
          std::stringstream ss;
          ss << "aggregating genomic regions failed -- regions were not "
             << "sorted: saw " << (*prev_region) << " before " << (*it_region);
          throw SMITHLABException(ss.str());
        }
        prev_region = it_region;
      }

      // only consider statistically significant bins -- we just skip over
      // non-sig bins. they will still be 'in' the cluster though, so 'f'
      // can do whatever it wants with them.
      if (*it_pval <= alpha) {
        //cerr << "looking at " << *it_region << endl;
        if (first) {
          clusterChrom = it_region->get_chrom();
          clusterStartCoord = it_region->get_start();
          clusterEndCoord = it_region->get_end();
          clusterStart_region = it_region;
          clusterEnd_region = it_region + 1;
          clusterStart_pval = it_pval;
          clusterEnd_pval = it_pval + 1;
          for (size_t i=0; i<c_starts.size(); ++i) {
            clusterStart_covars[i] = it_covars[i];
            clusterEnd_covars[i] = it_covars[i] + 1;
          }
          first = false;
        } else {
          if ((clusterChrom != it_region->get_chrom()) ||
              (clusterEndCoord + this->clusterDistance <= it_region->get_start())) {
            // we end the current cluster if this region is too far away from it
            // including if it's on a different chromosome
            f(clusterStart_region, clusterEnd_region,
              clusterStart_pval, clusterEnd_pval,
              clusterStart_covars, clusterEnd_covars);
            clusterChrom = it_region->get_chrom();
            clusterStartCoord = it_region->get_start();
            clusterEndCoord = it_region->get_end();
            clusterStart_region = it_region;
            clusterEnd_region = it_region + 1;
            clusterStart_pval = it_pval;
            clusterEnd_pval = it_pval + 1;
            for (size_t i=0; i<c_starts.size(); ++i) {
              clusterStart_covars[i] = it_covars[i];
              clusterEnd_covars[i] = it_covars[i] + 1;
            }
          } else {
            // we extend the current cluster if this region is within the
            // required distance
            if (it_region->get_end() > clusterEndCoord)
              clusterEndCoord = it_region->get_end();
            if (it_region->get_start() < clusterStartCoord)
              clusterStartCoord = it_region->get_start();
            clusterEnd_region = it_region + 1;
            clusterEnd_pval = it_pval + 1;
            for (size_t i=0; i<c_starts.size(); ++i) {
              clusterEnd_covars[i] = it_covars[i] + 1;
            }
          }
        }
      }

      ++it_pval;
      for (size_t i=0; i<c_starts.size(); ++i) ++(it_covars[i]);
    }
    // don't forget the last cluster
    f(clusterStart_region, clusterEnd_region,
      clusterStart_pval, clusterEnd_pval,
      clusterStart_covars, clusterEnd_covars);
  }

  /**
   * \brief TODO
   */
  template<typename Func>
  void aggregate(const std::vector<GenomicRegion>::const_iterator r_start,
      const std::vector<GenomicRegion>::const_iterator r_end,
      const std::vector<double>::const_iterator p_start,
      const std::vector<double>::const_iterator p_end,
      const double alpha,
      const Func f) const {
    const std::vector< vector<double>::const_iterator > dummy;
    aggregate(r_start, r_end, p_start, p_end, dummy, dummy, alpha, f);
  }

private:
  /** maximum distance that a region can be and still be considered part of the
   * cluster */
  size_t clusterDistance;
};


/**
 * \brief TODO
 */
class ClusterLimitsPrinter {
public:
  ClusterLimitsPrinter(std::ostream &ostrm) : outputStream(ostrm) {}
  void operator() (std::vector<GenomicRegion>::const_iterator s_region,
                   std::vector<GenomicRegion>::const_iterator e_region,
                   std::vector<double>::const_iterator s_pval,
                   std::vector<double>::const_iterator e_pval,
                   std::vector< std::vector<double>::const_iterator > s_covar,
                   std::vector< std::vector<double>::const_iterator > e_covar)
  const {
    // just to cut down on typing..
    typedef std::vector<GenomicRegion>::const_iterator gIter;
    typedef std::vector<double>::const_iterator dIter;

    // if name is consistent, use it -- else use "X"
    // use the average of the scores
    std::string name = "";
    double score = 0;
    size_t numBins = 0;
    for (gIter it = s_region; it != e_region; ++it) {
      numBins += 1;
      if (name == "") name = it->get_name();
      if (name != it->get_name()) {
        name = "X";
      }
      score += it->get_score();
    }
    score = score / numBins;

    // pick the strand that occurs most often
    size_t pos = 0, neg = 0;
    char strand = '+';
    for (gIter it = s_region; it != e_region; ++it) {
      if (it->get_strand() == '-') ++neg;
      else ++pos;
    }
    if (neg > pos) strand = '-';

    // average of p-values within the cluster
    double pval = std::accumulate(s_pval, e_pval, 0.0);
    pval = pval / numBins;

    // average of covariates within the cluster
    vector<double> cvars;
    for (size_t i=0; i<s_covar.size(); i++) {
      double sum = 0;
      for (dIter it = s_covar[i]; it != e_covar[i]; ++it) sum += (*it);
      cvars.push_back(sum / numBins);
    }

    GenomicRegion tmp;
    tmp.set_chrom(s_region->get_chrom());
    tmp.set_start(s_region->get_start());
    tmp.set_end((e_region-1)->get_end());
    tmp.set_name(name);
    tmp.set_score(score);
    tmp.set_strand(strand);
    this->outputStream << tmp << "\t" << pval;
    if (cvars.size() != 0) this->outputStream << "\t";
    for (size_t i=0; i<cvars.size(); ++i) {
      this->outputStream << cvars[i];
      if (i != cvars.size()-1) this->outputStream << "\t";
    }
    this->outputStream << std::endl;
  }
private:
  std::ostream &outputStream;
};


/**
 * \brief TODO
 */
class ClusterSummitPrinter {
public:
  ClusterSummitPrinter(std::ostream &ostrm, bool findMax) :
    fMax(findMax), outputStream(ostrm) {}
  void operator() (std::vector<GenomicRegion>::const_iterator s_region,
      std::vector<GenomicRegion>::const_iterator e_region,
      std::vector<double>::const_iterator s_pval,
      std::vector<double>::const_iterator e_pval,
      std::vector< std::vector<double>::const_iterator > s_covar,
      std::vector< std::vector<double>::const_iterator > e_covar) const {
    typedef std::vector<double>::const_iterator DIter;
    DIter mpval = (this->fMax) ? std::max_element(s_pval, e_pval) :
                                 std::min_element(s_pval, e_pval);

    outputStream << *(s_region + (mpval - s_pval)) << "\t" << *mpval;
    if (s_covar.size() != 0) this->outputStream << "\t";
    for (size_t i=0; i<s_covar.size(); ++i) {
      this->outputStream << *(s_covar[i] + (mpval - s_pval));
      if (i != s_covar.size()-1) this->outputStream << "\t";
    }
    this->outputStream << std::endl;
  }
  static const bool CLUSTER_MAX_SCORE = true;
  static const bool CLUSTER_MIN_SCORE = false;
private:
  // little nested class for comparing regions by score
  class CompScore {
  public :
    bool operator() (const GenomicRegion &lhs, const GenomicRegion &rhs) {
      return lhs.get_score() < rhs.get_score();
    }
  };

  bool fMax;
  std::ostream &outputStream;
};

#endif
