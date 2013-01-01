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
#define RREGAGG_HPP

#include <vector>
#include <string>

#include "GenomicRegion.hpp"

class GenomicRegionAggregator {
public :
  /**
   * \brief Construct a cluster aggregator with the given distance
   */
  GenomicRegionAggregator(size_t clusterDist) : clusterDistance(clusterDist) {};

  /**
   * \brief Find the region clusters and call the function f on each cluster.
   *        The regions must be sorted
   * \param start   iterator to start of regions
   * \param end     iterator to end of regions
   * \param f       function to call on each cluster, must take two const
   *                iterators to vector which will delineate the first and
   *                last+1 bins in the cluster
   * \throws SMITHLABException if regions aren't sorted
   */
  template<typename Func>
  void aggregate(std::vector<GenomicRegion>::const_iterator start,
                 std::vector<GenomicRegion>::const_iterator end,
                 Func f) const {
    // just to cut down on typing..
    typedef std::vector<GenomicRegion>::const_iterator gIter;

    std::string clusterChrom = "";
    size_t clusterStartCoord = 0, clusterEndCoord = 0;
    gIter clusterStart = start, clusterEnd = start;
    gIter prev;
    bool first = true;

    for (gIter it = start; it != end; ++it) {
      // make sure everything is sorted
      if (!first) {
        if (it < prev) {
          std::stringstream ss;
          ss << "aggregating genomic regions failed -- regions were not "
             << "sorted: saw " << (*prev) << " before " << (*it);
          throw SMITHLABException(ss.str());
        }
        prev = it;
      }

      if (first) {
        clusterChrom = it->get_chrom();
        clusterStartCoord = it->get_start();
        clusterEndCoord = it->get_end();
        first = false;
      } else {
        if ((clusterChrom != it->get_chrom()) ||
            (clusterEndCoord + this->clusterDistance < it->get_start())) {
          // we end the current cluster if this region is too far away from it
          // including if it's on a different chromosome
          f(clusterStart, clusterEnd);
          clusterChrom = it->get_chrom();
          clusterStartCoord = it->get_start();
          clusterEndCoord = it->get_end();
          clusterStart = it;
          clusterEnd = it + 1;
        } else {
          // we extend the current cluster if this region is within the
          // required distance
          if (it->get_end() > clusterEndCoord)
            clusterEndCoord = it->get_end();
          if (it->get_start() < clusterStartCoord)
            clusterStartCoord = it->get_start();
          clusterEnd = it + 1;
        }
      }
    }
    // don't forget the last cluster
    f(clusterStart, clusterEnd);
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
  void operator() (std::vector<GenomicRegion>::const_iterator s,
                   std::vector<GenomicRegion>::const_iterator e) {
    GenomicRegion tmp;
    tmp.set_chrom(s->get_chrom());
    tmp.set_start(s->get_start());
    tmp.set_end((e-1)->get_end());
    tmp.set_name((e-1)->get_name() == s->get_name() ? s->get_name() : "X");
    tmp.set_score(0);
    this->outputStream << tmp << std::endl;
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
  void operator() (std::vector<GenomicRegion>::const_iterator s,
                   std::vector<GenomicRegion>::const_iterator e) {
    if (this->fMax)
      this->outputStream << *std::max_element(s,e, CompScore()) << std::endl;
    else
      this->outputStream << *std::min_element(s,e, CompScore()) << std::endl;
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
