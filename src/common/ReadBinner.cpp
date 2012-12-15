/**
  \file  ReadBinner.cpp
  \brief Defines the ReadBinner class, which provides functionality
         for taking collections of sequencing reads and collecting them into
         bins

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
    14th Dec. 2012 -- Philip J. Uren -- rewrote to simplify logic; added option
                                        to provide pseudo-counts; added option
                                        to require certain bins in the output
**/

#include <vector>
#include <string>
#include <sstream>
#include <set>

// from smithlab common
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

// from Piranha common
#include "ReadBinner.hpp"

using std::set;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;

/**
 * \brief Count the number of reads falling into each equally sized bin
 *        from the start of every chromosome. Bins with zero reads are
 *        not included.
 * \param reads  vector of genome regions representing the reads
 * \param bins   vector to put the resultant bin counts into (will be cleared
 *               if there's anything in it
 */
void
ReadBinner::binReads(const vector<GenomicRegion> &reads,
                           vector<GenomicRegion> &bins,
                     const size_t pseudoCount) const {
  if (!check_sorted(reads)) {
    stringstream ss;
    ss << "binning reads failed. Reason: reads were not sorted in order of "
       << "chromosome and then end co-ordinate";
    throw SMITHLABException(ss.str());
  }

  bins.clear();
  vector<vector<GenomicRegion> > separated_by_chrom;
  separate_chromosomes(reads, separated_by_chrom);

  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    binChromosome(separated_by_chrom[i], bins, pseudoCount);
  }
}

/**
 * \brief Count the number of reads falling into each equally sized bin
 *        from the start of every chromosome. Bins with zero reads are
 *        not included (unless they're in requiredBins and pseudoCount is 0).
 * \param reads         vector of genome regions representing the reads.
 * \param bins          vector to put the resultant bin counts into (will be
 *                      cleared if there's anything in it.
 * \param requiredBins  a set of bins that must be included, regardless of
 *                      whether they have any reads or not
 * \param pseudoCount   a pseudoCount to add to all bins
 */
void
ReadBinner::binReads(const vector<GenomicRegion> &reads,
                           vector<GenomicRegion> &bins,
                     const vector<GenomicRegion> &requiredBins,
                     const size_t pseudoCount) const {
  if (!check_sorted(reads)) {
    stringstream ss;
    ss << "binning reads failed. Reason: reads were not sorted in order of "
       << "chromosome and then end coordinate";
    throw SMITHLABException(ss.str());
  }

  bins.clear();
  vector<vector<GenomicRegion> > separated_by_chrom, required_by_chrom;
  separate_chromosomes(reads, separated_by_chrom);
  separate_chromosomes(requiredBins, required_by_chrom);

  size_t readsChromIdx = 0, requiredChromIdx = 0;
  while ((readsChromIdx < separated_by_chrom.size()) ||
         (requiredChromIdx < required_by_chrom.size())) {
    string cnReads = "NONE", cnReqd = "NONE";
    if (readsChromIdx < separated_by_chrom.size())
      cnReads = separated_by_chrom[readsChromIdx].front().get_chrom();
    if (requiredChromIdx < required_by_chrom.size())
      cnReqd = required_by_chrom[requiredChromIdx].front().get_chrom();

    if (((readsChromIdx < separated_by_chrom.size()) &&
         (requiredChromIdx < required_by_chrom.size())) &&
         (cnReads == cnReqd)) {
      // have reads and required bins on this chrom
      binChromosome(separated_by_chrom[readsChromIdx], bins,
                    required_by_chrom[requiredChromIdx],
                    pseudoCount);
      readsChromIdx += 1;
      requiredChromIdx += 1;
    } else if ((readsChromIdx < separated_by_chrom.size()) &&
        ((requiredChromIdx == required_by_chrom.size()) ||
         (cnReads < cnReqd))) {
      // have reads on this chrom, but no required bins
      binChromosome(separated_by_chrom[readsChromIdx], bins, pseudoCount);
      readsChromIdx += 1;
    } else if ((requiredChromIdx < required_by_chrom.size()) &&
        ((readsChromIdx == separated_by_chrom.size()) ||
         (cnReads > cnReqd))) {
      // have required bins on this chrom, but no reads
      for (size_t i=0; i<required_by_chrom[requiredChromIdx].size(); ++i) {
        GenomicRegion tmp(required_by_chrom[requiredChromIdx][i]);
        tmp.set_score(pseudoCount);
        tmp.set_name("X");
        bins.push_back(tmp);
      }
      requiredChromIdx += 1;
    }
  }
}



/**
 * \brief                   Bin reads on a given chromosome
 * \param reads[IN]         Reads to be binned. Must all be on the same chrom
 * \param bins[OUT]         Vector to add new bins to. Won't be cleared, so
 *                          can be used to accumulate bins from multiple calls
 * \param pseudoCount[IN]   A pseudo count to add to each bin that is accepted
 *
 * TODO There is no check for correct sorting order!
 */
void
ReadBinner::binChromosome(const vector<GenomicRegion> &reads,
                                vector<GenomicRegion> &bins,
                          const size_t pseudoCount) const {
  if (reads.size() == 0) return;
  size_t lim = reads.back().get_end();
  size_t read_idx = 0;
  string chromName = reads.front().get_chrom();

  for (size_t j = 0; j < lim; j += this->binSize) {
    size_t counts = 0;
    while (read_idx < reads.size() &&
           reads[read_idx].get_start() < j + this->binSize) {
      if (reads[read_idx].get_start() >= j) ++counts;
      ++read_idx;

      if ((read_idx != reads.size()) &&
          (reads[read_idx].get_chrom() != chromName)) {
        stringstream ss;
        ss << "binning reads failed. Reason: binChromosome called with reads "
           << "from multiple chromosomes";
        throw SMITHLABException(ss.str());
      }
    }
    if (counts > 0) {
      bins.push_back(GenomicRegion(chromName, j, j+this->binSize, "X",
                                       counts + pseudoCount, '+'));
    }
  }
}


/**
 * \brief TODO
 *
 * TODO There is no check for correct sorting order!
 */
void
ReadBinner::binChromosome(const vector<GenomicRegion> &reads,
                                vector<GenomicRegion> &bins,
                          const vector<GenomicRegion> &requiredBins,
                          const size_t pseudoCount) const {
  if (reads.size() == 0) return;
  if (requiredBins.size() == 0) {
    binChromosome(reads, bins, pseudoCount);
  } else {
    string chromName = reads.front().get_chrom();
    size_t requiredIdx = 0;
    size_t readsIdx = 0;
    size_t requiredStart = requiredBins[requiredIdx].get_start();
    size_t requiredEnd = requiredBins[requiredIdx].get_end();
    if (requiredEnd - requiredStart != this->binSize) {
      stringstream ss;
      ss << "binning reads failed. Reason: size of required bins doesn't "
         << "match size of actual bins ";
      throw SMITHLABException(ss.str());
    }

    size_t lim = reads.back().get_end();
    for (size_t j = 0; j < lim; j += this->binSize) {
      size_t counts = 0;
      while (readsIdx < reads.size() &&
             reads[readsIdx].get_start() < j + this->binSize) {
        if (reads[readsIdx].get_start() >= j) ++counts;
        ++readsIdx;

        if ((readsIdx != reads.size()) &&
            (reads[readsIdx].get_chrom() != chromName)) {
          stringstream ss;
          ss << "binning reads failed. Reason: binChromosome called with "
             << "reads from multiple chromosomes";
          throw SMITHLABException(ss.str());
        }
      }

      // count warrants inclusion of this bin
      if (counts > 0) {
        bins.push_back(GenomicRegion(chromName, j, j+this->binSize, "X",
                                     counts + pseudoCount, '+'));
      }

      // this is a required bin
      if ((requiredIdx != requiredBins.size()) &&
          (requiredStart == j) && (requiredEnd == j + this->binSize)) {
        if (counts == 0)
          bins.push_back(GenomicRegion(chromName, j, j+this->binSize, "X",
                                       pseudoCount, '+'));
        requiredIdx += 1;
        if (requiredIdx != requiredBins.size()) {
          requiredStart = requiredBins[requiredIdx].get_start();
          requiredEnd = requiredBins[requiredIdx].get_end();
          if (requiredEnd - requiredStart != this->binSize) {
            stringstream ss;
            ss << "binning reads failed. Reason: size of required bins doesn't "
               << "match size of actual bins";
            throw SMITHLABException(ss.str());
          }
          if (requiredBins[requiredIdx].get_chrom() != chromName) {
            stringstream ss;
            ss << "binning reads failed. Reason: binChromosome called with "
               << "reads from multiple chromosomes";
            throw SMITHLABException(ss.str());
          }
        }
      }
    }

    // pick up any required bins that come after the last read
    for (size_t j = requiredIdx; j < requiredBins.size(); ++j) {
      GenomicRegion tmp(requiredBins[j]);
      tmp.set_score(pseudoCount);
      tmp.set_name("X");
      bins.push_back(tmp);
    }
  }
}


