/****
  \file ReadBinner.cpp
  \brief Defines the ReadBinner class, which takes provides functionality
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
**/

#include <vector>
#include <string>
#include <sstream>

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

#include "ReadBinner.hpp"

using std::string;
using std::vector;
using std::stringstream;
using std::cerr;
using std::endl;

/**
 * \brief TODO
 */
void
ReadBinner::binReads(const vector<GenomicRegion> &reads,
                           vector<GenomicRegion> &bins) const {
  if (!check_sorted(reads)) {
    stringstream ss;
    ss << "binning reads failed. Reason: reads were not sorted in order";
    throw SMITHLABException(ss.str());
  }

  bins.clear();
  vector<vector<GenomicRegion> > separated_by_chrom;
  separate_chromosomes(reads, separated_by_chrom);

  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    size_t read_idx = 0;
    string chromName = separated_by_chrom[i].front().get_chrom();
    size_t lim = separated_by_chrom[i].back().get_end();
    for (size_t j = 0; j < lim; j += this->binSize) {
      size_t counts = 0;
      while (read_idx < separated_by_chrom[i].size() &&
             separated_by_chrom[i][read_idx].get_start() < j + this->binSize) {
        if (separated_by_chrom[i][read_idx].get_start() >= j)
          ++counts;
        ++read_idx;
      }
      if (counts > 0) {
        bins.push_back(GenomicRegion(chromName, j, j+this->binSize, "X",
                                     counts, '+'));
      }
    }
  }
}


