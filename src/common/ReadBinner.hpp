/**
  \file  ReadBinner.hpp
  \brief Declares the ReadBinner class, which takes provides functionality
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

#ifndef RBINNER_HPP
#define RBINNER_HPP

#include <vector>
#include <string>

#include "GenomicRegion.hpp"

/**
 * \brief TODO
 */
class ReadBinner {
public:
  ReadBinner(size_t bsize) : binSize(bsize) {};
  void binReads(const std::vector<GenomicRegion> &reads,
                      std::vector<GenomicRegion> &bins,
                      const size_t pseudoCount=0) const;
  void binReads(const std::vector<GenomicRegion> &reads,
                      std::vector<GenomicRegion> &bins,
                const std::vector<GenomicRegion> &requiredBins,
                const size_t pseudoCount=0) const;
private:
  void binChromosome(const std::vector<GenomicRegion> &reads,
                     std::vector<GenomicRegion> &bins,
                     const size_t pseudoCount=0) const;
  void binChromosome(const std::vector<GenomicRegion> &reads,
                     std::vector<GenomicRegion> &bins,
                     const std::vector<GenomicRegion> &requiredBins,
                     const size_t pseudoCount=0) const;
  size_t binSize;
};

#endif
