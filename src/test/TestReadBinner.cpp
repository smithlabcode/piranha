/**
  \file TestReadBinner.cpp
  \brief This source file defines a set of unit tests for the ReadBinner class

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
****/

#include <vector>

#include "gtest/gtest.h"

#include "ReadBinner.hpp"
#include "GenomicRegion.hpp"

using std::vector;

/**
 * \brief test binning reads into single nucleotide bins
 */
TEST(ReadBinnerTest, testSingleNucBin) {
  vector<GenomicRegion> input;
  input.push_back(GenomicRegion("chr1", 2, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 3, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 6, 8,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 9, 11, "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 8,  "X", 0, '+'));

  vector<GenomicRegion> expect;
  expect.push_back(GenomicRegion("chr1", 2, 3,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 3, 4,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 5, 6,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 6, 7,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 9, 10,  "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 5, 6,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 6, 7,   "X", 2, '+'));

  vector<GenomicRegion> got;
  ReadBinner r(1);
  r.binReads(input, got);

  EXPECT_EQ(got.size(), expect.size());
  if (got.size(), expect.size()) {
    for (size_t i = 0; i < expect.size(); i++) {
      EXPECT_EQ(expect[i], got[i]);
    }
  }
}

