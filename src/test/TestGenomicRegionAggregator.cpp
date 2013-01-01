/**
  \file TestGenomicRegionAggregator
  \brief TODO

  \authors Philip J. Uren, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2012
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

#include "gtest/gtest.h"
#include "GenomicRegionAggregator.hpp"

#include <vector>
#include <sstream>

using std::endl;
using std::vector;
using std::stringstream;

/**
 * \brief simple test of GenomicRegionAggregator -- everything is on the same
 *        strand, names are unique so all output names will be "X" and all
 *        output strands will be "+"
 *
 */
TEST(GRegAggTest, testPrintSize10) {
  vector<GenomicRegion> regions, expect;
  regions.push_back(GenomicRegion("chr1", 10, 20, "R1",  0, '+'));
  regions.push_back(GenomicRegion("chr1", 13, 25, "R2",  0, '+'));
  regions.push_back(GenomicRegion("chr1", 29, 34, "R3",  0, '+'));
  regions.push_back(GenomicRegion("chr2", 10, 20, "R4",  0 , '+'));
  regions.push_back(GenomicRegion("chr2", 31, 32, "R5",  0, '+'));
  regions.push_back(GenomicRegion("chr3", 10, 20, "R6",  0, '+'));
  regions.push_back(GenomicRegion("chr3", 13, 21, "R7",  0, '+'));
  regions.push_back(GenomicRegion("chr3", 50, 53, "R8",  0, '+'));
  regions.push_back(GenomicRegion("chr3", 56, 59, "R9",  0, '+'));
  regions.push_back(GenomicRegion("chr4", 10, 20, "R10", 0, '+'));
  regions.push_back(GenomicRegion("chr4", 80, 88, "R11", 0, '+'));

  expect.push_back(GenomicRegion("chr1", 10, 34, "X",   0, '+'));
  expect.push_back(GenomicRegion("chr2", 10, 20, "R4",  0, '+'));
  expect.push_back(GenomicRegion("chr2", 31, 32, "R5",  0, '+'));
  expect.push_back(GenomicRegion("chr3", 10, 21, "X",   0, '+'));
  expect.push_back(GenomicRegion("chr3", 50, 59, "X",   0, '+'));
  expect.push_back(GenomicRegion("chr4", 10, 20, "R10", 0, '+'));
  expect.push_back(GenomicRegion("chr4", 80, 88, "R11", 0, '+'));
  stringstream ss_exp;
  for (size_t i=0; i<expect.size(); i++) ss_exp << expect[i] << endl;

  stringstream ss;
  GenomicRegionAggregator(10).aggregate(regions.begin(), regions.end(),
                                        ClusterLimitsPrinter(ss));

  EXPECT_EQ(ss_exp.str(), ss.str());
}

/**
 * \brief TODO
 */
TEST(GRegAggTest, testClusterSummitPrinter) {
  vector<GenomicRegion> regions, expect;
  regions.push_back(GenomicRegion("chr1", 10, 20, "R1",  0.01, '+'));
  regions.push_back(GenomicRegion("chr1", 13, 25, "R2",  0.02, '+'));
  regions.push_back(GenomicRegion("chr1", 29, 34, "R3",  0.05, '+'));
  regions.push_back(GenomicRegion("chr2", 10, 20, "R4",  0.05, '+'));
  regions.push_back(GenomicRegion("chr2", 31, 32, "R5",  0.03, '+'));
  regions.push_back(GenomicRegion("chr3", 10, 20, "R6",  0.03, '+'));
  regions.push_back(GenomicRegion("chr3", 13, 21, "R7",  0.02, '+'));
  regions.push_back(GenomicRegion("chr3", 50, 53, "R8",  0.01, '+'));
  regions.push_back(GenomicRegion("chr3", 56, 59, "R9",  0.03, '+'));
  regions.push_back(GenomicRegion("chr4", 10, 20, "R10", 0.03, '+'));
  regions.push_back(GenomicRegion("chr4", 80, 88, "R11", 0.02, '+'));

  expect.push_back(GenomicRegion("chr1", 10, 20, "R1",  0.01, '+'));
  expect.push_back(GenomicRegion("chr2", 10, 20, "R4",  0.05, '+'));
  expect.push_back(GenomicRegion("chr2", 31, 32, "R5",  0.03, '+'));
  expect.push_back(GenomicRegion("chr3", 13, 21, "R7",  0.02, '+'));
  expect.push_back(GenomicRegion("chr3", 50, 53, "R8",  0.01, '+'));
  expect.push_back(GenomicRegion("chr4", 10, 20, "R10", 0.03, '+'));
  expect.push_back(GenomicRegion("chr4", 80, 88, "R11", 0.02, '+'));
  stringstream ss_exp;
  for (size_t i=0; i<expect.size(); i++) ss_exp << expect[i] << endl;

  stringstream ss;
  GenomicRegionAggregator(10).aggregate(regions.begin(), regions.end(),
            ClusterSummitPrinter(ss, ClusterSummitPrinter::CLUSTER_MIN_SCORE));

  EXPECT_EQ(ss_exp.str(), ss.str());
}

