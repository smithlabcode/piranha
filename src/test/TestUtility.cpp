/**
  \file TestUtility.cpp
  \brief This source file defines a set of unit tests for utility functions
         defined in utility.hpp

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

#include "gtest/gtest.h"
#include "Utility.hpp"
#include <tr1/unordered_set>

#include <vector>
#include <string>

using std::vector;
using std::string;
using std::tr1::unordered_set;


/**
 * \brief TODO
 */
TEST(UtilTest, testCopyByIndex) {
  // simple basic test of functionality
  string vals_1[] = {"zero", "one", "two", "three", "four", "five", "six"};
  vector<string> vals_1_v (vals_1, vals_1+7), vals_2_v;
  size_t indices_1[] = {0, 2, 4, 5};
  unordered_set<size_t> indices (indices_1, indices_1+4);

  string expect_1[] = {"zero", "one", "two", "three", "four", "five", "six"};
  string expect_2[] = {"zero", "two", "four", "five"};

  copyByIndex<string>(vals_1_v, vals_2_v, indices);

  size_t e_size = 7;
  EXPECT_EQ(vals_1_v.size(), e_size);
  for (size_t i=0; i<e_size; i++) {
    EXPECT_EQ(expect_1[i], vals_1_v[i]);
  }
  e_size = 4;
  EXPECT_EQ(vals_2_v.size(), e_size);
  for (size_t i=0; i<e_size; i++) {
    EXPECT_EQ(expect_2[i], vals_2_v[i]);
  }
}

/**
 * \brief TODO
 */
TEST(UtilTest, testRemoveByIndex) {
  size_t indices_1[] = {0, 2, 4, 5};
  unordered_set<size_t> indices (indices_1, indices_1+4);
  string expect_1[] = {"one", "three", "six"};
  string vals_1[] = {"zero", "one", "two", "three", "four", "five", "six"};
  vector<string> vals (vals_1, vals_1+7);
  removeByIndex<string>(vals, indices);
  for (int i=0; i<3; i++) {
    EXPECT_EQ(expect_1[i], vals[i]);
  }
}

/**
 * \brief TODO
 */
TEST(UtilTest, testCopyRemoveByIndex) {
  size_t indices_1[] = {0, 2, 4, 5};
  unordered_set<size_t> indices (indices_1, indices_1+4);
  string expect_1[] = {"one", "three", "six"};
  string expect_2[] = {"zero", "two", "four", "five"};
  string vals_1[] = {"zero", "one", "two", "three", "four", "five", "six"};
  vector<string> vals_1_v (vals_1, vals_1+7);
  vector<string> vals_2_v;

  copyRemoveByIndex<string>(vals_1_v, vals_2_v, indices);
  size_t e_size = 3;
  EXPECT_EQ(vals_1_v.size(), e_size);
  for (size_t i=0; i<3; i++) {
    EXPECT_EQ(expect_1[i], vals_1_v[i]);
  }
  e_size = 4;
  EXPECT_EQ(vals_2_v.size(), e_size);
  for (size_t i=0; i<vals_2_v.size(); i++) {
    EXPECT_EQ(expect_2[i], vals_2_v[i]);
  }
}

/**
 * \brief TODO
 */
TEST(UtilTest, testCountItemsOrdered) {
  string items[] = {"apple", "orange", "banana", "banana", "orange", "banana"};
  vector<string> items_v (items, items + 6);
  string expect_keys[] = {"apple", "banana", "orange"};
  size_t expect_counts[] = {1, 3, 2};
  std::map<string, size_t> res = countItemsOrdered<string>(items_v);
  size_t indx = 0;
  typedef std::map<string, size_t>::iterator iter;
  for (iter it=res.begin(); it != res.end(); ++it) {
    EXPECT_EQ(it->first, expect_keys[indx]);
    EXPECT_EQ(it->second, expect_counts[indx]);
    indx += 1;
  }
}

/**
 * \brief TODO
 */
TEST(UtilTest, testMergeResponsesCovariatesPvals) {
  // 8 fg sites, 9 bg sites
  vector<GenomicRegion> sites_fg, sites_bg, expect_sites;
  sites_fg.push_back(GenomicRegion("chr1", 1,  2,  "C1R0",  1, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 2,  3,  "C1R1",  2, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 3,  4,  "C1R2",  5, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 4,  5,  "C1R3",  2, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 5,  6,  "C1R4",  4, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 6,  7,  "C1R5",  2, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 7,  8,  "C1R6",  4, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 8,  9,  "C1R7",  6, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 9,  10, "C1R8",  8, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 10, 11, "C1R9",  1, '+'));
  // ****
  sites_bg.push_back(GenomicRegion("chr2", 3, 4,  "C2R0",  1, '+'));
  sites_bg.push_back(GenomicRegion("chr2", 4, 5,  "C2R1",  4, '+'));
  sites_fg.push_back(GenomicRegion("chr2", 5, 6,  "C2R2",  1, '+'));
  sites_fg.push_back(GenomicRegion("chr2", 6, 7,  "C2R3",  7, '+'));
  sites_bg.push_back(GenomicRegion("chr2", 7, 8,  "C2R4",  2, '+'));
  sites_fg.push_back(GenomicRegion("chr2", 8, 9,  "C2R5",  5, '+'));
  sites_bg.push_back(GenomicRegion("chr2", 9, 10, "C2R6",  2, '+'));
  // ****
  expect_sites.push_back(GenomicRegion("chr1", 1,  2,  "C1R0",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 2,  3,  "C1R1",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 3,  4,  "C1R2",  5, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 4,  5,  "C1R3",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 5,  6,  "C1R4",  4, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 6,  7,  "C1R5",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 7,  8,  "C1R6",  4, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 8,  9,  "C1R7",  6, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 9,  10, "C1R8",  8, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 10, 11, "C1R9",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 3, 4,   "C2R0",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 4, 5,   "C2R1",  4, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 5, 6,   "C2R2",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 6, 7,   "C2R3",  7, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 7, 8,   "C2R4",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 8, 9,   "C2R5",  5, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 9, 10,  "C2R6",  2, '+'));

  // p-values
  const double pvs_fg[8] = {0.02,  0.0005, 0.001, 0.001,
                            0.002, 0.001,  0.05,  0.05};
  const double pvs_bg[9] = {0.05,  0.002,  0.001, 0.002,
                            0.2,   0.001,  0.2,   0.005, 0.01};
  vector<double> pvals_fg(pvs_fg, pvs_fg+8);
  vector<double> pvals_bg(pvs_bg, pvs_bg+9);
  const double pvs_exp[17] = {0.02,  0.05,  0.0005, 0.002, 0.001,
                              0.001, 0.002, 0.001,  0.002, 0.2,
                              0.001, 0.2,   0.001,  0.05,  0.005,
                              0.05,  0.01};
  vector<double> expect_pvals(pvs_exp, pvs_exp+17);

  // covariates
  const double c1_fg[8] = {1,2, 5,2,7,1,2,3};
  const double c1_bg[9] = {4,10,1,1,2,5,7,1,1};
  const double c2_fg[8] = {2,2, 7,1,5,2,2,1};
  const double c2_bg[9] = {9,11,2,1,2,7,8,7,3};
  const double cv1_exp[17] = {1, 4, 2, 10, 5,
                              1, 1, 2, 7,  2,
                              5, 7, 1, 2,  1,
                              3, 1};
  const double cv2_exp[17] = {2, 9, 2, 11, 7,
                              2, 1, 1, 5,  2,
                              7, 8, 2, 2,  7,
                              1, 3};
  vector<double> expect_cv1(cv1_exp, cv1_exp+17),
                 expect_cv2(cv2_exp, cv2_exp+17);

  vector< vector<double> > cvs_fg, cvs_bg;
  cvs_fg.push_back(vector<double>());
  cvs_fg.push_back(vector<double>());
  cvs_bg.push_back(vector<double>());
  cvs_bg.push_back(vector<double>());
  for (size_t i=0; i<8; i++) {
    cvs_fg[0].push_back(c1_fg[i]);
    cvs_fg[1].push_back(c2_fg[i]);
  }
  for (size_t i=0; i<9; i++) {
    cvs_bg[0].push_back(c1_bg[i]);
    cvs_bg[1].push_back(c2_bg[i]);
  }

  mergeResponsesCovariatesPvals(sites_fg, sites_bg, pvals_fg, pvals_bg,
                                cvs_fg, cvs_bg);

  EXPECT_EQ(expect_pvals, pvals_bg);
  EXPECT_EQ(expect_sites, sites_bg);
  EXPECT_EQ(expect_cv1, cvs_bg[0]);
  EXPECT_EQ(expect_cv2, cvs_bg[1]);
}


/**
 * \brief TODO
 */
TEST(UtilTest, testMergeResponsesPvals) {
  // 8 fg sites, 9 bg sites
  vector<GenomicRegion> sites_fg, sites_bg, expect_sites;
  sites_fg.push_back(GenomicRegion("chr1", 1,  2,  "C1R0",  1, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 2,  3,  "C1R1",  2, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 3,  4,  "C1R2",  5, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 4,  5,  "C1R3",  2, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 5,  6,  "C1R4",  4, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 6,  7,  "C1R5",  2, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 7,  8,  "C1R6",  4, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 8,  9,  "C1R7",  6, '+'));
  sites_fg.push_back(GenomicRegion("chr1", 9,  10, "C1R8",  8, '+'));
  sites_bg.push_back(GenomicRegion("chr1", 10, 11, "C1R9",  1, '+'));
  // ****
  sites_bg.push_back(GenomicRegion("chr2", 3, 4,  "C2R0",  1, '+'));
  sites_bg.push_back(GenomicRegion("chr2", 4, 5,  "C2R1",  4, '+'));
  sites_fg.push_back(GenomicRegion("chr2", 5, 6,  "C2R2",  1, '+'));
  sites_fg.push_back(GenomicRegion("chr2", 6, 7,  "C2R3",  7, '+'));
  sites_bg.push_back(GenomicRegion("chr2", 7, 8,  "C2R4",  2, '+'));
  sites_fg.push_back(GenomicRegion("chr2", 8, 9,  "C2R5",  5, '+'));
  sites_bg.push_back(GenomicRegion("chr2", 9, 10, "C2R6",  2, '+'));
  // ****
  expect_sites.push_back(GenomicRegion("chr1", 1,  2,  "C1R0",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 2,  3,  "C1R1",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 3,  4,  "C1R2",  5, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 4,  5,  "C1R3",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 5,  6,  "C1R4",  4, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 6,  7,  "C1R5",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 7,  8,  "C1R6",  4, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 8,  9,  "C1R7",  6, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 9,  10, "C1R8",  8, '+'));
  expect_sites.push_back(GenomicRegion("chr1", 10, 11, "C1R9",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 3, 4,   "C2R0",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 4, 5,   "C2R1",  4, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 5, 6,   "C2R2",  1, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 6, 7,   "C2R3",  7, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 7, 8,   "C2R4",  2, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 8, 9,   "C2R5",  5, '+'));
  expect_sites.push_back(GenomicRegion("chr2", 9, 10,  "C2R6",  2, '+'));

  // p-values
  const double pvs_fg[8] = {0.02,  0.0005, 0.001, 0.001,
                            0.002, 0.001,  0.05,  0.05};
  const double pvs_bg[9] = {0.05,  0.002,  0.001, 0.002,
                            0.2,   0.001,  0.2,   0.005, 0.01};
  vector<double> pvals_fg(pvs_fg, pvs_fg+8);
  vector<double> pvals_bg(pvs_bg, pvs_bg+9);
  const double pvs_exp[17] = {0.02,  0.05,  0.0005, 0.002, 0.001,
                              0.001, 0.002, 0.001,  0.002, 0.2,
                              0.001, 0.2,   0.001,  0.05,  0.005,
                              0.05,  0.01};
  vector<double> expect_pvals(pvs_exp, pvs_exp+17);

  mergeResponsesPvals(sites_fg, sites_bg, pvals_fg, pvals_bg);

  EXPECT_EQ(expect_pvals, pvals_bg);
  EXPECT_EQ(expect_sites, sites_bg);
}







