/**
  \file TestZTNB.cpp
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
