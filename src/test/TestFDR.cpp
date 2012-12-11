/**
  \file TestFDR.cpp
  \brief Defines tests for the FDR module

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
#include "FDR.hpp"

#include <vector>

using std::vector;

/**
 * \brief TODO
 */
TEST(FDRTest, testCorrectP1) {
  const double TOL = 0.0001;
  double p_1[] = {1.0000, 0.7590, 0.6528, 0.5719, 0.4262, 0.3240, 0.0459,
                  0.0344, 0.0298, 0.0278, 0.0201, 0.0095, 0.0019, 0.0004,
                  0.0001};
  double r_1[] = {1.00000000, 0.81321429, 0.75323077, 0.71487500, 0.58118182,
                  0.48600000, 0.07650000, 0.06450000, 0.06385714, 0.06385714,
                  0.06030000, 0.03562500, 0.00950000, 0.00300000, 0.00150000};
  vector<double> p_1_v (p_1, p_1 + 15);
  vector<double> expect_1_v (r_1, r_1 + 15);

  FDR::correctP(p_1_v);
  size_t e_size = 15;
  EXPECT_EQ(p_1_v.size(), e_size);
  for (size_t i=0; i<p_1_v.size(); i++) {
    EXPECT_NEAR(p_1_v[i], expect_1_v[i], TOL);
  }
}

/**
 * \brief TODO
 */
TEST(FDRTest, testCorrectP2) {
  const double TOL = 0.0001;
  double p_1[] = {0.40434733, 0.05974949, 0.21910558, 0.41934038,
                  0.60906895, 0.01364008, 0.87263729, 0.45343604,
                  0.76458524, 0.91082686};
  double r_1[] = {0.7557267, 0.2987475, 0.7303519, 0.7557267, 0.8700985,
                  0.1364008, 0.9108269, 0.7557267, 0.9108269, 0.9108269};
  vector<double> p_1_v (p_1, p_1 + 10);
  vector<double> expect_1_v (r_1, r_1 + 10);

  FDR::correctP(p_1_v);
  size_t e_size = 10;
  EXPECT_EQ(p_1_v.size(), e_size);
  for (size_t i=0; i<p_1_v.size(); i++) {
    EXPECT_NEAR(p_1_v[i], expect_1_v[i], TOL);
  }
}
