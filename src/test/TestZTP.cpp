/**
  \file TestZTP.cpp
  \brief This source file defines a set of unit tests for the ZTP distribution
         class.

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
#include "ZTP.hpp"

/**
 * \brief test the zero-truncated poisson PDF with a few sample values with
 *        known answers
 */
TEST(ZTPTest, testPDF) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {0, 0.5819767, 0.2909884, 0.09699612,
                           0.02424903, 0.004849806};
  ZeroTruncatedPoisson z(1);
  for (int i=0; i<6; i++) {
    EXPECT_NEAR(rightAnswers[i], z.pdf(inputs[i]), TOL);
  }
}

/**
 * \brief test the ZTP logLikelihood function with a few sample values with
 *        known answers
 */
TEST(ZTPTest, testLoglikelihood) {
  const double TOL = 0.0001;
  int inputs[] = {1, 2, 3, 4, 5};
  double rightAnswers[] = {-0.5413249, -1.234472, -2.333084,
                           -3.719379, -5.328817};
  ZeroTruncatedPoisson z(1);
  for (int i=0; i<5; i++) {
    EXPECT_NEAR(rightAnswers[i], z.loglikelihood(inputs[i]), TOL);
  }
}

/**
 * \brief test the ZTP CDF function with a few sample values with
 *        known answers
 */
TEST(ZTPTest, testCDF) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {0, 0.5819767, 0.872965, 0.9699612,
                           0.9942102, 0.99906};
  ZeroTruncatedPoisson z(1);
  for (int i=0; i<5; i++) {
    EXPECT_NEAR(rightAnswers[i], z.cdf(inputs[i]), TOL);
  }
}
