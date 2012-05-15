/**
  \file TestPoisson.cpp
  \brief This source file defines a set of unit tests for the Poisson
         distribution class.

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
#include "Poisson.hpp"

/**
 * \brief test the zero-truncated poisson PDF with a few sample values with
 *        known answers
 */
TEST(PoissonTest, testPDF) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {0.367879441, 0.367879441, 0.183939721, 0.061313240,
                           0.015328310, 0.003065662};
  Poisson z(1);
  for (int i=0; i<6; i++) {
    EXPECT_NEAR(rightAnswers[i], z.pdf(inputs[i]), TOL);
  }
}

/**
 * \brief test the ZTP logLikelihood function with a few sample values with
 *        known answers
 */
TEST(PoissonTest, testLoglikelihood) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {-1.000000, -1.000000, -1.693147, -2.791759,
                           -4.178054, -5.787492};
  Poisson z(1);
  for (int i=0; i<6; i++) {
    EXPECT_NEAR(rightAnswers[i], z.loglikelihood(inputs[i]), TOL);
  }
}

/**
 * \brief test the ZTP CDF function with a few sample values with
 *        known answers
 */
TEST(PoissonTest, testCDF) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {0.3678794, 0.7357589, 0.9196986, 0.9810118,
                           0.9963402, 0.9994058};
  Poisson z(1);
  for (int i=0; i<6; i++) {
    EXPECT_NEAR(rightAnswers[i], z.cdf(inputs[i]), TOL);
  }
}
