/**
  \file TestNB.cpp
  \brief This source file defines a set of unit tests for the negative binomial
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
#include "NegativeBinomial.hpp"

/**
 * \brief test the NB PDF with a few sample values with known answers
 */
TEST(NBTest, testPDF) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625};
  NegativeBinomial nb(1,1);
  for (int i=0; i<6; i++) {
    EXPECT_NEAR(rightAnswers[i], nb.pdf(inputs[i]), TOL);
  }
}

/**
 * \brief test the NB logLikelihood function with a few sample values with
 *        known answers
 */
TEST(NBTest, testLoglikelihood) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {log(0.5), log(0.25), log(0.125), log(0.0625),
                           log(0.03125), log(0.015625)};
  NegativeBinomial nb(1,1);
  for (int i=0; i<5; i++) {
    EXPECT_NEAR(rightAnswers[i], nb.loglikelihood(inputs[i]), TOL);
  }
}

/**
 * \brief test the NB CDF function with a few sample values with
 *        known answers
 */
TEST(NBTest, testCDF) {
  const double TOL = 0.0001;
  int inputs[] = {0, 1, 2, 3, 4, 5};
  double rightAnswers[] = {0.5, 0.75, 0.875, 0.9375, 0.96875, 0.984375};
  NegativeBinomial nb(1,1);
  for (int i=0; i<5; i++) {
    EXPECT_NEAR(rightAnswers[i], nb.cdf(inputs[i]), TOL);
  }
}
