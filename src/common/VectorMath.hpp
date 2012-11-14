/**
  \file VectorMath.hpp
  \brief This header file declares a set of helper functions for doing
         vector-level math operations.

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

#ifndef VMATH_HPP
#define VMATH_HPP

#include <vector>

double vectorDotProduct(const std::vector<double> a, 
                        const std::vector<double> b);
std::vector<double> vectorScalarProduct(const double s, 
                                        const std::vector<double> v);
std::vector<double> vectorSubtraction(const std::vector<double> a, 
                                      const std::vector<double> b);
bool vectorUnderThreshold(const std::vector<double> v, double t);
bool vectorUnderThresholdAbs(const std::vector<double> v, double t);


#endif

