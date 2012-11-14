/**
  \file VectorMath.cpp
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

#include <iostream>
#include <assert.h>
#include <cmath>
#include "VectorMath.hpp"


using std::vector;
using std::cout;
using std::endl;

/**
 * \brief TODO
 */
double
vectorDotProduct(const vector<double> a, const vector<double> b) {
  double res = 0;
  //assert(a.size() == b.size());
  for (size_t i=0; i<a.size(); i++) {
    res += a[i] * b[i]; 
  }
  return res;
}

/**
 * \brief TODO
 */
vector<double>
vectorScalarProduct(const double s, const vector<double> v) {
  vector<double> res;
  for (size_t i=0; i<v.size(); i++) res.push_back(v[i] * s);
  return res;
}

/**
 * \brief TODO
 */
vector<double>
vectorSubtraction(const vector<double> a, const vector<double> b) {
  vector<double> res;
  //assert(a.size() == b.size());
  for (size_t i=0; i<a.size(); i++) res.push_back(a[i] - b[i]);
  return res;
}

/**
 * \brief TODO
 */
bool 
vectorUnderThreshold(const vector<double> v, double t) {
  for (size_t i=0; i<v.size(); i++) {
    if (v[i] > t) return false;
  }
  return true;
}

/**
 * \brief TODO
 */
bool 
vectorUnderThresholdAbs(const vector<double> v, double t) {
  for (size_t i=0; i<v.size(); i++) {
    if (fabs(v[i]) > t) return false;
  }
  return true;
}

