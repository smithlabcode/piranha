/**
  \file Gaussian.cpp
  \brief This file implements the class for the Gaussian distribution

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


#define _USE_MATH_DEFINES

#include <vector>
#include <string> 
#include <sstream>
#include <iostream>
#include "Gaussian.hpp"

using std::istream;
using std::vector;
using std::string;
using std::stringstream;


/****
 * @summary: 
 */
void Gaussian::load(const SimpleXML& xml) {
  std::stringstream sm, sv; 
  vector<SimpleXML> mc = xml.getChildren("mean"), 
                    vc = xml.getChildren("variance");
  
  // sanity check
  assert(mc.size() == 1);
  assert(vc.size() == 1);
  assert(mc[0].isLeaf());
  assert(vc[0].isLeaf());
  
  sm << mc[0].getData();
  sv << vc[0].getData(); 
  sm >> this->mean;
  sv >> this->variance;
}


/****
 * @summary:
 */
double 
Gaussian::loglikelihood(const double x) const {
  double resqr = (x - this-> mean) * (x - this-> mean);
  return -0.5 * log(2*M_PI) - 0.5 * log(this->variance) - 
         (resqr / (2*this->variance));
}

/****
 * @summary:
 */
void 
Gaussian::estimateParams(const vector<double>& ys, const vector<double>& probs,
    const bool verbose) {
  assert(ys.size() == probs.size());
  size_t N = ys.size();
  
  // find weighted mean
  double sum = 0, count = 0;
  for (size_t i=0; i<N; i++) {
    sum += (ys[i] * probs[i]);
    count += probs[i];
  }
  this->mean = sum / count;
  
  // find weighted variance
  double sqrResidualSum = 0;
  for (size_t i=0; i<N; i++) { 
    double r = (ys[i] - mean) * (ys[i] - mean);
    sqrResidualSum += (r * probs[i]);
  }
  this->variance = sqrResidualSum / count;
}
