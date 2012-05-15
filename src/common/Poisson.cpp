/**
  \file Poisson.cpp
  \brief this source file defines the Poisson distribution class

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

#include <vector>
#include <string> 
#include <sstream>
#include <iostream>
#include "gsl/gsl_sf_gamma.h"
#include "Poisson.hpp"

using std::istream;
using std::vector;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::isfinite;

/******************************************************************************
 * Simple inspectors
 */

/**
 * \brief return an XML object which represents this distribution
 */
SimpleXML
Poisson::toXML() const {
  stringstream ss;
  ss << "<Poisson>"
     << "<mean>" << this->mean << "</mean>"
     << "</Poisson>";
  SimpleXML tmp(ss);
  return tmp;
}

/**
 * \brief return a string representation of this Poisson object in XML format
 */
string
Poisson::toXMLString() const {
  return this->toXML().toString();
}

/**
 * \brief return a string representation of this Poisson object.
 *
 * This is generally intended to be used for debugging and printing on the
 * console, rather than persistence.
 *
 */
string
Poisson::toString() const {
  stringstream s;
  s << "mean=" << this->mean;
  return s.str();
}

/**
 * \brief get the parameters for this Poisson distribtuion. This is just
 *        the lambda parameter, so actually this vector will have only a
 *        single element
 */
vector<double>
Poisson::getParams() const {
  vector<double> r;
  r.push_back(this->mean);
  return r;
}

/******************************************************************************
 *
 */

/**
 * \brief Set the parameters for this Poisson distribution
 *
 *        There is only one parameter, lambda (mean), so the vector must
 *        contain only a single element
 */
void
Poisson::setParams(const vector<double>& params) {
  if (params.size() != 2) {
    stringstream ss;
    ss << "Setting parameters for Poisson failed. Reason: "
       << "expected one parameter (lambda, the mean), but instead we got "
       << params.size() << " parameters";
    throw DistributionException(ss.str());
  }
  this->mean = params[0];
}

/**
 * \brief Randomize the parameters for this Poisson object
 */
void
Poisson::randomise() {
  this->mean = (rand() / double(RAND_MAX)) * MAX_RANDOM_LAMBDA;
}

/******************************************************************************
 * Poisson I/O
 */

/**
 * \brief save this Poisson distribution to the given stream; format is XML
 */
void
Poisson::save(std::ostream& ostrm) const {
  ostrm << this->toXMLString();
}

/**
 * \brief load this Poisson distribution from an XML description
 */
void Poisson::load(const SimpleXML& xml) {
  vector<SimpleXML> meanX = xml.getChildren("mean");
  if (meanX.size() != 1) {
    stringstream ss;
    ss << "loading Poisson from XML object " << endl
       << xml.toString() << endl << " failed. Reason: expected a single child "
       << "tag of type <mean> -- didn't find it";
    throw DistributionException(ss.str());
  }
  std::stringstream ss;
  ss << meanX[0].getData();
  ss >> this->mean;

  // sanity check on the values we loaded..
  if (!isfinite(this->mean)) {
    stringstream ss;
    ss << "loading Poisson from XML object " << endl
       << xml.toString() << endl << " failed. Reason: loaded value for "
       << "lambda was non-finite";
    throw DistributionException(ss.str());
  }
}

/******************************************************************************
 * Describing the distribution
 */

/**
 * \brief get the log-likelihood of the current set of parameters given the
 *        single data point x
 */
double 
Poisson::loglikelihood(const double x) const {
  double intpart;
  if (!(modf(x, &intpart) == 0.0)) {
    stringstream ss;
    ss << "evaluating Poisson log-likelihood value at " << x << " failed. "
       << "Reason: " << x << " is not an integer";
    throw DistributionException(ss.str());
  }
  // now we cast x to an int, but we're sure it was an int anyway..
  double res = (x * log(this->mean)) - this->mean - gsl_sf_lnfact(int(x));
  if (!isfinite(res)) {
    stringstream ss;
    ss << "evaluating Poisson log-likelihood function for lambda = "
       << this->mean << " at " << x << " failed. Reason: result was "
       << "non-finite";
    throw DistributionException(ss.str());
  }
  return res;
}

/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief Estimate the parameters for this Poisson distribution
 * \param ys vector of responses
 * \param probs vector of weights for each response
 * \param verbose print additional information to stderr during fitting
 */
void 
Poisson::estimateParams(const vector<double>& ys, const vector<double>& probs,
    const bool verbose) {
  if (ys.size() != probs.size()) {
    stringstream ss;
    ss << "Estimating parameters for Poisson distribution failed. "
       << "Reason: number of responses does not match number of probabilities";
    throw DistributionException(ss.str());
  }
  const size_t N = ys.size();
  
  // find weighted mean
  double sum = 0, count = 0;
  for (size_t i=0; i<N; i++) {
    sum += (ys[i] * probs[i]);
    count += probs[i];
  }
  this->mean = sum / count;

  // sanity check on the fitted value
  if (!isfinite(this->mean)) {
    stringstream ss;
    ss << "Estimating parameters for Poisson distribution failed. "
       << "Reason: fitted value of lambda is non-finite";
    throw DistributionException(ss.str());
  }

  if (verbose)
    cerr << "fitted Poisson distribution: " << this->toString() << endl;
}
