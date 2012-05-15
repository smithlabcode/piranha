/**
  \file ZTP.cpp
  \brief This source filed defines the zero-truncated Poisson distribution

  \authors Philip J. Uren, Timothy Daley, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,

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
#include <sstream>
#include <cmath>
#include <iostream>

// from GSL
#include <gsl/gsl_sf.h>

// from piranha common
#include "ZTP.hpp"

using std::stringstream;
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::isfinite;


/******************************************************************************
 * Static functions used in ZTP parameter estimation
 */

/**
 * \brief TODO
 */
static inline double
movement(const double a, const double b) {
  return fabs((a - b)/std::max(a, b)); //delta
}

/**
 * \brief TODO
 */
static inline double
lambda_score_funct(const double mean, const double lambda){
  return(lambda + mean*exp(-lambda) - mean);
}

/**
 * \brief TODO
 */
static double
trunc_estim_param_bisec(const double mean, const double tol){

  double lambda_low = mean-1;
  double lambda_high = mean;
  double lambda_mid = mean - 0.5;

  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();

  while(movement(lambda_high, lambda_low) > tol &&
         diff > tol){

    lambda_mid = (lambda_low + lambda_high)/2;

    const double mid_val = lambda_score_funct(mean, lambda_mid);

    if (mid_val < 0) lambda_low = lambda_mid;
    else lambda_high = lambda_mid;

    diff = fabs((prev_val - mid_val)/std::max(mid_val, prev_val));

    prev_val = mid_val;
  }

  if (!(isfinite(lambda_mid))) {
    stringstream ss;
    ss << "Zero-truncated Poisson parameter estimation failed. Reason: "
       << "got non-finite value for lambda";
    throw DistributionException(ss.str());
  }

  return lambda_mid;
}

/******************************************************************************
 * Defining the distribution
 */

/**
 * \brief ZTP probability mass function in log space
 */
double
ZeroTruncatedPoisson::loglikelihood(const double y) const {
  if (y==0) {
    stringstream ss;
    ss << "Zero-truncated Poisson log-likelihood calculation failed. Reason: "
       << "log-likelihood is undefined for response of 0";
    throw DistributionException(ss.str());
  }

  // we're going to cast y from double to size_t, but be careful doing it..
  double intpart;
  if ((y < 0) || (!(modf(y, &intpart) == 0.0))) {
    stringstream ss;
    ss << "evaluating ZTNB PDF value at " << y << " failed. "
       << "Reason: " << y << " is not a non-negative integer";
    throw DistributionException(ss.str());
  }

  double res = (-this->lambda + (y*log(this->lambda)) -
               gsl_sf_lnfact(size_t(y))) - log(1-exp(-this->lambda));
  if (!(isfinite(res))) {
    stringstream ss;
    ss << "Zero-truncated Poisson log-likelihood calculation failed. Reason: "
       << "got non-finite value for response of " << y << " with "
       << "distribution parameters of " << this->toString();
    throw DistributionException(ss.str());
  }
  return res;
}


/**
 * \brief The zero-truncated probability mass function -- we override this
 *        from the base class because we need to handle the zero carefully
 */
double
ZeroTruncatedPoisson::pdf(const double x) const {
  if (x == 0) return 0;
  return Distribution::pdf(x);
}


/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief estimate lambda (the untruncated mean) of this zero-truncated
 *        Poisson distribution.
 *
 * uses the bisection algorithm
 */
void
ZeroTruncatedPoisson::estimateParams(const vector<double>& obs,
                                     const vector<double>& probs,
                                     const bool verbose) {
  cout << verbose << endl;
  double mu = 0;
  for (size_t i=0; i<obs.size(); i++) {
    mu += (obs[i] * probs[i]);
  }
  mu = mu / obs.size();
  this->lambda = trunc_estim_param_bisec(mu, this->tolerance);
}

/******************************************************************************
 * Simple Inspectors
 */

/****
 * \brief return a vector of parameters for this distribution. In the case
 *        of the ZTP, it's just the value of lambda in a vector of size 1
 */
vector<double>
ZeroTruncatedPoisson::getParams() const {
  std::vector<double> r;
  r.push_back(this->lambda);
  return r;
};

/**
 * \brief get XML representation of this zero-truncated Poisson object
 */
SimpleXML
ZeroTruncatedPoisson::toXML() const {
  stringstream ss;
  ss << "<ZeroTruncatedPoisson>"
     << "<lambda>" << this->lambda << "</lambda>"
     << "</ZeroTruncatedPoisson>";
  SimpleXML tmp(ss);
  return tmp;
}

/**
 * \brief get string representation of this ZTP object in XML format
 */
string
ZeroTruncatedPoisson::toXMLString() const {
  return this->toXML().toString();
}

/**
 * \brief get string representation of this ZTP object
 */
string
ZeroTruncatedPoisson::toString() const {
  stringstream ss;
  ss << "lambda: " << this->lambda;
  return ss.str();
}

/******************************************************************************
 * Simple mutators
 */

/****
 * \brief  Set the parameters for this zero-truncated Poisson distribution
 * \param  params a vector with only a single value in it, which is the value
 *               of lambda
 * \throws DistributionException if more than one value is given in params
 */
void
ZeroTruncatedPoisson::setParams(const vector<double>& params) {
  if (params.size() != 1) {
    stringstream ss;
    ss << "Setting parameters for zero-truncated Poisson failed. Reason: "
       << "expected one parameters (lambda), but instead we got "
       << params.size() << " parameters";
    throw DistributionException(ss.str());
  }
  if (!(isfinite(params[0]))) {
    stringstream ss;
    ss << "Setting parameters for zero-truncated Poisson failed. Reason: "
       << "supplied parameter value for lambda was non-finite";
    throw DistributionException(ss.str());
  }
  this->lambda = params[0];
}

/**
 * \brief randomize the parameters for this distribution. In the case of the
 *        ZTP this just selects a random value for lambda between 0 and
 *        MAX_RANDOM_LAMBDA
 */
void
ZeroTruncatedPoisson::randomise() {
  this->lambda = (rand() / double(RAND_MAX)) * MAX_RANDOM_LAMBDA;
  if (!(isfinite(this->lambda))) {
    stringstream ss;
    ss << "Randomizing zero-truncated Poisson failed. Reason: "
       << "got non-finite value for lambda";
    throw DistributionException(ss.str());
  }
}

/******************************************************************************
 * Zero-truncated Poisson I/O
 */

/**
 * \brief  read this NegativeBinomial distribution from the provided XML object
 * \throws DistributionException if the XML object does not define the mean
 *         and alpha tags
 * \todo   throw exception if root tag is incorrectly named
 */
void
ZeroTruncatedPoisson::load(const SimpleXML& xml) {
  vector<SimpleXML> meanX = xml.getChildren("lambda");
  if (meanX.size() != 1) {
    stringstream ss;
    ss << "loading zero-truncated Poisson from XML object " << endl
       << xml.toString() << endl << " failed. Reason: expected a single child "
       << "tag of type <mean> -- didn't find it";
    throw DistributionException(ss.str());
  }

  std::stringstream ss1, ss2;
  ss1 << meanX[0].getData();
  ss1 >> this->lambda;

  if (!(isfinite(this->lambda))) {
    stringstream ss;
    ss << "Loading zero-truncated Poisson failed. Reason: "
       << "got non-finite value for lambda";
    throw DistributionException(ss.str());
  }
}

/****
 * \brief write this zero-truncated Poisson distribution to the provided
 *        ostream. Output format is XML
 */
void
ZeroTruncatedPoisson::save(std::ostream& ostrm) const {
  ostrm << this->toXMLString() << "\t";
}
