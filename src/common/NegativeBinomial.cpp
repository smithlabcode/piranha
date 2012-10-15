/**
  \file NegativeBinomial.cpp
  \brief This source file defines the negative binomial distribution

  \authors Philip J. Uren, Qiang Song, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Qiang Song, Andrew D. Smith

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

//#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <string> 
#include <limits>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <iostream>
#include <iterator>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_diff.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>

#include "SimpleXML.hpp"
#include "smithlab_utils.hpp"
#include "NegativeBinomial.hpp"

using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;
using std::string;
using std::bind2nd;
using std::divides;
using std::istream;
using std::make_pair;
using std::accumulate;
using std::max_element;
using std::stringstream;
using std::numeric_limits;

const double NegativeBinomial::max_allowed_alpha;
const double NegativeBinomial::min_allowed_alpha;
const double NegativeBinomial::alpha_allowed_error;

/******************************************************************************
 * Simple inspectors
 */

/****
 * @summary: TODO
 */
vector<double>
NegativeBinomial::getParams() const {
  std::vector<double> r;
  r.push_back(this->mean);
  r.push_back(this->alpha);
  return r;
};

/****
 * @summary: TODO
 */
SimpleXML
NegativeBinomial::toXML() const {
  stringstream ss;
  ss << "<NegativeBinomial>"
     << "<mean>" << this->mean << "</mean>"
     << "<alpha>" << this->alpha << "</alpha>"
     << "</NegativeBinomial>";
  cerr << "making XML object" << endl;
  SimpleXML tmp(ss);
  cerr << "returning object " << endl;
  return tmp;
}

/****
 * @summary: TODO
 */
string
NegativeBinomial::toXMLString() const {
  return this->toXML().toString();
}

/****
 * @summary: return a string representation of this Negative Binomial
 *           distribution.
 * TODO      change to XML representation
 */
string
NegativeBinomial::toString() const {
  stringstream s;
  s << "mean=" << this->mean << " alpha=" << this->alpha;
  return s.str();
}

/******************************************************************************
 * Simple mutators
 */

/****
 * @summary: Set the parameters for this NegativeBinomial distribution. We
 *           expect two parameters, mu and alpha, in that order.
 */
void
NegativeBinomial::setParams(const vector<double>& params) {
  if (params.size() != 2) {
    stringstream ss;
    ss << "Setting parameters for NegativeBinomial failed. Reason: "
       << "expected two parameters (mu and alpha), but instead we got "
       << params.size() << " parameters";
    throw DistributionException(ss.str());
  }
  this->mean = params[0];
  this->alpha = params[1];
  this->setHelpers();
}

/****
 * @summary: TODO
 */
void
NegativeBinomial::randomise() {
  // TODO fix magic number 4...
  this->mean = (rand() / double(RAND_MAX)) * 4;
  this->alpha = (rand() / double(RAND_MAX)) * 4;
  this->setHelpers();
}

/******************************************************************************
 * Specifying the distribution
 */

/****
 * @summary: TODO
 * TODO      throw exception for non-finite log-likelihood instead
 */
double 
NegativeBinomial::loglikelihood(const double x) const {
  const double P = (gsl_sf_lngamma(x + n_helper) - 
                   gsl_sf_lnfact(static_cast<size_t>(x))) +
                   n_log_p_minus_lngamma_n_helper + x*log_q_helper;
  if (!finite(P)) return -40;
  return P;
}

/******************************************************************************
 * Parameter estimation
 */

/****
 * @summary: Estimate mu and alpha for this NegativeBinomial distribution
 * @throws:  DistributionException if the size of <ys> does not match the size
 *           of <probs>
 */
void 
NegativeBinomial::estimateParams(const vector<double>& ys, 
                                 const vector<double>& probs,
                                 const bool verbose) {
  if (ys.size() != probs.size()) {
    stringstream ss;
    ss << "Estimating parameters for NegativeBinomial distribution failed. "
       << "Reason: number of responses does not match number of probabilities";
    throw DistributionException(ss.str());
  }
  this->estimateParamsML(ys, probs, verbose);
  if (verbose)
    cerr << "fitted NB distribution: " << this->toString() << endl;
}


/****
 * @summary: TODO
 */
void
NegativeBinomial::setHelpers() {
  n_helper = 1/this->alpha;
  p_helper = n_helper/(n_helper + this->mean);
  n_log_p_minus_lngamma_n_helper = n_helper*log(p_helper) - 
    gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
}


/****
 * @summary: TODO
 */
static inline double
score_fun_first_term(const vector<double> &vals_hist, 
             const double mu, const double alpha) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
    inner_sum += j/(1 + alpha*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}


/****
 * @summary: TODO
 */
static inline double
alpha_score_function(const vector<double> &vals_hist, const double mu, 
             const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_fun_first_term(vals_hist, mu, alpha)/vals_count + 
      (log(one_plus_alpha_mu)/alpha - mu)/alpha);
}


/****
 * @summary: TODO
 */
void
NegativeBinomial::estimateParamsML(const vector<double> &vals,
                                   const vector<double> &probs,
                                   const bool verbose) {
  vector<double> workspace_vals, workspace_probs;
    
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);// - centering_value;
    workspace_vals[i] = log(vals[i]) + log(probs[i]);// - centering_value;
  }
  
  const double vals_count =
      exp(smithlab::log_sum_log_vec(workspace_probs, lim));
  this->mean =
      exp(smithlab::log_sum_log_vec(workspace_vals, lim))/vals_count;
  
  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.begin() + lim);
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];
  
  const double mu = this->mean;
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  
  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error && 
     fabs((a_high - a_low)/max(a_high, a_low)) > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val =
        alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/std::max(mid_val, prev_val));
    prev_val = mid_val;
  }
  this->alpha = a_mid;

  setHelpers();
}

/******************************************************************************
 * Negative Binomial I/O
 */

/****
 * @summary: read this NegativeBinomial distribution from the provided XML
 *           object.
 * @throws:  DistributionException if the XML object does not define the mean
 *           and alpha tags
 * TODO      throw exception if root tag is incorrectly named
 */
void
NegativeBinomial::load(const SimpleXML& xml) {
  vector<SimpleXML> meanX = xml.getChildren("mean");
  if (meanX.size() != 1) {
    stringstream ss;
    ss << "loading Negative binomial from XML object " << endl
       << xml.toString() << endl << " failed. Reason: expected a single child "
       << "tag of type <mean> -- didn't find it";
    throw DistributionException(ss.str());
  }
  vector<SimpleXML> alphaX = xml.getChildren("alpha");
  if (alphaX.size() != 1) {
    stringstream ss;
    ss << "loading Negative binomial from XML object " << endl
       << xml.toString() << endl << " failed. Reason: expected a single child "
       << "tag of type <alpha> -- didn't find it";
    throw DistributionException(ss.str());
  }

  std::stringstream ss1, ss2;
  ss1 << meanX[0].getData();
  ss2 << alphaX[0].getData();
  ss1 >> this->mean;
  ss2 >> this->alpha;
  this->setHelpers();
}

/****
 * @summary: write this NegativeBinomial distribution to the provided ostream.
 *           output format is XML
 */
void
NegativeBinomial::save(std::ostream& ostrm) const {
  ostrm << this->toXMLString() << "\t";
}

