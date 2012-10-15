/**
  \file ZTNB.cpp
  \brief This source file defines the zero-truncated negative binomial
         distribution

  \authors Philip J. Uren, Song Qiang, Timothy Daley, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Song Qiang, Timothy Daley, Andrew D. Smith

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

#include <sstream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "ZTNB.hpp"
#include "smithlab_utils.hpp"

using std::max;
using std::cin;
using std::setw;
using std::fabs;
using std::ceil;
using std::cout;
using std::sort;
using std::endl;
using std::cerr;
using std::pair;
using std::string;
using std::vector;
using std::ostream;
using std::greater;
using std::isfinite;
using std::make_pair;
using std::stringstream;
using std::numeric_limits;

const double ZTNB::max_allowed_alpha;
const double ZTNB::min_allowed_alpha;
const double ZTNB::tolerance;
const double ZTNB::MAX_RANDOM_MU;
const double ZTNB::MAX_RANDOM_ALPHA;

/******************************************************************************
 * Static inline functions
 *****************************************************************************/
static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}
static inline double
compute_mean(const vector<size_t> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
           vals_hist.end(), 0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}
static inline double
compute_mean(const vector<double> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
           vals_hist.end(), 0.0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}

/******************************************************************************
 * ZTNB class
 *****************************************************************************/


/******************************************************************************
 * Simple Inspectors
 */

/**
 * \brief Return the parameters for this zero-truncated negative binomial
 *        distribution as a vector of double. Order is mu, alpha
 */
vector<double>
ZTNB::getParams() const {
  vector<double> r;
  r.push_back(this->mu);
  r.push_back(this->alpha);
  return r;
};

/**
 * \brief Get an XML object representation of this zero-truncated negative
 *        binomial distribution.
 */
SimpleXML
ZTNB::toXML() const {
  stringstream ss;
  ss << "<ZTNB>"
     << "<mean>" << this->mu << "</mean>"
     << "<alpha>" << this->alpha << "</alpha>"
     << "</ZTNB>";
  return SimpleXML(ss);
}

/**
 * \brief Get a string representation fo this zero-truncated negative binomial
 */
string
ZTNB::toString() const {
  stringstream s;
  s << "mean=" << this->mu << " alpha=" << this->alpha;
  return s.str();
}

/**
 * \brief Get a string representation of this zero-truncated negative
 *        binomial distribution in XML format.
 */
string
ZTNB::toXMLString() const {
  return this->toXML().toString();
}


/******************************************************************************
 * Simple mutators
 */

/**
 * \brief Set the parameters for this zero-truncated negative binomial
 *        distribution. We expect two parameters, mu and alpha, in that
 *        order
 */
void
ZTNB::setParams(const std::vector<double>& params) {
  if (params.size() != 2) {
    stringstream ss;
    ss << "Setting parameters for zero-truncated NegativeBinomial failed. "
       << "Reason: expected two parameters (mu and alpha), but instead we got "
       << params.size() << " parameters";
    throw DistributionException(ss.str());
  }
  this->mu = params[0];
  this->alpha = params[1];
  this->set_helpers();
}

/**
 * \brief Randomize the parameters for this zero-truncated negative binomial
 *        distribution.
 */
void
ZTNB::randomise() {
  this->mu = (rand() / double(RAND_MAX)) * MAX_RANDOM_MU;
  this->alpha = (rand() / double(RAND_MAX)) * MAX_RANDOM_ALPHA;
  this->set_helpers();
}


/******************************************************************************
 * Public members for describing the distribution
 */

/**
 * \brief The zero-truncated negative binomial probability mass function
 *
 * We override this here from the base class because we're going to handle
 * zero carefully, and we need to convert the parameter to size_t at this
 * point to make use of the existing code.
 */
double
ZTNB::pdf(const double x) const {
  if (x == 0) return 0;

  // we're going to cast x from double to size_t, but be careful doing it..
  double intpart;
  if ((x < 0) || (!(modf(x, &intpart) == 0.0))) {
    stringstream ss;
    ss << "evaluating ZTNB PDF value at " << x << " failed. "
       << "Reason: " << x << " is not a non-negative integer";
    throw DistributionException(ss.str());
  }

  return exp(this->trunc_log_pdf(size_t(x)));
}

/**
 * \brief The zero-truncated negative binomial log-likelihood function
 *        evaluated at a single point as defined by the given parameter
 * \param x the response to evaluate the function at; although this is a
 *          double, the value it contains must be a non-zero, positive integer
 * \throws DistributionException if x is not >= 1 or x is not an integer
 */
double
ZTNB::loglikelihood(const double x) const {
  // can't handle zero..
  if (x==0) {
    stringstream ss;
    ss << "Zero-truncated negative binomial log-likelihood calculation failed. "
       << "Reason: log-likelihood is undefined for response of 0";
    throw DistributionException(ss.str());
  }

  // we're going to cast x from double to size_t, but be careful doing it..
  double intpart;
  if ((x < 0) || (!(modf(x, &intpart) == 0.0))) {
    stringstream ss;
    ss << "evaluating ZTNB PDF value at " << x << " failed. "
       << "Reason: " << x << " is not a non-negative integer";
    throw DistributionException(ss.str());
  }

  return this->trunc_log_pdf(size_t(x));
}

/**
 * \brief Calculate a p-value for the given response and the current values
 *        of the distribution parameters
 * \throws DistributionException if x is not an integer, or is less than 0,
 *                               or if the resultant p-value is > 1, < 0 or
 *                               non-finite.
 */
double
ZTNB::pvalue(const double x) const {
  // handle zero carefully
  if (x==0) return 1;

  // we're going to cast x from double to size_t, but be careful doing it..
  double intpart;
  if ((x < 0) || (!(modf(x, &intpart) == 0.0))) {
    stringstream ss;
    ss << "evaluating ZTNB PDF value at " << x << " failed. "
       << "Reason: " << x << " is not a non-negative integer";
    throw DistributionException(ss.str());
  }

  // get the answer and do a couple of sanity checks on it
  double res = this->trunc_pval(size_t(x));
  if (!isfinite(res)) {
    stringstream ss;
    ss << "Calculating p-value for zero-truncated negative binomial with "
       << "response of " << x << " failed."
       << "Reason: p-value is non-finite";
    throw DistributionException(ss.str());
  }
  if ((res > 1) || (res < 0)) {
    stringstream ss;
    ss << "Calculating p-value for zero-truncated negative binomial with "
       << "response of " << x << " failed."
       << "Reason: p-value is less than 0 or greater than 1";
    throw DistributionException(ss.str());
  }

  return res;
}

/******************************************************************************
 * Public members for parameter estimation
 */

/**
 * \brief  estimate the parameters of this distribution given the
 *         responses <ys> and the probabilities for each response <probs>
 * \bug probs is not currently used
 */
void
ZTNB::estimateParams(const std::vector<double>& ys,
                    const std::vector<double>& probs,
                    const bool verbose) {
  for (size_t i=0; i<probs.size(); i++) {
    if (probs[i] != 1)
      throw DistributionException("ZTNB fitting with probs not implemented");
  }

  //const vector<double>& probst = (probs.size() > 0) ?
  //    probs : vector<double>(ys.size(), 1);
  //assert(ys.size() == probst.size());  // TODO: throw exception here

  size_t max = static_cast<size_t>(*std::max_element(ys.begin(), ys.end()));
  vector<size_t> hist(max+1,0);
  //vector<double> hprobs(max+1,1);

  for (size_t i=0; i<ys.size(); i++) {
    //TODO: messy cast to int here..
    hist[int(ys[i])] += 1;
    //hprobs[ys[i]] *= probst[i];
  }
  EM_estim_params (0.001, 100, hist); // TODO: fix magic numbers..

  // sanity check on the fitted values
  if ((!isfinite(this->mu)) || (!isfinite(this->alpha))) {
    stringstream ss;
    ss << "Fitting zero-truncated negative binomial failed "
       << "Reason: fitted value for mu or alpha was non-finite";
    throw DistributionException(ss.str());
  }
}


/******************************************************************************
 * ZTNB I/O
 */

/**
 * \brief load this ZTNB object from an XML description
 */
void
ZTNB::load(const SimpleXML& xml) {
  vector<SimpleXML> meanX = xml.getChildren("mean");
  if (meanX.size() != 1) {
    stringstream ss;
    ss << "loading zero-truncated negative binomial from XML object " << endl
       << xml.toString() << endl << " failed. Reason: expected a single child "
       << "tag of type <mean> -- didn't find it";
    throw DistributionException(ss.str());
  }
  vector<SimpleXML> alphaX = xml.getChildren("alpha");
  if (alphaX.size() != 1) {
    stringstream ss;
    ss << "loading zero-truncated negative binomial from XML object " << endl
       << xml.toString() << endl << " failed. Reason: expected a single child "
       << "tag of type <alpha> -- didn't find it";
    throw DistributionException(ss.str());
  }

  std::stringstream ss1, ss2;
  ss1 << meanX[0].getData();
  ss2 << alphaX[0].getData();
  ss1 >> this->mu;
  ss2 >> this->alpha;

  // sanity check on the values we loaded..
  if ((!isfinite(this->mu)) || (!isfinite(this->alpha))) {
    stringstream ss;
    ss << "loading zero-truncated negative binomial from XML object " << endl
       << xml.toString() << endl << " failed. Reason: loaded value for "
       << "mu or alpha was non-finite";
    throw DistributionException(ss.str());
  }

  this->set_helpers();
}

/****
 * @summary: write this ZTNB object to the given ostream. Format is XML
 */
void
ZTNB::save(std::ostream& ostrm) const {
  ostrm << this->toXML().toString();
}


/******************************************************************************
 * Private Member functions for distribution specification
 */

/**
 * \brief TODO
 */
double
ZTNB::trunc_log_L(const vector<size_t> &vals_hist) const {
  const double log_1_plus_a_mu = log(1 + alpha*mu);
  const double log_1_minus_prob_zero = log(1.0 - exp(-log(1+alpha*mu)/alpha));
  const double alpha_inv = 1.0/alpha;
  const size_t lim = vals_hist.size();
  double holding_val = 0.0;
  double log_L = 0.0;
  for (size_t i = 1; i < lim; i++){
    holding_val += log(1+alpha*(i-1));
    log_L += vals_hist[i]*
      (holding_val - gsl_sf_lngamma(i) + 
       i*log(mu)-(i+alpha_inv)*log_1_plus_a_mu - log_1_minus_prob_zero);
  }
  
  return(log_L);
}

/**
 * \brief TODO
 */
double
ZTNB::trunc_pval(const size_t val) const {
  double pval = 1.0;
  for(size_t i = 1; i < val; i++) {
    pval -= exp(trunc_log_pdf(i));
    if (i == 20000) break;
    if (pval <= 0) {
      pval = 0;
      break;
    }
  }
  return(pval);
}

/**
 * \brief Computes the expected number of summands needed to reach a sum of
 *        "sum"
 */
double
ZTNB::expected_inverse_sum(const double mean,
                            const size_t sample_size,
                            const size_t sum) const {
  const double alph = get_alpha();
  const double m = get_mu();
  const double prob_zero = exp(-log(1+alph*m)/alph);
  return(exp(log(sample_size) - log(1-prob_zero) +
       log(1-exp(-log(1+exp(log(alph)+log(sum)-log(sample_size)
          +log(m) - log(mean)))/alph))));
}


/**
 * \brief TODO
 */
double
ZTNB::expected_zeros(const double pseudo_size) const {
  const double alpha = get_alpha();
  const double mu = get_mu();
  const double prob_zero = pow(1+alpha*mu, -1/alpha);
  const double expected_zeros = pseudo_size*(prob_zero/(1-prob_zero));
  return expected_zeros;
}

/**
 * \brief TODO
 */
double 
ZTNB::trunc_log_pdf(const size_t val) const {
  double return_val = 0.0;
  double holding_val = 0.0;
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
  if(val > 0){
    for(size_t j =0; j < val; ++j){
      holding_val += log(1+alpha*j);
    }
    holding_val -= gsl_sf_lngamma(val+1);
  }
  if(val > 0)
    return_val  = holding_val +val*log(mu) - (val+1/alpha)*log(1+alpha*mu)
    -log(1-prob_zero);
  else
    return_val  = -log(1+alpha*mu)/alpha;
  
  return return_val;
}

/**
 * \brief TODO
 */
double
ZTNB::log_pdf(const size_t val) const {
  double holding_val = 0.0;
  if(val > 0){
    for(size_t j =0; j < val; ++j){
      holding_val += log(1+alpha*j);
    }
    holding_val -= gsl_sf_lngamma(val+1);
  }
  return (holding_val + val*log(mu) - (val+1/alpha)*log(1+alpha*mu));
}

/**
 * \brief TODO
 */
void
ZTNB::set_helpers() {
  n_helper = 1/alpha;
  p_helper = n_helper/(n_helper + mu);
  n_log_p_minus_lngamma_n_helper = 
    n_helper*log(p_helper) - gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
}

/**
 * \brief TODO
 */
double 
ZTNB::operator()(const int val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) - 
        gsl_sf_lnfact(static_cast<size_t>(val))) + 
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!finite(P)) return -40;
  return P;
}

/**
 * \brief TODO
 */
double
ZTNB::score_fun_first_term(const vector<size_t> &vals_hist,
          const double a_mid){
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
  inner_sum += j/(1 + a_mid*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}

/**
 * \brief TODO
 */
double
ZTNB::alpha_score_function(const vector<size_t> &vals_hist,
          const double mean,
          const double a_mid,
          const double vals_count){
  const double one_plus_alpha_mu = 1 + a_mid*mean;
  return (score_fun_first_term(vals_hist, a_mid)/vals_count + 
    (log(one_plus_alpha_mu)/a_mid - mean)/a_mid); 
}

/**
 * \brief TODO
 */
double
ZTNB::score_fun_first_term(const vector<double> &vals_hist,
        const double a_mid) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
  inner_sum += j/(1 + a_mid*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum; 
}

/**
 * \brief TODO
 */
double
ZTNB::alpha_score_function(const vector<double> &vals_hist, 
        const double mean,
        const double a_mid,
        const double vals_count){
  const double one_plus_alpha_mu = 1 + a_mid*mean;
  return (score_fun_first_term(vals_hist, a_mid)/vals_count + 
    (log(one_plus_alpha_mu)/a_mid - mean)/a_mid); 
}


/******************************************************************************
 * Private parameter estimation members
 */

/**
 * \brief TODO
 * \note  this seems to be fitting just the normal NB?
 */
void 
ZTNB::estim_params(const vector<size_t> &vals_hist,
                    const vector<double> &probs) {
  // TODO: throw exception if sizes don't match on inputs
  vector<double> pseudo_hist(vals_hist.size(), 0.0);
  for(size_t i = 0; i < vals_hist.size(); i++){
    pseudo_hist[i] = vals_hist[i]*probs[i];
  }
  mu = compute_mean(pseudo_hist);
  
  const double pseudo_size = 
    accumulate(pseudo_hist.begin(), pseudo_hist.end(), 0.0);
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(pseudo_hist,mu, a_mid, pseudo_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  // bisection, but what happened to the terms involving
                  // the gamma func? See Zhang et al. top of page 7
  set_helpers();
} 


/**
 * \brief TODO
 * \note  this seems to be fitting just the normal NB?
 */
void 
ZTNB::estim_params(const vector<size_t> &vals_hist){
  mu = compute_mean(vals_hist);
  //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
  
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), 
           vals_hist.end(), 0));
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  // bisection, but what happened to the terms involving 
                  // the gamma func? See Zhang et al. top of page 7
  set_helpers();
}

/**
 * \brief TODO
 */
double
ZTNB::EM_estim_params(const double tol, const size_t max_iter,
          vector<size_t> &vals_hist){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(),0));
  double error = numeric_limits<double>::max();
  double prev_score = numeric_limits<double>::max();
  double score = 0.0;
  for(size_t i = 0; i < max_iter; i++){
    // TODO fix messy cast..
    vals_hist[0] = size_t(expected_zeros(vals_size));
    estim_params(vals_hist);
    vals_hist[0] = 0;
    score = trunc_log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol)
      break;
    prev_score = score;
  }
  return(trunc_log_L(vals_hist));
}




