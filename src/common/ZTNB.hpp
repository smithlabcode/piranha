/**
  \file ZTNB.hpp
  \brief This header file declares the zero-truncated negative binomial
         distribution

  \authors Philip J. Uren, Timothy Daley, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Timothy Daley, Andrew D. Smith

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

#ifndef ZTNBINOM_HPP
#define ZTNBINOM_HPP

#include <fstream>
#include <iomanip>
#include <numeric>
#include <limits>
#include <vector>
#include <string> 

// from GSL
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

// from smithlab common
#include "smithlab_utils.hpp"

// from piranha common
#include "Distribution.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

using std::max;
using std::setw;
using std::fabs;
using std::ceil;
using std::greater;
using std::numeric_limits;
using std::cin;

/**
 * \brief The zero-truncated negative binomial distribution
 *
 * This class is something of a wrapper around a set of already written ZTNB
 * methods; the older code has been hidden away in private members and the
 * public interface has been made consistent with the Distribution base class
 */
class ZTNB : public Distribution {
public:
  /*** constructors, destructors, object initialisation ***/
  ZTNB(const double m_, const double a_): mu(m_), alpha(a_) { set_helpers(); }
  ZTNB() : mu(1), alpha(1) { set_helpers(); }
  ~ZTNB() {;}
  
  /*** simple inspectors ***/
  double get_mu() const {return this->mu;}
  double get_alpha() const {return this->alpha;}
  std::vector<double> getParams() const;
  std::string toString() const;
  SimpleXML toXML() const;
  std::string toXMLString() const;

  /*** describing the distribution ***/
  double operator()(int val) const;
  double loglikelihood(const double x) const;
  double pdf(const double x) const;
  double pvalue(const double x) const;
  
  /*** ZTNB I/O ***/
  void save(std::ostream& ostrm) const;
  void load(const SimpleXML& xml);

  /*** simple mutators ***/
  void set_mu(const double m_) {mu = m_;}
  void set_alpha(const double a_) {alpha = a_;}
  void setParams(const std::vector<double>& params);
  void randomise();

  /*** Parameter estimation ***/
  void estimateParams(const std::vector<double>& ys,
                      const std::vector<double>& probs = std::vector<double>(),
                      const bool verbose=false);
  
private:
  /*** private constants ***/
  static const double max_allowed_alpha = 1000;
  static const double min_allowed_alpha = 1e-20;
  static const double tolerance = 1e-10;
  static const double MAX_RANDOM_MU = 100;
  static const double MAX_RANDOM_ALPHA = 10;

  /*** private member functions ***/
  double score_fun_first_term(const vector<size_t> &vals_hist,
			      const double a_mid);
  double alpha_score_function(const vector<size_t> &vals_hist,
			      const double mean, const double a_mid,
			      const double vals_count);
  double score_fun_first_term(const vector<double> &pseudo_hist,
			      const double a_mid);
  double alpha_score_function(const vector<double> &pseudo_hist,
			      const double mean, const double a_mid,
			      const double vals_count);
  double EM_estim_params(const double tol, const size_t max_iter,
          vector<size_t> &vals_hist); // returns log_like
  void estim_params(const vector<size_t> &val_hist);
  void estim_params(const vector<size_t> &vals_hist,
                    const vector<double> &probs);
  void set_helpers();
  
  /** private member functions for distribution specification **/
  double trunc_log_pdf(const size_t val) const;
  double log_pdf(const size_t val) const;
  double trunc_log_L(const vector<size_t> &vals_hist) const;
  double expected_zeros(const double pseudo_size) const;
  double trunc_pval(const size_t val) const;
  double expected_inverse_sum(const double mean, const size_t sample_size,
                              const size_t sum) const;

  /*** private instance variables ***/
  /** \brief mu, the mean of the (untruncated) distribution **/
  double mu;
  /** \brief alpha, the NB dispersion parameter **/
  double alpha;
  // the rest are helpers to avoid recomputing things
  double n_helper;
  double p_helper;
  double n_log_p_minus_lngamma_n_helper;
  double log_q_helper;
};

#endif

