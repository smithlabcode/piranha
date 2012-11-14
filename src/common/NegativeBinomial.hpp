/**
  \file NegativeBinomial.hpp
  \brief This header file declares the negative binomial distribution class

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

#ifndef NBINOM_HPP
#define NBINOM_HPP

#include <vector>
#include <string>
#include "SimpleXML.hpp"
#include "Distribution.hpp"

class NegativeBinomial : public Distribution {
public:
  /*** constructors ***/
  NegativeBinomial() : mean(0), alpha(1) {this->setHelpers();}
  NegativeBinomial(std::vector<double> params) { this->setParams(params); }
  NegativeBinomial(double mean, double alpha) : mean(mean), alpha(alpha) {
    this->setHelpers();
  }
  
  /*** Simple inspectors ***/
  std::vector<double> getParams() const;
  SimpleXML toXML() const;
  std::string toXMLString() const;
  std::string toString() const;

  /*** Specifying the distribution ***/
  double loglikelihood(const double x) const;

  /*** Negative Binomial I/O ***/
  void save(std::ostream& ostrm) const;
  void load(const SimpleXML& xml);
  
  /*** Simple mutators ***/
  void setParams(const std::vector<double>& params);
  void randomise();

  /*** Parameter estimation ***/
  void estimateParams(const std::vector<double>& ys,
                      const std::vector<double>& probs = std::vector<double>(),
                      const bool verbose=false);
  /*void estimateParams(const std::vector<double>& ys, const bool verbose) {
    this->estimateParams(ys, std::vector<double>(ys.size(),1), verbose);
  }*/

private:
  /*** private constants ***/
  const static double max_allowed_alpha = 100;
  const static double min_allowed_alpha = 1e-20;
  const static double alpha_allowed_error = 1e-10;

  /*** private member functions -- mutators ***/
  void setHelpers();
  void estimateParamsML(const std::vector<double> &vals,
                        const std::vector<double> &probs,
                        const bool verbose=false);

  /*** private instance variables ***/
  double mean;
  double alpha;
  // the rest are just helpers to avoid recomputing things
  double n_helper;
  double p_helper;
  double n_log_p_minus_lngamma_n_helper;
  double log_q_helper;
};

#endif

