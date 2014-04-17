/**
  \file ZTNBRegression.hpp
  \brief This header declares the zero-truncated version of the negative
         binomial regression.

  This is not a GLM, since it isn't an exponential family distribution, but
  it's still a regression model.

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

#ifndef ZTNBR_HPP
#define ZTNBR_HPP

#include <vector>
#include <string>
#include "Matrix.hpp"
#include "RegressionModel.hpp"

/** \brief The zero-truncated version of the negative binomial regression. **/
class ZTNBRegression : public RegressionModel {

public:
  /*** constructors, destructors and object initialization ***/
  ZTNBRegression(const size_t P) :
    RegressionModel(P), alpha(0.00001) {;}

  /*** describing the distribution ***/
  double loglikelihood(const double response,
                      const std::vector<double>& covariates) const;
  using RegressionModel::loglikelihood; // else the above hides these
  double zeroProbability(size_t i,
                         const std::vector<std::vector<double> > &c) const;
  std::vector<double> gradient(const std::vector<double>& response,
                      const std::vector< std::vector<double> >& covariates,
                      const std::vector<double>& probs) const;
  Matrix Hessian(const std::vector<double>& response,
                      const std::vector< std::vector<double> >& covariates,
                      const std::vector<double>& probs) const;
  double deviance(const std::vector<double>& responses,
                      const std::vector<double>& probs,
                      const std::vector< std::vector<double> >& covars) const;
  double expectedValue(const size_t i,
                      const std::vector< std::vector<double> >& covars) const;
  double variance(const size_t i,
                      const std::vector< std::vector<double> >& covars) const;
  double chiSquaredStatistic(const std::vector<double>& response,
                             const std::vector< std::vector<double> >& covars,
                             const std::vector<double>& probs) const;
  double pdf(const double response, const std::vector<double>& covars) const;
  
  /*** Simple inspectors ***/
  double getAlpha() const;
  std::string toString() const;
  SimpleXML toXML() const;
  std::string toXMLString() const;
  std::vector<double> getParams() const;

  /*** ZTNB model I/O ***/
  void save(std::ostream& ostrm) const;
  void load(std::istream& ostrm);
  void load(SimpleXML& xml);

  /*** simple mutators ***/
  void setAlpha(const double alpha);
  void setParams(const std::vector<double>& params);

  /*** Parameter estimation ***/
  void estimateParams(const std::vector<double>& response,
                      const std::vector< std::vector<double> >& covariates,
                      const std::vector<double>& probs,
                      const std::vector<double>& startingPoint,
                      const bool verbose=false);
  // TODO move up to Regression Model?
  void estimateParams(const std::vector<double>& response,
                      const std::vector< std::vector<double> >& covariates,
                      const bool verbose);
  
  /*** goodness of fit statistics ***/
  // TODO these should be moved up the hierarchy, and checked!
  double BayesianInformationCriterion(const std::vector<double>& y,
                   const std::vector< std::vector<double> >& x) const;
  double AkaikeInformationCriterion(const std::vector<double>& y,
                   const std::vector< std::vector<double> >& x) const;

private:
  /*** private constants ***/
  static const double dispersonThreshold;
  static const double maxAlpha;
  static const double minAlpha;
  static const double maxZeroProb;
  static const double minZeroProb;
  static const double minVariance;
  static const double maxExpectedValue;
  /** \brief floor for those calculations where zero cannot be allowed **/
  static const double minZeroTolerance;
  /** \brief ceiling for those calculations where one cannot be allowed **/
  static const double maxOneTolerance;
  static const size_t maxDispersionDampeningIters;

  /*** private instance variables ***/
  double alpha;

  /*** private member functions ***/
  double mu(const size_t i,
            const std::vector< std::vector<double> >& covariates) const;
};

#endif

