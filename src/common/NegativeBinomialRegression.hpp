/**
  \file NegativeBinomialRegression.cpp
  \brief This header declares the negative binomial regression model class

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
****/

#ifndef NBR_HPP
#define NBR_HPP

#include <vector>
#include <string>
#include "Matrix.hpp"
#include "RegressionModel.hpp"
#include "GeneralisedLinearModel.hpp"

/** \brief The negative binomial regression model
 *
 *         This class is an implementation of the GeneralisedLinearModel
 *         abstract interface. We're using the Template Method Pattern so the
 *         fitting algorithms are given in the superclass and this class just
 *         provides the implementation of those components specific to the
 *         NegativeBinomial regression.
 */
class NegativeBinomialRegression : public GeneralisedLinearModel {
public:
  /*** constructors, destructors, object initialization ***/
  // TODO fix magic numbers
  NegativeBinomialRegression() :
    GeneralisedLinearModel(), alpha(0.00001) {;}
  NegativeBinomialRegression(const size_t P) :
    GeneralisedLinearModel(P), alpha(0.00001) {;}
  NegativeBinomialRegression(const size_t P, FittingMethod fitmthd) :
    GeneralisedLinearModel(P, fitmthd), alpha(0.00001){;}
  
  /*** simple inspectors ***/
  double getAlpha() const;
  std::string toString() const;
  std::string toXML() const;
  std::vector<double> getRegressionCoefficients() const {
    return this->coefficients;
  }
  std::vector<double> getParams() const;

  /*** specifying the distribution ***/
  double loglikelihood(const double response, 
                   const std::vector<double>& covariates) const;
  using GeneralisedLinearModel::loglikelihood; // else the above hides these
  double deviance(const std::vector<double>& responses,
                   const std::vector<double>& probs,
                   const std::vector< std::vector<double> >& covariates) const;
  double linkFirstDerivative(const size_t i,
                   const std::vector< std::vector<double> >& covariates) const;
  double variance(const size_t i,
                   const std::vector< std::vector<double> >& covariates) const;
  // first and second derivatives with respect to beta -- needed for NR
  std::vector<double> gradient(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs) const;
  Matrix Hessian(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs) const;


  /*** Negative Binomial regression I/O ***/
  void save(std::ostream& ostrm) const { ostrm << this->toXML() << std::endl; }
  // TODO move to cpp,
  void load(std::istream& istrm) {
      SimpleXML xml(istrm);
      this->load(xml);
  }
  void load(SimpleXML&);

  /*** goodness of fit metrics ***/
  // TODO force from superclass?
  double AkaikeInformationCriterion(const std::vector<double>& y,
                   const std::vector< std::vector<double> >& x);
  double BaysianInformationCriterion(const std::vector<double>& y,
                   const std::vector< std::vector<double> >& x);
  
  /*** mutators ***/
  void setAlpha(const double alpha);
  void randomise(size_t numCovariates);
  void setRegressionCoefficients(const std::vector<double> & c) {
    this->coefficients = c;
  }
  void setParams(const std::vector<double>& params); 

  /*** Parameter Estimation ***/
  void estimateParams(const std::vector<double>& response, 
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs,
                   const std::vector<double>& startingPoint,
                   const bool verbose=false);
  void estimateParams(const std::vector<double>& response, 
                   const std::vector< std::vector<double> >& covariates,
                   const bool verbose); 
  // TODO make this private?
  // TODO change this so don't need to pass alpha, but rather uses member
  /*void estimateBetaNR(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates,
                   const double alpha);*/
  
private:
  /*** private constants ***/
  /** \brief TODO **/
  static const double threshold = 0.001; 
  
  /*** private instance variables ***/
  /** \brief the negative binomial dispersion parameter, alpha **/
  double alpha;
  
  /*** static private functions ***/
  // TODO make non-static
  static double mu(size_t i,
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& coefficients);
};


#endif

