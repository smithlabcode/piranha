/**
  \file PoissonRegression.hpp
  \brief this header file declares the PoissonRegression class, which is an
         implementation of the GeneralisedLinearModel abstract interface.
         We're using the Template Method Pattern so the fitting algorithms
         are given in the superclass and this class just provides the
         implementation of those components specific to the Poisson regression.

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

#ifndef POISR_HPP
#define POISR_HPP

#include <vector>
#include <string>
#include "Matrix.hpp"
#include "RegressionModel.hpp"
#include "GeneralisedLinearModel.hpp"

class PoissonRegression : public GeneralisedLinearModel {
public:
  /* constructors, destructors, object initialization */
  PoissonRegression() {;}
  PoissonRegression(const size_t P) : GeneralisedLinearModel(P) {;}
  PoissonRegression(const size_t P, FittingMethod fitmthd) :
    GeneralisedLinearModel(P, fitmthd) {;}
  
  /* inspectors */
  std::string toString() const;
  std::string toXML() const;
  std::vector<double> getRegressionCoefficients() const {
    return this->coefficients;
  }
  std::vector<double> getParams() const {
    return this->getRegressionCoefficients();
  }

  /* specifying the distribution */
  double loglikelihood(const double response, 
                  const std::vector<double>& covariates) const;

  double linkFirstDerivative(const size_t i,
                  const std::vector< std::vector<double> >& covariates) const;
  double variance(const size_t i,
                  const std::vector< std::vector<double> >& covariates) const;
  double deviance(const std::vector<double>& responses,
                  const std::vector<double>& probs,
                  const std::vector< std::vector<double> >& covariates) const;
  
  /* Poisson Regression IO */
  void load(std::istream& istrm) {
      SimpleXML xml(istrm);
      this->load(xml);
  }
  void load(SimpleXML&);
  void save(std::ostream& ostrm) const { ostrm << this->toXML() << std::endl; }

  /* Simple mutators */
  void setRegressionCoefficients(const std::vector<double, 
      std::allocator<double> >& c) {
    this->coefficients = c;
  }
  void setParams(const std::vector<double>& params) {
    this->setRegressionCoefficients(params);
  }

  /* Parameter estimation */
  void estimateParams(const std::vector<double>& response, 
               const std::vector< std::vector<double> >& covariates,
               const std::vector<double>& probs,
               const std::vector<double>& startingPoint,
               const bool verbose=false);
  void estimateParams(const std::vector<double>& response, 
               const std::vector< std::vector<double> >& covariates,
               const bool verbose);
private:
  /* private constants */
  static const double threshold;
  static const double maxIter;
  static const double tiny;
  
  /* static private functions */
  // TODO remove static
  static long double mu(size_t i,
               const std::vector< std::vector<double> >& covariates,
               const std::vector<double>& coefficients);

  
  /* private member functions for parameter estimation */
  void estimateParamsNR(const std::vector<double>& response, 
               const std::vector< std::vector<double> >& covariates,
               const bool verbose=false);
  void estimateParamsNR(const std::vector<double>& response, 
               const std::vector< std::vector<double> >& covariates,
               const std::vector<double>& probs,
               const bool verbose=false);
  void estimateParamsIRLS(const std::vector<double>& response, 
               const std::vector< std::vector<double> >& covariates,
               const std::vector<double>& probs = std::vector<double>(),
               const std::vector<double>& startingPoint = std::vector<double>(),
               const bool verbose=false);

  /* private member functions defining the distribution */
  std::vector<double> gradient(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs) const;
  std::vector<double> gradient(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates) const;
  Matrix Hessian(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs) const;
  Matrix Hessian(const std::vector<double>& response,
                   const std::vector< std::vector<double> >& covariates) const;
};

#endif

