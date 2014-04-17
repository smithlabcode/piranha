/**
  \file LinearRegression.hpp
  \brief This header file declares the class for the LinearRegression model

  \author Philip J. Uren

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren

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

#ifndef LINR_HPP
#define LINR_HPP

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "Matrix.hpp"
#include "RegressionModel.hpp"
#include "Gaussian.hpp"

class LinearRegression : public RegressionModel {
private:
  // TODO disable copy constructor for now, implement later
  LinearRegression(LinearRegression& lr);
  
  // private instance variables
  std::vector<double> coefficients;
  double helperVar;
  double helperLogVar;
  Gaussian helperG;
  
  // private constants
  const static double minVar;
  
  // private member functions
  void estimateParamsOLS(const std::vector<double>& response, 
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs = std::vector<double>(),
                   const bool verbose=false);
public:
  // constructor -- initialise all coefficients to 1
  LinearRegression(const size_t P) : 
    coefficients(std::vector<double>(P,1)),
    helperVar(1),
    helperLogVar(0)
  {;}
  
  /**** inspectors ****/
  std::vector<double> getRegressionCoefficients() const {
    return this->coefficients;
  };
  std::vector<double> getParams() const { 
    return this->getRegressionCoefficients(); 
  }
  double loglikelihood(const double response, 
                   const std::vector<double>& covariates) const;
  std::string toString() const;
  double fittedValue(const std::vector<double>& covariates) const;
  double getVariance() const { return this->helperVar; };
  std::string toXML() const {
    std::stringstream ss;
    ss << "<linearRegression>" << std::endl;
    for (size_t i=0; i<this->coefficients.size(); i++) 
      ss << "<coefficient> " << this->coefficients[i] << " </coefficient>" << std::endl;
    ss << this->helperG.toXML() << std::endl;
    ss << "</linearRegression>";
    return ss.str();
  }
  void save(std::ostream& ostrm) const {
    ostrm << this->toXML() << std::endl;
  }
  
  /**** mutators ****/
  void load(std::istream& istrm) {
    SimpleXML xml(istrm);
    this->load(xml);
  }
  void load(SimpleXML& xml);
  void setVariance(double v) {
    this->helperVar = v;
    this->helperLogVar = log(v);
  }
  void randomise(size_t numCovariates);
  void setParams(const std::vector<double>& params) {
    this->setRegressionCoefficients(params);
  }
  void setRegressionCoefficients(const std::vector<double>& params) {
    this->coefficients = params;
  }
  void estimateParams(const std::vector<double>& response, 
                   const std::vector< std::vector<double> >& covariates,
                   const std::vector<double>& probs = std::vector<double>(),
                   const std::vector<double>& startingPoint = std::vector<double>(),
                   const bool verbose=false) {
    this->estimateParamsOLS(response, covariates, probs);
  };
  void estimateParams(const std::vector<double>& response, 
                     const std::vector< std::vector<double> >& covariates,
                     const bool verbose=false) {
    this->estimateParams(response, covariates, std::vector<double>(), 
                         std::vector<double>(), verbose);
  };



  // TODO fix these later..
  double deviance(const std::vector<double>& responses,
            const std::vector<double>& probs,
            const std::vector< std::vector<double> >& covariates) const {
    throw RegressionModelException("Not implemented");
  }

  /*** Specifying the distribution ***/
  std::vector<double> gradient(const std::vector<double>& response,
            const std::vector< std::vector<double> >& covariates,
            const std::vector<double>& probs) const {
    throw RegressionModelException("Not implemented");
  }
  Matrix Hessian(const std::vector<double>& response,
            const std::vector< std::vector<double> >& covariates,
            const std::vector<double>& probs) const {
    throw RegressionModelException("Not implemented");
  }

};

#endif

