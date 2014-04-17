/**
  \file RegressionModel.hpp
  \brief This is an abstract base class for all regression models. It defines
         certain common members that must be implemented by all regression
         model subclasses and gives implementations for those operations that
         will be unchanged amongst subclasses

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

#ifndef REGMDL_HPP
#define REGMDL_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "smithlab_utils.hpp"
#include "Matrix.hpp"
#include "SimpleXML.hpp"

class RegressionModelException : public SMITHLABException {
public:
  RegressionModelException(std::string s = std::string()) :
    SMITHLABException(s) {}
};

class RegressionModel {
public:
  /*** Constructors, destructors and object initialization ***/
  RegressionModel() {;}
  RegressionModel(const size_t P) :
      coefficients(std::vector<double>(P,1)) {;}
  // force destruction via derived class
  virtual ~RegressionModel() {;}
  
  /*** Specifying the distribution ***/
  virtual double loglikelihood(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covars) const;
  virtual double loglikelihood(const double response, 
                  const std::vector<double>& covariates) const = 0;
  virtual double loglikelihood(const std::vector<double>& responses,
                  const std::vector<double>& probs,
                  const std::vector< std::vector<double> >& covars) const;
  virtual std::vector<double> gradient(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs) const = 0;
  virtual Matrix Hessian(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs) const = 0;
  virtual double deviance(const std::vector<double>& responses,
                  const std::vector<double>& probs,
                  const std::vector< std::vector<double> >& covars) const = 0;
  virtual double deviance(const std::vector<double>& responses,
                  const std::vector< std::vector<double> >& covars) const;
  virtual double pvalue(const double response,
                  const std::vector<double>& covariates) const;
  virtual double pdf (const double response,
                  const std::vector<double>& covariates) const;
  virtual double cdf(const double response,
                  const std::vector<double>& covariates) const;

  /*** Simple inspectors ***/
  virtual std::vector<double> getRegressionCoefficients() const;
  virtual size_t getNumCoefficients() const;
  virtual std::string toString() const = 0;
  // get all parameters, including underlying distribution parameters
  // (if any) and regression parameters. Order determined by derived classes
  virtual std::vector<double> getParams() const = 0;

  /*** Regression model I/O ***/
  virtual void save(std::ostream& ostrm) const = 0;
  virtual void save(std::string fn) const;
  virtual void load(std::istream& ostrm) = 0;
  virtual void load(SimpleXML& xml) = 0;
  virtual void load(std::string fn);

  /*** simple mutators ***/
  virtual void randomise(const size_t numCovariates);
  virtual void setRegressionCoefficients(const std::vector<double>& params);
  // set all parameters, including underlying distribution parameters
  // (if any) and regression parameters. Order determined by derived classes
  virtual void setParams(const std::vector<double>& params) = 0;
  virtual void estimateParams(const std::vector<double>& response, 
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs,
                  const std::vector<double>& startingPoint,
                  const bool verbose=false) = 0;
  virtual void estimateParams(const std::vector<double>& response, 
                  const std::vector< std::vector<double> >& covariates,
                  const bool verbose=false) = 0;

  /*** Parameter estimation ***/
  virtual void estimateBeta(const std::vector<double>& response,
      const std::vector< std::vector<double> >& covariates,
      const std::vector<double>& probs,
      const std::vector<double>& startingPoint,
      const bool verbose=false);
  virtual void estimateBetaNR(const std::vector<double>& response,
      const std::vector< std::vector<double> >& covariates,
      const std::vector<double>& probs,
      const std::vector<double>& startingPoint,
      const bool verbose=false);

protected:
  /*** protected instance variables ***/
  std::vector<double> coefficients;

  /*** protected constants ***/
  static const double devianceThreshold;
  static const double maxLinearPredictorValue;
  static const double minLinearPredictorValue;

private:
  /*** private constants ***/
  static const size_t maxNRIters = 100;
};

#endif

