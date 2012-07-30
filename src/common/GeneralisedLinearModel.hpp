/**
  \file GeneralisedLinearModel.hpp
  \brief This is abstract base class defining an interface for GLMs; it
         inherits a lot from RegressionModel and is basically just a place to
         put the IRLS algorithm and its accoutrement so that it needn't be
         re-written for each GLM implementation.

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

#ifndef GLM_HPP
#define GLM_HPP

#include <vector>
#include <string>
#include "RegressionModel.hpp"
#include "FittingMethod.hpp"
#include "Matrix.hpp"

class GeneralisedLinearModel : public RegressionModel {
public:
  // constructors 
  GeneralisedLinearModel() :
      fittingMethod("DEFAULT") {;}
  GeneralisedLinearModel(FittingMethod fitmthd) :
      fittingMethod(fitmthd) {;}
  GeneralisedLinearModel(const size_t P) :
      RegressionModel(P),
      fittingMethod("DEFAULT")  {;}
  GeneralisedLinearModel(const size_t P, FittingMethod fitmthd) :
      RegressionModel(P),
      fittingMethod(fitmthd)  {;}
  
  // inspectors
  virtual FittingMethod getFittingMethod() { return this->fittingMethod; }


  /* specifying the distribution */
  virtual double linkFirstDerivative(const size_t i,
        const std::vector< std::vector<double> >& covariates) const = 0;
  virtual double variance(const size_t i,
        const std::vector< std::vector<double> >& covariates) const = 0;
  
  // mutators
  virtual void estimateBeta(const std::vector<double>& response,
      const std::vector< std::vector<double> >& covariates,
      const std::vector<double>& probs,
      const std::vector<double>& startingPoint,
      const bool verbose=false);
  virtual void setFittingMethod(FittingMethod f) { this->fittingMethod = f; }

protected:
  // protected methods
  virtual void estimateBetaIRLS(const std::vector<double>& response,
      const std::vector< std::vector<double> >& covariates,
      const std::vector<double>& probs = std::vector<double>(),
      const std::vector<double>& startingPoint = std::vector<double>(),
      const bool verbose=false);

  // protected constants 
  static const size_t maxIRLSIters = 100;
  
  // protected instance variables
  FittingMethod fittingMethod;
};

#endif

