/**
  \file MixtureModel.hpp
  \brief This header file declares an abstract base class for mixture models
         in the package. This defines a consistent interface to all mixtures.

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

#ifndef MIXMOD_HPP
#define MIXMOD_HPP

#include <string>
#include <vector>
#include <exception>
#include "smithlab_utils.hpp"
#include "FittingMethod.hpp"

class MixtureModelException : public SMITHLABException {
public:
  MixtureModelException(std::string s = std::string()) : SMITHLABException(s) {}
};

class MixtureModel {
public:
  // force derived to implement destructor; have to provide empty 
  // implm here because base destructor is always called
  virtual ~MixtureModel() {;}
  
  // Factory pattern 
  static MixtureModel* create(const std::string& type, 
                           const size_t numComponents,
                           const size_t numCovariates=0,
                           const FittingMethod ftmthd=FittingMethod("DEFAULT"));
  static MixtureModel* create(const std::string& fn,
                           const std::string& type,
                           const size_t numComponents,
                           const size_t numCovariates=0,
                           const FittingMethod ftmthd=FittingMethod("DEFAULT"));
  
  // mutators
  virtual void estimateParams(const std::vector<double>& responses,
                              const std::vector<std::vector<double> >& covars, 
                              bool verbose=false) = 0;
  virtual void estimateParams(const std::vector<double>& obs, 
                              bool verbose=false) { 
    throw MixtureModelException("Not Implemented"); 
  }
  virtual void load(const std::string& fn) = 0;
  
  // inspectors, const methods
  // Note: we don't force derived classes to implement the likelihood methods
  //       as some don't make sense in certain situations -- but we have them 
  //       here to provide a consistent interface.
  virtual double loglikelihood(const std::vector<double> &responses,
                               const std::vector<std::vector<double> > &covars, 
                               bool debug=false) const {
    throw MixtureModelException("Not Implemented");
  }
  virtual double loglikelihood(const std::vector<double> &obs) const {
    throw MixtureModelException("Not Implemented");
  }
  virtual double getComponentProb(const std::vector<double>& instance, 
      const size_t component) const = 0;
  virtual std::vector<double> getComponentProbs(
      const std::vector<double>& instance) const = 0;
  virtual std::vector<double> getComponentParameters(const int k) const = 0;
  virtual std::string getComponentAsString(size_t k) const = 0;
  virtual size_t numberOfComponents() const = 0;
  virtual void save(std::string& fn) const = 0;
  virtual void save(std::ostream& fh) const = 0;
  virtual std::string toString() const = 0; 
};

#endif
