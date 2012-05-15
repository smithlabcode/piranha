/**
  \file Distribution.hpp
  \brief This source file declares the abstract Distribution class, which
         serves as a base class for all simple distributions

  \authors Philip J. Uren, Qiang Song, Andrew D. Smith

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

#ifndef DISTRO_HPP
#define DISTRO_HPP

#include <vector>
#include <string>
#include <iostream> 
#include <fstream>
#include <math.h>
#include "smithlab_utils.hpp"
#include "SimpleXML.hpp"

/**
 * \brief An exception class to be used for distribution-based exceptions
 */
class DistributionException : public SMITHLABException {
public:
  DistributionException(std::string s = std::string()) : SMITHLABException(s) {}
};

/**
 * \brief The abstract base class for distributions.
 */
class Distribution {
public:  
  /*** constructors, destructors and object initialisation ***/
  // force derived to implement destructor; have to provide empty 
  // implementation here because base destructor is always called
  virtual ~Distribution() {;}
  
  /*** Factory pattern ***/
  static Distribution* create(const std::string& type);
  static Distribution* create(const std::string& fn, const std::string& type);
  
  /*** Simple Inspectors ***/
  virtual std::vector<double> getParams() const = 0;
  virtual std::string toString() const = 0;

  /*** Simple Mutators ***/
  virtual void setParams(const std::vector<double>& params) = 0;
  virtual void randomise() = 0;

  /*** Describing the distribution ***/
  virtual double loglikelihood(const double x) const = 0;
  virtual double pdf(const double x) const { 
    return exp(this->loglikelihood(x)); 
  }
  virtual double cdf(const double x) const {
    double res=0;
    for (size_t i=0; i<=x; i++) res += pdf(i);
    return res;
  }
  virtual double pvalue(const double x) const {
    return 1 - this->cdf(x);
  }

  /*** Distribution I/O ***/
  virtual void save(std::ostream& ostrm) const = 0;
  virtual void load(const SimpleXML& xml) = 0;
  virtual void save(std::string fn) {
    std::ofstream ostrm(fn.c_str());
    this->save(ostrm);
  }
  virtual void load(std::istream& istrm) {
    SimpleXML xml(istrm);
    this->load(xml);
  }
  virtual void load(const std::string& fn) {
    std::ifstream istrm(fn.c_str());
    this->load(istrm);
  }

  /*** Parameter estimation ***/
  virtual void estimateParams(const std::vector<double>& ys,
                              const std::vector<double>& probs,
                              const bool verbose=false) = 0;
  virtual void estimateParams(const std::vector<double>& ys,
                   const bool verbose) {
    this->estimateParams(ys, std::vector<double>(ys.size(), 1), verbose);
  }
  virtual void estimateParams(const std::vector<double>& ys) {
    this->estimateParams(ys, std::vector<double>(ys.size(), 1), false);
  }
};

#endif

