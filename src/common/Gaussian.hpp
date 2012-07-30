/**
  \file Gaussian.hpp
  \brief This header file defines the class for the Gaussian distribution

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

#ifndef GUASS_HPP
#define GUASS_HPP

#include <vector>
#include <string>
#include <sstream>
#include "Distribution.hpp"
#include "SimpleXML.hpp"

class Gaussian : public Distribution {
private:
  // prviate constants
  static const double threshold = 0.000001; 
  static const double maxIter=100;
  static const double tiny = 1e-254;
  
  // private instance variables
  double mean;
  double variance;
  
public:
  // constructors 
  Gaussian() : mean(0), variance(1) {;}
  Gaussian(std::vector<double> params) { this->setParams(params); }
  Gaussian(double mean, double variance) : mean(mean), variance(variance) {;}
  Gaussian(SimpleXML& xml) { this->load(xml); }
  
  // inspectors, const methods
  std::vector<double> getParams() const {
    std::vector<double> r;
    r.push_back(this->mean);
    r.push_back(this->variance);
    return r;
  };
  double loglikelihood(const double x) const;
  std::string toString() const {
    std::stringstream s;
    s << "mean=" << this->mean << " variance=" << this->variance;
    return s.str();
  }
  std::string toXML() const {
    std::stringstream ss;
    ss << "<gaussian>" << "<mean>" << this->mean << "</mean>"
       << "<variance>" << this->variance << "</variance>"
       << "</gaussian>";
    return ss.str();
  }
  void save(std::ostream& ostrm) const {
    ostrm << this->toXML();
  }
  
  // mutators
  void load(const SimpleXML& xml);
  void setParams(const std::vector<double>& params) {
    assert(params.size() == 2);
    this->mean = params[0];
    this->variance = params[1];
  }
  void randomise() {
    // TODO magic number 4... 
    this->mean = (rand() / double(RAND_MAX)) * 4;
    this->variance = (rand() / double(RAND_MAX)) * 4;
  }
  void estimateParams(const std::vector<double>& ys,
      const std::vector<double>& probs = std::vector<double>(),
      const bool verbose=false);
};

#endif

