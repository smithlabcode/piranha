/**
  \file Poisson.hpp
  \brief this header file declares the Poisson distribution class

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

#ifndef POISS_HPP
#define POISS_HPP

#include <vector>
#include <string>
#include "Distribution.hpp"

class Poisson : public Distribution {
public:
  /*** constructors, destructors and object initialization ***/
  Poisson() : mean(1) {;}
  Poisson(std::vector<double> params) { this->setParams(params); }
  Poisson(double mean) : mean(mean) {;}
  
  /*** simple inspectors ***/
  std::vector<double> getParams() const;
  SimpleXML toXML() const;
  std::string toXMLString() const;
  std::string toString() const;
  
  /*** simple mutators ***/
  void setParams(const std::vector<double>& params);
  void randomise();

  /*** Specifying the distribtuion ***/
  double loglikelihood(const double x) const;

  /*** Poisson I/O ***/
  void save(std::ostream& ostrm) const;
  void load(const SimpleXML& xml);

  /*** parameter estimation ***/
  void estimateParams(const std::vector<double>& ys,
      const std::vector<double>& probs = std::vector<double>(),
      const bool verbose=false);
private:
  /*** private instance variables ***/
  /** \brief the mean of this Poisson; lambda **/
  double mean;

  /*** private static constants ***/
  /** \brief lambda won't exceed this value when randomise() is called **/
  static const double MAX_RANDOM_LAMBDA = 10;
};

#endif

