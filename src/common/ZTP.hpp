/**
  \file ZTP.hpp This header declares the zero-truncated Poisson distribution

  \authors Philip J. Uren, Timothy Daley, Andrew D. Smith

  \section copyright Copyright Details

  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Timothy Daley, Andrew D. Smith

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

  \section revisions Revision History
**/

#ifndef ZTP_HPP
#define ZTP_HPP

#include <string>
#include <vector>
#include <exception>
#include <limits>

#include "SimpleXML.hpp"
#include "Distribution.hpp"

/**
 * \brief the zero-truncated Poisson distribution
 */
class ZeroTruncatedPoisson : public Distribution {
public:
  /*** constructors, destructors and object initialisation ***/
  ZeroTruncatedPoisson() : lambda(1) {;}
  ZeroTruncatedPoisson(const double lambda) : lambda(lambda) {;}

  /*** simple inspectors ***/
  std::string toString() const;
  std::vector<double> getParams() const;
  SimpleXML toXML() const;
  std::string toXMLString() const;

  /*** describing the distribution ***/
  double loglikelihood(const double x) const;
  double pdf(const double x) const;

  /*** ZTP I/O ***/
  void save(std::ostream& ostrm) const;
  void load(const SimpleXML& xml);

  /*** ZTP mutators ***/
  void setParams(const std::vector<double>& params);
  void randomise();

  /*** Parameter estimation ***/
  void estimateParams(const std::vector<double>& ys,
                      const std::vector<double>& probs,
                      const bool verbose=false);
private:
  /** \brief Lambda, the (untruncated) mean of the distribution and its
             only parameter **/
  double lambda;
  /** \brief tolerance for estimation of lambda by bisection algorithm **/
  static const double tolerance = 1e-20;
  /** \brief max value that can be given to lambda when randomise is called **/
  static const double MAX_RANDOM_LAMBDA = 100;
};

#endif

