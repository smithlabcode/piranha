/****
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Andrew D. Smith

  Authors: Philip J. Uren, Andrew D. Smith

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

  --------------------

  Known Bugs:    None

  Revision
  History:       None

  TODO:          None
****/

#ifndef TRUNC_POISR_HPP
#define TRUNC_POISR_HPP

#include <vector>
#include <string>
#include "Matrix.hpp"
#include "RegressionModel.hpp"

class TruncatedPoissonRegression : public RegressionModel {
public:
  /*** Constructors, destructors and object initialization ***/
  TruncatedPoissonRegression(size_t truncAt, size_t P) :
    RegressionModel(P),
    truncateAt(truncAt) {;}

  /*** Specifying the distribution ***/
  double loglikelihood(const double response,
                  const std::vector<double>& covariates) const;
  std::vector<double> gradient(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs) const;
  Matrix Hessian(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs) const;
  double deviance(const std::vector<double>& responses,
                  const std::vector<double>& probs,
                  const std::vector< std::vector<double> >& covars) const;
  double pdf(const double response, const std::vector<double>& covars) const;

  /*** Simple inspectors ***/
  std::string toString() const;
  SimpleXML toXML() const;
  std::string toXMLString() const;
  std::vector<double> getParams() const;

  /*** Regression model I/O ***/
  void save(std::ostream& ostrm) const;
  void load(std::istream& ostrm);
  void load(SimpleXML& xml);

  /*** simple mutators ***/
  void setParams(const std::vector<double>& params);

  /*** Parameter estimation ***/
  void estimateParams(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs,
                  const std::vector<double>& startingPoint,
                  const bool verbose=false);
  void estimateParams(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const bool verbose=false);

private:
  /*** private instance variables ***/
  size_t truncateAt;

  /*** Private member functions ***/
  double mu(const size_t i,
               const std::vector< std::vector<double> >& covars) const;
  double delta(const size_t i, const std::vector<double>& response,
               const std::vector< std::vector<double> >& covariates) const;

  /*** Private constants ***/
  static const double tiny = 0.00001; // just used in a few places to stop
                                      // things going to zero
};

#endif

