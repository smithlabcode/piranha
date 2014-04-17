/**
  \file PoissonRegression.cpp
  \brief this source file defines the PoissonRegression class, which is an
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

#include "PoissonRegression.hpp"
#include "VectorMath.hpp"
#include "Matrix.hpp"
#include <sstream>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <gsl/gsl_sf_gamma.h>

using std::vector;
using std::string;
using std::stringstream;
using std::cout;
using std::cin;
using std::endl;

// initialise static data members
const double PoissonRegression::threshold = 0.000001;
const double PoissonRegression::maxIter = 100;
const double PoissonRegression::tiny = 1e-254;

/******************************************************************************
 * Parameter estimation
 */

/*****
 * @summary: estimate parameters for the Poisson regression weighting each 
 *           input data point with the provided weights and starting with an
 *           initial guess at the regression coefficients. This is pretty easy, 
 *           we just use the IRLS or NR algorithm to estimate beta and we're
 *           done -- both are implemented in the parent class for GLMs.
 */
void 
PoissonRegression::estimateParams(const vector<double>& response, 
                       const vector< vector<double> >& covariates,
                       const vector<double>& probs,
                       const vector<double>& startingPoint,
                       const bool verbose) {
  this->estimateBeta(response, covariates, probs, startingPoint, verbose);
}

/*****
 * @summary: estimating parameters for the Poisson regression without weighting
 *           the data points (i.e. all weights set to 1) and without giving the
 *           algorithm a starting set of parameters.
 */
void 
PoissonRegression::estimateParams(const vector<double>& response, 
                       const vector< vector<double> >& covariates,
                       const bool verbose) {
  const int N = response.size();
  RegressionModel::randomise(this->coefficients.size());
  this->estimateParams(response, covariates, vector<double>(N,1), 
                       this->coefficients, verbose);
}


/******************************************************************************
 * Specifying the distribution
 */

/*****
 * @summary: calculate mu for the ith observation given a set of
 *           regression coefficients
 */
long double
PoissonRegression::mu(size_t i, const vector< vector<double> >& covariates,
                      const vector<double>& coefficients) {
  return exp(vectorDotProduct(covariates[i], coefficients));
}

/******
 * @summary: compute the log-likelihood that a single observation (with
 *           associated covariates) comes from this Poisson regression 
 *           model with it's current set of coefficients. 
 */
double 
PoissonRegression::loglikelihood(const double response, 
                   const vector<double>& covariates) const {
  if (this->coefficients.size() != covariates.size()) {
    stringstream ss;
    ss << "failed calculating log likelihood for Poisson "
       << "regression model. Reason: "
       << "number of coefficients does not match number of covariates";
    throw RegressionModelException(ss.str());
  }

  double eta = 0;
  for (size_t i=0; i<covariates.size(); i++) {
    eta += (covariates[i] * this->coefficients[i]);
  }
  double mu = exp(eta);
  return response * eta - mu - gsl_sf_lnfact(int(response));
}

/****
 * @summary: TODO
 * @note:    Poisson deviance is given in Hilbe, pg83 
 */
double 
PoissonRegression::deviance(const vector<double>& responses,
        const vector<double>& probs,
        const vector< vector<double> >& covariates) const {
  double dev = 0;
  const size_t N = responses.size();
  for (size_t i=0; i<N; i++) {
    double mui = mu(i, covariates, this->coefficients) +\
                 PoissonRegression::tiny;
    double yi = responses[i] + PoissonRegression::tiny;
    dev += (probs[i] * (yi * (log(yi) - log(mui)) - yi + mui));
  }
  dev *= 2;
  return dev;
}

/****
 * @summary: TODO
 */
double
PoissonRegression::linkFirstDerivative(const size_t i,
                    const vector< vector<double> >& covariates) const {
  return 1 / this->mu(i, covariates, this->coefficients);
}

/****
 * @summary: TODO
 */
double
PoissonRegression::variance(const size_t i,
                    const vector< vector<double> >& covariates) const {
  return this->mu(i, covariates, this->coefficients);
}

/*****
 * @summary: calculate the first derivative at the point specified by
 *           the current regression coefficients
 * @return:  vector of length p; return by value but p is generally small (<10)
 * TODO we do not check that response and covariates are of the same length
 * TODO we do not check that all covaraite vectors are the same size
 */
vector<double> 
PoissonRegression::gradient(const vector<double>& response,
                            const vector< vector<double> >& covariates,
                            const vector<double>& probs) const {
  size_t n = response.size();      // num observations
  size_t p = covariates[0].size(); // num covariates
  
  Matrix gradient(1,p);
  for (size_t i=0; i<n; i++) {
    Matrix xi(covariates[i]);
    gradient = gradient + 
        (xi * probs[i] * (response[i] - mu(i, covariates, this->coefficients)));
  }
  
  return gradient.asVector();
}
vector<double>
PoissonRegression::gradient(const vector<double>& response,
                            const vector< vector<double> >& covariates) const {
  vector<double> probs(response.size(),1);
  return gradient(response, covariates, probs);
}

/*****
 * @summary: calculate the Hessian with respect to beta at the point specified
 *           by the current regression coefficients
 * @return:  Hessian matrix; returns by value, but this matrix is p*p, and
 *           p (number of covariates) is generally small (less than 10)
 */
Matrix
PoissonRegression::Hessian(const vector<double>& response,
                           const vector< vector<double> >& covariates,
                           const vector<double>& probs) const {
  size_t n = response.size();      // num observations
  Matrix hessian(this->coefficients.size(), this->coefficients.size());
  for (size_t i=0; i<n; i++) {
    vector< vector<double> > tmp;
    tmp.push_back(covariates[i]);
    

    // we make a row vector of the covariates using a matrix
    // the column vector is then just the transpose of the matrix
    Matrix cov_i_rowVec (tmp);
    Matrix cov_i_colVec = cov_i_rowVec.transpose();

    // add to the sum
    hessian = hessian + 
              ((cov_i_colVec * cov_i_rowVec) * 
               (-1 * probs[i] * mu(i, covariates, this->coefficients)));
  }
  
  // done, return the Hessian matrix
  return hessian;
}
Matrix
PoissonRegression::Hessian(const vector<double>& response,
                           const vector< vector<double> >& covariates) const {
  vector<double> probs(response.size(),1);
  return Hessian(response, covariates, probs);
}


/******************************************************************************
 * Inspectors
 */

/*****
 * @summary: return a string representation of this object
 */
string
PoissonRegression::toString() const {
  stringstream ss;
  ss << "beta: ";
  for (size_t i=0; i<this->coefficients.size(); i++) {
    ss << this->coefficients[i];
    if (i <this->coefficients.size()-1) ss << ", ";
  }
  return ss.str();
}

/*****
 * @summary: return an XML string representation of this object
 */
string
PoissonRegression::toXML() const {
  std::stringstream ss;
  ss << "<PoissonRegression>" << std::endl;
  for (size_t i=0; i<this->coefficients.size(); i++) {
    ss << "<coefficient> " << this->coefficients[i]
       << " </coefficient>" << std::endl;
  }
  ss << "</PoissonRegression>";
  return ss.str();
}


/******************************************************************************
 * Poisson Regression I/O
 */

/****
 * @summary: load this PoissonRegression model from an XML description
 */
void
PoissonRegression::load(SimpleXML& xml) {
  if (xml.getTagName() != "PoissonRegression") {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString()
       << " as Poisson regression model. Reason: "
       << "root tag is not PoissonRegression" << endl;
    throw RegressionModelException(ss.str());
  }

  vector<SimpleXML> coeffs = xml.getChildren("coefficient");
  this->coefficients = vector<double>(coeffs.size(),0);
  for (size_t i=0; i<coeffs.size(); i++) {
    if (!coeffs[i].isLeaf()) {
      stringstream ss;
      ss << "failed parsing " << endl << xml.toString()
         << " as Poisson regression model. Reason: "
         << "coefficient tag is not leaf" << endl;
      throw RegressionModelException(ss.str());
    }
    std::stringstream ss;
    ss << coeffs[i].getData();
    ss >> this->coefficients[i];
  }
}

