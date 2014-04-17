/**
  \file LinearRegression.hpp
  \brief This source file defines the class for the LinearRegression model

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

#define _USE_MATH_DEFINES

#include "LinearRegression.hpp"
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
using std::cerr;
using std::endl;

// initialise static data members
const double LinearRegression::minVar = 1e-100;

/****
 * @summary: load this linear regression model from an XML description 
 */
void LinearRegression::load(SimpleXML& xml) {
  assert(xml.getTagName() == "linearRegression");
  vector<SimpleXML> coeffs = xml.getChildren("coefficient");
  vector<SimpleXML> errDist = xml.getChildren("gaussian");
  assert(errDist.size() == 1);
  this->helperG = Gaussian(errDist[0]);
  this->coefficients = vector<double>(coeffs.size(),0);
  for (size_t i=0; i<coeffs.size(); i++) {
    assert(coeffs[i].isLeaf());
    std::stringstream ss;
    ss << coeffs[i].getData();
    ss >> this->coefficients[i]; 
  }
}


/*****
 * @summary: randomly pick a set of parameters for this model
 */
void
LinearRegression::randomise(size_t numCovariates) {
  double rmax = (double) RAND_MAX;
  this->coefficients.clear();
  for (size_t i=0; i<numCovariates; i++) {
      this->coefficients.push_back((rand() / rmax) * 4);
  }
  this->helperVar = (rand() / rmax) * 100;
  this->helperLogVar = log(this->helperVar);
}


/****
 * @summary: 
 */
double 
LinearRegression::fittedValue(const vector<double>& covariates) const {
  double res = 0;
  for (size_t i=0; i<covariates.size(); i++) {
    res += (covariates[i] * this->coefficients[i]);
  }
  return res;
}


/*****
 * @summary:
 */
void 
LinearRegression::estimateParamsOLS(const vector<double>& y, 
                                    const vector< vector<double> >& covariates,
                                    const vector<double>& probsIn,
                                    const bool verbose) {
  // sanity check
  size_t n = y.size();
  assert(n > 0);
  assert(n == covariates.size());
  size_t p = covariates[0].size();
  assert(p > 0);
  
  // construct covariate matrices  
  Matrix x(covariates);
  Matrix xt = x.transpose();
  
  // weight inputs?
  vector<double> probs;
  if (probsIn.size() != n) probs = vector<double>(n,1);
  else probs = probsIn;
  
  // do regression 
  Matrix zt(y);
  Matrix z = Matrix(zt).transpose();
  DiagonalMatrix w(n,n);
  for (size_t i=0; i<n; i++) w(i,i) = probs[i];
  this->coefficients = ((((xt * w) * x).inverse()) * ((xt * w) * z)).asVector();
  
  vector<double> errs;
  for (size_t i=0; i<n; i++) errs.push_back(y[i] - fittedValue(covariates[i]));
  this->helperG.estimateParams(errs, probs);
  // todo magic number
  if (this->helperG.getParams()[1] <= 1e-100) {
    vector<double> p;
    p.push_back(this->helperG.getParams()[0]);
    p.push_back(1e-100);
    this->helperG.setParams(p);
  }
  
  /*// remember the variance for doing log-likelihood later
  this->helperVar = 0;
  for (size_t i=0; i<n; i++) {
    double err = y[i] - fittedValue(covariates[i]);
    this->helperVar += (probs[i] * err * err);
  }
  this->helperVar /= n;
  // don't allow 0
  if (this->helperVar <= this->minVar) this->helperVar = this->minVar; 
  this->helperLogVar = log(this->helperVar);
  
  //cout << "finished fitting, setting var to " << this->helperVar << " and logvar to " << this->helperLogVar << endl;*/
}



/******
 * @summary: compute the log-likelihood that a single observation (with
 *           associated covariates) comes from this linear regression 
 *           model with it's current set of coefficients. 
 */
double 
LinearRegression::loglikelihood(const double response, 
                   const vector<double>& covariates) const {
  double err = response - fittedValue(covariates);
  return this->helperG.loglikelihood(err);
  /*cout << "response is " << response << " covariates are: ";
  for (size_t i=0; i<covariates.size(); i++) cout << covariates[i] << ", ";
  cout << "var is " << this->helperVar << " logvar is " << this->helperLogVar;
  cout << "loglike is " << (-0.5 * (log(2*M_PI) + this->helperLogVar)) - 
      ((err * err) / (2*this->helperVar)) << endl;*/
  /*return (-0.5 * (log(2*M_PI) + this->helperLogVar)) - 
         ((err * err) / (2*this->helperVar));*/
}


/*****
 * @summary: return a string representation of this object 
 */
string
LinearRegression::toString() const {
  stringstream ss;
  ss << "beta=";
  for (size_t i=0; i<this->coefficients.size(); i++) { 
    ss << this->coefficients[i];
    if (i <this->coefficients.size()-1) ss << ", ";
  }
  ss << " " << this->helperG.toString(); 
  return ss.str();
}

