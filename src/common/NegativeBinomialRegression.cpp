/**
  \file NegativeBinomialRegression.cpp
  \brief This source file defines the negative binomial regression model class

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

#include <cmath>
#include <sstream>
#include <numeric>
#include <cassert>
#include <iostream>

#include <gsl/gsl_sf_gamma.h>

#include "NegativeBinomialRegression.hpp"
#include "PoissonRegression.hpp"
#include "VectorMath.hpp"
#include "Matrix.hpp"

using std::stringstream;
using std::accumulate;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

const double NegativeBinomialRegression::threshold;

/******************************************************************************
 * Simple inspectors
 */

/**
 * \brief Get alpha, the dispersion parameter, for this negative binomial
 *        regression model.
 */
double
NegativeBinomialRegression::getAlpha() const {
  return this->alpha;
}

/**
 * \brief return a string representation of this object
 */
string
NegativeBinomialRegression::toString() const {
  stringstream ss;
  ss << "beta: ";
  for (size_t i=0; i<this->coefficients.size(); i++)
    ss << this->coefficients[i] << ", ";
  ss << " --- alpha: " << this->alpha;
  return ss.str();
}

/**
 * \brief return an XML string representation of this object
 */
string
NegativeBinomialRegression::toXML() const {
  stringstream ss;
  ss << "<NegativeBinomialRegression>" << endl;
  for (size_t i=0; i<this->coefficients.size(); i++) {
    ss << "<coefficient> " << this->coefficients[i]
       << " </coefficient>" << endl;
  }
  ss << "<alpha>" << this->alpha << "</alpha>" << endl;
  ss << "</NegativeBinomialRegression>";
  return ss.str();
}

/**
 * \brief Get the parameters for this negative binomial regression model.
 *        alpha is the first value in the returned vector, the rest are the
 *        covariates
 */
vector<double>
NegativeBinomialRegression::getParams() const {
  vector<double> res;
  res.push_back(alpha);
  for (size_t i=0; i<this->coefficients.size(); i++) {
    res.push_back(this->coefficients[i]);
  }
  return res;
}

/******************************************************************************
 * Specifying the distribution
 */

/**
 * \brief calculate mu for the ith observation given a set of regression
 *        coefficients
 */
double
NegativeBinomialRegression::mu(size_t i, const vector< vector<double> >& covar,
                               const vector<double>& coefficients) {
  double res = exp(vectorDotProduct(covar[i], coefficients));
  // bound mu
  //if (res > 1000) res = 1000;
  return res;
}

/**
 * \brief TODO
 * \todo  check the math here
 * \todo  shouldn't this be the weighted version and the non-weighted is in
 *        the superclass?
 */
double
NegativeBinomialRegression::loglikelihood(const double response,
                                      const vector<double>& covariates) const {
  if (this->coefficients.size() != covariates.size()) {
    stringstream ss;
    ss << "failed calculating log likelihood for negative binomial "
       << "regression model. Reason: "
       << "number of coefficients does not match number of covariates";
    throw RegressionModelException(ss.str());
  }

  double eta = 0;
  for (size_t i=0; i<covariates.size(); i++) {
    eta += (covariates[i] * this->coefficients[i]);
  }
  double mu = exp(eta);
  double delta = 1 + alpha * mu;

  return ((response * log((alpha*mu) / (delta))) -
          (1/alpha) * log(delta) +
          gsl_sf_lngamma(response + 1/alpha) -
          gsl_sf_lngamma(response + 1) -
          gsl_sf_lngamma(1/alpha));
}

/**
 * \brief TODO
 * \todo give reference for this
 */
double
NegativeBinomialRegression::variance(const size_t i,
                    const vector< vector<double> >& covariates) const {
  double mui = mu(i, covariates, this->coefficients);
  return mui + (this->alpha * mui * mui);
}

/**
 * \brief TODO
 */
double
NegativeBinomialRegression::linkFirstDerivative(const size_t i,
                    const vector< vector<double> >& covariates) const {
  return 1 / this->mu(i, covariates, this->coefficients);
}

/**
 * \brief TODO
 *
 *        NB2 deviance is given in Hilbe on page 211, but excludes prior
 *        weights. See the following web-site for the weighted version:
 *        See http://www.sfu.ca/sasdoc/sashtml/stat/chap29/sect27.htm
 * \todo  fix magic number;
 * \todo  check discrepancy with alpha vs 1/alpha in t1 and t2
 */
double
NegativeBinomialRegression::deviance(const vector<double>& responses,
        const vector<double>& probs,
        const vector< vector<double> >& covariates) const {
  const size_t N = responses.size();
  double dev = 0;
  for (size_t i=0; i<N; i++) {
    double mui = mu(i, covariates, this->coefficients) + 0.00001;
    double yi = responses[i] + 0.00001;
    double oneOnAlpha = 1 / (this->alpha + 0.00001);
    double t0 = (yi + oneOnAlpha);
    double t1 = log(yi + (alpha));
    double t2 = log(mui + (alpha));
    dev += (probs[i] * ((yi * (log(yi) - log(mui))) -
            (t0 * (t1 - t2))));
  }
  dev *= 2;
  return dev;
}

/**
 * \brief calculate the first derivative at the point specified by the
 *        current regression coefficients
 * \note  this returns the gradient by value, which requires a copy, but as
 *        the vector is only of length P (number of covariates), and P
 *        is generally very small (less than 10), this is not a significant
 *        performance issue.
 * \todo  check that response and covariates are of the same length
 * \todo  check that all covariate vectors are the same size
 */
vector<double>
NegativeBinomialRegression::gradient(const vector<double>& response,
                            const vector< vector<double> >& covariates,
                            const vector<double>& probs) const {
  const size_t P = covariates[0].size(); // number of covariates
  const size_t N = response.size();      // number of observations

  Matrix gradient(1,P);
  for (size_t i=0; i<N; i++) {
    Matrix xi(covariates[i]);
    double mui = mu(i, covariates, coefficients);
    double frac = ((response[i] - mui) / (1+alpha*mui));
    gradient = gradient + (xi * frac * probs[i]);
  }

  return gradient.asVector();
}

/**
 * \brief calculate the Hessian with respect to beta
 *        at the point specified by the given regression coefficients
 * \note  this returns the Hessian by value, which requires a copy, but as
 *        the matrix is only of size P*P (P = number of covariates), and P
 *        is generally very small (less than 10), this is not a significant
 *        performance issue.
 * \todo  check that response and covariates are of the same length
 * \todo  check that all covariate vectors are the same size
 */
Matrix
NegativeBinomialRegression::Hessian(const vector<double>& response,
                           const vector< vector<double> >& covariates,
                           const vector<double>& probs) const {
  const size_t N = response.size();      // number of observations

  // Initialize result
  Matrix hessian(coefficients.size(), coefficients.size());

  // calculate Hessian
  for (size_t i=0; i<N; i++) {
    double mui = mu(i, covariates, coefficients);
    double denom = (1 + alpha*mui) * (1 + alpha*mui);
    double num = -1 * mui * (response[i] * alpha + 1);
    double frac = num / denom;
    Matrix xi(covariates[i]);
    Matrix xt = xi.transpose();
    Matrix tmp = ((xt * xi) * frac);
    hessian = hessian + ((xt * xi) * frac * probs[i]);
  }

  return hessian;
}

/******************************************************************************
 * Mutators
 */

/**
 * \brief Set alpha, the dispersion parameter, for this negative binomial
 *        regression model.
 */
void
NegativeBinomialRegression::setAlpha(const double alpha) {
  this->alpha = alpha;
}

/**
 * \brief Set the parameters for this Negative Binomial Regression.
 *
 * The first parameter must be alpha, then the rest are taken as beta (the
 * regression coefficients)
 */
void 
NegativeBinomialRegression::setParams(const std::vector<double>& params) {
  this->alpha = params.front();
  this->coefficients.clear();
  this->setRegressionCoefficients(vector<double>(params.begin() + 1, 
                                                 params.end()));
}

/**
 * \brief randomly pick a set of parameters for this model
 */
void
NegativeBinomialRegression::randomise(size_t numCovariates) {
  double rmax = (double) RAND_MAX;
  this->coefficients.clear();
  for (size_t i=0; i<numCovariates; i++) {
      this->coefficients.push_back((rand() / rmax) * 4);
  }
  this->alpha = (rand() / rmax) * 4;
}

/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief TODO
 */
void 
NegativeBinomialRegression::estimateParams(const vector<double>& response, 
                       const vector< vector<double> >& covariates,
                       const vector<double>& probs,
                       const vector<double>& startingPoint,
                       const bool verbose) {
  double n = response.size();
  double p = covariates[0].size();
  double df = accumulate(probs.begin(), probs.end(), 0.0) - p;

  // begin by fitting a Poisson regression model
  if (verbose) {
    cerr << "initial parameter estimate using Poisson Regression" << endl;
  }
  PoissonRegression reg(this->coefficients.size(), this->fittingMethod);
  reg.estimateParams(response, covariates, probs, startingPoint, verbose);
  this->setRegressionCoefficients(reg.getRegressionCoefficients());
  if (verbose) {
    cerr << "finished initial Poisson estimate" << endl;
  }

  // calculate the dispersion from poisson model
  double chi2 = 0;
  for (size_t i=0; i<n; i++) {
    double mui = mu(i, covariates, this->getRegressionCoefficients());
    chi2 += (probs[i] * (((response[i] - mui) * (response[i] - mui)) / mui));
  }
  double dispersion = chi2 / df;
  double phi = 1 / dispersion;
  double oldDisp = 0;
  double deltaDisp = NegativeBinomialRegression::threshold + 1;

  // iteratively estimate alpha using dispersion dampening, then fit
  // beta with fixed alpha.
  while (fabs(deltaDisp) > NegativeBinomialRegression::threshold) {
    oldDisp = dispersion;
    this->alpha = phi;
    this->estimateBeta(response, covariates, probs,
                       this->getRegressionCoefficients(), verbose);
    chi2 = 0;
    for (size_t i=0; i<n; i++) {
      double mui = mu(i, covariates, this->coefficients);
      double sigma2 = (response[i] - mui) * (response[i] - mui);
      chi2 += (probs[i] * (sigma2 / (mui + this->alpha*mui*mui)));
    }
    dispersion = chi2 / df;
    phi = dispersion * phi;
    deltaDisp = dispersion - oldDisp;

    if (verbose) {
      cerr << "finished dispersion dampening iteration; current model "
           << "estimate: " << this->toString() << " - " << "old dispersion "
           << "was: " << oldDisp << " - new dispersion is: " << dispersion
           << " - delta dispersion is: " << deltaDisp << endl;
    }
  }

  if (verbose) {
    cerr << "AIC: " << this->AkaikeInformationCriterion(response, covariates)
         << endl;
         //<< "BIC: " << this->BayesianInformationCriterion(response, covariates)
         //<< endl;
  }
}

/**
 * \brief TODO
 */
void 
NegativeBinomialRegression::estimateParams(const std::vector<double>& response, 
                 const std::vector< std::vector<double> >& covariates,
                 const bool verbose) {
  const int N = response.size();
  this->randomise(this->coefficients.size());
  this->estimateParams(response, covariates, vector<double>(N,1), 
                       this->coefficients, verbose);
}

/******************************************************************************
 * Negative Binomial Regression I/O
 */

/**
 * \brief load this NegativeBinomialRegression model from an XML description
 */
void 
NegativeBinomialRegression::load(SimpleXML& xml) {
  if (xml.getTagName() != "NegativeBinomialRegression") {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString()
       << " as negative binomial regression model. Reason: "
       << "root tag is not NegativeBinomialRegression" << endl;
    throw RegressionModelException(ss.str());
  }
  
  // get alpha
  vector<SimpleXML> alphaTags = xml.getChildren("alpha");
  if ((alphaTags.size() != 1) || (!alphaTags[0].isLeaf())) {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString()
       << " as negative binomial regression model. Reason: "
       << "missing or incorrectly formed alpha tag" << endl;
    throw RegressionModelException(ss.str());
  }
  std::stringstream ssTmp;
  ssTmp << alphaTags[0].getData();
  ssTmp >> this->alpha; 
  
  // get beta
  vector<SimpleXML> coeffs = xml.getChildren("coefficient");
  this->coefficients = vector<double>(coeffs.size(),0);
  for (size_t i=0; i<coeffs.size(); i++) {
    if (!coeffs[i].isLeaf()) {
      stringstream ss;
      ss << "failed parsing " << endl << xml.toString()
         << " as negative binomial regression model. Reason: "
         << "coefficient tag is not leaf" << endl;
      throw RegressionModelException(ss.str());
    }
    std::stringstream ss;
    ss << coeffs[i].getData();
    ss >> this->coefficients[i]; 
  }
}

/******************************************************************************
 * Goodness of fit measures
 */

/**
 * \brief AkaikeInformationCriterion (AIC) for the model with its current
 *        parameters given the responses <y> and covaraites <x>
 * \todo  weighted version of this?
 * \todo  this is duplicated in several classes, consolidate and check
 */
double
NegativeBinomialRegression::AkaikeInformationCriterion(const vector<double>& y,
                              const vector< vector<double> >& x) {
  const double P = this->coefficients.size();
  const double N = y.size();
  return (-2 * (this->loglikelihood(y,x) - P)) / N;
}

/*****
 * \brief BaysianInformationCriterion (BIC) for the model with its current
 *        parameters given the responses <y> and covariates <x>
 * \todo  weighted version of this?
 * \todo  this is duplicated in several classes, consolidate and check
 */
double 
NegativeBinomialRegression::BaysianInformationCriterion(const vector<double>& y,
                              const vector< vector<double> >& x) {
  //const double P = this->coefficients.size();
  //const double N = y.size();
  string msg = "BIC not yet implemented for Negative Binomial Regression.";
  throw RegressionModelException(msg);
}

