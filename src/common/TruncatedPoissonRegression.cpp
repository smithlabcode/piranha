/**
  \file TruncatedPoissonRegression.cpp
  \brief This source file defines the truncated Poisson regression class, which
         is an implementation of the RegressionModel abstract interface.
         We're using the Template Method Pattern so the fitting algorithms are
         given in the superclass and this class just provides the
         implementation of those components specific to the Truncated Poisson
         regression.

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

#include "TruncatedPoissonRegression.hpp"
#include "PoissonRegression.hpp"
#include "VectorMath.hpp"
#include "Matrix.hpp"

#include <sstream>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <float.h>
#include <string>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::isfinite;
using std::stringstream;

const double TruncatedPoissonRegression::tiny;

/******************************************************************************
 *  Specifying the distribution
 */

/**
 * \brief calculate mu for the ith observation given a set of regression
 *        coefficients
 * \param covariates matrix of covariate values
 */
double
TruncatedPoissonRegression::mu(const size_t i,
                           const vector< vector<double> >& covariates) const {
  double res = exp(vectorDotProduct(covariates[i], this->coefficients));
  if (!(isfinite(res))) {
    stringstream ss;
    ss << "Calculating inverse link function for truncated Poisson regression "
       << "failed. Reason: resultant mu is non-finite. Regression "
       << "coefficients were ";
    for (size_t j=0; j<this->coefficients.size(); j++) {
      ss << this->coefficients[j];
      if (i<this->coefficients.size()-1) ss << ",";
    }
    ss << " covariate values were ";
    for (size_t j=0; j<covariates[i].size(); j++) {
      ss << covariates[i][j];
      if (j<covariates[i].size()-1) ss << ",";
    }
    throw RegressionModelException(ss.str());
  }
  return res;
}

/**
 * \brief Delta is a helper function used in the calculation of the
 *        truncated Poisson likelihood function derivatives (gradient and
 *        Hessian).
 */
double
TruncatedPoissonRegression::delta(const size_t i,
                  const vector<double>& response,
                  const vector< vector<double> >& covariates) const {
  double mui = this->mu(i, covariates);
  double pdf = gsl_ran_poisson_pdf(this->truncateAt, mui);
  double cdf = gsl_cdf_poisson_P(this->truncateAt, mui);
  if (cdf == 1) cdf = 0.99; // save us from DBZ below
  double res = (mui * pdf) / (1 - cdf);

  // sanity check
  if (!(isfinite(res))) {
    stringstream ss;
    ss << "Calculating delta function for truncated Poisson regression "
       << "failed. Reason: resultant value is non-finite";
    throw RegressionModelException(ss.str());
  }

  return res;
}

/**
 * \brief Evaluate the truncated Poisson log-likelihood function at the
 *        point specified by the given response and associated covariates
 */
double
TruncatedPoissonRegression::loglikelihood(const double response,
                                          const vector<double>& covars) const {
  if (this->coefficients.size() != covars.size()) {
    stringstream ss;
    ss << "failed calculating log likelihood for Poisson "
       << "regression model. Reason: "
       << "number of coefficients does not match number of covariates";
    throw RegressionModelException(ss.str());
  }

  // although we allow the response to be a double, it must contain
  // a positive, non-zero integer
  double intpart;
  if ((!(modf(response, &intpart) == 0.0)) || (response < 1)) {
    stringstream ss;
    ss << "evaluating truncated Poisson regression "
       << "log-likelihood function with response " << response << " failed. "
       << "Reason: " << response << " is not a positive non-zero integer";
    throw RegressionModelException(ss.str());
  }

  double eta = 0;
  for (size_t i=0; i<covars.size(); i++) {
    eta += (covars[i] * this->coefficients[i]);
  }
  double mu = exp(eta);
  double pdf = gsl_ran_poisson_pdf(size_t(response), mu) +\
                   TruncatedPoissonRegression::tiny;
  double cdf = gsl_cdf_poisson_P(this->truncateAt,mu);
  double res = log(pdf) - log(1 - cdf);

  // sanity check on result
  if (!isfinite(res)) {
    stringstream ss;
    ss << "evaluating truncated Poisson regression "
       << "log-likelihood function with response " << response << " and "
       << "distribution parameters " << this->toString() << " failed. "
       << "Reason: result was non-finite";
    throw RegressionModelException(ss.str());
  }

  return res;
}

/**
 * \brief The probability mass function for the zero-truncated Poisson
 *        regression model.
 * \param response TODO
 * \param covariates TODO
 *
 * We override this here because we need to handle zero carefully -- otherwise
 * we'll get undefined behaviour from the log-likelihood function at zero.
 */
double
TruncatedPoissonRegression::pdf(const double response,
                     const std::vector<double>& covariates) const {
  if (response == 0) return 0;
  return exp(this->loglikelihood(response, covariates));
}

/**
 * \brief Calculate the first derivative of the log-likelihood function at
 *        the point specified by the given response and covariates weighted
 *        by the given probabilities
 * \todo  check that all covariate vectors are the same size
 */
vector<double>
TruncatedPoissonRegression::gradient(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs) const {
  if (response.size() != covariates.size()) {
    stringstream ss;
    ss << "failed evaluating truncated Poisson regression gradient. "
       << "Reason: size of response vector does not match size of covariate "
       << "matrix";
    throw RegressionModelException(ss.str());
  }

  const size_t N = response.size();      // number of observations
  const size_t P = covariates[0].size(); // number of covariates

  // TODO can initialize matrix directly...
  vector< vector<double> > tmp(1,vector<double>(P,0));
  Matrix res(tmp);

  for (size_t i=0; i<N; i++) {
    double mui = this->mu(i, covariates);
    double deltai = this->delta(i,response,covariates);
    Matrix xi(covariates[i]);
    res = res + (xi * (response[i] - mui - deltai) * probs[i]);
  }

  return res.asVector();
}

/**
 * \brief Calculate the second derivative of the log-likelihood function at
 *        the point specified by the given response and covariates weighted
 *        by the given probabilities
 * \todo  check that all covariates are the same size
 * \todo  check that response and covariates are the same length
 */
Matrix
TruncatedPoissonRegression::Hessian(const std::vector<double>& response,
                  const std::vector< std::vector<double> >& covariates,
                  const std::vector<double>& probs) const {
  size_t n = response.size();      // num observations
  Matrix hessian(coefficients.size(), coefficients.size());
  for (size_t i=0; i<n; i++) {
    vector< vector<double> > tmp;
    tmp.push_back(covariates[i]);

    // we make a row vector of the covariates using a matrix
    // the column vector is then just the transpose of the matrix
    Matrix cov_i_rowVec (tmp);
    Matrix cov_i_colVec = cov_i_rowVec.transpose();

    // evaluate mu and delta at i
    double mui = this->mu(i, covariates);
    double deltai = delta(i, response, covariates);

    // add to the sum
    double frac = -1 * probs[i] * (mui +\
                  (deltai * (this->truncateAt - mui - deltai + 1)));
    hessian = hessian + ((cov_i_colVec * cov_i_rowVec) * frac);
  }

  // done, return the Hessian matrix
  return hessian;
}

/**
 * \brief Calculate the deviance for this truncated Poisson regression
 *        model. This is a log-likelihood ratio between the saturated
 *        model and the current model.
 * \note  The truncated deviance is just the Poisson deviance plus an
 *        adjustment for truncation.
 */
double
TruncatedPoissonRegression::deviance(const std::vector<double>& responses,
                  const std::vector<double>& probs,
                  const std::vector< std::vector<double> >& covars) const {
  double adjustment = 0;
  const size_t N = responses.size();
  for (size_t i=0; i<N; i++) {
    double mui = this->mu(i, covars);
    double cdfMu = gsl_cdf_poisson_P(this->truncateAt, mui);
    double cdfY = gsl_cdf_poisson_P(this->truncateAt, responses[i]);
    adjustment += (log(1 - cdfMu) - log(1 - cdfY));
  }
  adjustment *= 2;

  PoissonRegression p(this->getNumCoefficients());
  p.setRegressionCoefficients(this->getRegressionCoefficients());
  return p.deviance(responses, probs, covars) + adjustment;
}


/******************************************************************************
 *  Simple inspectors
 */

/**
 * \brief  Return a string representation of this object
 */
string
TruncatedPoissonRegression::toString() const {
  stringstream ss;
  ss << "beta: ";
  for (size_t i=0; i<this->coefficients.size(); i++)
    ss << this->coefficients[i] << ", ";
  return ss.str();
}

/**
 * \brief Return an XML representation of this object
 */
SimpleXML
TruncatedPoissonRegression::toXML() const {
  return SimpleXML(this->toXMLString());
}

/**
 * \brief Return an XML representation of this object as a string
 */
string
TruncatedPoissonRegression::toXMLString() const {
  std::stringstream ss;
  ss << "<TruncatedPoissonRegression>" << std::endl;
  ss << "<truncateAt> " << this->truncateAt << " </truncateAt>" << endl;
  for (size_t i=0; i<this->coefficients.size(); i++) {
    ss << "<coefficient> " << this->coefficients[i]
       << " </coefficient>" << std::endl;
  }
  ss << "</TruncatedPoissonRegression>";
  return ss.str();
}

/**
 * \brief  get the set of parameters for this truncated Poisson regression
 *         model
 * \return vector of double, first value is the truncation point, rest are
 *         the regression parameters
 */
vector<double>
TruncatedPoissonRegression::getParams() const {
  vector<double> res;
  res.push_back(this->truncateAt);
  for (size_t i=0; i<this->coefficients.size(); i++)
    res.push_back(this->coefficients[i]);
  return res;
}

/******************************************************************************
 * Regression model I/O
 */

/**
 * \brief write this truncated Poisson regression model to the given output
 *        stream
 */
void
TruncatedPoissonRegression::save(std::ostream& ostrm) const {
  ostrm << this->toXMLString() << std::endl;
}

/**
 * \brief read the values for this truncated Poisson regression model from the
 *        given input stream
 */
void
TruncatedPoissonRegression::load(std::istream& istrm) {
  SimpleXML xml(istrm);
  this->load(xml);
}

/**
 * \brief load the member variable values for this truncated Poisson
 *        regression model from an XML representation
 */
void
TruncatedPoissonRegression::load(SimpleXML& xml) {
  if (xml.getTagName() != "TruncatedPoissonRegression") {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString() << endl
       << " as truncated Poisson regression model. Reason: "
       << "root tag is not TruncatedPoissonRegression" << endl;
    throw RegressionModelException(ss.str());
  }

  // get the truncation point
  vector<SimpleXML> truncTags = xml.getChildren("truncateAt");
  if ((truncTags.size() != 1) || (!truncTags[0].isLeaf())) {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString() << endl
       << " as truncated Poisson regression model. Reason: "
       << "missing or incorrectly formed truncation tag" << endl;
    throw RegressionModelException(ss.str());
  }
  std::stringstream ssTmp;
  ssTmp << truncTags[0].getData();
  ssTmp >> this->truncateAt;

  // get beta
  vector<SimpleXML> coeffs = xml.getChildren("coefficient");
  this->coefficients = vector<double>(coeffs.size(),0);
  for (size_t i=0; i<coeffs.size(); i++) {
    if (!coeffs[i].isLeaf()) {
      stringstream ss;
      ss << "failed parsing " << endl << xml.toString() << endl
         << " as truncated Poisson regression model. Reason: "
         << "coefficient tag is not leaf" << endl;
      throw RegressionModelException(ss.str());
    }
    std::stringstream ss;
    ss << coeffs[i].getData();
    ss >> this->coefficients[i];
  }
}

/******************************************************************************
 * Simple mutators
 */

/**
 * \brief Set the parameters for this truncated Poisson model.
 * \param params First parameter in the vector is assumed to be the value for
 *                the truncation point, the rest are the regression parameters.
 * \todo throw exception if truncation point is invalid (must be non-negative
 *       integer)
 */
void
TruncatedPoissonRegression::setParams(const vector<double>& params) {
  if (params.size() != this->getNumCoefficients() + 1) {
    stringstream ss;
    ss << "Failed to set parameters for truncated Poisson regression model. "
       << "Reason: wrong number of parameters. Found " << params.size()
       << " but expected " << (this->getNumCoefficients() + 1);
    throw RegressionModelException(ss.str());
  }
  this->truncateAt = size_t(params[0]);
  vector<double> tmp(params.begin() + 1, params.end());
  this->setRegressionCoefficients(tmp);
}

/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief Estimate parameters for this truncated Poisson model, given a
 *        set of weights (probs) for the observed data and a set of
 *        parameter values (startingPoint) to start searching from.
 */
void
TruncatedPoissonRegression::estimateParams(const vector<double>& response,
                                    const vector< vector<double> >& covariates,
                                    const vector<double>& probs,
                                    const vector<double>& startingPoint,
                                    const bool verbose) {
  this->estimateBetaNR(response, covariates, probs, startingPoint, verbose);
}

/**
 * \brief Estimate parameters of this truncated Poisson model given a set
 *        of observed data in the form of responses and their associated
 *        covariate values.
 * \todo should move this up to RegressionModel
 */
void
TruncatedPoissonRegression::estimateParams(const std::vector<double>& response,
                                    const vector< vector<double> >& covariates,
                                    const bool verbose) {
  const size_t P = this->getNumCoefficients();
  this->randomise(P);
  this->estimateParams(response, covariates, vector<double>(P,1),
                       this->getRegressionCoefficients(), verbose);
}



