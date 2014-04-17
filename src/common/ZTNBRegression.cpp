/**
  \file ZTNBRegression.cpp
  \brief This source file defines the zero-truncated version of the negative
         binomial regression.

  This is not a GLM, since it isn't an exponential family distribution, but
  it's still a regression model.

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

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <numeric>
#include <limits>

#include "ZTNBRegression.hpp"
#include "TruncatedPoissonRegression.hpp"
#include "NegativeBinomialRegression.hpp"
#include "PoissonRegression.hpp"
#include "VectorMath.hpp"
#include "Matrix.hpp"
#include "SimpleXML.hpp"

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::isfinite;
using std::accumulate;
using std::stringstream;

const double ZTNBRegression::dispersonThreshold = 0.00000001;
const double ZTNBRegression::maxAlpha = 100;
const double ZTNBRegression::minAlpha = 1e-3;
const double ZTNBRegression::maxZeroProb = 0.9999999999;
const double ZTNBRegression::minZeroProb = 0.0000000001;
const double ZTNBRegression::minVariance = 0.0000000001;
const double ZTNBRegression::maxExpectedValue = 10000;
const double ZTNBRegression::minZeroTolerance = 0.01;
const double ZTNBRegression::maxOneTolerance = 0.99;
const size_t ZTNBRegression::maxDispersionDampeningIters = 100;

/******************************************************************************
 * Simple Inspectors
 */

/**
 * \brief Get alpha, the dispersion parameter, for this zero-truncated
 *        negative binomial regression model.
 */
double
ZTNBRegression::getAlpha() const {
  return this->alpha;
}

/**
 * \brief Create and return a string representation of this object
 */
string
ZTNBRegression::toString() const {
  stringstream ss;
  ss << "beta: ";
  for (size_t i=0; i<this->coefficients.size(); i++)
    ss << this->coefficients[i] << ", ";
  ss << " --- alpha: " << this->alpha;
  return ss.str();
}

/**
 * \brief Create and return a SimpleXML object which represents this
 *        ZTNBRegression
 */
SimpleXML
ZTNBRegression::toXML() const {
  return SimpleXML(this->toXMLString());
}

/**
 * \brief Create and return a string representation, in XML format, of this
 */
string
ZTNBRegression::toXMLString() const {
  stringstream ss;
  ss << "<ZeroTruncatedNegativeBinomialRegression>" << endl;
  for (size_t i=0; i<this->coefficients.size(); i++) {
    ss << "<coefficient> " << this->coefficients[i]
       << " </coefficient>" << endl;
  }
  ss << "<alpha>" << this->alpha << "</alpha>" << endl;
  ss << "</ZeroTruncatedNegativeBinomialRegression>";
  return ss.str();
}

/**
 * \brief Get the parameters for this zero-truncated negative binomial
 *        regression model.
 *
 * Alpha is the first value in the returned vector, the rest are coefficients
 */
vector<double>
ZTNBRegression::getParams() const {
  vector<double> res;
  res.push_back(alpha);
  for (size_t i=0; i<this->coefficients.size(); i++) {
    res.push_back(this->coefficients[i]);
  }
  return res;
}

/******************************************************************************
 * Simple mutators
 */

/**
 * \brief Set the parameters for this zero-truncated negative binomial
 *        regression model.
 *
 * The first value in the vector is alpha, the rest are coefficients.
 */
void
ZTNBRegression::setParams(const vector<double>& ps) {
  this->alpha = ps[0];
  this->setRegressionCoefficients(vector<double>(ps.begin()+1, ps.end()));
}

/**
 * \brief Set alpha, the dispersion parameter, for this zero-truncated
 *        negative binomial regression model.
 */
void
ZTNBRegression::setAlpha(const double alpha) {
  this->alpha = alpha;
}

/******************************************************************************
 * Describing the distribution
 */

/**
 * \brief Calculate mu for the ith observation given a set of regression
 *        coefficients
 */
double
ZTNBRegression::mu(const size_t i,
                   const vector< vector<double> >& covariates) const {
  double res = exp(vectorDotProduct(covariates[i], this->coefficients));

  // make sure mu hasn't become infinite
  // in theory this can't happen, since we bound beta in the base class,
  // but we check just in case...
  if (!isfinite(res)) {
    stringstream cvars;
    for (size_t j=0; j<covariates[i].size(); j++) {
      cvars << covariates[i][j];
      if (i != covariates.size() - 1) cvars << ",";
    }
    stringstream ss;
    ss << "evaluating inverse log-link of zero-truncated negative binomial "
       << "regression with covariate values " << cvars.str() << " and "
       << "parameters " << this->toString() << " failed. "
       << "Reason: result was non-finite";
    throw RegressionModelException(ss.str());
  }

  return res;
}

/**
 * \brief what is the probability of seeing a zero in the non-truncated
 *        parent distribution?
 */
double
ZTNBRegression::zeroProbability(size_t i, const vector<vector<double> > &cvrs) const {
  double mui = this->mu(i, cvrs);
  double zeroProb = pow(1+alpha*mui,-1/alpha);
  // fix this so it won't go to zero or one -- we'll have issues with logs
  // of zeroProb and 1-zeroProb if we don't do this...
  if (zeroProb > maxZeroProb) zeroProb = ZTNBRegression::maxZeroProb;
  if (zeroProb < minZeroProb) zeroProb = ZTNBRegression::minZeroProb;
  return zeroProb;
}

/**
 * \brief Calculate the expected value for y_i, given its covariates
 * \param covars the full set of covariate values (as a reference). Only the
 *               the ith vector will be used.
 *
 * The expected value for the zero-truncated negative binomial is given on page
 * 119 in Cameron and Trivedi
 */
double
ZTNBRegression::expectedValue(const size_t i,
                              const vector< vector<double> > &covars) const {
  double mui = mu(i, covars);
  double denom = 1 - zeroProbability(i, covars);
  // no need to check for denom == 0, since zeroProb is never equal to 1
  double res = mui / denom;
  if (res > ZTNBRegression::maxExpectedValue)
    res = ZTNBRegression::maxExpectedValue;
  return res;
}

/**
 * \brief Return the variance of the zero-truncated negative binomial
 *        for the given response with the given set of covariates
 *
 *        The variance for the zero-truncated negative binomial is given on
 *        page 119 in Cameron and Trivedi, but seems to be wrong. See
 *        http://vosesoftware.com/ModelRiskHelp/index.htm#Distributions/
 *        Discrete_distributions/Zero-truncated_Negative_Binomial_equations.htm
 *        instead.
 * \todo  throw exception if NaN / INF
 */
double
ZTNBRegression::variance(const size_t i,
                         const vector< vector<double> >& covariates) const {
  double mui = mu(i, covariates);
  double pi = 1 / (1-zeroProbability(i, covariates));
  double varNB = mui + mui * mui * alpha;
  double p = mui * alpha / (1 + alpha * mui);
  double r = 1 / alpha;
  double res = pi * varNB * (1 + p*r) - pi*pi * varNB * p * r;

  if (!isfinite(res)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "variance function with "
       << "parameters " << this->toString() << " failed. "
       << "Reason: result was non-finite";
    cerr << "covar one " << covariates[i][0] << endl;
    cerr << "mu " << mui << endl;
    cerr << "pi " << pi << endl;
    cerr << "varNB " << varNB << endl;
    cerr << "p " << p << endl;
    cerr << "r " << r << endl;
    cerr << "varZTNB " << res << endl;
    throw RegressionModelException(ss.str());
  }

  return res;
}

/**
 * \brief Calculate the Pearson Chi^2 statistic for this zero-truncated
 *        negative binomial regression model
 *
 *        Uses the current values for alpha and regression coefficients and the
 *        data (response and associated covariates) given as parameters. The
 *        Chi^2 statistic is defined on page 52 of Hilbe
 *
 * \todo  could be moved up in the class hierarchy if everything is forced to
 *        define variance?
 */
double
ZTNBRegression::chiSquaredStatistic(const vector<double>& response,
                                    const vector< vector<double> >& covars,
                                    const vector<double>& probs) const {
  const size_t N = response.size();
  double chi2 = 0;
  for (size_t i=0; i<N; i++) {
    double expect = this->expectedValue(i,covars);
    double sigma2 = (response[i] - expect) * (response[i] - expect);

    // can't let variance go to zero here, or we'll get DBZ error below
    double var = variance(i,covars);
    if (var < ZTNBRegression::minVariance) var = ZTNBRegression::minVariance;

    chi2 += (probs[i] * (sigma2 / var));
    /*cerr << "chi2 " << chi2 << endl;
    cerr << "response " << response[i] << endl;
    cerr << "expect" << expect << endl;
    cerr << probs[i] << endl;
    cerr << sigma2 << endl;
    cerr << var << endl;*/
  }

  // sanity check
  if (!isfinite(chi2)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "chi2 function with "
       << "parameters " << this->toString() << " failed. "
       << "Reason: result was non-finite";
    throw RegressionModelException(ss.str());
  }

  return chi2;
}

/**
 * \brief The probability mass function for the zero-truncated Negative
 *        Binomial regression model.
 * \param response TODO
 * \param covariates TODO
 *
 * We override this here because we need to handle zero carefully -- otherwise
 * we'll get undefined behaviour from the log-likelihood function at zero.
 */
double
ZTNBRegression::pdf(const double response,
                     const std::vector<double>& covariates) const {
  if (response == 0) return 0;
  return exp(this->loglikelihood(response, covariates));
}

/**
 * \brief Evaluate the log-likelihood function for this zero-truncated negative
 *        binomial regression model
 * \todo remove coupling with NB class from this function
 */
double
ZTNBRegression::loglikelihood(const double response,
                              const vector<double>& covariates) const {
  // although we allow the response to be a double, it must contain
  // a positive, non-zero integer
  double intpart;
  if ((!(modf(response, &intpart) == 0.0)) || (response < 1)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "log-likelihood function with response " << response << " failed. "
       << "Reason: " << response << " is not a positive non-zero integer";
    throw RegressionModelException(ss.str());
  }

  double mu = exp(vectorDotProduct(covariates, this->coefficients));

  // do this by hand, can't call the function because we don't know i
  // and we don't have the full matrix of covariates here
  double zeroProb = pow(1+alpha*mu,-1/alpha);
  if (zeroProb > maxZeroProb) zeroProb = ZTNBRegression::maxZeroProb;
  if (zeroProb < minZeroProb) zeroProb = ZTNBRegression::minZeroProb;

  double adjustment = log(1 - zeroProb);
  NegativeBinomialRegression nb2 (coefficients.size());
  nb2.setAlpha(this->getAlpha());
  nb2.setRegressionCoefficients(this->getRegressionCoefficients());
  double res = nb2.loglikelihood(response, covariates) - adjustment;

  // sanity check on result
  if (!isfinite(res)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "log-likelihood function with response " << response << " and "
       << "distribution parameters " << this->toString() << " failed. "
       << "Reason: result was non-finite";
    throw RegressionModelException(ss.str());
  }

  return res;
}

/**
 * \brief calculate the weighted deviance for this zero-truncated negative
 *        binomial regression model given a set of weighted responses and
 *        their associated covariates
 * \todo  provide reference for this
 */
double
ZTNBRegression::deviance(const vector<double>& responses,
                         const vector<double>& probs,
                         const vector< vector<double> >& covars) const {
  double adjustment = 0;
  const size_t N = responses.size();
  for (size_t i=0; i<N; i++) {
    double mui = this->mu(i, covars);
    double oneAlphaMu = pow(1 + this->alpha * mui, -1/this->alpha);
    double oneAlphaY = pow(1 + this->alpha * responses[i], -1/this->alpha);

    // make sure we don't take the log of zero below
    if (oneAlphaMu > ZTNBRegression::maxOneTolerance)
      oneAlphaMu = ZTNBRegression::maxOneTolerance;
    if (oneAlphaY > ZTNBRegression::maxOneTolerance)
      oneAlphaY = ZTNBRegression::maxOneTolerance;

    adjustment += (probs[i] * (log(1 - oneAlphaMu) - log(1 - oneAlphaY)));
  }
  adjustment *= 2;
  if (!isfinite(adjustment)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "deviance function with "
       << "parameters " << this->toString() << " failed. "
       << "Reason: truncation adjustment was non-finite";
    throw RegressionModelException(ss.str());
  }

  NegativeBinomialRegression nb2;
  nb2.setAlpha(this->getAlpha());
  nb2.setRegressionCoefficients(this->getRegressionCoefficients());

  double nb2dev = nb2.deviance(responses, probs, covars);
  // sanity check on nb2 deviance
  if (!isfinite(nb2dev)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "deviance function with "
       << "distribution parameters " << this->toString() << " failed. "
       << "Reason: nb2 variance was non-finite";
    throw RegressionModelException(ss.str());
  }

  double res = nb2dev + adjustment;

  // sanity check on result
  if (!isfinite(res)) {
    stringstream ss;
    ss << "evaluating zero-truncated negative binomial regression "
       << "deviance function with "
       << "distribution parameters " << this->toString() << " failed. "
       << "Reason: result was non-finite";
    throw RegressionModelException(ss.str());
  }

  return res;
}

/**
 * \brief calculate the first derivative at the point specified by
 *        the given regression coefficients weighted by the given response
 *        weights
 */
vector<double> 
ZTNBRegression::gradient(const vector<double>& response,
                         const vector< vector<double> >& covariates,
                         const vector<double>& probs) const {
  if (response.size() != covariates.size()) {
    stringstream ss;
    ss << "Calculating gradient vector for zero-truncated negative binomial "
       << "regression model failed. Reason: number of responses does not "
       << "match number of covariate vectors";
    throw RegressionModelException(ss.str());
  }
  if (covariates.size() == 0) {
    stringstream ss;
    ss << "Calculating gradient vector for zero-truncated negative binomial "
       << "regression model failed. Reason: empty response vector and "
       << "covariate matrix";
    throw RegressionModelException(ss.str());
  }

  const size_t P = covariates[0].size(); // number of covariates
  const size_t N = response.size();      // number of observations
  
  if (P != this->coefficients.size()) {
    stringstream ss;
    ss << "Calculating gradient vector for zero-truncated negative binomial "
       << "regression model failed. Reason: covariate matrix dimension "
       << "doesn't match coefficient vector dimension";
    throw RegressionModelException(ss.str());
  }

  // initialize our result
  vector<double> res;
  for (size_t j=0; j<P; j++) res.push_back(0);
  Matrix adjustmentVec(res);
  
  // calculate the adjustment to the NB2 gradient
  for (size_t i=0; i<N; i++) {
    if (covariates[i].size() != P) {
      stringstream ss;
      ss << "Calculating gradient vector for zero-truncated negative binomial "
         << "regression model failed. Reason: covariates matrix is ragged";
      throw RegressionModelException(ss.str());
    }

    double mui = mu(i, covariates);
    double oneAlphaMui = 1 + alpha*mui;
    double denom = pow((oneAlphaMui), 1 + 1/alpha) - oneAlphaMui;
    double num = mui;

    // be careful with the division below..don't want frac == INF
    if (denom < ZTNBRegression::minZeroTolerance)
      denom = ZTNBRegression::minZeroTolerance;

    double frac = num / denom;
    Matrix xi(covariates[i]);
    adjustmentVec = adjustmentVec + (xi * frac * probs[i]);
  }
  
  // calculate nb2 case and subtract correction for zero-truncation
  NegativeBinomialRegression nb2 (coefficients.size());
  nb2.setAlpha(this->getAlpha());
  nb2.setRegressionCoefficients(this->getRegressionCoefficients());
  Matrix gradNBVec(nb2.gradient(response, covariates, probs));
  return (gradNBVec - adjustmentVec).asVector();
}

/**
 * \brief calculate the Hessian with respect to beta
 *        at the point specified by the given regression coefficients
 */
Matrix
ZTNBRegression::Hessian(const vector<double>& response,
                        const vector< vector<double> >& covariates,
                        const vector<double>& probs) const {
  // num observations
  const size_t N = response.size();
  
  // init result matrix
  Matrix adjustmentMatrix(coefficients.size(), coefficients.size());
  
  // calculate Hessian
  for (size_t i=0; i<N; i++) {
    double mui = mu(i, covariates);
    double delta = (1 + alpha * mui);
    double t = pow(delta, 1/alpha);
    double to = t * delta;
    double num = mui * ((to - delta) - ((mui + alpha*mui) * t - alpha*mui));
    double denom = (to - delta) * (to - delta);

    // t may overflow for small alpha, if it does, it implies that the zero
    // probability is essentially zero and so the adjustment can safely be
    // set to zero, same if the denominator of the fraction overflows to inf
    if (isfinite(t) && isfinite(denom)) {
      // be careful with the division below..don't want frac == INF
      if (denom < ZTNBRegression::minZeroTolerance)
        denom = ZTNBRegression::minZeroTolerance;

      double frac = num / denom;
      Matrix xi(covariates[i]);
      Matrix xt = xi.transpose();
      adjustmentMatrix = adjustmentMatrix + ((xt * xi) * frac * probs[i]);
    }
  }
  
  NegativeBinomialRegression nb2 (coefficients.size());
  nb2.setAlpha(this->getAlpha());
  nb2.setRegressionCoefficients(this->getRegressionCoefficients());
  Matrix hessianNB(nb2.Hessian(response, covariates, probs));

  return (hessianNB - adjustmentMatrix);
}

/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief  Estimate the regression coefficients and the NB scale parameter.
 *
 * uses a dispersion dampening algorithm for the estimation of alpha
 */
void 
ZTNBRegression::estimateParams(const vector<double>& response,
                               const vector< vector<double> >& covariates,
                               const vector<double>& probs,
                               const vector<double>& startingPoint,
                               const bool verbose) {
  const double N = response.size();
  const double P = covariates[0].size();
  double df = accumulate(probs.begin(), probs.end(), 0.0) - P;

  // begin by fitting a truncated Poisson regression model
  if (verbose) {
    cerr << "Initial parameter estimate using truncated Poisson regression"
         << endl;
  }
  TruncatedPoissonRegression reg(0, this->coefficients.size());
  reg.estimateParams(response, covariates, probs, startingPoint, verbose);
  this->setRegressionCoefficients(reg.getRegressionCoefficients());
  if (verbose) {
    cerr << "Finished initial truncated Poisson estimate" << endl;
  }

  // calculate the dispersion from poisson model
  double chi2 = 0;
  for (size_t i=0; i<N; i++) {
    double mui = mu(i, covariates);
    chi2 += (probs[i] * (((response[i] - mui) * (response[i] - mui)) / mui));
  }
  double dispersion = chi2 / df;
  double phi = 1 / dispersion;
  double oldDisp = 0;
  double deltaDisp = ZTNBRegression::dispersonThreshold + 1;

  // iteratively estimate alpha using dispersion dampening, then fit
  // beta with fixed alpha.
  size_t numIters = 0;
  while (fabs(deltaDisp) > ZTNBRegression::dispersonThreshold) {
    oldDisp = dispersion;

    // set alpha, but don't allow it to get too big or too small
    this->alpha = phi;
    if (this->alpha > ZTNBRegression::maxAlpha)
      this->alpha = ZTNBRegression::maxAlpha;
    if (this->alpha < ZTNBRegression::minAlpha)
      this->alpha = ZTNBRegression::minAlpha;

    this->estimateBeta(response, covariates, probs,
                       this->getRegressionCoefficients(), verbose);
    chi2 = this->chiSquaredStatistic(response,covariates,probs);
    dispersion = chi2 / df;
    phi = dispersion * phi;
    deltaDisp = dispersion - oldDisp;

    if (verbose) {
      cerr << "finished dispersion dampening iteration; current model "
           << "estimate: " << this->toString() << " - " << "old dispersion "
           << "was: " << oldDisp << " - new dispersion is: " << dispersion
           << " - delta dispersion is: " << deltaDisp << endl;
    }

    numIters += 1;
    if (numIters > ZTNBRegression::maxDispersionDampeningIters) break;
  }

  if (verbose) {
  cerr << "AIC: " << this->AkaikeInformationCriterion(response, covariates)
       << endl
       << "BIC: " << this->BayesianInformationCriterion(response, covariates)
       << endl;
  }
}

/**
 * \brief estimate parameters for this zero-truncated negative binomial
 *        regression model given the unweighted responses and their
 *        associated covariates
 * TODO Move up to RegressionModel?
 */
void
ZTNBRegression::estimateParams(const std::vector<double>& response,
                 const std::vector< std::vector<double> >& covariates,
                 const bool verbose) {
  const int N = response.size();
  this->randomise(this->coefficients.size());
  this->estimateParams(response, covariates, vector<double>(N,1),
                       this->coefficients, verbose);
}

/******************************************************************************
 * Goodness of fit metrics
 */

/**
 * \brief TODO
 * \todo check this!
 * \todo move up to RegressionModel
 * \todo should this be weighted by probs?
 */
double
ZTNBRegression::BayesianInformationCriterion(const vector<double>& y,
                              const vector< vector<double> >& x) const {
  const size_t P = this->getNumCoefficients();
  const size_t N = y.size();
  return (-2 * (this->loglikelihood(y,x) - P*log(P))) / N;
}

/**
 * \brief TODO
 * \todo check this!
 * \todo move up to RegressionModel
 * \todo should this be weighted by probs?
 */
double
ZTNBRegression::AkaikeInformationCriterion(const vector<double>& y,
                              const vector< vector<double> >& x) const {
  const size_t P = this->getNumCoefficients();
  const size_t N = y.size();
  return (-2 * (ZTNBRegression::loglikelihood(y,x) - P)) / N;
}

/******************************************************************************
 * Zero-truncated negative binomial I/O
 */

/**
 * \brief Save this Zero-truncated negative binomial regression model to
 *        the given ostream. Saves in XML format.
 */
void
ZTNBRegression::save(std::ostream& ostrm) const {
  ostrm << this->toXMLString() << std::endl;
}

/**
 * \brief Load this zero-truncated negative binomial regression model from
 *        the given istream. Expects XML format.
 */
void
ZTNBRegression::load(std::istream& istrm) {
  SimpleXML xml(istrm);
  this->load(xml);
}

/**
 * \brief Load this zero-truncated negative binomial regression model from
 *        an XML object
 */
void
ZTNBRegression::load(SimpleXML& xml) {
  if (xml.getTagName() != "ZeroTruncatedNegativeBinomialRegression") {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString()
       << " as zero-truncated negative binomial regression model. Reason: "
       << "root tag is not ZeroTruncatedNegativeBinomialRegression" << endl;
    throw RegressionModelException(ss.str());
  }

  // get alpha
  vector<SimpleXML> alphaTags = xml.getChildren("alpha");
  if ((alphaTags.size() != 1) || (!alphaTags[0].isLeaf())) {
    stringstream ss;
    ss << "failed parsing " << endl << xml.toString()
       << " as zero-truncated negative binomial regression model. Reason: "
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

