/**
  \file RegressionModel.cpp
  \brief This is an abstract base class for all regression models. It defines
  certain common members that must be implemented by all regression model
  subclasses and gives implementations for those operations that will
  be unchanged amongst subclasses

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

#include <vector>
#include <numeric>
#include <sstream>

// from Piranha common
#include "RegressionModel.hpp"
#include "VectorMath.hpp"

using std::vector;
using std::string;
using std::endl;
using std::cerr;
using std::accumulate;
using std::stringstream;

const double RegressionModel::devianceThreshold;
const double RegressionModel::maxLinearPredictorValue;
const double RegressionModel::minLinearPredictorValue;

/******************************************************************************
 * Simple mutators
 */

/**
 * \brief Set the regression coefficients, beta, for this regression model
 * \todo  check this isn't unduly overridden in NegativeBinomialRegression
 * \todo  check this isn't unduly overridden in ZTNBRegression
 * \todo  check this isn't unduly overridden in PoissonRegression
 * \todo  check this isn't unduly overridden in TruncatedPoissonRegression
 */
void
RegressionModel::setRegressionCoefficients(const vector<double>& params) {
  this->coefficients = params;
}

/**
 * \brief TODO
 * \todo  check this isn't unduly overridden in NegativeBinomialRegression
 * \todo  check this isn't unduly overridden in ZTNBRegression
 * \todo  check this isn't unduly overridden in PoissonRegression
 * \todo  check this isn't unduly overridden in TruncatedPoissonRegression
 */
void
RegressionModel::randomise(const size_t numCovariates) {
  double rmax = (double) RAND_MAX;
  this->coefficients.clear();
  for (size_t i=0; i<numCovariates; i++) {
    this->coefficients.push_back((rand() / rmax) * 4);
  }
}


/******************************************************************************
 * Simple inspectors
 */

/**
 * \brief Get a copy of the regression coefficients for this regression model
 * \todo  check this isn't unduly overridden in NegativeBinomialRegression
 * \todo  check this isn't unduly overridden in ZTNBRegression
 * \todo  check this isn't unduly overridden in PoissonRegression
 * \todo  check this isn't unduly overridden in TruncatedPoissonRegression
 */
vector<double>
RegressionModel::getRegressionCoefficients() const {
  return this->coefficients;
}

/**
 * \brief return the number of regression coefficients this regression model
 *        currently has
 */
size_t
RegressionModel::getNumCoefficients() const {
  return this->coefficients.size();
}

/******************************************************************************
 * Specifying the distribution
 */

/**
 * \brief Calculate the deviance function for this regression model given
 *        a set of responses and covariates.
 */
double
RegressionModel::deviance(const vector<double>& responses,
                          const vector< vector<double> >& covariates) const {
  const size_t N = responses.size();
  return this->deviance(responses, vector<double>(N,1), covariates);
}

/**
 * \brief Evaluate the log-likelihood function for this regression model as
 *        the sum over the responses and covariate values given
 * \param response the count responses
 * \param covars   the covariate values associated with the responses given
 */
double
RegressionModel::loglikelihood(const vector<double>& response,
                               const vector< vector<double> >& covars) const {
  const size_t N = response.size();
  return this->loglikelihood(response, std::vector<double>(N,1), covars);
}

/**
 * \brief TODO
 */
double
RegressionModel::loglikelihood(const vector<double>& responses,
                               const vector<double>& probs,
                               const vector< vector<double> >& covars) const {
  double ll = 0;
  for (size_t i=0; i<responses.size(); i++) {
    ll += (this->loglikelihood(responses[i], covars[i]) + log(probs[i]));
  }
  return ll - log(accumulate(probs.begin(), probs.end(), 0.0));
}

/**
 * \brief calculate the p-value for this regression model given the
 *        current set of parameters for the given response and it's
 *        associated covariates
 */
double
RegressionModel::pvalue(const double response,
                        const std::vector<double>& covariates) const {
  double res = 1 - this->cdf(response, covariates);
  if (res < 0) {
    stringstream cvars;
    for (size_t i=0; i<covariates.size(); i++) {
      cvars << covariates[i];
      if (i != covariates.size() - 1) cvars << ",";
    }
    stringstream ss;
    ss << "evaluating p-value for response" << response << " and covariates "
       << cvars.str() << " in regression "
       << "with parameters " << this->toString() << " failed. "
       << "Reason: result less than zero";
    throw RegressionModelException(ss.str());
  }
  return res;
}

/**
 * \brief todo
 */
double
RegressionModel::pdf(const double response,
                     const std::vector<double>& covariates) const {
  return exp(this->loglikelihood(response, covariates));
}

/**
 * \brief calculate the cumulative distribution function for this regression
 *        model given it's current
 */
double
RegressionModel::cdf(const double response,
                     const std::vector<double>& covariates) const {
  double res=0;
  for (size_t i=0; i<=response; i++) res += this->pdf(i, covariates);

  // this might happen because we don't allow -INF from the log-likelihood
  // function and so for values very far down the tail the error can become
  // large, but by the time we hit that, really the value is 1, so..
  if (res > 1) res = 1;

  return res;
}

/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief Given a set of responses, covariates, response probabilities and
 *        an initial estimate, estimate the correct value of the regression
 *        coefficients, beta, for this regression model
 *
 * This class has an implementation for this which uses the NR algorithm, but
 * derived classes can override this behavior
 */
void
RegressionModel::estimateBeta(const vector<double>& response,
                              const vector< vector<double> >& covariates,
                              const vector<double>& probs,
                              const vector<double>& startingPoint,
                              const bool verbose) {
  this->estimateBetaNR(response, covariates, probs, startingPoint, verbose);
}

/**
 * \brief TODO
 */
void
RegressionModel::estimateBetaNR(const vector<double>& y,
                                const vector< vector<double> >& covariates,
                                const vector<double>& probsIn,
                                const vector<double>& startingPoint,
                                const bool verbose) {
  Matrix x(covariates);
  Matrix xt = x.transpose();
  vector<double> b;
  const size_t N = y.size();
  const size_t P = x.numCols();
  double curDeviance = 0;
  double deltaDev = RegressionModel::devianceThreshold + 1;
  double oldDeviance = 0;

  if (verbose) cerr << "estimating beta in GLM by Newton-Raphson" << endl;



  // weight inputs?
  // TODO should be able to remove this..
  vector<double> probs;
  if (probsIn.size() != N) probs = vector<double>(N,1);
  else probs = probsIn;

  // start with a random guess of the regression coefficients, or did the
  // caller give us one?
  // TODO should be able to remove this, assume all params set properly
  if (startingPoint.size() == P) b = startingPoint;
  else {
    for (size_t i=0; i<P; i++) {
      double r = (static_cast<double> (rand() % 100)) / 100;
      b.push_back(r);
    }
    this->coefficients = b;
  }

  if (verbose) {
    cerr << "initial parameter estimates: " << this->toString()
         << " -- log-likelihood: " << this->loglikelihood(y, probs, covariates)
         << endl;
  }

  // find the max covariate values -- we need to know these to set an
  // appropriate upper bound for beta. Also, figure out what is the most
  // each covariate can contribute to the linear predictor
  double maxExpMuContrib =
      RegressionModel::maxLinearPredictorValue / double(P);
  double minExpMuContrib =
      RegressionModel::minLinearPredictorValue / double(P);
  vector<double> maxCovars(P, 0);
  vector<double> minCovars(P, 0);
  for (size_t i=0; i<P; i++) {
    maxCovars[i] = x.maxInCol(i);
    minCovars[i] = x.minInCol(i);
  }

  bool first = true;            // is this the first iteration?
  size_t iterCount = 0;         // counter for number of iterations done
  size_t boundStraight = 0;     // number of consecutive iterations we've
                                // had to adjust beta to avoid over/underflow
  vector<double> prevCoefficients;
  while ((fabs(deltaDev) > RegressionModel::devianceThreshold) &&
         (iterCount < RegressionModel::maxNRIters)) {
    // get the Hessian and invert it
    Matrix inverseHessian = Hessian(y, covariates, probs).inverse();

    // Get the gradient and make a column vector out of it by using a matrix
    // transpose
    vector<double> del = gradient(y, covariates, probs);
    vector< vector<double> > tmp;
    tmp.push_back(del);
    Matrix gradientColumnVec(tmp);
    gradientColumnVec = gradientColumnVec.transpose();

    // multiply the inverted hessian and the gradient to get the update vector
    vector<double> updateVector =\
        (inverseHessian * gradientColumnVec).asVector();

    // do the update
    // TODO - replace with matrix operation? Could remove vector operations
    this->coefficients = vectorSubtraction(this->coefficients, updateVector);

    // tell the user how we're doing?
    if (verbose) {
      cerr << "current parameter estimates: " << this->toString() << " ";
    }

    // bound beta
    for (size_t i=0; i<P; i++) {
      if (this->coefficients[i] * maxCovars[i] > maxExpMuContrib) {
        this->coefficients[i] = maxExpMuContrib / maxCovars[i];
        if (verbose) {
          cerr << "warning: set coefficient " << i << " to "
               << this->coefficients[i];
        }
      }
      else if (this->coefficients[i] * maxCovars[i] < minExpMuContrib) {
        this->coefficients[i] = minExpMuContrib / maxCovars[i];
        if (verbose) {
          cerr << "warning: set coefficient " << i << " to "
               << this->coefficients[i];
        }
      }
    }

    // stopping criteria
    if (first) {
      first = false;
      curDeviance = this->deviance(y, probs, covariates);
      oldDeviance = curDeviance;
    }
    else {
      oldDeviance = curDeviance;
      curDeviance = this->deviance(y, probs, covariates);
      deltaDev = curDeviance - oldDeviance;
    }

    // stop if we're just bouncing between extreme values
    bool boundAllThisIter = true;
    bool boundAllPrevIter = true;
    for (size_t i=0; i<P; i++) {
      if (!((this->coefficients[i] == maxExpMuContrib / maxCovars[i]) ||
          (this->coefficients[i] == minExpMuContrib / maxCovars[i]))) {
        boundAllThisIter = false;
      }
      if (prevCoefficients.size() != this->coefficients.size()) {
        boundAllPrevIter = false;
      }
      else if (!((prevCoefficients[i] == maxExpMuContrib / maxCovars[i]) ||
                 (prevCoefficients[i] == minExpMuContrib / maxCovars[i]))) {
        boundAllPrevIter = false;
      }
    }
    if (boundAllThisIter && boundAllPrevIter) boundStraight += 1;
    if (!boundAllThisIter && boundAllPrevIter) boundStraight = 0;
    if (boundStraight > 3) {
      double currll = this->loglikelihood(y, probs, covariates);
      this->coefficients.swap(prevCoefficients);
      double prevll = this->loglikelihood(y, probs, covariates);
      if (prevll < currll) this->coefficients.swap(prevCoefficients);
      iterCount = RegressionModel::maxNRIters;
    }

    // tell the user how we're doing?
    if (verbose) {
      cerr << " -- delta deviance: " << deltaDev << " old deviance: "
           << oldDeviance << " new deviance: " << curDeviance << endl;
      if (boundStraight > 3) cerr << "reached oscillation max" << endl;
    }

    prevCoefficients = this->coefficients;
    iterCount += 1;
  }
}

/******************************************************************************
 * Regression Model I/O
 */

/**
 * \brief TODO
 */
void
RegressionModel::save(std::string fn) const {
  std::ofstream fs(fn.c_str());
  this->save(fs);
}

/**
 * \brief TODO
 */
void
RegressionModel::load(std::string fn) {
  std::ifstream fs(fn.c_str());
  if (!fs.good()) {
    std::stringstream ss;
    ss << "failed to open " << fn;
    throw RegressionModelException(ss.str());
  }
  this->load(fs);
}
