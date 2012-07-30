/**
  \file RegressionMixtureModel.cpp
  \brief This source file defines the RegressionMixtureModel class, which is
         an implementation of the MixtureModel abstract interface. It defines
         a kind of mixture where the components are regression models.

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

#include <cmath>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "LinearRegression.hpp"
#include "PoissonRegression.hpp"
#include "RegressionMixtureModel.hpp"
#include "NegativeBinomialRegression.hpp"
#include "RegressionBuilder.hpp"

using std::max;
using std::endl;
using std::cout;
using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;

/******************************************************************************
 * CONSTRUCTORS, DESTRUCTORS AND OBJECT INITIALIZATION
 */

/*****
 * @summary: Regression Mixture constructor
 */
RegressionMixtureModel::RegressionMixtureModel(string type, size_t K, size_t P,
                                               FittingMethod fitMthd,
                                               string fn) {
  this->init(type, K, P, fitMthd);
  this->load(fn);
}

/*****
 * @summary: Regression Mixture constructor
 */
RegressionMixtureModel::RegressionMixtureModel(string type, size_t K, size_t P,
                                               FittingMethod fitMthd) {
  this->init(type, K, P, fitMthd);
}

/*****
 * @summary: object initialization
 */
void RegressionMixtureModel::init(string type, size_t K,
                                  size_t P, FittingMethod fitMthd) {
  for (size_t k = 0; k<K; k++) {
    this->ds.push_back(RegressionBuilder::buildRegression(type,P,fitMthd));
    this->mixing.push_back(1/double(K));
  }
}

/*****
 * @summary: RegressionMixtureModel destructor
 */
RegressionMixtureModel::~RegressionMixtureModel() {
  // components belong to us; we were good, 
  // we didn't leak any pointers or reference, can safely 86 them
  for (size_t k = 0; k < this->ds.size(); k++) delete this->ds[k];
}

/******************************************************************************
 * SIMPLE SETTERS
 */

/****
 * @summary: set the parameters for this regression mixture model.
 * @param params: vector of vector, outer vector must be K+1 in size
 *                <params[0]> are the mixing params and there must be K of them
 *                <params[k]> 1<k<=K are the parameters for the k-1 component
 */
void
RegressionMixtureModel::setParams(const vector< vector<double> >& params) {
  if (params.size() != this->ds.size() + 1) {
    stringstream ss;
    ss << "failed to set parameters for regression mixture. "
       << "Reason: expected " << this->ds.size() << " sets of parameters "
       << "(in addition to mixing parameters), but found only "
       << (params.size() - 1);
    throw MixtureModelException(ss.str());
  }

  this->setMixing(params[0]);
  for (size_t i=1; i<=this->ds.size(); i++) {
    this->ds[i-1]->setParams(params[i]);
  }
}

/****
 * @summary: set the mixing parameters for this regression mixture model
 */
void
RegressionMixtureModel::setMixing(const std::vector<double> mix) {
  if (mix.size() != this->ds.size()) {
    stringstream ss;
    ss << "failed to set mixing parameter for regression mixture. "
       << "Reason: expected " << this->ds.size() << " mixing values, but found "
       << "only " << mix.size();
    throw MixtureModelException(ss.str());
  }
  // TODO throw exception if mixing doesn't sum to approximately 1
  this->mixing = mix;
}

/******************************************************************************
 * INSPECTORS FOR INTEROGATING REGRESSION MIXTURE MODEL DETAILS
 */

/*****
 * @summary: Return a string representation of component <k>
 */
string
RegressionMixtureModel::getComponentAsString(size_t k) const {
  stringstream res;
  res << "component " << k << " mixing: " << this->mixing[k] << " details: ";
  res << this->ds[k]->toString() << endl;
  return res.str();
}

/*****
 * @summary: Return the parameters for component <k> as vector<double>
 */
vector<double>
RegressionMixtureModel::getComponentParameters(const int k) const {
  return this->ds[k]->getParams();
}

/*****
 * @summary: return a string representation of this object.
 * @note:    this is not suitable for saving/loading the object; it's intended
 *           for display on the terminal window
 */
string
RegressionMixtureModel::toString() const {
  stringstream ss;
  const size_t K = this->ds.size();
  for (size_t k=0; k<K; k++) {
    ss << "d[" << k << "] = " << this->ds[k]->toString();
    ss << ", mixing = " << this->mixing[k] << endl;
  }
  return ss.str();
}

/******************************************************************************
 * LOG LIKELIHOOD FUNCTIONS
 */

/****
 * @summary: log-likelihood that a given set of observations was produced
 *           by this mixture
 */
double
RegressionMixtureModel::loglikelihood(const vector<double>& response,
                        const vector<vector<double> >& probs,
                        const vector<vector<double> >& covariates) const {
  const size_t N = response.size();
  const size_t K = this->ds.size();

  vector<double> zs(K,0);
  double sum = 0;
  for (size_t i=0; i<N; i++) {
    for (size_t k=0; k<K; k++) {
      zs[k] = log(mixing[k]) + log(probs[k][i]) +\
              this->ds[k]->loglikelihood(response[i], covariates[i]);
    }
    double m = *std::max_element(zs.begin(), zs.end());
    double d = 0;
    for (size_t k=0; k<K; k++) d += exp(zs[k] - m);
    sum += (m + log(d));
  }
  return sum;
}

/****
 * @summary: TODO
 */
double
RegressionMixtureModel::loglikelihood(const vector<double>& response,
                        const vector<vector<double> >& covariates) const {
  const size_t N = response.size();
  const size_t K = this->ds.size();
  vector< vector<double> > probs = \
      vector< vector<double> >(K,vector<double>(N,1));
  return this->loglikelihood(response,probs,covariates);
}

/******************************************************************************
 * INSPECTORS FOR QUERYING COMPONENT PROBABILITIES
 */

/*****
 * @summary: get the probability that <instance> was generated from
 *           component number <component>
 * @note:    Honor base class contract to implement this, 0th element of
 *           instance is assumed to be response, rest are covariates
 */
double
RegressionMixtureModel::getComponentProb(const vector<double>& instance,
                                         const size_t component) const {
  if (instance.size() < 2) {
    stringstream ss;
    ss << "failed to get component probabilities for regression mixture. "
       << "Reason: instance should be defined by response and at least one "
       << "covariate; found no covariates";
    throw MixtureModelException(ss.str());
  }

  double tmpRes = instance[0];
  vector<double> covars;
  covars.resize(instance.size()-1);
  copy(instance.begin()+1, instance.end(), covars.begin());
  return this->getComponentProb(tmpRes, covars, component);
}

/*****
 * @summary: get the k-dimensional vector of probabilities that <inst> was
 *           generated from each of the k components of this mixture.
 * @note:    Honor base class contract to implement this, 0th element of
 *           instance is assumed to be response, rest are covariates
 */
vector<double>
RegressionMixtureModel::getComponentProbs(const vector<double>& inst) const {
  vector<double> res;
  for (size_t k=0; k<this->ds.size(); k++)
    res.push_back(this->getComponentProb(inst,k));
  return res;
}

/****
 * @summary: given a collection of responses and associated covariates,
 *           determine the probability that each response was produced
 *           by component <component> of this mixture.
 */
vector<double>
RegressionMixtureModel::getComponentProbs(const vector<double>& responses,
                                        const vector< vector<double> >& covars,
                                        const size_t component) const {
  // make sure the caller gave us sensible input
  if (responses.size() != covars.size()) {
    stringstream ss;
    ss << "failed to get component probabilities for regression mixture. "
       << "Reason: number of coefficients does not match number of covariates";
    throw MixtureModelException(ss.str());
  }
  if (covars.size() <= 0) {
    stringstream ss;
    ss << "failed to get component probabilities for regression mixture. "
       << "Reason: empty covariate set";
    throw MixtureModelException(ss.str());
  }
  std::vector<double> r;
  for (size_t i=0; i<responses.size(); i++)
    r.push_back(this->getComponentProb(responses[i], covars[i], component));
  return r;
}

/****
 * @summary: given a collection of responses and associated covariates,
 *           determine the probability that each response was produced
 *           by each component of this mixture.
 * @return:  n * p matrix (vec of vec) where entry [i][j] is the prob that
 *           response[i] with covariates[i] was produced by component [j] of
 *           this mixture.
 */
vector< vector<double> >
RegressionMixtureModel::getComponentProbs (const vector<double>& response,
                             const vector<vector<double> >& covariates) const {
  // make sure the caller gave us sensible input
  if (response.size() != covariates.size()) {
    stringstream ss;
    ss << "failed to get component probabilities for regression mixture. "
       << "Reason: number of coefficients does not match number of covariates";
    throw MixtureModelException(ss.str());
  }
  if (covariates.size() <= 0) {
    stringstream ss;
    ss << "failed to get component probabilities for regression mixture. "
       << "Reason: empty covariate set";
    throw MixtureModelException(ss.str());
  }

  const size_t N = response.size();
  const size_t K = this->ds.size();

  // compute probs given current params
  // probs[k][i] = prob of ith response in kth component
  vector< vector<double> > probs(K, vector<double>(N,0));
  vector<double> zs(K,0);
  for (size_t i=0; i<N; i++) {
    for (size_t k=0; k<K; k++) {
      zs[k] = log(mixing[k]) +
              this->ds[k]->loglikelihood(response[i], covariates[i]);
    }
    double m = *std::max_element(zs.begin(), zs.end());
    double sum = 0;
    for (size_t k=0; k<K; k++) sum += exp(zs[k] - m);
    for (size_t k=0; k<K; k++) {
      probs[k][i] = exp(zs[k] - m) / sum;
    }
  }
  return probs;
}

/****
 * @summary: what is the probability that the <component>'th component
 *           produced the response <response> given the covariates <covars>?
 */
double
RegressionMixtureModel::getComponentProb(const double response,
                                         const std::vector<double> covars,
                                         const size_t component) const {
  const size_t K = this->ds.size();

  // sanity check, can't get more components than we have!
  if (component >= K) {
    stringstream ss;
    ss << "failed retrieving component " << component << " probability "
       << "from mixture. Reason: mixture only has " << K << " components!";
    throw MixtureModelException(ss.str());
  }

  // compute probabilities given current parameters
  vector<double> zs(K,0);
  vector<double> probs(K,0);

  for (size_t k=0; k<K; k++)
    zs[k] = log(mixing[k]) + this->ds[k]->loglikelihood(response, covars);
  double m = *std::max_element(zs.begin(), zs.end());
  double sum = 0;
  for (size_t k=0; k<K; k++) sum += exp(zs[k] - m);
  for (size_t k=0; k<K; k++) {
    probs[k] = exp(zs[k] - m) / sum;
  }

  return probs[component];
}

/******************************************************************************
 * PARAMETER ESTIMATION
 */

/****
 * @summary: Estimate the parameters for this RegressionMixtureModel using
 *           the EM algorithm
 */
void
RegressionMixtureModel::estimateParamsEM(const vector<double>& response, 
                            const vector<vector<double> >& covariates,
                            const bool verbose) {

  const size_t N = response.size();
  const size_t P = covariates[0].size();
  const size_t K = this->ds.size();

  // it's easy if there's only 1 component!
  if (K == 1) {
    this->ds[0]->estimateParams(response, covariates, verbose);
    this->mixing[0] = 1.0;
    return;
  }

  // we need to randomly initialize the components
  for (size_t k=0; k<K; k++) this->ds[k]->randomise(P);

  // do the EM to fit the mixture
  if (verbose) cout << "Estimate parameters in regression mixture "
                    << "model by EM" << endl;
  double oll = 0;
  bool first = true, done = false;
  while(!done) {
    // expectation step -- compute probabilities given current parameters
    vector<double> zs(K,0), probSums(K,0);
    // probs[k][i] = probability of ith response in kth component
    vector< vector<double> > probs(K, vector<double>(N,0));
    for (size_t i=0; i<N; i++) {
      for (size_t k=0; k<K; k++) {
        zs[k] = log(mixing[k]) +\
                this->ds[k]->loglikelihood(response[i], covariates[i]);
      }
      double m = *std::max_element(zs.begin(), zs.end());
      double sum = 0;
      for (size_t k=0; k<K; k++) sum += exp(zs[k] - m);
      for (size_t k=0; k<K; k++) {
        probs[k][i] = exp(zs[k] - m) / sum;
        probSums[k] += probs[k][i];
      }
    }
    
    // maximisation step -- update model parameters
    for (size_t k=0; k<K; k++) {
      this->mixing[k] = (probSums[k] / (double) response.size());
      cout << "estimate params for component " << k << endl;
      cout << this->ds[k]->toString() << endl;
      this->ds[k]->estimateParams(response, covariates, probs[k],
                                  this->ds[k]->getRegressionCoefficients(),
                                  verbose);
    }
    
    if (verbose) {
      cout << "model " << endl;
      cout << this->toString() << endl;
      cout << "------" << endl;
    }
    
    // finished?
    double nll = this->loglikelihood(response, covariates);
    if (verbose)
      cout << "log likelihood: " << nll << " -- " << oll
           << " -- " << fabs(nll - oll) << " -- "
           << STOPPING_CRITERIA << endl;
    if ((!first) && (fabs(nll - oll) < STOPPING_CRITERIA)) done = true;
    first = false;
    oll = nll;
  }
}

/******************************************************************************
 * REGRESSION MIXTURE MODEL I/O
 */

/*****
 * @summary: Load this regression mixture model from a string representation.
 *           We expect the string to be a concrete XML description.
 */
void 
RegressionMixtureModel::load(const std::string& fn) {
  ifstream istrm(fn.c_str());
  if (!istrm.good()) {
    stringstream ss;
    ss << "failed to open " << fn;
    throw MixtureModelException(ss.str()); 
  }
  SimpleXML xml(istrm);
  this->load(xml);
}

/*****
 * @summary: Load this regression mixture model from an abstract XML object
 */
void RegressionMixtureModel::load(SimpleXML& xml) {
  vector<SimpleXML> comps = xml.getChildren("component");
  for (size_t k=0; k<comps.size(); k++) {
    vector<SimpleXML> cs = comps[k].getChildren();
    if (cs.size() != 2) {
      stringstream ss;
      ss << "failed parsing " << endl << xml.toString() << " as regression "
         << "mixture model. Reason: component tag should have two sub-tags:"
         << "one for <mixing> and one for component regression. Found "
         << cs.size() << " instead";
      throw MixtureModelException(ss.str());
    }

    for (size_t j=0; j<2; j++) {
      if (cs[j].getTagName() == "mixing") {
        // this is the mixing parameter
        std::stringstream ss;
        ss << cs[j].getData();
        ss >> this->mixing[k]; 
      } else {
        // must be the regression XML object
        this->ds[k]->load(cs[j]);
      }
    }
  }
}

/*****
 * @summary: Save the RegressionMixtureModel to an output stream
 * @note:    Concrete syntax is XML
 */
void
RegressionMixtureModel::save(std::ostream& ostrm) const {
  size_t K = this->numberOfComponents();
  ostrm << "<RegressionMixtureModel>" << endl;
  for (size_t k=0; k<K; k++) {
    ostrm << "<component> " << endl;
    this->ds[k]->save(ostrm);
    ostrm << "<mixing>" << this->mixing[k] << "</mixing>"; 
    if (k != K-1) ostrm << endl;
    ostrm << "</component> " << endl;
  }
  ostrm << "</RegressionMixtureModel>" << endl;
}

/****
 * @summary: Save the RegressionMixtureModel to a file with the given filename
 * @note:    Concrete syntax is XML
 */
void
RegressionMixtureModel::save(std::string& fn) const {
  std::ofstream fs(fn.c_str());
  save(fs);
}


