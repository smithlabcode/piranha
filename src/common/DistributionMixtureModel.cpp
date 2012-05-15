/**
  \file DistributionMixtureModel.hpp
  \brief This source file defines the DistributionMixtureModel which allows
         modeling and fitting of mixtures of simple distributions

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

#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <numeric>
#include "Gaussian.hpp"
#include "Poisson.hpp"
#include "ZTP.hpp"
#include "ZTNB.hpp"
#include "NegativeBinomial.hpp"
#include "DistributionMixtureModel.hpp"

#include "gsl/gsl_randist.h"
 
using std::stringstream;
using std::max_element;
using std::ofstream;
using std::ifstream;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::set;


/******************************************************************************
 * Constructors, destructors and object initialization
 */

/**
 * \brief TODO
 */
DistributionMixtureModel::DistributionMixtureModel(string type, 
    size_t K, string fn) {
  this->init(type, K);
  this->load(fn);
}

/**
 * \brief TODO
 */
DistributionMixtureModel::DistributionMixtureModel(string type, size_t K) {
  this->init(type, K);
}

/**
 * \brief TODO
 * \todo Move construction to separate static factory method
 */
void 
DistributionMixtureModel::init(string type, size_t K) {
  for (size_t k = 0; k<K; k++) {
    if (type == "Gaussian") this->ds.push_back(new Gaussian());
    else if (type == "Poisson") this->ds.push_back(new Poisson());
    else if (type == "ZeroTruncatedPoisson") this->ds.push_back(new Poisson());
    else if (type == "NegativeBinomial") this->ds.push_back(new NegativeBinomial());
    else if (type == "ZeroTruncatedNegativeBinomial") this->ds.push_back(new ZTNB());
    else {
      stringstream msg;
      msg << "DistributionMixtureModel - unknown distribution -> '"
          << type << "'";
      throw MixtureModelException(msg.str()); 
    }
    
    this->mixing.push_back(1/double(K));
  }
}
   
/**
 * \brief DistirbutionMixtureModel destructor
 */
DistributionMixtureModel::~DistributionMixtureModel() {
  // components belong to us; we were good, we didn't leak any pointers or
  // references, we can safely 86 them.
  for (size_t k = 0; k < this->ds.size(); k++) delete this->ds[k];
}


/******************************************************************************
 * Describing the mixture
 */

/**
 * \brief log-likelihood that a given vector of observations was produced
 *        by this mixture
 * \note  we scale values in log-space to avoid an underflow
 */
double
DistributionMixtureModel::loglikelihood(const std::vector<double> &obs) const {
  size_t N = obs.size(), K = this->numberOfComponents();
  double sum = 0;
  for (size_t i=0; i<N; i++) {
    // calculate z's for this i
    vector<double> z;
    for (size_t k=0; k<K; k++) {
      double zi = log(this->mixing[k]) + this->ds[k]->loglikelihood(obs[i]);
      z.push_back(zi);
    }

    // find max z and get sum over k
    double m = *max_element(z.begin(), z.end());
    double d = 0;
    for (size_t k=0; k<K; k++) d += exp(z[k] - m);

    // add contribution from ith observation
    sum += (m + log(d));
  }

  return sum;
}
    

/******************************************************************************
 * Simple inspectors
 */

/**
 * \brief TODO
 */
string
DistributionMixtureModel::toString() const {
  stringstream res;
  for (size_t k=0; k<this->numberOfComponents(); k++) {
    res << "component " << k << " mixing: " << mixing[k] << " details: "
        << this->ds[k]->toString() << endl;
  }
  return res.str();
}

/**
 * \brief TODO
 */
SimpleXML
DistributionMixtureModel::toXML() const {
  stringstream ss;
  size_t K = this->numberOfComponents();
  ss << "<DistributionMixtureModel>" << endl;
  for (size_t k=0; k<K; k++) {
    ss << "<component> ";
    this->ds[k]->save(ss);
    ss << "<mixing>" << this->mixing[k] << "</mixing>";
    ss << "</component> ";
  }
  ss << "</DistributionMixtureModel>";
  return SimpleXML(ss);
}

/**
 * \brief TODO
 */
string
DistributionMixtureModel::getComponentAsString(size_t k) const {
  stringstream res;
  res << "component " << k << " mixing: " << this->mixing[k] << " details: ";
  res << this->ds[k]->toString() << endl;
  return res.str();
}

/**
 * \brief TODO
 */
size_t
DistributionMixtureModel::numberOfComponents() const {
  return this->ds.size();
}


/******************************************************************************
 * Inspectors for querying component probabilities
 */

/**
 * \brief TODO
 */
double
DistributionMixtureModel::getComponentProb(const std::vector<double>& instance,
        const size_t component) const {
  size_t K = this->numberOfComponents();
  if (instance.size() != 1) {
    stringstream msg;
    msg << "DistributionMixtureModel - instance must be "
        << "described by single value, but found " << instance.size();
    throw MixtureModelException(msg.str());
  }
  double val = instance[0];

  // calculate z's
  vector<double> z;
  for (size_t k=0; k<K; k++) {
    z.push_back(log(this->mixing[k]) + this->ds[k]->loglikelihood(val));
  }

  // find max z and get sum over k
  double m = *max_element(z.begin(), z.end());
  double d = 0;
  for (size_t k=0; k<K; k++) d += exp(z[k] - m);

  return exp(z[component] - m) / d;
}

/**
 * \brief TODO
 */
vector<double>
DistributionMixtureModel::getComponentProbs(const vector<double>& instance) const {
  vector<double> res;
  for (size_t k=0; k<this->numberOfComponents(); k++) {
    res.push_back(getComponentProb(instance, k));
  }
  return res;
}

/**
 * \brief TODO
 */
vector<double>
DistributionMixtureModel::getComponentParameters(const int k) const {
  return this->ds[k]->getParams();
}


/******************************************************************************
 * Parameter estimation
 */

/**
 * \brief Fit this DistributionMixtureModel using the given responses
 */
void 
DistributionMixtureModel::estimateParams(const vector<double>& obs, 
                                         bool verbose) {
  this->estimateParamsEM(obs, verbose);
}

/**
 * \brief Fit this DistributionMixtureModel using the given responses and
 *        covariates.
 * \note If covariates are non-empty, this fails because it's not possible
 *       to use them for fitting simple distirbutions. This is just here to
 *       provide uniformity to the interface between DistributionMixtureModel
 *       and RegressionMixtureModel classes.
 */
void 
DistributionMixtureModel::estimateParams(const vector<double>& responses,
                                         const vector< vector<double> >& covs, 
                                         bool verbose) {
  if ((covs.size() > 0) && (covs[0].size() > 0)) {
    stringstream ss;
    ss << "Distribution mixture models can't be fit using covariates, "
       << "but a non-empty set of covariates was provided to "
       << "DistributionMixture::estimateParams()";
    throw MixtureModelException(ss.str());
  }
  this->estimateParams(responses, verbose);
}


/**
 * \brief Given a set of observations, fit this mixture using
 *        expectation maximization
 */
void DistributionMixtureModel::estimateParamsEM(const vector<double>& obs,
                                                bool verbose) {
  bool debug = false;
  double oldLogLike = 0;
  bool first = true;
  bool done = false;
  size_t numIter = 0;
  
  size_t K = this->numberOfComponents();
  
  // only one component? easy peasy..
  if (K == 1) {
    this->ds[0]->estimateParams(obs, verbose);
    this->mixing[0] = 1;
    return;
  }

  // fail if there are more components in this mixture than there are 
  // distinct observation values -- we can't fit that!
  set<double> unique(obs.begin(), obs.end());
  if (obs.size() < this->ds.size()) {
    stringstream s;
    s << "cannot fit mixture model with " << this->ds.size()
        << " using only " << obs.size() << " data points!";
    throw MixtureModelException(s.str().c_str());
  }

  // make an initial random guess 
  for (size_t k=0; k<K; k++) this->ds[k]->randomise();

  if (verbose)
    cout << endl << "initial mixture starting point: " << this->toString()
         << endl;

  while (!done) {
    // expectation step -- compute weights given current params
    vector<vector<double> > w;
    for (size_t i = 0; i < obs.size(); i++) {
      // calculate z's for this i
      vector<double> z;
      for (size_t k = 0; k<K; k++) {
        z.push_back(log(this->mixing[k]) + this->ds[k]->loglikelihood(obs[i]));
      }

      // find max z and denominator for w
      double m = *max_element(z.begin(), z.end());
      double d = 0;
      for (size_t k=0; k<K; k++) {
        d += exp(z[k] - m);
      }
      
      // calc weights for each component
      vector<double> wi;
      for (size_t k = 0; k<K; k++) {
        wi.push_back(exp(z[k] - m) / d);
      }
      w.push_back(wi);
    }
    if (debug) {
      for (size_t k = 0; k<K; k++) {
        for (size_t i = 0; i < obs.size(); i++) {
          cout << w[i][k] << ", "
              << this->ds[k]->loglikelihood(obs[i]) << " (" << obs[i]
              << "," << this->ds[k]->getParams()[0] << ","
              << this->ds[k]->getParams()[1] << ") " << ", ";
        }
        cout << endl;
      }
    }

    // maximisation step -- update model parameters
    mixing.clear();
    for (size_t k=0; k<K; k++) {
      // update kth component 
      vector<double> wk;
      for (size_t i = 0; i < obs.size(); i++)
        wk.push_back(w[i][k]);
      this->ds[k]->estimateParams(obs, wk);
      // update mixing parameter for kth component
      double sum = 0;
      for (size_t i = 0; i < obs.size(); i++)
        sum += w[i][k];
      mixing.push_back(sum / obs.size());
    }
    if (debug)
      cout << "after maximisation step, model is: " << this->toString() << endl;

    // set a lower bound on mixing params
    for (size_t k=0; k<K; k++) {
      if (this->mixing[k] < this->mixingLowerBound)
        this->mixing[k] = this->mixingLowerBound;
    }
    double total = 0;
    for (size_t i = 0; i < this->mixing.size(); i++)
      total += this->mixing[i];
    for (size_t k = 0; k<K; k++)
      this->mixing[k] = this->mixing[k] / total;

    // finished because likelihood is stable?
    double newLogLike = this->loglikelihood(obs);
    if (verbose) {
      cout << this->toString() << endl;
      cout << "new loglike: " << newLogLike << " old loglike: " << oldLogLike
           << endl;
    }
    if ((!first) && (fabs(newLogLike - oldLogLike) < this->stoppingEM))
      done = true;
    first = false;
    oldLogLike = newLogLike;

    // finished because we're giving up?
    numIter++;
    if (numIter >= this->emMaxIter)
      done = true;
  }
}


/******************************************************************************
 * Distribution mixture model I/O
 */

/**
 * \brief load this distributionMixtureModel from the indicated file
 * \bug this is currently not implemented!
 */
void
DistributionMixtureModel::load(const std::string& fn) {
  /*ifstream istrm(fn.c_str());
  size_t K = this->numberOfComponents();
  for (size_t k=0; k<K; k++) {
    this->ds[k]->load(istrm);
    istrm >> this->mixing[k];
  }*/
  // TODO fix this
  throw DistributionException("DistributionMixtureModel::load not implm.");
}

/**
 * \brief open a stream to the indicated file and write this
 *        DistributionMixtureModel to the stream.
 */
void
DistributionMixtureModel::save(std::string& fn) const {
  ofstream fs(fn.c_str());
  this->save(fs);
}

/**
 * \brief write this DistributionMixtureModel to the given ostream; format
 *        is XML
 */
void
DistributionMixtureModel::save(std::ostream& ostrm) const {
  ostrm << this->toXML().toString();
}
