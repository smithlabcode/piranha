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
          
  TODO NR fitting code needs to move up to RgeressionModel, since this is
       used in ZTP and ZTNB which are not GLMs. That means coefficients will
       also need to be kicked up the hierarchy
****/

#include <sstream>
#include <vector>
#include "math.h"

#include "GeneralisedLinearModel.hpp"
#include "PoissonRegression.hpp"
#include "VectorMath.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::stringstream;





/****
 * @summary: TODO
 */
void 
GeneralisedLinearModel::estimateBeta(const std::vector<double>& response,
                                     const vector< vector<double> >& covariates,
                                     const vector<double>& probs,
                                     const vector<double>& startingPoint,
                                     const bool verbose) {
  if (this->fittingMethod.is("IRLS") || this->fittingMethod.is("DEFAULT")) {
    this->estimateBetaIRLS(response, covariates, probs, startingPoint,
                           verbose);
  } else if (this->fittingMethod.is("NR")) {
    this->estimateBetaNR(response, covariates, probs, startingPoint,
                         verbose);
  } else {
    stringstream ss;
    ss << "GLM fitting method not recognized: "
       << this->fittingMethod.toString();
    throw RegressionModelException(ss.str());
  }
}



/****
 * @summary: TODO
 */
void
GeneralisedLinearModel::estimateBetaIRLS(const vector<double>& y, 
                                    const vector< vector<double> >& covariates,
                                    const vector<double>& probsIn,
                                    const vector<double>& startingPoint,
                                    const bool verbose) {
  Matrix x(covariates);
  Matrix xt = x.transpose();
  size_t n = y.size();
  size_t p = x.numCols();
  double curDeviance = 0; 
  double deltaDev = RegressionModel::devianceThreshold + 1;
  double oldDeviance = 0;
  
  if (verbose) cerr << "estimating beta in GLM by IRLS" << endl;

  // weight inputs?
  vector<double> probs;
  if (probsIn.size() != n) probs = vector<double>(n,1);
  else probs = probsIn;
  
  // start with a random guess of the regression coefficients, or did the
  // caller give us one?
  if (startingPoint.size() == p) this->coefficients = startingPoint;
  else {
    vector<double> b;
    for (size_t i=0; i<p; i++) {
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
      
  size_t iterCount = 0;
  while ((fabs(deltaDev) > RegressionModel::devianceThreshold) &&
         (iterCount < GeneralisedLinearModel::maxIRLSIters)) {
    // weights
    DiagonalMatrix w(n,n);
    for (size_t i=0; i<n; i++) {
      // e.g. for Poisson: w(i,i) = mu(i, covariates, b) * probs[i];
      double v = this->variance(i, covariates);
      double linkDeriv = this->linkFirstDerivative(i, covariates);
      double linkDerivSqr = linkDeriv * linkDeriv;
      w(i,i) = probs[i] / (v * linkDerivSqr);
    }
    
    // working response
    vector<double> zt;
    for (size_t i=0; i<n; i++) {
      double eta = vectorDotProduct(x(i), this->coefficients);
      double mui = exp(eta);
      zt.push_back(eta + ((y[i] - mui) * this->linkFirstDerivative(i, covariates)));
    }
    Matrix z = Matrix(zt).transpose();
    
    // update
    this->coefficients = ((((xt * w) * x).inverse()) * ((xt * w) * z)).asVector();
    
    // stopping criteria
    oldDeviance = curDeviance;
    curDeviance = this->deviance(y, probs, covariates);
    deltaDev = log(curDeviance) - log(oldDeviance);
    if (verbose)
      cerr << "current parameter estimates: " << this->toString()
           << " -- delta deviance: " << deltaDev << " old dev: " 
           << oldDeviance << " new dev: " << curDeviance
           << " -- log-likelihood: " << this->loglikelihood(y,covariates) << endl;
    
    iterCount += 1;
  } 
}
