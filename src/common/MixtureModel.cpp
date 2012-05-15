/****   
 
  Description: TODO

  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Andrew D. Smith
  
  Authors: Philip J. Uren
  
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
****/

#include <string>
#include <vector>
#include <sstream>
#include <assert.h>
#include "MixtureModel.hpp"
#include "RegressionMixtureModel.hpp"
#include "DistributionMixtureModel.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;

/****
 * @summary: implementation of the factory pattern for MixtureModels
 * 
 */  
MixtureModel* 
MixtureModel::create(const std::string& type,  
                     const size_t numComponents, 
                     const size_t numCovariates,
                     const FittingMethod fitMthd) {
  MixtureModel* res = NULL;
  if ((type == "LinearRegression") ||
      (type == "NegativeBinomialRegression") ||
      (type == "PoissonRegression") ||
      (type == "ZeroTruncatedPoissonRegression") ||
      (type == "ZeroTruncatedNegativeBinomialRegression")) {
    res = new RegressionMixtureModel(type, numComponents,
                                     numCovariates, fitMthd);
  }
  else if ((type == "NegativeBinomial") ||
           (type == "Poisson") ||
           (type == "ZeroTruncatedNegativeBinomial") ||
           (type == "ZeroTruncatedPoisson") ||
           (type == "Gaussian")) {
    if (numCovariates != 0) {
      stringstream ss; 
      ss << type << " is a mixture of plain distributions, " << 
          "cannot create a mixture model with covariates"; 
      throw MixtureModelException(ss.str());
    }
    if (!fitMthd.is("DEFAULT")) {
      // TODO this check doesn't belong here.. we should just be passing
      // the fit method to the constructor and letting it make this decision
      stringstream ss;
      ss << type << " cannot be fit using " << fitMthd.toString();
      throw MixtureModelException(ss.str());
    }
    res = new DistributionMixtureModel(type, numComponents);
  } else {
    stringstream ss; 
    ss << type << " is an unknown mixture type"; 
    throw MixtureModelException(ss.str());
  }

  if (res == NULL) {
    // we screwed up if we get to here..
    stringstream ss;
    ss << "failed to build mixture model of type " << type
       << " reason: unknown";
    throw MixtureModelException(ss.str());
  } else return res;
}

/****
 * @summary: implementation of the factory pattern for MixtureModels
 * 
 */ 
MixtureModel* 
MixtureModel::create(const string& fn,
                     const string& type,  
                     const size_t numComponents, 
                     const size_t numCovariates,
                     const FittingMethod fitMthd) {
  MixtureModel* r = create(type, numComponents, numCovariates, fitMthd);
  r->load(fn);
  return r;
}

