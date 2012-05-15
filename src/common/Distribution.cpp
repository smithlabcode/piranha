/**
  \file Distribution.cpp
  \brief This source file defines the abstract Distribution class, which
         serves as a base class for all simple distributions

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
#include <assert.h>
#include "Distribution.hpp"
#include "Poisson.hpp"
#include "NegativeBinomial.hpp"
#include "ZTNB.hpp"
#include "ZTP.hpp"

using std::string;
using std::vector;
using std::stringstream;

/**
 * \brief implementation of the factory pattern for Distribution objects
 * \returns a pointer to a distribution of the type requested initialized to
 *          its default state
 */  
Distribution* 
Distribution::create(const std::string& type) {
  Distribution* res = NULL;
  if (type == "Poisson") res = new Poisson();
  else if (type == "ZeroTruncatedPoisson") res = new ZeroTruncatedPoisson();
  else if (type == "ZeroTruncatedNegativeBinomial") res = new ZTNB();
  else if (type == "NegativeBinomial") res = new NegativeBinomial();
  else {
    stringstream ss; 
    ss << type << " is an unknown distribution type";
    throw DistributionException(ss.str());
  }
  assert(res != NULL);
  return res;
}


/**
 * \brief implementation of the factory pattern for Distribution objects
 * \returns a pointer to a distribution of the type requested initialized to
 *          using the given file.
 * \param fn filename for a file that contains a representation of the
 *           distribution.
 */ 
Distribution* 
Distribution::create(const string& fn, const string& type) {
  Distribution* r = create(type);
  r->load(fn);
  return r;
}
