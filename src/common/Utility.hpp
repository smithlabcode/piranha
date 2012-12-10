/**
  \file Utility.hpp
  \brief This header file declares utility functions and templates for use in
         the Piranha package.

  \authors Philip J. Uren, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2012
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

#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

#include "smithlab_utils.hpp"
#include "SimpleXML.hpp"

/**
 * \brief An exception class to be used for utility-based exceptions
 */
class UtilException : public SMITHLABException {
public:
  UtilException(std::string s = std::string()) : SMITHLABException(s) {}
};

/**
 * \brief TODO
 */
struct Contained {
  Contained(const std::tr1::unordered_set<size_t> &x) : theSet(x) {}
  int operator()(size_t y) { return theSet.count(y); }
private:
  const std::tr1::unordered_set<size_t> &theSet;
};

/**
 * \brief TODO
 */
template <typename T>
void removeByIndex(std::vector<T>& vals,
                   const std::tr1::unordered_set<size_t> &indices) {
  Contained contained(indices);
  int last = 0;
  for (size_t i=0; i<vals.size(); ++i, ++last) {
    while (contained(i)) ++i;
    if (i >= vals.size()) break;
    vals[last] = vals[i];
  }
  vals.resize(last);
}

/**
 * \brief TODO
 */
template <typename T>
void copyByIndex(const std::vector<T>& vals,
                 std::vector<T>& newVals,
                 const std::tr1::unordered_set<size_t> &indices) {
  Contained contained(indices);
  for (size_t i=0; i<vals.size(); ++i)
    if (contained(i)) newVals.push_back(vals[i]);
}

/**
 * \brief TODO
 */
template <typename T>
void copyRemoveByIndex(std::vector<T>& vals,
                 std::vector<T>& newVals,
                 const std::tr1::unordered_set<size_t> &indices) {
  copyByIndex(vals, newVals, indices);
  removeByIndex(vals, indices);
}


/**
 * \brief TODO
 */
template <typename T>
std::map<T, size_t> countItemsOrdered (std::vector<T> &v) {
  std::map<T, size_t> res;
  for (size_t i = 0; i< v.size(); i++) {
    if (res.count(v[i]) == 0) res[v[i]] = 0;
    res[v[i]] += 1;
  }
  return res;
}

#endif
