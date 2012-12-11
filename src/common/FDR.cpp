/**
  \file FDR.hpp
  \brief This source file defines the functions for correcting p-values for
         multiple hypothesis testing

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

#include "FDR.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

using std::vector;

struct CompP : std::binary_function<double, double, bool>
{
    CompP(const std::vector<double>& pvals) : m_pvals(pvals) {}
    bool operator()(size_t Lhs, size_t Rhs) const {
        return m_pvals[Lhs] < m_pvals[Rhs];
    }
    const std::vector<double>& m_pvals;
};

void
FDR::correctP(vector<double> &pvals) {
  // make vector of indices
  std::vector<size_t> pos(pvals.size());
  for (size_t i = 0; i != pos.size(); ++i) pos[i] = i;

  // sort indices
  std::sort(pos.begin(), pos.end(), CompP(pvals));

  const double N = pvals.size();
  for (int i=N-1; i>=0; i--) {
    const double tmp = (N / (i + 1)) * pvals[pos[i]];
    if (i + 1 == N) pvals[pos[i]] = tmp;
    else pvals[pos[i]] = std::min(tmp, pvals[pos[i+1]]);
  }
}
