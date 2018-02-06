/**
  \file Utility.cpp
  \brief TODO

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
****/

// stl includes
#include <algorithm>
#include <vector>
#include <sstream>

// piranha includes
#include "Utility.hpp"

using std::vector;
using std::stringstream;

/**
 * \brief TODO
 */
void
mergeResponsesPvals(vector<GenomicRegion> &sites_fg,
                    vector<GenomicRegion> &sites_bg,
                    vector<double> &pvals_fg,
                    vector<double> &pvals_bg) {
  vector< vector<double> > dummy;
  mergeResponsesCovariatesPvals(sites_fg, sites_bg, pvals_fg,
                                pvals_bg, dummy, dummy);
}

/**
 * \brief in-place merge of sorted response, covariates and p-value vectors
 *        results are placed into the bg vectors, fg ones are left empty
 * \param TODO
 */
void
mergeResponsesCovariatesPvals(vector<GenomicRegion> &sites_fg,
                              vector<GenomicRegion> &sites_bg,
                              vector<double> &pvals_fg,
                              vector<double> &pvals_bg,
                              vector< vector<double> > &covars_fg,
                              vector< vector<double> > &covars_bg) {
  // step 0 -- sanity check
  if (covars_fg.size() != covars_bg.size()) {
    stringstream ss;
    ss << "failed merging covariates -- number of covariates in "
       << "foreground and background didn't match. Foreground had "
       << covars_fg.size() << " while background had " << covars_bg.size();
    throw SMITHLABException(ss.str());
  }

  // step 1 -- concatenate vectors and free the fg vectors, but remember
  //           where partition point is
  size_t lastSite_bg = sites_bg.size();
  size_t lastPval_bg = pvals_bg.size();
  vector<size_t> lastCvar_bg;
  copy(sites_fg.begin(), sites_fg.end(), back_inserter(sites_bg));
  copy(pvals_fg.begin(), pvals_fg.end(), back_inserter(pvals_bg));
  for (size_t i = 0; i < covars_fg.size(); i++) {
    lastCvar_bg.push_back(covars_bg[i].size());
    copy(covars_fg[i].begin(), covars_fg[i].end(), back_inserter(covars_bg[i]));
    covars_fg[i].clear();
  }
  sites_fg.clear();
  pvals_fg.clear();

  // step 2 -- merge pvals based on sort order defined by sites
  typedef vector<double>::iterator dIt;
  vector< vector<double>::iterator > pvals_idx;
  for (dIt it = pvals_bg.begin(); it != pvals_bg.end(); ++it)
    pvals_idx.push_back(it);
  inplace_merge(pvals_idx.begin(), pvals_idx.begin() + lastPval_bg,
                pvals_idx.end(), SiteComp(sites_bg, pvals_bg));
  vector<double> pvals_res;
  for (size_t i = 0; i < pvals_idx.size(); i++)
    pvals_res.push_back(*(pvals_idx[i]));
  pvals_bg.swap(pvals_res);
  pvals_res.clear();

  // step 3 -- merge covariates based on sort order defined by sites
  for (size_t i = 0; i < covars_fg.size(); i++) {
    vector< vector<double>::iterator > cvar_idx;
    for (dIt it = covars_bg[i].begin(); it != covars_bg[i].end(); ++it)
      cvar_idx.push_back(it);
      inplace_merge(cvar_idx.begin(), cvar_idx.begin() + lastCvar_bg[i],
                    cvar_idx.end(), SiteComp(sites_bg, covars_bg[i]));
      vector<double> cvars_res;
      for (size_t j = 0; j < cvar_idx.size(); j++)
        cvars_res.push_back(*(cvar_idx[j]));
      covars_bg[i].swap(cvars_res);
      cvars_res.clear();
  }

  // step 4 -- merge the sites. tah-da, we're done.
  inplace_merge(sites_bg.begin(), sites_bg.begin() + lastSite_bg,
                  sites_bg.end());
}
