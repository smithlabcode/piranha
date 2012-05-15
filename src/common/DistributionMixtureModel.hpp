/**
  \file DistributionMixtureModel.hpp
  \brief This header declares the DistributionMixtureModel which allows
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

#ifndef DMIXMOD_HPP
#define DMIXMOD_HPP

#include <string>
#include <vector>
#include <exception>
#include "smithlab_utils.hpp"
#include "Distribution.hpp"
#include "MixtureModel.hpp"

/**
 * \brief TODO
 */
class DistributionMixtureModel : public MixtureModel {
public:
  /*** constructors, destructors and object initialization ***/
  DistributionMixtureModel(std::string type, size_t K, std::string fn); 
  DistributionMixtureModel(std::string type, size_t K);
  ~DistributionMixtureModel();
  void init(std::string type, size_t K);
    
  /*** parameter estimation ***/
  void estimateParams(const std::vector<double>& obs, bool verbose=false);
  void estimateParams(const std::vector<double>& responses,
                      const std::vector<std::vector<double> >& covars, 
                      bool verbose=false);

  /*** Distribution mixture model I/O ***/
  void load(const std::string& fn);
  void save(std::string& fn) const;
  void save(std::ostream& ostrm) const;

  /*** Inspectors ***/
  // simple inspectors
  SimpleXML toXML() const;
  std::string toString() const;
  std::string getComponentAsString(size_t k) const;
  // inspectors for describing the mixture
  double loglikelihood(const std::vector<double> &obs) const;
  // inspectors for querying components
  size_t numberOfComponents() const;
  double getComponentProb(const std::vector<double>& instance,
                          const size_t component) const;
  std::vector<double> getComponentProbs(const std::vector<double>& inst) const;
  std::vector<double> getComponentParameters(const int k) const;
private:
  /*** private instance variables ***/
  std::vector<Distribution*> ds;
  std::vector<double> mixing;

  /*** private constants ***/
  const static double stoppingEM = 0.01;
  const static size_t emMaxIter = 1000;
  const static double mixingLowerBound = 0.0001;

  /*** private member functions ***/
  void estimateParamsEM(const std::vector<double>& obs, bool debug=false);
};

#endif
