/****
  @summary: This class is an implementation of the MixtureModel abstract
            interface. It defines a kind of mixture where the components
            are regression models

  --------------------

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

  TODO:          None

****/

#ifndef REGMIX_HPP
#define REGMIX_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "SimpleXML.hpp"
#include "MixtureModel.hpp"
#include "RegressionModel.hpp"

class RegressionMixtureModel : public MixtureModel {
public :
  /* constructors, destructors and object initialization */
  RegressionMixtureModel(std::string type, size_t K, size_t P,
                         FittingMethod fitMthd, std::string fn);
  RegressionMixtureModel(std::string type, size_t K, size_t P,
                         FittingMethod fitMthd);
  ~RegressionMixtureModel();
  
  /* simple mutators */
  void setParams(const std::vector< std::vector<double> >& params);
  void setMixing(const std::vector<double> mix);

  /*parameter estimation */
  void estimateParams(const std::vector<double>& response, 
                      const std::vector<std::vector<double> >& covariates,
                      const bool verbose=false) {
    estimateParamsEM(response, covariates, verbose);
  }
  
  /* inspectors -- interrogating component probabilities */
  std::vector< std::vector<double> > getComponentProbs (
                    const std::vector<double>& response,
                    const std::vector<std::vector<double> >& covariates) const;
  std::vector<double> getComponentProbs(const std::vector<double>& responses,
                    const std::vector< std::vector<double> >& covars,
                    const size_t component) const;
  double getComponentProb(const double response, 
                    const std::vector<double> covars,
                    const size_t component) const;
  double getComponentProb(const std::vector<double>& instance,
                    const size_t component) const;
  std::vector<double> getComponentProbs(const std::vector<double>& i) const;
  
  /* inspectors */
  std::vector<double> getMixing() const { return this->mixing; } 
  std::string toString() const;
  std::string getComponentAsString(size_t k) const;
  size_t numberOfComponents() const { return this->ds.size(); }
  std::vector<double> getComponentParameters(const int k) const;
  double loglikelihood(const std::vector<double>& response,
      const std::vector<std::vector<double> >& covariates) const;
  double loglikelihood(const std::vector<double>& response,
      const std::vector< std::vector<double> >& probs,
      const std::vector< std::vector<double> >& covariates) const;
  std::vector<std::vector<double> > getParams() const {
    std::vector<std::vector<double> > res;
    for (size_t i=0; i<this->ds.size(); i++) 
      res.push_back(this->ds[i]->getRegressionCoefficients());
    return res;
  }
  
  /* RegressionMixtureModel I/O */
  void load(const std::string& fn);
  void load(SimpleXML& xml);
  void save(std::string& fn) const;
  void save(std::ostream& fh) const;
  
private:
  /* constructors, destructors and object initialization */
  void init(std::string type, size_t K, size_t P, FittingMethod fitMthd);
    
  /* parameter estimation */
  void estimateParamsEM(const std::vector<double>& response, 
                        const std::vector<std::vector<double> >& covariates,
                        const bool verbose=false);
  
  /* private instance variables */
  std::vector<RegressionModel*> ds;
  std::vector<double> mixing; 
  
  /* private constants */
  const static double STOPPING_CRITERIA = 0.01;
};

#endif
