/**
  \file RegressionBuilder.hpp
  \brief This header file declares and defines the RegressionBuilder class,
         which is a static factory for constructing regression models

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


#ifndef RBL_HPP
#define RBL_HPP

#include <string>
#include <sstream>
#include "RegressionModel.hpp"
#include "PoissonRegression.hpp"
#include "NegativeBinomialRegression.hpp"
#include "TruncatedPoissonRegression.hpp"
#include "LinearRegression.hpp"
#include "ZTNBRegression.hpp"
#include "FittingMethod.hpp"

class RegressionBuilder {
public:
  static RegressionModel* buildRegression(const std::string& type,
                                          const size_t P,
                                          const FittingMethod& fitMthd) {
    if (type == "LinearRegression") {
      // todo - pass fitMthd and check in constructor
      return new LinearRegression(P);
    } else if (type == "PoissonRegression") {
      return new PoissonRegression(P, fitMthd);
    } else if (type == "NegativeBinomialRegression") {
      return new NegativeBinomialRegression(P, fitMthd);
    } else if (type == "ZeroTruncatedPoissonRegression") {
      // TODO warn ZTP can only be fit by NR and so we've discarded fitMthd?
      return new TruncatedPoissonRegression(0,P);
    } else if (type == "ZeroTruncatedNegativeBinomialRegression") {
      // TODO warn ZTNB can only be fit by NR and so we've discarded fitMthd?
      return new ZTNBRegression(P);
    }
    else {
      std::stringstream msg;
      msg << "Unkown regression type -> '" << type << "'";
      throw MixtureModelException(msg.str());
    }
  }
};


#endif





