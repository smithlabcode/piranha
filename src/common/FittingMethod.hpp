/****
  @summary: This class is a representation of fitting method for models.
            Right now it's really just a simple wrapper around string, but
            this might change later.

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


#ifndef FMT_HPP
#define FMT_HPP

#include <string>

class FittingMethod {
public:
  explicit FittingMethod(std::string s) : type(s) {;}
  bool is(std::string s) const { return this->type == s; }
  std::string toString() const { return this->type; }
private:
  std::string type;
};

static FittingMethod defaultFittingMethod("DEFAULT");

#endif
