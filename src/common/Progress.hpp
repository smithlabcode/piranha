/****   
 
	Description:   Simple class for tracking and displaying progress towards
								 completing a task 

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

#ifndef PROG_HPP
#define PROG_HPP

#include <string>
#include <sys/time.h>

/******************************************************************************
 * ProgressIndicator class definition
 *****************************************************************************/

class ProgressIndicator {
public:
	ProgressIndicator(const long total, const std::string messagePrefix, const std::string messageSuffix);
	void showProgress();
	void setDone(const long ndone);
	void incrementDone();
private:
	long total;
	long done;
	std::string prefix;
	std::string suffix;
	
	std::string previous;
	std::string previousTimeMsg;
	
	int previousPercentage;
	
	int charsWritten;
	
	bool finished;
	time_t startTime;
};

#endif 
