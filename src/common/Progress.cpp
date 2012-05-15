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

#include "Progress.hpp"
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <ctime>

using std::string;
using std::ostringstream;
using std::cerr;
using std::endl;

/******************************************************************************
 * Progress Indicator class implementation 
 *****************************************************************************/

/****
 * @sumarry: constructor for ProgressIndicator 
 */
ProgressIndicator::ProgressIndicator(const long total, const string messagePrefix = "", const string messageSuffix = "") {
	this->total = total;
	this->done = 0;
	this->prefix = messagePrefix;
	this->suffix = messageSuffix;
	this->previous = "";
	this->finished = false;
	time (&(this->startTime));
	this->previousTimeMsg = " (estimated time remaining: <calculating>)";
	this->previousPercentage = -1;
	this->charsWritten = 0;
}

/****
 * @summary: update how much of the task has been completed
 */
void
ProgressIndicator::setDone(const long ndone){
	this->done = ndone;
}

/****
 * @summary: update how much of the task has been completed by one
 */
void
ProgressIndicator::incrementDone(){
	this->done++;
}

/****
 * @summary: display current progress message on stderr
 */
void
ProgressIndicator::showProgress() {
  if (!this->finished) {
	// figure out the percentage done
	int percent = (int) ceil (100 * this->done / (double) this->total);
	
	// if percentage done has changed, we'll issue a new message
	if (this->previousPercentage != percent) {
	  ostringstream s,t,b;
	  s << "\r" << this->prefix << " " << percent << "% "<< this->suffix;
		  
      // first, clear what was there before...
      b << "\r";
      for (int i=0; i<=this->charsWritten; i++) {
                      b << " ";
      }
      cerr << b.str();
      
      // work out time elapsed since start
      time_t now;
      time (&now);
      double elapsedSeconds = difftime (now,this->startTime);
		  
      if (percent == 100) t << " (estimated time remaining: done)"; 
      else {
        // time has very low resolution.. if we're getting zero for the 
        // amount done, don't do it this time, wait longer before we check again.
        long secondsLeft = (long) (((100-percent) / (double) percent) * elapsedSeconds);
        if (secondsLeft > 0.0) {
          if (secondsLeft < 60) {
            t << " (estimated time remaining: " << secondsLeft << " seconds)";
          } else {
            if (secondsLeft < 3600) {
              long minsLeft = secondsLeft / 60;
              secondsLeft = secondsLeft % 60;
              t << " (estimated time remaining: " << minsLeft << " minutes " << secondsLeft << " seconds)";
            } else {
              long hoursLeft = secondsLeft / 3600;
              secondsLeft = secondsLeft - (hoursLeft * 3600);
              long minsLeft = secondsLeft / 60;
              secondsLeft = secondsLeft - (minsLeft * 60);
              t << " (estimated time remaining: " << hoursLeft << " hours " << minsLeft << " minutes " << secondsLeft << " seconds)";
            }
          }
          this->previousTimeMsg = t.str();
        } else {
          t << this->previousTimeMsg;
        }
      }
			
      // output the message
      cerr << s.str() << t.str();
      cerr.flush();
      this->previous = s.str();
      this->previousPercentage = percent;
      this->charsWritten = s.str().size() + t.str().size();
	}
	if (percent == 100) {
	  cerr << endl;
	  this->finished = true;
	}
  }
}
