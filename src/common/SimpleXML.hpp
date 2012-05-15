/**
  \file SimpleXML.hpp
  \brief This header declares the XMLToken, XMLNode classes, the SimpleXML
         type and the XMLException class. These are all used for parsing
         simple XML files

  \author Philip J. Uren

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren

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

#ifndef SIMXML_HPP
#define SIMXML_HPP

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

// from smithlab common
#include "smithlab_utils.hpp"

/**
 * \brief An exception class for XML exceptions
 */
class XMLException : public SMITHLABException {
public:
  XMLException(std::string s = std::string()) : SMITHLABException(s) {}
};

/**
 * \brief An XMLToken is either a single start tag, a single end tag or a
 *        data string
 *
 * This class if mainly used as a helper for tokenizing an XML string/stream
 */
class XMLToken { 
public:
	XMLToken(std::string s) : name(s), istag(false), isendtag(false) {;}
	XMLToken(std::string s, bool end) : name(s), istag(true), isendtag(end) {;}

	bool isTag() const { return this->istag; }
	bool isData() const { return !this->isTag(); }
	bool isEndTag() const { return this->isendtag; }
	std::string getName() const { return this->name; }
	std::string toString() const;

	/*** public static functions ***/
	static std::vector<XMLToken> tokenize(std::istream& istrm);

private:
	/** \brief name of this token (tag name if it's a tag, else data value) **/
	std::string name;

	/** \brief true if this is a tag, false if it's just data **/
	bool istag;

	/** \brief true if this is an end tag (e.g. </someTag>), else false **/
	bool isendtag;
};

/**
 * \brief An XML Node tag pair with the same name and has either a data body
 *        or other nested XMLNodes
 */
class XMLNode {
public:
  /*** constructors, destructors and object initialisation ***/
  XMLNode(std::string fn);
  XMLNode(std::istream& istrm);
  XMLNode(const XMLNode& x);
  XMLNode(std::vector<XMLToken>::iterator tstart, 
      std::vector<XMLToken>::iterator tend) {
  	this->build(tstart, tend);
  }
  ~XMLNode() {
    for (size_t i=0; i<this->children.size(); i++) 
      delete this->children[i];
  }
 
  /*** building the tree ***/
  void build(std::vector<XMLToken>& tokens) {
    this->build(tokens.begin(), tokens.end()-1);
  }
  void build(std::vector<XMLToken>::iterator tstart, 
      std::vector<XMLToken>::iterator tend);

  /*** inspectors ***/
  std::string toString() const;
  std::string toStringPretty(size_t tab) const;
  size_t numChrildren() const { return this->children.size(); }
  std::string getTagName () const { return this->tagName; }
  bool isLeaf() const { return this->numChildren() == 0; }
  std::string getData() const { return data; }
  size_t numChildren() const { return this->children.size(); }
  
  /*** inspectors -- finding sub trees ***/
  std::vector<XMLNode> getChildren() const;
  std::vector<XMLNode> getChildren(const std::string& tagname) const;
  
private:
  /** \brief name of tag forming this XMLNode (e.g. <name> ... </name>) **/
  std::string tagName;
  /** \brief data in this node; node will have data or children, not both **/
  std::string data;
  /** \brief child nodes nested within this node; note that the node will
   *         contain either data or children, but not both
   */
  std::vector<XMLNode*> children;
};

typedef XMLNode SimpleXML;

#endif

