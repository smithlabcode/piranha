/**
  \file SimpleXML.cpp
  \brief This source file defines the XMLToken, XMLNode classes, the SimpleXML
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

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include "SimpleXML.hpp"

using std::vector;
using std::string;
using std::istream;
using std::ifstream;
using std::stringstream;
using std::cerr;
using std::cout;
using std::endl;

/******************************************************************************
 * Static helper functions
 */

/**
 * \brief Return true if <ch> is a whitespace character (tab, space or newline)
 */
static bool
isWhiteSpace (char ch) {
  if ((ch == ' ') ||  
      (ch == '\t') ||  
      (ch == '\n')) 
    return true; 
  return false; 
}

/**
 * \brief Remove whitespace from stream until the next char is a non-whitespace
 *        char (but don't remove it).
 */
static void
chomp (std::istream& istrm) {
  while (istrm.good() && isWhiteSpace(istrm.peek())) istrm.get();
};

/******************************************************************************
 * XMLToken Class
 */

/**
 * \brief Tokenize a stream of XML by chopping it up into XMLTokens and
 *        returning a vector of these tokens.
 *
 * This function consumes all data in the stream
 */
vector<XMLToken>
XMLToken::tokenize(istream& istrm) {
  vector<XMLToken> tokens;
  while (istrm.good()) {
    chomp(istrm);
    if (!istrm.good()) break;

    // parse tag
    if (istrm.peek() == '<') {
      string tagname;
      bool isEndTag = false;

      istrm.get(); // skip <
      chomp(istrm);
      if (!istrm.good()) break;
      if (istrm.peek() == '/') {
          istrm.get(); // skip /
          isEndTag = true;
      }

      while (istrm.good() && (istrm.peek() != '>'))
          tagname += istrm.get();
      istrm.get(); // skip >
      tokens.push_back(XMLToken(tagname,isEndTag));
    } else {
      // parse string
      string buf = "";
      while (istrm.good() && (istrm.peek() != '<'))
          buf += istrm.get();
      tokens.push_back(XMLToken(buf));
    }
  }
  return tokens;
}

/**
 * \brief Return a string representation of this XMLToken object
 */
string
XMLToken::toString() const {
  stringstream ss;
  if (this->isTag()) ss << "<";
  if (this->isEndTag()) ss << "/";
  ss << this->getName();
  if (this->isTag()) ss << ">";
  return ss.str();
}

/******************************************************************************
 * XMLNode Class
 */

/**
 * \brief  Constructor, builds the XMLNode object by reading it from a
 *         file
 * \throws XMLException if the file couldn't be opened
 */
XMLNode::XMLNode(std::string fn) {
  ifstream strm(fn.c_str());
  if (!strm) {
    stringstream ss;
    ss << "failed to open file " << fn << " for reading XML object";
    throw XMLException(ss.str());
  }
  vector<XMLToken> tokens = XMLToken::tokenize(strm);
  this->build(tokens);
}

/**
 * \brief Constructor, builds the XMLNode object by reading it from an
 *        input stream (istream)
 *
 * This constructor will consume the whole stream
 */
XMLNode::XMLNode(std::istream& istrm) {
  vector<XMLToken> tokens = XMLToken::tokenize(istrm);
  this->build(tokens);
}


/**
 * \brief Copy constructor
 */
XMLNode::XMLNode(const XMLNode& x) {
  this->tagName = x.tagName;
  this->data = x.data;
  // handle with care..
  for (size_t i=0; i<x.children.size(); i++) {
    this->children.push_back(new XMLNode(*(x.children[i])));
  }
}


/**
* \brief  build a tree from a vector of tokens.
*
* We assume the tokens have a single outer tag for the root; this is checked.
*
* \param  tstart iterator pointing to first element in XMLToken sequence
* \param  tend iterator pointing to one beyond the last element in the XMLToken
*         sequence
* \throws XMLException if the first and last XMLTokens are not tags
*                      or if the first and last tags don't have matching names
*                      or it the last tag is not a closing tag (i.e. </tag>)
*                      or if there is only one Token between the start and end
*                            tags and it's a tag (i.e. orphaned tag)
*                      or if a tag nests both child tags and a data element
*/
void 
XMLNode::build(vector<XMLToken>::iterator tstart,
               vector<XMLToken>::iterator tend) {
  if (!tstart->isTag() || !tend->isTag()) {
    stringstream ss;
    ss << "failed to build XML node from XMLToken vector. Reason: first or"
       << "last token is not XML tag";
    throw XMLException(ss.str());
  }
	
  // first and last tokens are known to be tags now
  if (tstart->getName() != tend->getName()) {
    stringstream ss;
    ss << "failed to build XML node from XMLToken vector."
       << "Reason: first and last tags have different names";
    throw XMLException(ss.str());
  }
  if (!tend->isEndTag()) {
    stringstream ss;
    ss << "failed to build XML node from XMLToken vector."
       << "Reason: final tag is not a closing tag";
    throw XMLException(ss.str());
  }

  // if we get this far, first and last tags are present, there's nothing
  // else before/after them and they match each other, and the final 
  // tag is a valid closing tag
  this->tagName = tstart->getName();

  // this might be a leaf with no data, in which case we're done
  if (tstart+1 == tend) return;

  // or it might be a leaf with data
  if (tstart+2 == tend) {
    if ((tstart+1)->isTag()) {
      stringstream ss;
      ss << "failed to build XML node from XMLToken vector."
         << "Reason: XMLToken appears to contain an orphan tag";
      throw XMLException(ss.str());
    }
    this->data = (tstart+1)->toString(); 
  }
  // okay, it's definitely other tag pairs
  else {
    vector<XMLToken>::iterator childStart = tstart + 1;
    while (childStart < tend) {
      if (!childStart->isTag()) {
        stringstream ss;
        ss << "failed to build XML node from XMLToken vector."
           << "Reason: XMLToken has mixed data and child tags";
        throw XMLException(ss.str());
      }

      // scan for the closing tag
      vector<XMLToken>::iterator childEnd = childStart + 1;
      while ((!childEnd->isTag()) ||
             (childEnd->getName() != childStart->getName())) childEnd++;
			
      // build child tree
      this->children.push_back(new XMLNode(childStart, childEnd));

      childStart = childEnd + 1;	
    }
  }
}


/**
 * \brief get a vector of child nodes from this one which have a given tag
 *        name.
 *
 * The resulting vector might be empty. This method returns copies of the
 * nodes, not references to them, so this cannot be used to mess around with
 * the underlying data
 */
std::vector<XMLNode> 
XMLNode::getChildren(const std::string& tagname) const {
  std::vector<XMLNode> res;
  for (size_t i=0; i<this->numChildren(); i++) {
    if (this->children[i]->getTagName() == tagname) 
      res.push_back(*(this->children[i]));
  }
  return res;
}

/**
 * \brief get a vector of all child nodes nested in this XMLNode
 *
 * The resulting vector might be empty. This method returns copies of the
 * nodes, not references to them, so this cannot be used to mess around with
 * the underlying data
 */
std::vector<XMLNode> 
XMLNode::getChildren() const {
  std::vector<XMLNode> r;
  size_t n = this->numChildren();
  for (size_t i=0; i<n; i++) r.push_back(*(this->children[i]));
  return r;
}

/**
 * \brief Return a 'pretty' string representation of this XMLNode object
 * \param tab This is used for formatting. All lines will be tabbed in by
 *            this many tabs. Nested XMLNodes will be tabbed in by tab+1, etc.
 */
string 
XMLNode::toStringPretty(size_t tab) const {
  string tabs(tab,'\t');
  stringstream ss;
  ss << tabs << "<" << tagName << ">";
  if (this->children.size() == 0) {
      ss << " " <<  this->data << " ";
  }
  else {
      for (size_t i=0; i<this->children.size(); i++) 
        ss << std::endl << this->children[i]->toStringPretty(tab+1);
  }
  ss << tabs << "</" << this->tagName << ">" << endl;
  return ss.str();
}

/**
 * \brief Return a string representation of this XMLNode object
 */
string
XMLNode::toString() const {
  return this->toStringPretty(0);
}

