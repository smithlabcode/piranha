/**
  \file Matrix.cpp
  \brief This source file defines the matrix and diagonal matrix classes for
         matrix math

  \authors Philip J. Uren

  \copyright Copyright Details
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
  
  \section bug Known Bugs
  
  \section history Revision History
**/

#include "Matrix.hpp"

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <assert.h>
#include <cmath>

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::isfinite;
using std::stringstream;

// initialise static data members
const double Matrix::threshold = 0.001;

/******************************************************************************
 *  MATRIX CLASS
 */

/**
 * \brief copy constructor
 */
Matrix::Matrix(const Matrix& m) : rows(m.numRows()), cols(m.numCols()) {
  this->self = gsl_matrix_alloc(this->numRows(),this->numCols());
  gsl_matrix_memcpy(this->self, m.self);
}

/**
 * \brief construct a matrix of the given dimensions and set all of the
 *        elements to be equal to zero
 */
Matrix::Matrix(size_t r, size_t c) : 
     self(gsl_matrix_calloc(r,c)), rows(r), cols(c) {;}

/**
 * \brief construct a matrix object from a vector of vector of double
 * \todo fix parameters -- should be by reference
 */
Matrix::Matrix(const vector< vector<double> > a) : self(NULL) {
  init(a);
}

/**
 * \brief construct a matrix with a single row from a vector of double
 * \todo fix parameters -- should be by reference
 */
Matrix::Matrix(const vector<double> a) : self(NULL) {
  vector< vector<double> > t;
  t.push_back(a);
  init(t);
}

/**
 * \brief construct a matrix from a string representation
 */
Matrix::Matrix(string m) : self(NULL) {
  vector<string> lines = smithlab::squash(smithlab::split(m,"\n"));
  vector< vector<double> > data;
  for (size_t i=0; i<lines.size(); i++) {
    vector<string> tmp = smithlab::split(lines[i], ",");
    vector<double> row;
    for (size_t j=0; j<tmp.size(); j++) row.push_back(atof(tmp[j].c_str()));
    data.push_back(row);
  }
  init(data);
}

/**
 * \brief construct a matrix from a GSL matrix;
 *
 * this takes possession of the GSL matrix, which should not be modified
 * externally after passing to this constructor and will be freed by this
 * object when it is destroyed.
 */
Matrix::Matrix(gsl_matrix* s, size_t r, size_t c) : 
    self(s), rows(r), cols(c) {;}

/**
 * \brief matrix destructor, frees underlying memory
 */
Matrix::~Matrix() {
  if (this->self != NULL) gsl_matrix_free(this->self);
}

/**
 * \brief initialize the object with the given data.
 *
 * This is a private member and is intended for use only internally.
 * External initialization should be done through constructors.
 *
 * \param data A vector of data, must have at least 1 row and 1 column. Treated
 *             as row-major
 * \throws MatrixException if data doesn't have at least 1 row and 1 column,
 *                         or if the data is ragged (uneven row lengths)
 */
void
Matrix::init(const vector< vector<double> >& data) {
  // matrices can't be empty
  if (data.size() == 0) {
    stringstream ss;
    ss << "Failed to create matrix. Reason: provided data vector (rows) was "
       << "empty";
    throw MatrixException(ss.str());
  }
  if (data[0].size() == 0) {
    stringstream ss;
    ss << "Failed to create matrix. Reason: provided data vector (cols) was "
       << "empty";
    throw MatrixException(ss.str());
  }
  
  // matrices can't be ragged
  size_t rows = data.size(), cols = data[0].size();
  for (size_t row = 0; row < rows; row++) {
    if (data[row].size() != cols) {
      stringstream ss;
      ss << "Failed to create matrix. Reason: provided data has ragged rows. "
         << "Row 0 has " << cols << " but row "
         << row << " has " << data[row].size();
      throw MatrixException(ss.str());
    }
  }
  
  // all good, remember our instance variables.. 
  this->rows = rows;
  this->cols = cols;  
  this->self = gsl_matrix_calloc(this->rows,this->cols);

  for (size_t i = 0; i < this->rows; i++) {
    for (size_t j = 0; j < this->cols; j++) {
      if (!isfinite(data[i][j])) {
        stringstream ss;
        ss << "Failed to create matrix. Reason: element (" << i << "," << j
           << ")" << " is non-finite --> " << data[i][j];
        throw MatrixException(ss.str());
      }
      gsl_matrix_set (this->self, i, j, data[i][j]);
    }
  }
}

/**
 * \brief swap the contents of this matrix object with another
 */
void
Matrix::swap(Matrix& rhs) {
  std::swap(this->rows, rhs.rows);
  std::swap(this->cols, rhs.cols);
  std::swap(this->self, rhs.self);
}

/**
 * \brief get the number of rows in this matrix
 */
size_t Matrix::numRows() const { return this->rows; }

/**
 * \brief get the number of columns in this matrix
 */
size_t Matrix::numCols() const { return this->cols; }

/**
 * \brief return the maximum value in the ith row
 * \todo add unit test for this
 */
double
Matrix::maxInRow(size_t i) const {
  double maxVal = 0;
  bool first = true;
  for (size_t j=0; j<this->numCols(); j++) {
    double v = (*this)(i,j);
    if ((first) || (v > maxVal)) maxVal = v;
    first = false;
  }
  return maxVal;
}

/**
 * \brief return the maximum value in the ith col
 * \todo add unit test for this
 */
double
Matrix::maxInCol(size_t j) const {
  double maxVal = 0;
  bool first = true;
  for (size_t i=0; i<this->numRows(); i++) {
    double v = (*this)(i,j);
    if ((first) || (v > maxVal)) maxVal = v;
    first = false;
  }
  return maxVal;
}

/**
 * \brief return the maximum value in the ith row
 * \todo add unit test for this
 */
double
Matrix::minInRow(size_t i) const {
  double minVal = 0;
  bool first = true;
  for (size_t j=0; j<this->numCols(); j++) {
    double v = (*this)(i,j);
    if ((first) || (v < minVal)) minVal = v;
    first = false;
  }
  return minVal;
}

/**
 * \brief return the maximum value in the ith col
 * \todo add unit test for this
 */
double
Matrix::minInCol(size_t j) const {
  double minVal = 0;
  bool first = true;
  for (size_t i=0; i<this->numRows(); i++) {
    double v = (*this)(i,j);
    if ((first) || (v < minVal)) minVal = v;
    first = false;
  }
  return minVal;
}

/**
 * \brief multiply this matrix by another matrix
 */
Matrix 
Matrix::operator* (const Matrix &b) const {
  gsl_matrix* C = gsl_matrix_alloc(this->rows,b.cols);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, this->self, b.self, 0.0, C);
  return Matrix (C, this->rows, b.cols);
}

/**
 * \brief multiply this matrix by a diagonal matrix
 */
Matrix 
Matrix::operator* (const DiagonalMatrix& b) const {
  Matrix res(this->numRows(), b.numCols());
  for (size_t i=0; i<res.numRows(); i++) {
    for (size_t j=0; j<res.numCols(); j++) {
      double r = ((*this)(i,j) * b(j,j));
      if (!isfinite(r)) {
        stringstream ss;
        ss << "Matrix multiplication failed. Reason: resultant matrix is "
           << "non-finite at element (" << i << "," << j << ")" << " --> "
           << r;
        throw MatrixException(ss.str());
      }
      res(i,j) = r;
    }
  }
  return res;
}

/**
 * \brief multiply this matrix by a scalar
 */
Matrix
Matrix::operator* (const double s) const {  
  gsl_matrix* C = gsl_matrix_alloc(this->rows,this->cols);
  Matrix res(C,this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) {
      double r = (*this)(i,j) * s;
      if (!isfinite(r)) {
        stringstream ss;
        ss << "Matrix multiplication failed. Reason: resultant matrix is "
           << "non-finite at element (" << i << "," << j << ")" << " --> "
           << r;
        throw MatrixException(ss.str());
      }
      res(i,j) = r;
    }
  }
  return res;
}

/**
 * \brief add this matrix to another matrix
 */
Matrix
Matrix::operator+ (const Matrix &b) const {
  gsl_matrix* C = gsl_matrix_alloc(this->rows,this->cols);
  Matrix res(C,this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) {
      double r = (*this)(i,j) + b(i,j);
      if (!isfinite(r)) {
        stringstream ss;
        ss << "Matrix addition failed. Reason: resultant matrix is "
           << "non-finite at element (" << i << "," << j << ")" << " --> "
           << r;
        throw MatrixException(ss.str());
      }
      res(i,j) = r;
    }
  }
  return res;
}

/**
 * \brief add a scaler to this matrix
 */
Matrix
Matrix::operator+ (const double b) const {
  gsl_matrix* C = gsl_matrix_alloc(this->rows,this->cols);
  Matrix res(C,this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) {
      double r = (*this)(i,j) + b;
      if (!isfinite(r)) {
        stringstream ss;
        ss << "Matrix addition failed. Reason: resultant matrix is "
           << "non-finite at element (" << i << "," << j << ")" << " --> "
           << r;
        throw MatrixException(ss.str());
      }
      res(i,j) = r;
    }
  }
  return res;
}

/**
 * \brief subtract a matrix from this matrix
 */
Matrix
Matrix::operator- (const Matrix &b) const {
  gsl_matrix* C = gsl_matrix_alloc(this->rows,this->cols);
  Matrix res(C,this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) {
      double r = (*this)(i,j) - b(i,j);
      if (!isfinite(r)) {
        stringstream ss;
        ss << "Matrix subtraction failed. Reason: resultant matrix is "
           << "non-finite at element (" << i << "," << j << ")" << " --> "
           << r;
        throw MatrixException(ss.str());
      }
      res(i,j) = r;
    }
  }
  return res;
}

/**
 * \brief subtract a scalar from this matrix
 */
Matrix
Matrix::operator- (const double b) const {
  gsl_matrix* C = gsl_matrix_alloc(this->rows,this->cols);
  Matrix res(C,this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) {
      double r = (*this)(i,j) - b;
      if (!isfinite(r)) {
        stringstream ss;
        ss << "Matrix subtraction failed. Reason: resultant matrix is "
           << "non-finite at element (" << i << "," << j << ")" << " --> "
           << r;
        throw MatrixException(ss.str());
      }
      res(i,j) = r;
    }
  }
  return res;
}

/**
 * \brief compare two matrices for equality
 */
bool 
Matrix::operator== (const Matrix &b) const {
  if (this->rows != b.rows) return false;
  if (this->cols != b.cols) return false;
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) {
      if (fabs((*this)(i,j) - b(i,j)) > Matrix::threshold) return false;  
    }
  }
  return true; 
}

/**
 * \brief calculate the transposition of this matrix
 */
Matrix
Matrix::transpose() const {
  size_t newRows = this->cols;
  size_t newCols = this->rows;
  
  // make an empty matrix..
  vector< vector<double> > res;
  for (size_t row = 0; row < newRows; row++) {
    res.push_back(vector<double> (newCols, 0));
  }
  
  // fill it in
  for (size_t row = 0; row < newRows; row++) {
    for (size_t col = 0; col < newCols; col++) {
      res[row][col] = (*this)(col,row);
    }
  }
  
  // done
  return Matrix(res);
}

/**
 * \brief calculate the inverse of this matrix.
 *
 * This is just a wrapper for the corresponding method in gsl
 */
Matrix
Matrix::inverse() const {
  // make a gsl matrix first..
  gsl_matrix* m = gsl_matrix_alloc(this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) 
    for (size_t j=0; j<this->cols; j++) {
      gsl_matrix_set(m,i,j,(*this)(i,j));
    }
  
  // decompose
  gsl_permutation* p = gsl_permutation_alloc(this->rows);
  int sign;
  gsl_linalg_LU_decomp(m, p, &sign);  
  
  // get the inverse as a gsl matrix
  gsl_matrix* inv = gsl_matrix_alloc(this->rows, this->cols);
  gsl_linalg_LU_invert(m, p, inv);
    
  // now make a Matrix object from the inverse gsl matrix
  Matrix res(this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) 
      for (size_t j=0; j<this->cols; j++) {
        res(i,j) = gsl_matrix_get(inv,i,j);
      }
  
  // done..
  gsl_matrix_free(m);
  return res;
}

/**
 * \brief calculate the determinant of this matrix.
 *
 * This is just a wrapper for the corresponding method in gsl
 */
double
Matrix::determinant() const {
  // make a gsl matrix first..
  gsl_matrix* m = gsl_matrix_alloc(this->rows, this->cols);
  for (size_t i=0; i<this->rows; i++) 
    for (size_t j=0; j<this->cols; j++) {
      gsl_matrix_set(m,i,j,(*this)(i,j));
    }
  
  // decomp..
  gsl_permutation* p = gsl_permutation_alloc(this->rows);
  int sign;
  gsl_linalg_LU_decomp(m, p, &sign);
  
  // return determinant..
  double res = exp(gsl_linalg_LU_lndet(m));
  gsl_matrix_free(m);
  return res;
}

/**
 * \brief row-indexing operator, returns a vector representing the
 *        data in the requested row
 *
 * This returns a copy of the data, not a reference to it, so it cannot be
 * used on the LHS of assignment
 */
const vector<double>
Matrix::operator() (size_t row) const {
  vector<double> res;
  for (size_t j=0; j<this->cols; j++) res.push_back((*this)(row,j));
  return res;
}

/**
 * \brief return a vector representation of this matrix. Only works if the
 *        matrix has only a single column or a single row.
 */
vector<double> 
Matrix::asVector() const {
  // one of the dimensions of this matrix must be 1 for this to work.
  if ((this->rows != 1) && (this->cols != 1)) {
    stringstream ss;
    ss << "converting " << this->numCols() << " x " << this->numRows() << " "
       << "failed. Reason: at least one dimension must be equal to 1";
    throw MatrixException(ss.str());
  }
  
  vector<double> res;
  if (this->rows == 1) {
    res = (*this)(0);
  } else if (this->cols == 1) {
    for (size_t i=0; i<this->rows; i++) res.push_back((*this)(i,0));
  }  
  return res;
}

/**
 * \brief return a copy of this matrix as a vector of vectors, in row-major
 */
vector< vector<double> >
Matrix::asVectorOfVector() const {
  vector< vector<double> > res;
  for (size_t i = 0; i < this->numRows(); i++) {
    vector<double> tmp = (*this)(i);
    res.push_back(tmp);
  }
  return res;
}

/**
 * \brief return a string representation of this matrix object
 * \todo pad cells so columns line up
 */
string
Matrix::toString() const {
  stringstream ss;
  for (size_t i=0; i<this->rows; i++) {
    for (size_t j=0; j<this->cols; j++) 
      ss << (*this)(i)[j] << '\t';
    ss << '\n';
  }
  return ss.str();
}


/******************************************************************************
 * DIAGONAL MATRIX CLASS
 */

/**
 * \brief construct a diagonal matrix from a gsl matrix
 * \bug   not implemented yet
 * \todo  parse the gsl_matrix, throw exception if not diagonal
 */
DiagonalMatrix::DiagonalMatrix(gsl_matrix* s, size_t r, size_t c) {
  throw MatrixException("not implemented");
}

/**
 * \brief construct a matrix of the given dimensions and set all of the
 *           elements to be equal to zero.
 *
 * We allow this so as to duplicate the normal Matrix construtors, but we
 * force square dimensions.
 *
 * \todo include dimensions found in error
 */
DiagonalMatrix::DiagonalMatrix(size_t rows, size_t cols) :
    self(vector<double>(rows,0)) {
  if (rows != cols)
    throw MatrixException("Diagonal matrix should be square");
}

/**
 * \brief override this ctor and cry if somebody tries to use it; we
 *        don't know how to handle this yet..
 * \todo  parse matrix, throw exception if not diagonal
 * \todo  check that this is marked explicit
 */
DiagonalMatrix::DiagonalMatrix(vector< vector<double> > a) {
  throw MatrixException("not implemented");
}

/**
 * \brief behaves differently to equiv. ctor in base - here we treat the
 *        vector as being the diagonal of the matrix, not the first row
 */
DiagonalMatrix::DiagonalMatrix(vector<double> a) : self(a) {;}

/**
 * \brief construct a matrix from a string representation -- not
 *        implemented yet
 * \todo  check that this is marked explicit
 * \todo  parse the string, build matrix, throw exception if not diagonal
 */
DiagonalMatrix::DiagonalMatrix(string m) {
  throw MatrixException("not implemented");
}

/**
 * \brief initialise the object from a vector of vectors, checking for
 *        violation of constraints.
 * \todo  how to handle empty vectors?
 * \todo  not implemented yet, should check for diagonal data, ragged 2d array.
 */
void
DiagonalMatrix::init(vector< vector<double> > data) {
  throw MatrixException("not implemented");
}

/**
 * \brief multiply this diagonal matrix by a normal matrix
 *
 * In theory this is the same as the normal matrix class class, but that uses
 * CBLAS_DGEMM, which we don't want to do because it requires us to explicitly
 * build the full diagonal matrix, including zeros -> waste of space and time.
 */
Matrix
DiagonalMatrix::operator* (const Matrix &b) const {
  if (this->numCols() != b.numRows()) {
    stringstream ss;
    ss << "Diagonal matrix multiplication failed. Reason: columns in first "
       << "matrix (" << this->numCols() << ") does not match rows in second "
       << "matrix (" << b.numRows() << ")";
    throw MatrixException(ss.str());
  }

  Matrix res(this->numRows(), b.numCols());
  for (size_t i=0; i<res.numRows(); i++) {
    for (size_t j=0; j<res.numCols(); j++) {
      for (size_t k=0; k<this->numCols(); k++) {
        res(i,j) += ((*this)(i,k) * b(k,j));
        if (!isfinite(res(i,j))) {
          stringstream ss;
          ss << "Diagonal matrix multiplication failed. Reason: resultant "
             << "matrix is non-finite at element (" << i << "," << j << ")"
             << " --> " << res(i,j);
          throw MatrixException(ss.str());
        }
      }
    }
  }
  return res;
}

/**
 * \brief find the inverse of this matrix.
 * \todo can this be done without needing the full matrix? Even if it's slower
 *       than gsl, we'd save a lot overall by not having to build the full
 *       matrix
 */
Matrix
DiagonalMatrix::inverse() const {
  throw MatrixException("not implemented");
}

/**
 * \brief find the determinant of this matrix.
 * \todo  can this be done without needing the full matrix? Even if it's
 *        slower than gsl, we'd save a lot overall by not having to build
 *        the full matrix
 */
double
DiagonalMatrix::determinant() const {
  throw MatrixException("not implemented");
}

