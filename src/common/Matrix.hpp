/**
  \file Matrix.hpp
  \brief This header file declares the matrix and diagonal matrix classes for
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

#ifndef MAT_HPP
#define MAT_HPP

#include <vector>
#include <string>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

// from smithlab common
#include "smithlab_utils.hpp"

// from Piranha common
#include "Matrix.hpp"

// forward declare the diagonal matrix class
class DiagonalMatrix;

/**
 * \brief Exception class for matrix-related exceptions
 */
class MatrixException : public SMITHLABException {
public:
  MatrixException(std::string s = std::string()) : SMITHLABException(s) {}
};

/**
 * \brief Used for representing matrix data and performing matrix math
 *        operations
 */
class Matrix {
public :  
  /*** constructors, destructors ***/
  Matrix(size_t r, size_t c);
  Matrix(const std::vector< std::vector<double> > a);
  Matrix(const std::vector<double> a);
  explicit Matrix(std::string);
  Matrix(const Matrix& m);
  ~Matrix();
  
  /*** mutators ***/
  void swap(Matrix& rhs);

  /*** non-const operators ***/
  Matrix& operator=(const Matrix& rhs) {
    Matrix tmp(rhs);
    this->swap(tmp);
    return *this;
  }
  double& operator() (size_t row, size_t col) {
    return *gsl_matrix_ptr(this->self, row, col);
  }

  /*** const operators ***/
  Matrix operator* (const Matrix &b) const;
  Matrix operator* (const double s) const;
  Matrix operator* (const DiagonalMatrix& s) const;
  Matrix operator+ (const Matrix &b) const;
  Matrix operator+ (const double b) const;
  Matrix operator- (const Matrix &b) const;
  Matrix operator- (const double b) const;
  bool operator== (const Matrix &b) const;
  const std::vector<double> operator() (size_t row) const;
  const double& operator() (size_t row, size_t col) const {
    return *gsl_matrix_ptr(this->self, row, col) ;
  }
  
  /*** inspectors ***/
  size_t numRows() const;
  size_t numCols() const;
  double maxInRow(size_t i) const;
  double maxInCol(size_t i) const;
  double minInRow(size_t i) const;
  double minInCol(size_t i) const;
  std::string toString() const;
  std::vector<double> asVector() const;
  std::vector< std::vector<double> > asVectorOfVector() const;
  double determinant() const;
  Matrix transpose() const;
  Matrix inverse() const;

  /*** Matrix I/O ***/
  friend std::ostream& operator<<(std::ostream& s, const Matrix &m) {
    return s << m.toString();
  }

private:
  /*** private constants ***/
  static const double threshold;

  /*** private instance variables ***/
  /** \brief contains the actual data for the matrix. **/
  gsl_matrix* self;
  /** \brief the number of rows in the matrix **/
  size_t rows;
  /** \brief the number of columns in the matrix **/
  size_t cols;

  /*** private constructor ***/
  Matrix(gsl_matrix* s, size_t rows, size_t cols);

  /*** private methods ***/
  void init(const std::vector< std::vector<double> >& d);
};

/**
 * \brief A diagonal matrix contains only zeros off the left to right, top to
 *        bottom diagonal. It behaves exactly like a normal matrix except we
 *        don't store the 0s that are off the diagonal to save space.
 */
class DiagonalMatrix {
public :
  /*** constructors, destructors ***/
  DiagonalMatrix(gsl_matrix* s, size_t r, size_t c);
  DiagonalMatrix(size_t rows, size_t cols);
  DiagonalMatrix(std::vector< std::vector<double> > a);
  DiagonalMatrix(std::vector<double> a);
  DiagonalMatrix(std::string m);
  ~DiagonalMatrix() {;}

  /*** operators, mutators ***/
  const double& operator() (size_t row, size_t col) const {
    if (row != col)
      throw MatrixException("DiagonalMatrix - out of bounds");
    return this->self[row];
  }
  double& operator() (size_t row, size_t col) {
    if (row != col)
      throw MatrixException("DiagonalMatrix - out of bounds");
    return this->self[row];
  }
  Matrix operator* (const Matrix &b) const;

  /*** inspectors, const methods ***/
  size_t numRows() const { return this->self.size(); }
  size_t numCols() const { return this->self.size(); }
  Matrix inverse() const;
  double determinant() const;

private:
  /*** private methods ***/
  void init(std::vector< std::vector<double> > data);

  /*** private data members ***/
  /** \brief the actual underlying data for the matrix **/
  std::vector<double> self;
};

#endif
