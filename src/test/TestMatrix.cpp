/**
  \file TestMatrix.cpp
  \brief This source file defines a set of unit tests for the Matrix class

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
****/

#include <vector>

#include "Matrix.hpp"
#include "gtest/gtest.h"

using std::vector;
using std::cerr;
using std::endl;

class MatrixTest : public testing::Test {
protected:
  virtual void SetUp() {
    // 0 1 2 3
    // 1 2 3 4
    // 2 3 4 5
    // 3 4 5 6
    for (size_t i=0; i<4; i++) {
      vector<double> tmp;
      for (size_t j =0; j<4; j++) {
        tmp.push_back(i+j);
      }
      m1v.push_back(tmp);
    }
    m1 = new Matrix(m1v);

    // 0 0 0 0
    // 0 1 2 3
    // 0 2 4 6
    // 0 3 6 9
    for (size_t i=0; i<4; i++) {
      vector<double> tmp;
      for (size_t j=0; j<4; j++) {
        tmp.push_back(i*j);
      }
      m2v.push_back(tmp);
    }
    m2 = new Matrix(m2v);

    // 15 16
    // 27 36
    // 4  88
    vector<double> tmp1, tmp2, tmp3;
    tmp1.push_back(15);
    tmp1.push_back(16);
    tmp2.push_back(27);
    tmp2.push_back(36);
    tmp3.push_back(4);
    tmp3.push_back(88);
    m3v.push_back(tmp1);
    m3v.push_back(tmp2);
    m3v.push_back(tmp3);
    m3 = new Matrix(m3v);
  }

  virtual void TearDown() {
    delete m1;
    delete m2;
    delete m3;
  }


  // Declares the variables your tests want to use.
  vector< vector<double> > m1v;
  vector< vector<double> > m2v;
  vector< vector<double> > m3v;
  Matrix* m1;
  Matrix* m2;
  Matrix* m3;
};

/**
 * \brief test getting max in column
 */
TEST_F(MatrixTest, testMaxInCol) {
  EXPECT_EQ(3, m1->maxInCol(0));
  EXPECT_EQ(4, m1->maxInCol(1));
  EXPECT_EQ(5, m1->maxInCol(2));
  EXPECT_EQ(6, m1->maxInCol(3));

  EXPECT_EQ(0, m2->maxInCol(0));
  EXPECT_EQ(3, m2->maxInCol(1));
  EXPECT_EQ(6, m2->maxInCol(2));
  EXPECT_EQ(9, m2->maxInCol(3));

  EXPECT_EQ(27, m3->maxInCol(0));
  EXPECT_EQ(88, m3->maxInCol(1));
}

/**
 * \brief test getting max in row
 */
TEST_F(MatrixTest, testMaxInRow) {
  EXPECT_EQ(3, m1->maxInRow(0));
  EXPECT_EQ(4, m1->maxInRow(1));
  EXPECT_EQ(5, m1->maxInRow(2));
  EXPECT_EQ(6, m1->maxInRow(3));

  EXPECT_EQ(0, m2->maxInRow(0));
  EXPECT_EQ(3, m2->maxInRow(1));
  EXPECT_EQ(6, m2->maxInRow(2));
  EXPECT_EQ(9, m2->maxInRow(3));

  EXPECT_EQ(16, m3->maxInRow(0));
  EXPECT_EQ(36, m3->maxInRow(1));
  EXPECT_EQ(88, m3->maxInRow(2));
}

/**
 * \brief test adding two matrices together
 */
TEST_F(MatrixTest, testAddition) {
  vector< vector<double> > ansn;
  for (size_t i=0; i<4; i++) {
    vector<double> tmp;
    for (size_t j=0; j<4; j++) {
      tmp.push_back(i * j + i + j);
    }
    ansn.push_back(tmp);
  }
  EXPECT_EQ(Matrix(ansn), (*m1) + (*m2));
}

/**
 * \brief test subtracting two matrices
 */
TEST_F(MatrixTest, testSubtraction) {
  vector< vector<double> > ansn;
  for (size_t i=0; i<4; i++) {
    vector<double> tmp;
    for (size_t j=0; j<4; j++) {
      double v = i + j - double(i*j);
      tmp.push_back(v);
    }
    ansn.push_back(tmp);
  }
  EXPECT_EQ(Matrix(ansn), (*m1) - (*m2));
}

/**
 * \brief test multiplying two matrices
 */
TEST_F(MatrixTest, testMultiplication) {
  int r1[] = {0, 14, 28, 42};
  int r2[] = {0, 20, 40, 60};
  int r3[] = {0, 26, 52, 78};
  int r4[] = {0, 32, 64, 96};
  vector< vector<double> > ansn;
  ansn.push_back(vector<double>(r1, r1+4));
  ansn.push_back(vector<double>(r2, r2+4));
  ansn.push_back(vector<double>(r3, r3+4));
  ansn.push_back(vector<double>(r4, r4+4));
  EXPECT_EQ(Matrix(ansn), (*m1) * (*m2));
}

/**
 * \brief test the conversion from Matrix to vector< vector<double> >
 */
TEST_F(MatrixTest, testToVectorOfVector) {
  EXPECT_EQ(m1->asVectorOfVector(), m1v);
  EXPECT_EQ(m2->asVectorOfVector(), m2v);
}

/**
 * \brief test transposing matrix columns and rows
 */
TEST_F(MatrixTest, testTranspose) {
  // 0 14 28 42       0  0  0  0
  // 0 20 40 60  -->  14 20 26 32
  // 0 26 52 78  -->  28 40 52 64
  // 0 32 64 96       42 60 78 96
  int r1[] = {0,  0,  0,  0};
  int r2[] = {14, 20, 26, 32};
  int r3[] = {28, 40, 52, 64};
  int r4[] = {42, 60, 78, 96};
  vector< vector<double> > ansn;
  ansn.push_back(vector<double>(r1, r1+4));
  ansn.push_back(vector<double>(r2, r2+4));
  ansn.push_back(vector<double>(r3, r3+4));
  ansn.push_back(vector<double>(r4, r4+4));
  EXPECT_EQ(Matrix(ansn), ((*m1) * (*m2)).transpose());
}
