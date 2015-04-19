/**
  \file TestReadBinner.cpp
  \brief This source file defines a set of unit tests for the ReadBinner class

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

#include "gtest/gtest.h"

#include "ReadBinner.hpp"
#include "GenomicRegion.hpp"

using std::vector;

/**
 * \brief test binning reads into single nucleotide bins
 */
TEST(ReadBinnerTest, testSingleNucBin) {
  vector<GenomicRegion> input;
  input.push_back(GenomicRegion("chr1", 2, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 3, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 6, 8,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 9, 11, "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 8,  "X", 0, '+'));

  vector<GenomicRegion> expect;
  expect.push_back(GenomicRegion("chr1", 2, 3,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 3, 4,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 5, 6,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 6, 7,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 9, 10,  "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 5, 6,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 6, 7,   "X", 2, '+'));

  vector<GenomicRegion> got;
  ReadBinner r(1);
  const bool UNSTRANDED = true;
  r.binReads(input, got, UNSTRANDED, 0);

  EXPECT_EQ(got.size(), expect.size());
  if (got.size(), expect.size()) {
    for (size_t i = 0; i < expect.size(); i++) {
      EXPECT_EQ(expect[i], got[i]);
    }
  }
}

TEST(ReadBinnerTestStrand, testPreservingStrand) {
  vector<GenomicRegion> input;
  input.push_back(GenomicRegion("chr1", 2, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 3, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 6, 8,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 9, 11, "X", 0, '-'));
  input.push_back(GenomicRegion("chr2", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 8,  "X", 0, '+'));

  vector<GenomicRegion> expect;
  expect.push_back(GenomicRegion("chr1", 2, 3,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 3, 4,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 5, 6,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 6, 7,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 9, 10,  "X", 1, '-'));
  expect.push_back(GenomicRegion("chr2", 5, 6,   "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 6, 7,   "X", 2, '+'));

  vector<GenomicRegion> got;
  ReadBinner r(1);
  const bool UNSTRANDED = false;
  r.binReads(input, got, UNSTRANDED, 0);

  EXPECT_EQ(got.size(), expect.size());
  if (got.size(), expect.size()) {
    for (size_t i = 0; i < expect.size(); i++) {
      EXPECT_EQ(expect[i], got[i]);
    }
  }
}

/**
 * \brief test binning reads into two nucleotide bins where we require
 *        some bins to exist even if they have zero read-counts, and
 *        adding pseudo-count of 1 to all accepted bins
 */
TEST(ReadBinnerTest, testTwoNucBin_require) {
  vector<GenomicRegion> input;
  input.push_back(GenomicRegion("chr1", 2, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 3, 5,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 6, 8,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr1", 9, 11, "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 5, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 7,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr2", 6, 8,  "X", 0, '+'));
  input.push_back(GenomicRegion("chr4", 0, 5,  "X", 0, '+'));

  vector<GenomicRegion> requiredBins;
  requiredBins.push_back(GenomicRegion("chr1", 0,  2,  "X", 0, '+'));
  requiredBins.push_back(GenomicRegion("chr1", 6,  8,  "X", 0, '+'));
  requiredBins.push_back(GenomicRegion("chr1", 10, 12, "X", 0, '+'));
  requiredBins.push_back(GenomicRegion("chr2", 0,  2,  "X", 0, '+'));
  requiredBins.push_back(GenomicRegion("chr2", 6,  8,  "X", 0, '+'));
  requiredBins.push_back(GenomicRegion("chr2", 8,  10, "X", 0, '+'));
  requiredBins.push_back(GenomicRegion("chr3", 8,  10, "X", 0, '+'));

  vector<GenomicRegion> expect;
  expect.push_back(GenomicRegion("chr1", 0,  2,  "X", 1, '+'));
  expect.push_back(GenomicRegion("chr1", 2,  4,  "X", 3, '+'));
  expect.push_back(GenomicRegion("chr1", 4,  6,  "X", 2, '+'));
  expect.push_back(GenomicRegion("chr1", 6,  8,  "X", 2, '+'));
  expect.push_back(GenomicRegion("chr1", 8,  10, "X", 2, '+'));
  expect.push_back(GenomicRegion("chr1", 10, 12, "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 0,  2,  "X", 1, '+'));
  expect.push_back(GenomicRegion("chr2", 4,  6,  "X", 2, '+'));
  expect.push_back(GenomicRegion("chr2", 6,  8,  "X", 3, '+'));
  expect.push_back(GenomicRegion("chr2", 8,  10, "X", 1, '+'));
  expect.push_back(GenomicRegion("chr3", 8,  10, "X", 1, '+'));
  expect.push_back(GenomicRegion("chr4", 0,  2,  "X", 2, '+'));

  vector<GenomicRegion> got;
  ReadBinner r(2);
  const bool UNSTRANDED = true;
  r.binReads(input, got, requiredBins, UNSTRANDED, 1);

  EXPECT_EQ(got.size(), expect.size());
  if (got.size(), expect.size()) {
    for (size_t i = 0; i < expect.size(); i++) {
      EXPECT_EQ(expect[i], got[i]);
    }
  }
}

/**
 * \brief Another test for the binning of reads -- this one exercises the
 *        splitting of chromosomes, which had a bug in an earlier version.
 */
TEST(ReadBinnerTest, testSplitByChromBug) {

  // responses -- define the required bins. Binsize is 50
  vector<GenomicRegion> input;
  input.push_back(GenomicRegion("chr2" "\t" "61719585"  "\t" "61719616"  "\t"
                                "X" "\t" "0" "\t" "+"));
  input.push_back(GenomicRegion("chr2" "\t" "61719701"  "\t" "61719719"  "\t"
                                "X" "\t" "0" "\t" "+"));
  input.push_back(GenomicRegion("chr5" "\t" "131927086" "\t" "131927098" "\t"
                                "X" "\t" "0" "\t" "-"));
  input.push_back(GenomicRegion("chr5" "\t" "131927568" "\t" "131927605" "\t"
                                "X" "\t" "0" "\t" "-"));
  input.push_back(GenomicRegion("chr6" "\t" "31603521"  "\t" "31603526"  "\t"
                                "X" "\t" "0" "\t" "-"));
  input.push_back(GenomicRegion("chr6" "\t" "31603743"  "\t" "31603787"  "\t"
                                "X" "\t" "0" "\t" "-"));
  input.push_back(GenomicRegion("chr3" "\t" "197677833" "\t" "197677849" "\t"
                                "X" "\t" "0" "\t" "+"));
  input.push_back(GenomicRegion("chr3" "\t" "197678029" "\t" "197678062" "\t"
                                "X" "\t" "0" "\t" "+"));
  input.push_back(GenomicRegion("chr1" "\t" "150939919" "\t" "150939960" "\t"
                                "X" "\t" "0" "\t" "-"));
  input.push_back(GenomicRegion("chr1" "\t" "150940139" "\t" "150940147" "\t"
                                "X" "\t" "0" "\t" "-"));

  // covars -- are the things we're binning
  vector<GenomicRegion> covar;
  covar.push_back(GenomicRegion("chr3"  "\t" "52730293"  "\t" "52730334"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chr3"  "\t" "52730557"  "\t" "52730565"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chr3"  "\t" "52562926"  "\t" "52562951"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chr3"  "\t" "52563165"  "\t" "52563189"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chr1"  "\t" "223933103" "\t" "223933141" "\t"
                                "X" "\t" "0" "\t" "-"));
  covar.push_back(GenomicRegion("chr1"  "\t" "223934698" "\t" "223934709" "\t"
                                "X" "\t" "0" "\t" "-"));
  covar.push_back(GenomicRegion("chr11" "\t" "17332788"  "\t" "17332830"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chr11" "\t" "17333418"  "\t" "17333425"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chrX"  "\t" "47085420"  "\t" "47085426"  "\t"
                                "X" "\t" "0" "\t" "+"));
  covar.push_back(GenomicRegion("chrX"  "\t" "47085677"  "\t" "47085720"  "\t"
                                "X" "\t" "0" "\t" "+"));

  // Expected required bins
  vector<GenomicRegion> expectReq;
  expectReq.push_back(GenomicRegion("chr1"  "\t" "150939900"  "\t" "150939950"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr1"  "\t" "150940100"  "\t" "150940150"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr2"  "\t" "61719550"   "\t" "61719600"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr2"  "\t" "61719700"   "\t" "61719750"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr3"  "\t" "197677800"  "\t" "197677850"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr3"  "\t" "197678000"  "\t" "197678050"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr5"  "\t" "131927050"  "\t" "131927100"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr5"  "\t" "131927550"  "\t" "131927600"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr6"  "\t" "31603500"   "\t" "31603550"
                                    "\t" "X"  "\t" "1"  "\t" "+"));
  expectReq.push_back(GenomicRegion("chr6"  "\t" "31603700"   "\t" "31603750"
                                    "\t" "X"  "\t" "1"  "\t" "+"));

  // Expected covariate binning results
  vector<GenomicRegion> expectCovarBinned;
  expectCovarBinned.push_back(GenomicRegion("chr1" "\t" "150939900" "\t"
                                     "150939950" "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr1" "\t" "150940100" "\t"
                                     "150940150" "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr1" "\t" "223933100" "\t"
                                     "223933150" "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr1" "\t" "223934650" "\t"
                                     "223934700" "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr11" "\t" "17332750" "\t"
                                     "17332800"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr11" "\t" "17333400" "\t"
                                     "17333450"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr2" "\t" "61719550"  "\t"
                                     "61719600"  "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr2" "\t" "61719700"  "\t"
                                     "61719750"  "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr3" "\t" "52562900"  "\t"
                                     "52562950"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr3" "\t" "52563150"  "\t"
                                     "52563200"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr3" "\t" "52730250"  "\t"
                                     "52730300"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr3" "\t" "52730550"  "\t"
                                     "52730600"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr3" "\t" "197677800" "\t"
                                     "197677850" "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr3" "\t" "197678000" "\t"
                                     "197678050" "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr5" "\t" "131927050" "\t"
                                     "131927100" "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr5" "\t" "131927550" "\t"
                                     "131927600" "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr6" "\t" "31603500"  "\t"
                                     "31603550"  "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chr6" "\t" "31603700"  "\t"
                                     "31603750"  "\t" "X"  "\t" "1"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chrX" "\t" "47085400"  "\t"
                                     "47085450"  "\t" "X"  "\t" "2"  "\t" "+"));
  expectCovarBinned.push_back(GenomicRegion("chrX" "\t" "47085650"  "\t"
                                     "47085700"  "\t" "X"  "\t" "2"  "\t" "+"));

  std::sort(input.begin(), input.end());
  std::sort(covar.begin(), covar.end());

  vector<GenomicRegion> gotBinnedResponses;
  ReadBinner r(50);
  const bool UNSTRANDED = true;
  r.binReads(input, gotBinnedResponses, UNSTRANDED, 0);
  EXPECT_EQ(gotBinnedResponses, expectReq);

  vector<GenomicRegion> gotCovarBinned;
  r.binReads(covar, gotCovarBinned, gotBinnedResponses, UNSTRANDED, 1);
  EXPECT_EQ(gotCovarBinned, expectCovarBinned);
}
