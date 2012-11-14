#!/usr/bin/python

""" 
  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Andrew D. Smith
  
  Authors: Philip J. Uren, Andrew D. Smith
  
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
  
  Description:   Given a file of mapped reads and a set of chromosome files in
                 fasta format, for each location in the genome with at least
                 one read starting, determine the number of reads in the 
                 cluster that have T -> C conversions
"""

import os, sys, unittest
import collections
from optparse import OptionParser

class MappedRead :
  def __init__(self, s):
    chrom, start, end, name, score, strand, sequence, quality = s.split("\t")
    self._chrom = chrom.strip()
    self._start = int(start)
    self._end = int(end)
    self._name = name.strip()
    self._score = float(score)
    self._strand = strand.strip()
    self._sequence = sequence.strip()
    self._quality = quality.strip()
    self._containsConversion = False
  
  def __str__(self):
    return self._chrom + "\t" + str(self._start) + "\t" + str(self._end) + "\t" +\
           self._name + "\t" + str(self._score) + "\t" + str(self._strand) + "\t" +\
           str(self._sequence) + "\t" + str(self._quality)
  
  def __getitem__(self, index):
    return self._sequence[index]
  
  def getReadLength(self):
    return self._end - self._start
  
  def getChrom(self):
    return self._chrom
  
  def getStart(self):
    return self._start
  
  def getEnd(self):
    return self._end
  
  def containsConversion(self):
    return self._containsConversion
  
  def setContainsConversion(self, v):
    self._containsConversion = v
    

class MappedReadIter :
  def __init__(self, mrIter):
    self.finished = False
    self.nextRead = None
    self.mrIter = mrIter
    try:
      self.nextRead = MappedRead(self.mrIter.next())
    except StopIteration:
      self.finished = True

  def __iter__(self):
    return self
  
  def peek(self):
    return self.nextRead
  
  def exhausted(self):
    return self.finished

  def next(self):
    if self.finished:
      raise StopIteration()
    res = self.nextRead
    try:
      self.nextRead = MappedRead(self.mrIter.next())
    except StopIteration:
      self.nextRead = None
      self.finished = True
    return res
  
def chromIter(chromfh):
  foundHeader = False
  chromName = None 
  pos = 0
  for line in chromfh :
    line = line.strip()
    if line == "" : continue 
    if line[0] == ">" : 
      if foundHeader == True : 
        sys.stderr.write("found multiple sequences, need one per chromosome")
      foundHeader = True
      chromName = line[1:].strip()
    else :
      for nuc in line : 
        pos += 1
        yield chromName, nuc, pos
  

def processChrom(cIter, readIter, ofh, readLength):
  debug = False
  openAnswers = collections.deque()
  openReads = collections.deque()
  for chrom, nuc, pos in cIter :
    # output any answers that are ready...
    while len(openAnswers) > 0 and openAnswers[0] + readLength <= pos :
      score = sum([1 for r in openReads if r.containsConversion()]) / float(len(openReads))
      ofh.write(chrom + "\t" + str(openAnswers[0]) + "\t" +\
                str(openAnswers[0]+1) + "\t" + "X" + "\t" +\
                str(score) + "\t" + "+" + "\n")
      openAnswers.popleft()
    
    # add any new reads that start at this position to the open list
    
    while (not readIter.exhausted()) and \
          (readIter.peek().getStart() == pos) and \
          (readIter.peek().getChrom() == chrom) :
      prevStart = readIter.peek().getStart()
      openReads.append(readIter.next())
      if (not readIter.exhausted()) and \
         (readIter.peek().getChrom() == chrom) and \
         (prevStart > readIter.peek().getStart()) :
        sys.stderr.write("ERROR: input reads not sorted by start\n")
        sys.exit()
      if len(openAnswers) == 0 or openAnswers[-1] != pos : openAnswers.append(pos)
      
    # update any open reads that have T -> C conversion
    if nuc == "t" or nuc == "T" :
      for r in openReads :
        if pos >= r.getEnd() : continue 
        if r[pos - r.getStart()] == "c" or r[pos  - r.getStart()] == "C" : 
          r.setContainsConversion(True)
    
    # drop reads that are no longer affecting any open answers
    while (len(openReads) > 0) and \
          (len(openAnswers) == 0 or openReads[0].getEnd() < openAnswers[0]): 
      openReads.popleft()
    
    if debug :
      sys.stderr.write(str(nuc) + "\t" + str(pos) + "\n")
      sys.stderr.write("\ttop read: " + str(readIter.peek()) + "\n")
      sys.stderr.write("\topen reads\n") 
      for r in openReads : sys.stderr.write("\t\t" + str(r) + "\n")
      sys.stderr.write("\topen answers\n") 
      for r in openAnswers : sys.stderr.write("\t\t" + str(r) + "\n")
      
    # we can stop if there are no reads in the open list and no more on this 
    # chromosome (or no more at all)
    if (len(openReads) == 0) and \
       (readIter.exhausted() or readIter.peek().getChrom() != chrom) : 
      break
      
      
def processChroms(chromIterators, mappedReadIterator, ofh):
  rl = mappedReadIterator.peek().getReadLength()
  prevChrom = None
  for cIter in chromIterators :
    processChrom(cIter, mappedReadIterator, ofh, rl)
    if not mappedReadIterator.exhausted() :
      if prevChrom != None and prevChrom > mappedReadIterator.peek().getChrom() :
        sys.stderr.write("ERROR: input reads not sorted by chromosome\n") 
      prevChrom = mappedReadIterator.peek().getChrom() 
      
      
def loadChroms(chromsDir): 
  return [open(x) for x in os.listdir(chromsDir)]
      
      
def main():
  parser = OptionParser()
  parser.add_option("-c", "--chroms", dest="chromsDir",
                  help="load chromosome fasta files from DIR", metavar="DIR")
  parser.add_option("-t", "--test", dest="test", default=False, action="store_true",
                  help="run unit tests")

  options, args = parser.parse_args()
  if options.test :
    unittest.main(argv = [sys.argv[0]])
    sys.exit()
  if len(args) == 0:
    parser.print_help()
    sys.exit()
  if options.chromsDir == None : 
    sys.stderr.write("must provide chroms directory\n")
    sys.exit()
  
  chromsIters = loadChroms(options.chromsDir)
  readIter = MappedReadIter(open(args[1]))
  
  processChroms(chromsIters, readIter, sys.stdout)
  
  

    
###############################################################################
#                               UNIT TESTING                                  #
###############################################################################
class DummyInputStream :
  def __init__(self, lines, name = None):
    if type(lines).__name__ == "str" :
      lines = lines.split("\n")
    lines = [x + "\n" for x in lines]
    self.lines = lines
    self.current = 0
    self.length = len(self.lines)
    
    self.name = "none"
    if name != None : self.filename = name
  
  def readline(self):
    if self.current >= self.length : 
      return ""
    self.current += 1 
    return self.lines[self.current-1]
  
  def reset(self):
    self.current = 0
    
  def tell(self):
    return self.current
  
  def seek(self, s):
    self.current = s
  
  def __iter__(self):
    #return self.lines.__iter__()
    return self
  
  def next(self):
    self.current += 1
    if self.current > (len(self.lines) - 1) : raise StopIteration 
    return self.lines[self.current-1]
  
  def __eq__(self, o):
    if o == None : return False
    if o == sys.stdin: return True
    return False
  def __ne__(self, o):
    if o == None : return True
    if o == sys.stdin: return False
    return True
  
  def close(self):
    pass
  
  
class DummyOutputStream :
  def __init__(self):
    self.stored = []
    self.prev = ""
    
  def write(self, sth):
    sth = self.prev + sth
    if sth.find("\n") != -1 :
      lines = str(sth).split("\n")
      for line in lines :
        if line == "" : continue 
        self.stored.append(line + "\n")
      self.prev = ""
    else :
      self.prev = sth
      
  def close(self):
    if self.prev != "" :
      self.stored.append(self.prev)
    
  def itemsWritten(self):
    return self.stored
  
  def numItemsWritten(self):
    return len(self.itemsWritten())  
  
  def __str__(self):
    return "\n".join(self.itemsWritten())


class PseudoCountTests(unittest.TestCase):
  """
    Unit tests for paired end filter 
  """
  
  def test_processChroms(self):
    debug = False
    chrom1 = ">chr1\n" +\
             "AAAAAATATAAAAAAAAAAAAAA\n" 
    chrom2 = ">chr2\n" +\
             "AAAAAATATAAAAAAAAAAAAAA\n"
    mappedReads =\
      "chr1"+"\t"+"3"+"\t"+"9" +"\t"+"R1"+"\t"+"0"+"\t"+"+"+"\t"+"AAAACA"+"\t"+"BBBBBB\n" +\
      "chr1"+"\t"+"4"+"\t"+"10"+"\t"+"R2"+"\t"+"0"+"\t"+"+"+"\t"+"AAACAC"+"\t"+"BBBBBB\n" +\
      "chr1"+"\t"+"6"+"\t"+"12"+"\t"+"R3"+"\t"+"0"+"\t"+"+"+"\t"+"ACATAA"+"\t"+"BBBBBB\n" +\
      "chr1"+"\t"+"9"+"\t"+"15"+"\t"+"R4"+"\t"+"0"+"\t"+"+"+"\t"+"CAAAAA"+"\t"+"BBBBBB\n" +\
      "chr1"+"\t"+"9"+"\t"+"15"+"\t"+"R5"+"\t"+"0"+"\t"+"+"+"\t"+"TAAATA"+"\t"+"BBBBBB\n" +\
      "chr2"+"\t"+"3"+"\t"+"9" +"\t"+"R21"+"\t"+"0"+"\t"+"+"+"\t"+"AAAACA"+"\t"+"BBBBBB\n" +\
      "chr2"+"\t"+"4"+"\t"+"10"+"\t"+"R22"+"\t"+"0"+"\t"+"+"+"\t"+"AAACAC"+"\t"+"BBBBBB\n" +\
      "chr2"+"\t"+"6"+"\t"+"12"+"\t"+"R23"+"\t"+"0"+"\t"+"+"+"\t"+"ACATAA"+"\t"+"BBBBBB\n" +\
      "chr2"+"\t"+"9"+"\t"+"15"+"\t"+"R24"+"\t"+"0"+"\t"+"+"+"\t"+"CAAAAA"+"\t"+"BBBBBB\n" +\
      "chr2"+"\t"+"9"+"\t"+"15"+"\t"+"R25"+"\t"+"0"+"\t"+"+"+"\t"+"TAAATA"+"\t"+"BBBBBB\n" 

    answer =\
      ["chr1"+"\t"+"3"+"\t"+"4" +"\t"+"X"+"\t"+"1.0"+"\t"+"+",
       "chr1"+"\t"+"4"+"\t"+"5" +"\t"+"X"+"\t"+"0.8"+"\t"+"+",
       "chr1"+"\t"+"6"+"\t"+"7" +"\t"+"X"+"\t"+"0.8"+"\t"+"+",
       "chr1"+"\t"+"9"+"\t"+"10"+"\t"+"X"+"\t"+"0.8"+"\t"+"+",
       "chr2"+"\t"+"3"+"\t"+"4" +"\t"+"X"+"\t"+"1.0"+"\t"+"+",
       "chr2"+"\t"+"4"+"\t"+"5" +"\t"+"X"+"\t"+"0.8"+"\t"+"+",
       "chr2"+"\t"+"6"+"\t"+"7" +"\t"+"X"+"\t"+"0.8"+"\t"+"+",
       "chr2"+"\t"+"9"+"\t"+"10"+"\t"+"X"+"\t"+"0.8"+"\t"+"+"] 
      
    chrom1 = chromIter(DummyInputStream(chrom1))
    chrom2 = chromIter(DummyInputStream(chrom2))
    readsfh = DummyInputStream(mappedReads)
    ofs = DummyOutputStream()
    processChroms([chrom1, chrom2], MappedReadIter(readsfh), ofs)
    gotOutput = ofs.itemsWritten()
    
    self.assertTrue(len(gotOutput) == len(answer))
    for i in range(0, len(answer)) :
      if debug : sys.stderr.write(answer[i] + "\t-->\t" + gotOutput[i] + "\n")
      self.assertTrue(answer[i].strip() == gotOutput[i].strip())

if __name__ == "__main__":
  main()