#!/usr/bin/python

""" 
  Date of Creation: 25rd Jun 2011     
                      
  Description:    Simulate input data for the peak caller using a 
                  distribution chosen by the user.

  Copyright (C) 2010  
  University of Southern California,
  Philip J. Uren,
  Andrew D. Smith
  
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
  
  TODO:          
                 None
"""

import random, sys, os, copy
from math import exp
from rpy2.robjects import r
from interface.cli import CLI, Option
from mapping.bed import BEDElement

DEFAULT_VERBOSITY = False
MAX_COV_VAL = 2
MAX_COEF_VAL = 1.5
MAX_MU_DISTRIBUTION = 1000
MAX_ALPHA = 10
DEFAULT_CHROMS = ["chr" + str(i) for i in range(0,22)] + ["chrX", "chrY"]
MAX_RANDOM_DISTANCE = 1000000
MAX_RANDOM_SIZE = 1000000


def getUI():
  programName = os.path.basename(sys.argv[0])
  longDescription =   "..."
  shortDescription =  "..."
  
  ui = CLI(programName, shortDescription, longDescription)
  ui.minArgs = 0
  ui.maxArgs = 1
  ui.addOption(Option(short = "r", long = "resOutput", argName = "filename",
                      description = "output resultant read counts in BED " +\
                                    "format to this file",
                      required = True, type = str))
  ui.addOption(Option(short = "c", long = "covOutput", argName = "filename",
                      description = "output resultant covariate values in " +\
                                    "BED format to these files; wrap " +\
                                    "multiple filenames in quotes",
                      required = False, type = str))
  ui.addOption(Option(short = "k", long="components", argName="int",
                      description="number of components to simulate",
                      required=False, type=int))
  ui.addOption(Option(short = "v", long = "verbose",
                      description = "output additional messages to stderr " +\
                                    "about run (default: " +\
                                    str(DEFAULT_VERBOSITY) + ")",
                      default = DEFAULT_VERBOSITY, required = False))
  ui.addOption(Option(short = "d", long = "distribution", argName = "type",
                      description = "what kind of distribution to use for " +\
                                    "generating read counts? Currently " +\
                                    "supports " +\
                                    "Poisson, " +\
                                    "NegativeBinomial, " +\
                                    "ZeroTruncatedPoisson, " +\
                                    "ZeroTruncatedNegativeBinomial, " +\
                                    "PoissonRegression, " +\
                                    "NegativeBinomialRegression, " +\
                                    "ZeroTruncatedPoissonRegression, " +\
                                    "ZeroTruncatedNegativeBinomialRegression",
                      required = False, type = str))
  ui.addOption(Option(short = "n", long = "num", argName = "int",
                      description = "how many data points to generate?",
                      required = True, type = int))
  ui.addOption(Option(short="h", long="help", 
                      description="show this help message ", special=True))
  ui.addOption(Option(short="u", long="test", 
                      description="run unit tests ", special=True))
  
  ui.parseCommandLine(sys.argv[1:])
  return ui


def _main():
  # get options and arguments 
  ui = getUI()
  
  # just run unit tests
  if ui.optionIsSet("test") :
    unittest.main(argv = [sys.argv[0]])
    sys.exit()
  
  # just show help
  if ui.optionIsSet("help") :
    ui.usage()
    sys.exit()
  verbose = (ui.optionIsSet("verbose") == True) or DEFAULT_VERBOSITY
  
  # get output file handles. We know the result output filename must be 
  # provided because we marked it as required above
  resOutfh = open(ui.getValue("resOutput"), "w")
  covOutfhs = None
  p = 0
  if ui.optionIsSet("covOutput") :
    covOutfhs = [open(fn, "w") for fn in ui.getValue("covOutput").split()]
    p = len(covOutfhs)
  
  # how much data to make? Again, we know this will be set, since we made it
  # a required option
  n = ui.getValue("num")
  if n <= 0 :
    sys.stderr.write("must provide a positive number of elements to " +\
                     "produce. Found " + str(n))
    sys.exit()
  
  # what type of distribution to use?
  type = "Poisson"
  if ui.optionIsSet("distribution") : type = ui.getValue("distribution")
  
  # how many components to simulate?
  numComponents = 1
  if ui.optionIsSet("components") : numComponents = ui.getValue("components")
  if numComponents <= 0 :
    sys.stderr.write("must provide a positive number of components to " +\
                     "produce. Found " + str(numComponents))
    sys.exit()
  perComponent = sumRandom(numComponents, n)
  if verbose :
    sys.stderr.write("mixing: " +\
                     str([i/float(n) for i in perComponent]) + "\n")
  
  # make the data..
  assert(numComponents >= 1)
  coefficients = [[random.random() * MAX_COEF_VAL for i in range(0,p)] 
                  for j in range(0,numComponents)]
  if coefficients[0] == [] : coefficients = None
  regions = simulateRandomRegions(n)
  resElements = regions
  covElements = [[copy.copy(regions[i]) for j in range(0,p)] 
                                        for i in range(0,n)]
  if covElements[0] == [] : covElements = None
  for i in range(0,numComponents) :
    s,e = sum(perComponent[0:i]), sum(perComponent[0:i+1])
    cElems, coeffs = None, None
    if covElements != None : cElems = covElements[s:e]
    if coefficients != None : coeffs = coefficients[i] 
    d = simulateScores(resElements[s:e], type, cElems, 
                          coeffs, verbose)
  
  # write the data
  for i in range(0, n) :
    resOutfh.write(str(resElements[i]) + "\n")
    for j in range(0, p) :
      covOutfhs[j].write(str(covElements[i][j]) + "\n")
      
  # tell the user what we did, if they want.. 
  if verbose : 
    sys.stderr.write("Generated " + str(n) + " simulated data points \n")
    if d > 0 : sys.stderr.write("truncated " + str(d) + " elements \n")



def sumRandom(n, total):
  """
    @summary: Return a randomly chosen list of n positive integers 
              summing to total.
  """
  dividers = sorted(random.sample(xrange(1, total), n - 1))
  return [a - b for a, b in zip(dividers + [total], [0] + dividers)]


def simulateRandomRegions(num, chroms=DEFAULT_CHROMS):
  """
    @summary: generate <num> BED elements, no overlap. 
  """
  perChrom = [divmod(num,len(chroms))[0] if i < len(chroms) - 1 else 
              divmod(num,len(chroms))[0] + divmod(num,len(chroms))[1] 
              for i in range(0,len(chroms))]
  res = []
  for i in range(0, len(chroms)) :
    prev = 0
    for j in range(0, perChrom[i]) :
      s = prev + 1 + random.random() * (MAX_RANDOM_DISTANCE - 1)
      e = s + 1 + random.random() * (MAX_RANDOM_SIZE - 1)
      res.append(BEDElement(chroms[i], s, e, "X", 0, "+"))
      prev = e
  return res
  
    
def simulateScores(responseElements, distributionType, 
                   covariateElements=None, coefficients=None, verbose=False):
  """
    @summary: TODO
    @param covariateElements: an n x p matrix of BED elements
    @param responseElements: a vector of n BED elements
    @param distributionType: 
  """
  n = len(responseElements)
  assert((n > 0) and \
         (covariateElements == None or len(covariateElements) == n))
  p = 0 if covariateElements == None else len(covariateElements[0])
  alpha = None
  totalDropped = 0
  
  # complain if we're asked to build something that can't use covariates, but
  # we've been given covariates in some fashion. 
  if distributionType.lower().strip() in ["poisson", "negativebinomial",
                                  "zerotruncatedpoisson", 
                                  "zerotruncatednegativebinomial"] and \
     (covariateElements != None or coefficients != None):
    #TODO throw exception here instead, so we can cleanup in main
    sys.stderr.write("Distribution type " + str(distributionType) + " does " +\
                     "not support creation with covariates. Did you forget " +\
                     "to specifying that it is a regression?")
    sys.exit()
  
  # complain if we're asked to build something that needs covariates, but we 
  # haven't been given any
  if distributionType.lower().strip() in \
                          ["poissonregression", 
                           "negativebinomialregression",
                           "zerotruncatedpoissonregression", 
                           "zerotruncatednegativebinomialregression"] and \
                          (covariateElements == None or coefficients == None) :
    #TODO throw exception here instead, so we can cleanup in main
    sys.stderr.write("Distribution type " + str(distributionType) + " needs" +\
                     " covariates for creation, but not are provided. " +\
                     "Did you forget to list covariate output filenames?")
    sys.exit()
    
  # tell the user what regression coefficients we're using?
  if verbose and coefficients != None :
    sys.stderr.write("Regression coefficients used: " +\
                     ", ".join([str(x) for x in coefficients]) + "\n")
    
  # we might need to make up mu and alpha, we do it here regardless and just
  # won't use them below if they're not needed
  mu = random.random() * MAX_MU_DISTRIBUTION
  alphad = random.random() * MAX_ALPHA
  if distributionType.lower().strip() in ["poisson", "negativebinomial",
                                          "zerotruncatedpoisson", 
                                          "zerotruncatednegativebinomial"] :
    sys.stdout.write("mu = " + str(mu) + "\n")
  if distributionType.lower().strip() in ["negativebinomial",
                                          "zerotruncatednegativebinomial"] :
    sys.stdout.write("alpha = " + str(alphad) + "\n")
  
    
  for i in range(0, n) :
    # get the ith set of covariate values, if we have them
    covars = None
    if covariateElements != None :
      for e in covariateElements[i] : e.score = random.random() * MAX_COV_VAL
      covars = [e.score for e in covariateElements[i]]
      
    ### Poisson regression 
    if distributionType.lower().strip() == "poissonregression" :
      s = simulateReadCountPoissonReg(covars, coefficients, verbose=verbose)
      
    ### NB regression 
    elif distributionType.lower().strip() == "negativebinomialregression" :
      if alpha == None : 
        alpha = random.random() * MAX_ALPHA  
        sys.stderr.write("alpha = " + str(alpha) + "\n")
      s = simulateReadCountNBReg(covars, coefficients, alpha, verbose=verbose)
    
    ### ZTP regression  
    elif distributionType.lower().strip() == "zerotruncatedpoissonregression" :
      s, dropped = simulateReadCountTruncatedPoissonReg(covars, coefficients, 
                                                        verbose=verbose)
      totalDropped += dropped
    
    ### ZTNB regression
    elif distributionType.lower().strip() == \
                                    "zerotruncatednegativebinomialregression" :
      if alpha == None : 
        alpha = random.random() * MAX_ALPHA
        sys.stderr.write("alpha = " + str(alpha) + "\n")
      s, dropped = simulateReadCountTruncatedNBReg(covars, coefficients, 
                                                   alpha, verbose=verbose)
      totalDropped += dropped
    
    ### regular poisson distribution ###
    elif distributionType.lower().strip() == "poisson" :
      s = simulateReadCountPoisson(mu)
    
    ### regular NB distribution 
    elif distributionType.lower().strip() == "negativebinomial" :
      s = simulateReadCountNB(mu, alphad)
      
    ### regular zero-truncated poisson distribution ###
    elif distributionType.lower().strip() == "zerotruncatedpoisson" :
      s, dropped = simulateReadCountTruncatedPoisson(mu)
      totalDropped += dropped
    
    ### regular zero-truncated negative binomial ###
    elif distributionType.lower().strip() == "zerotruncatednegativebinomial" :
      s, dropped = simulateReadCountTruncatedNB(mu, alphad)
      totalDropped += dropped 
      
    ### Oops, don't know this one.. 
    else :
      sys.stderr.write("unknown type: " + str(distributionType))
      sys.stderr.write("; quitting\n")
      sys.exit()
    responseElements[i].score = s
  return totalDropped
    
def simulateReadCountPoissonReg(covariates, regressionParams, verbose=False):
  """
    @summary: TODO  
  """
  assert(len(covariates) == len(regressionParams))
  p = len(covariates)
  mu = exp(sum([covariates[j] * regressionParams[j] for j in range(0,p)]))
  return r("rpois(1, lambda=" + str(mu) + ")")[0]

def simulateReadCountNBReg(covariates, regressionParams, alpha, 
                           verbose=False):
  """
    @summary: TODO 
  """ 
  assert(len(covariates) == len(regressionParams))
  p = len(covariates)
  mu = exp(sum([covariates[j] * regressionParams[j] for j in range(0,p)]))
  return r("rnbinom(1, mu=" + str(mu) + ", size=" + str(1/alpha) + ")")[0]    

def simulateReadCountTruncatedPoissonReg(covariates, regressionParams, 
                                         truncateAt=0, verbose=False):
  """
    @summary: TODO 
  """ 
  res = None
  dropped = -1
  while res == None or res <= truncateAt :
    dropped += 1 
    res = simulateReadCountPoissonReg(covariates, regressionParams, verbose)
  return res, dropped

def simulateReadCountTruncatedNBReg(covariates, regressionParams, alpha,
                                    truncateAt=0, verbose=False):
  """
    @summary: TODO 
  """ 
  res = None
  dropped = -1
  while res == None or res <= truncateAt :
    dropped += 1 
    res = simulateReadCountNBReg(covariates, regressionParams, alpha, verbose)
  return res, dropped

def simulateReadCountPoisson(mu, verbose=False):
  """
    @summary: TODO 
  """ 
  return r("rpois(1, lambda=" + str(mu) + ")")[0]

def simulateReadCountNB(mu, alpha, verbose=False):
  """
    @summary: TODO
  """
  return r("rnbinom(1, mu=" + str(mu) + ", size=" + str(1/alpha) + ")")[0]

def simulateReadCountTruncatedPoisson(mu, truncateAt=0, verbose=False):
  """
    @summary: TODO
    
  """
  res = None
  dropped = -1
  while res == None or res <= truncateAt :
    dropped += 1 
    res = simulateReadCountPoisson(mu)
  return res, dropped

def simulateReadCountTruncatedNB(mu, alpha, truncateAt=0, verbose=False):
  """
    @summary: TODO
  """
  res = None
  dropped = -1
  while res == None or res <= truncateAt :
    dropped += 1 
    res = simulateReadCountNB(mu, alpha, verbose)
  return res, dropped


if __name__ == "__main__":
  _main()
  
  