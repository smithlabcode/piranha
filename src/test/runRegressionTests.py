import sys, os, shlex, re, math
import subprocess
import tempfile
import StringIO
from xml.dom.minidom import parseString

class TextColours:
  OK_GREEN = '\033[92m'
  WARNING_YELLOW = '\033[93m'
  FAIL_RED = '\033[91m'
  ENDC = '\033[0m'

class Test :
  def __init__(self, name, command, expected) :
    self.name = str(name).strip()
    self.command = str(command).strip()
    self.expected = str(expected).strip()
  def __str__(self) :
    return "name=" + self.name + " command=" + self.command + " expect=" + self.expected
  
def compareStrings (s1, s2, absTol=None, relTol=None, logfile=None): 
  numpat = re.compile(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
  eps = 0.000001
  
  ## first we're going to compare the numbers ##
  numbers_s1 = numpat.findall(s1)
  numbers_s2 = numpat.findall(s2)
  if len(numbers_s1) != len(numbers_s2) :
    logfile.write("count of numbers doesn't match\n")
    logfile.write("expected numbers: " + ", ".join([str(x) for x in numbers_s2]) + "\n")
    logfile.write("got numbers: " + ", ".join([str(x) for x in numbers_s1]) + "\n")
    return False
  for i in range(0, len(numbers_s1)) :
    # for now we'll just treat everything as floats
    n1, n2 = float(numbers_s1[i]), float(numbers_s2[i])
    
    # do we have an absolute allowed tolerance? 
    if absTol != None :
      if math.fabs(n1 - n2) > absTol : 
        logfile.write("comparing " + str(n1) + " with " + str(n2) +\
                      " failed to satisfy absolute tolerance limit of " +\
                      str(absTol) + "\n")
        return False
    # what about a relative tolerance?
    elif relTol != None : 
      if math.fabs(n1 - n2) / (n2 + eps) > relTol :
        logfile.write("comparing " + str(n1) + " with " + str(n2) +\
                      " failed to satisfy relative tolerance limit of " +\
                      str(relTol) + "\n")
        return False
    # if nothing else, require exact match
    elif n1 != n2 : 
      return False
      logfile.write("comparing " + str(n1) + " with " + str(n2) +\
                    " failed to satisfy absolute equality requirement \n")
    
  ## next we're going to compare the strings
  strings_s1 = numpat.split(s1)
  strings_s2 = numpat.split(s2)
  if len(strings_s1) != len(strings_s2) :
    logfile.write("count of strings doesn't match\n")
    logfile.write("expected strings: " + ", ".join([str(x) for x in numbers_s2]) + "\n")
    logfile.write("got strings: " + ", ".join([str(x) for x in numbers_s2]) + "\n")
    return False
  for i in range(0, len(strings_s1)) :
    str1, str2 = strings_s1[i], strings_s2[i]
    if str1.strip() != str2.strip() :
      logfile.write("comparing " + str(str1.strip()) + " with " +\
                    str(str2.strip()) + " failed.\n") 
      return False
  
  # everything matched up
  return True
    

def loadTests(fn) : 
  data = open(fn).read() 
  dom = parseString(data)

  res = []
  tests=dom.getElementsByTagName("regressionTests")[0].getElementsByTagName('test')
  for test in tests :
    name = test.getElementsByTagName("name")[0].firstChild.data
    command = test.getElementsByTagName("command")[0].firstChild.data
    expected = test.getElementsByTagName("expected")[0].firstChild.data
    res.append(Test(name, command, expected))
  return res

def runTest(test, colWidth, logfile=None) :
  print "Running '" + test.name + "'... ",
  for i in range(0, colWidth - len(test.name)): print "",

  # we make a temporary file to store STDOUT, STDERR from the test
  tmpf, tmpnm = tempfile.mkstemp(dir=os.getcwd())
  tmpf2, tmpnm2 = tempfile.mkstemp(dir=os.getcwd());

  try :
    # run the test..
    print TextColours.WARNING_YELLOW + "[exec..]" + TextColours.ENDC,
    sys.stdout.flush()
    args = shlex.split(test.command)
    subprocess.Popen(args, stdout=tmpf, stderr=tmpf2)
    os.wait()

    actualOutput = open(tmpnm).read()
    expectedOutput = open(test.expected).read()
    logbuffer = StringIO.StringIO()
    print "\b\b\b\b\b\b\b\b\b\b",
    if not (compareStrings(actualOutput, expectedOutput, relTol=0.05, logfile=logbuffer)) :
      print TextColours.FAIL_RED + "[failed]" + TextColours.ENDC
      if logfile != None :
        logfile.write("******************* FAILED TEST *******************\n") 
        logfile.write("testname = "+ str(test.name) + "\n\n")
        logfile.write("------------------ ACTUAL STDOUT ------------------\n")
        logfile.write(open(tmpnm).read())
        logfile.write("\n")
        logfile.write("------------------ ACTUAL STDERR ------------------\n")
        logfile.write(open(tmpnm2).read())
        logfile.write("\n")
        logfile.write("----------------- EXPECTED STDOUT -----------------\n")
        logfile.write(open(test.expected).read())
        logfile.write("\n")
        logfile.write("-------------- DETAILED FAILURE NOTE --------------\n")
        logfile.write(logbuffer.getvalue())
        logfile.write("\n")
    else :
      print TextColours.OK_GREEN + "[passed]" + TextColours.ENDC
  finally :
    # cleanup our temporary files 
    os.close(tmpf)
    os.remove(tmpnm)
    os.close(tmpf2)
    os.remove(tmpnm2)

def runTests(testfn):
  tests = loadTests(testfn)
  colWidth = max([len(t.name) for t in tests])
  logbuffer = StringIO.StringIO()
  for test in tests :
    runTest(test, colWidth, logbuffer)
  logoutput = logbuffer.getvalue()
  if logoutput.strip() != "" :
    logfile = open("testlog.txt", "w")
    sys.stderr.write("\n")
    sys.stderr.write("Some regression tests failed.\n")
    sys.stderr.write("A log was created in src/test/testlog.txt\n")
    logfile.write(logoutput)

runTests(sys.argv[1])