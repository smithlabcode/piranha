
                    _____  _                 _            
                   |  __ \(_)               | |           
                   | |__) |_ _ __ __ _ _ __ | |__   __ _  
                   |  ___/| | '__/ _` | '_ \| '_ \ / _` | 
                   | |    | | | | (_| | | | | | | | (_| | 
                   |_|    |_|_|  \__,_|_| |_|_| |_|\__,_|
                   **************************************
                   *               V1.2.1               *
                   ************************************** 
                                        

*********************************
Copyright and License Information
*******************************************************************************
Copyright (C) 2012
University of Southern California,
Philip J. Uren, Andrew D. Smith
  
Authors: Philip J. Uren, Emad Bahrami-Samani, Andrew D. Smith
  
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

This software package contains Google Test and BAMTools -- see their 
respective directories for copyright and license information for those 
packages. 


*****************
Table of contents
*******************************************************************************
1. Building and Installing Piranha 
2. Basic usage of Piranha 
3. Command-line options of Piranha 
4. Input and output formats of Piranha
5. Simulation of input data for Piranha
6. Contacts and bug reports


**********************************
1. Building and Installing Piranha
*******************************************************************************
Piranha has been designed to operate in a UNIX-like environment.
It has been tested on MacOS X Snow Leopard, USCLinux, Ubuntu and SUSE. We've
made every effort to ensure it's compatible with a wide range of systems, but
we can't test on all possible platforms and configurations. If you find an
incompatibility when building or using it, please let us know!

Step 0 -- dependencies and requirements
---------------------------------------
  64-bit machine and GCC version > 4.1 (to support TR1). 
  
  Piranha additionally requires a functioning installation of the GNU 
  Scientific Library (GSL). If you don't already have this installed, you will
  need to download and install it from http://www.gnu.org/software/gsl/
  This release has been tested with GSL version 1.14
  Your installation of GSL must be built for a 64 bit architecture.
  
  [OPTIONAL] Tests
    The regression tests (which you don't need to run, but can to test that 
    Piranha was built correctly) require Python. You can download Python 
    from http://www.python.org/
    This release has been tested with Python 2.6.1 
  
  [OPTIONAL] BAM support
    If you want BAM support, you must download and build BAMTools.
    If Piranha cannot find the BAMTools development headers and library,
    it will be compiled without support for BAM. See below for details on
    how to specify the location of BAMTools during the build process. 
    BAMTools is available from https://github.com/pezmaster31/bamtools
    This release has been tested with BAMTools version 2.1.1   

Step 1 -- configuring
---------------------
  First configure the installation. To do this, where '>' is your prompt and 
  the CWD is the root of the distribution, type:
  
  > ./configure 
  
  [OPTIONAL] BAM support
    configure will warn you if BAM support is not available. If you don't want 
    BAM support, this is fine. If you do want it however, and the configure 
    script did not find BAMTools, you will have to specify where the BAMTools 
    development headers and libraries are located. You can do that as follows
  
    > ./configure --with-bam_tools_headers="/some/path/BAMTools/include/" \
      --with-bam_tools_library="/some/path/BAMTools/lib/"
  
    the command is split into two lines here for aesthetic reasons, but it 
    need not be so. There are a few caveats:
      * If you subsequently move the BAMTools library, Piranha will not work.
        (this should be obvious really).
      * The BAMTools installer currently doesn't update the linker info 
        after installation. That means you might install it, but Piranha will
        still not find it when running ./configure -- in this case, I suggest
        manually specifying the location as shown above.

Step 2 -- building
------------------
  
  To build the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution  
  
  > make all 
  
Step 3 -- installing
--------------------
  To install the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution
  
  > make install
  
  This will place the binaries in the bin directory under the package root.
  They can be used directly from there without any additional steps. You can
  add that directory to your PATH environment variable to avoid having to 
  specify their full paths, or you can copy the binaries to another directory
  of your choice in your PATH.
  
Step 4 -- testing [OPTIONAL]
----------------------------
  You can verify that Piranha was built correctly by running the included
  unit and regression tests. At the command prompt (assuming your prompt is 
  '>') type the following:
  
  > make test 
  
  
*************************
2. Basic usage of Piranha 
*******************************************************************************
Piranha has two main modes of operation. In the first, a single regular 
distribution is fit to the input data and each region is assigned a p-value
based on this distribution. The default distribution is the zero-truncated 
negative binomial. To run Piranha in this way, use the following, where > is 
your prompt  

> ./Piranha input.bed

The second mode of operation fits a regression model relating counts to 
covariates. The default is a zero-truncated negative binomial regression. To
run Piranha in this way, do the following, where > is your prompt 

> ./Piranha input.bed covariate1.bed covariate2.bed

Note that the first file is always assumed to be the counts, and the remaining
files are assumed to be covariates. 

For further options, run Piranha without any arguments, or see the section 
'Command Line Options' below.


**********************************
3. Command-line options of Piranha 
*******************************************************************************
Usage: Piranha [OPTIONS] *.bed

Options:
  -o, -output             Name of output file, STDOUT if omitted 
  -s, -sort               indicates that input is unsorted and Piranha should 
                          sort it for you
  -p, -p_threshold        significance threshold for sites 
  -a, -background_thresh  indicates that this proportion of the lowest scores 
                          should be considered the background. Default is 0.99
  -b, -bin_size           indicates that input is raw reads and should be binned 
                          into bins of this size 
  -i, -bin_size_covars    indicates that the covariates (all except first 
                          file) are raw reads and should be binned into bins of 
                          this size 
  -z, -bin_size_both      synonymous with -b x -i x for any x 
  -u, -cluster_dist       merge significant bins within this distance. 
                          Setting to 0 disables merging, default is 1 (merge 
                          adjacent) 
  -r, -suppress_covars    don't print covariate values in output
  -f, -fit                Fit only, output model to file 
  -d, -dist               Distribution type. Currently supports Poisson, 
                          NegativeBinomial, ZeroTruncatedPoisson, 
                          ZeroTruncatedNegativeBinomial, 
                          PoissonRegression, NegativeBinomialRegression, 
                          ZeroTruncatedPoissonRegression, 
                          ZeroTruncatedNegativeBinomialRegression 
  -t, -fitMethod          component fitting method 
  -m, -model              Use the specified model file instead of fitting to 
                          input data 
  -v, -VERBOSE            output additional messages about run to stderr if set 
  -x, -unstranded         Don't preserve strand (puts all the peaks in positive 
                          strand) 
  -n, -no_normalisation   don't normalise covariates 
  -l, -log_covars         convert covariates to log scale 

Help options:
  -?, -help               print this help message 
      -about              print about message 
 


**************************************
4. Input and output formats of Piranha
*******************************************************************************
    NOTE: BAM support is only available if Piranha was linked with BAMTools
    =======================================================================
 
Piranha takes three possible inputs:
  1.) The response file   (BED or BAM format) [REQUIRED]
  2.) The covariate files (BED format)        [OPTIONAL]      
  3.) The model file      (XML format)        [OPTIONAL]
  
The response file (always the first argument, and required) may contain:
  1.) Binned read counts. In this case the input can be provided in BED 
      format, but not BAM format. The Score field of the BED format file will 
      be treated as the count of the number of reads mapping into that bin. 
      This is the default and is what is expected if no contrary options are 
      specified.
  2.) Raw read locations, in which case Piranha will bin the reads for you. 
      If you provide your input like this, you MUST set the bin size option, 
      or Piranha will treat your input as being already binned. If Piranha is 
      creating bins from raw reads, it will always start at the first index of 
      each chromosome (i.e. index 1) and move at a step size equal to the bin 
      size. Bins with no reads in them will not be retained. 
The file format will be determined by looking at the file extension (.bam or 
.bed). Case is not important; unknown extensions are treated as .bed files. 

The covariates files are optional. You can provide one or more of them.
If none are provided, the program fits a regular distribution (zero-truncated 
negative binomial by default) to the read counts. If covariates are provided,
the program performs regression (zero truncated negative binomial regression
by default) using the provided covariates. Covariates must be provided in
BED format, where the score field (the fifth) is taken as value of the 
covariate for the genomic region defined by the other fields. Every covariate
file must contain genomic loci that match the response file.

Piranha has a special option for handling the case where a single covariate
is present and this covariate is the number of reads from another sequencing
run (for example, a control IP in the case of RIP or an RNA-seq experiment
in the case of CLIP). In this case, the covariate can be provided as a BED
file that contains the raw read locations; the -i option should be used to
indicate that the covariate file is raw reads -- Piranha will then bin these
reads into bins that match the response.

Piranha's output is a tab delimited file with the input regions from the 
response (in BED format) and a single additional column giving the p-value 
for each bin. If covariates were provided, additional columns are added
containing the value of each of the covariates for the given locus.

If the -f option is specified, Piranha will output the model instead of the
scored bins. Models are given in XML format, and can be loaded back into
the program later for scoring a new (or the same) input using the -m option. 


***************************************
5. Simulation of input data for Piranha
*******************************************************************************
We've also included a program for simulating input data - Simulate. After 
building and installing, it can be found in the bin directory. This program 
can produce simulated data from a range of distributions, both regular and 
regression based. Run it without any arguments to see its command line 
options. As an example, to generate 1000 simulated regions where the response 
is dependent on two covariates from a Zero-truncated negative binomial 
regression distribution, you would execute the program as follows, where > is
your prompt (here the command is split over two lines for formatting reasons,
but this need not be the case). 

> ./Simulate -d ZeroTruncatedNegativeBinomialRegression -n 1000 \
> -r response.bed -c "cov1.bed cov2.bed" 


*****************************
6. Frequently asked questions
*******************************************************************************

Q. I get the error: "Failed to split responses, smallest response accounts for 
   more than 99% of all responses. Try increasing the threshold". What does
   this mean? How do I fix it?
A. Piranha assumes most bins in your input are background noise. The program
   models this background and then tests each read count to see if it 
   significantly exceeds the background. This doesn't work if all or most of
   your bins contain the same number of reads though. That is what this 
   message is telling you: more than 99% of the bins in your input contain
   the same number of reads as the smallest read count (usually 1 read, but 
   not neccessarily). Most often this is because you forgot to use the -b
   option to bin your reads and so Piranha thinks each read is its own bin.
   
Q. I am running Piranha with one or more covariates and I got the following
   error (or similiar): "ERROR: evaluating zero-truncated negative binomial 
   regression log-likelihood function with response 1 and distribution 
   parameters beta: 6000,  --- alpha: 100 failed. Reason: result was 
   non-finite"
A. The model fitting algorithm failed to converge. This often happens when
   the covariate(s) contain large values, which is most often the case when
   you provide sequencing data like RIP control IP or RNA-seq data as a
   covariate. Try running the program again with the -l option, which
   converts the read counts from your covariate into log-space; this 
   generally solves the problem.

***************************
7. Contacts and bug reports
*******************************************************************************
Philip J. Uren
uren@usc.edu

Andrew D. Smith
andrewds@usc.edu

If you found a bug in Piranha, we'd like to know about it. Before contacting us
though, please check the following list:

  1.) Are you using the latest version? The bug you found may already have 
      been fixed.
  2.) Check that your input is in the correct format and you have selected
      the correct options.
  3.) Please reduce your input to the smallest possible size that still 
      produces the bug; we will need your input data to reproduce the 
      problem.

