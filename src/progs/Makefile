#    Copyright (C) 2011 University of Southern California,
#                       Philip J. Uren, Andrew D. Smith
#
#    Authors: Philip J. Uren, Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

PROGS = Piranha Simulate

CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS = -g 
COMMON_DIR = ../common/
BAMTOOLS_ROOT = ../BAMTools
BAMTOOLS_INCLUDE = ../BAMTools/include
BAMTOOLS_LIB = $(BAMTOOLS_ROOT)/lib 
LIBS = -lgsl -lgslcblas -L$(BAMTOOLS_LIB) -lbamtools

INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR) $(BAMTOOLS_INCLUDE)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif

all:    $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< $(INCLUDEARGS)
	
%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c -o $@ $< $(INCLUDEARGS)	

%: %.cpp
	$(CXX) $(CFLAGS) -Wl,-rpath,`cd $(BAMTOOLS_LIB); pwd` -o $@ \
	$(filter-out $(BAMTOOLS_LIB), $^)  $(INCLUDEARGS) $(LIBS)
	
$(BAMTOOLS_LIB) :
	cd $(BAMTOOLS_ROOT); mkdir build; cd build; cmake ..; make

Piranha: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o MappedRead.o)
Piranha : $(addprefix $(COMMON_DIR), \
	Gaussian.o NegativeBinomial.o \
	PoissonRegression.o Poisson.o SimpleXML.o \
	NegativeBinomialRegression.o MixtureModel.o \
	ZTNBRegression.o RegressionMixtureModel.o \
	LinearRegression.o Distribution.o \
	TruncatedPoissonRegression.o VectorMath.o Matrix.o \
	RegressionModel.o ZTP.o \
	DistributionMixtureModel.o ZTNB.o GeneralisedLinearModel.o ReadBinner.o)
Piranha: $(BAMTOOLS_LIB)

Simulate : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o)
Simulate: $(BAMTOOLS_LIB)

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
	
distclean:
	@-rm -f $(PROGS) *.o *.so *.a *~
	@rm -rf $(BAMTOOLS_ROOT)/build $(BAMTOOLS_ROOT)/include 
	@rm -rf $(BAMTOOLS_ROOT)/bin $(BAMTOOLS_ROOT)/lib
	
test:
	make -C ../test SMITHLAB_CPP=$(SMITHLAB_CPP) test; \

.PHONY: clean
