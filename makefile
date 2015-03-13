# This makefile assumes a POSIX compliant system (linux)
#
.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h
#
#

VPATH = .
MARCH = $(shell if [ `uname -m` = "x86_64" ]; then echo "athlon64"; else echo "i686"; fi)

all: estimator statestimator exactgroupby sample generate1drandomdata

exactgroupby: exactgroupby.cpp
	g++  -Wall -O2 -g2  -march=$(MARCH) -o exactgroupby exactgroupby.cpp

statestimator: statestimator.cpp
	g++  -Wall -O2 -g2  -march=$(MARCH) -o statestimator statestimator.cpp

estimator: estimator.cpp randomhasher.h trailingzeros.h estimator.h stogibbons.h
	g++  -Wall -O2 -g2  -march=$(MARCH) -o estimator estimator.cpp

sample: sample.cpp
	g++  -Wall -O2 -g2  -march=$(MARCH) -o sample sample.cpp

generate1drandomdata: generate1drandomdata.cpp
	g++  -Wall -O2 -g2  -march=$(MARCH) -o generate1drandomdata generate1drandomdata.cpp

data: USCensus1990.data.txt

unit: view data
	./unittesting.py

test: view data
	./estimator --groupby dAge,iClass,iMarital --8bits USCensus1990.data.txt
	echo
	./estimator --groupby dAge,iClass,iMarital --32bits USCensus1990.data.txt
	echo
	./estimator --groupby dAge,iClass,iMarital  USCensus1990.data.txt

USCensus1990.data.txt:
	wget http://kdd.ics.uci.edu/databases/census1990/USCensus1990.data.txt

clean :
	rm -f *.o
	rm -f estimator
	rm -f *.a
	rm -f tags
