# This makefile assumes a POSIX compliant system (linux)
#
.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h
#
#

all: estimator statestimator exactgroupby sample generate1drandomdata

exactgroupby: src/exactgroupby.cpp
	c++  -std=c++17 -Wall -O3  -o exactgroupby src/exactgroupby.cpp -Iinclude

statestimator: src/statestimator.cpp
	c++ -std=c++17 -Wall -O3  -o statestimator src/statestimator.cpp -Iinclude

estimator: src/estimator.cpp include/randomhasher.h include/trailingzeros.h include/estimator.h include/stogibbons.h
	c++ -std=c++17  -Wall -O3  -o src/estimator src/estimator.cpp -Iinclude 

sample: src/sample.cpp
	c++ -std=c++17  -Wall -O3  -o sample src/sample.cpp -Iinclude

generate1drandomdata: src/generate1drandomdata.cpp
	c++ -std=c++17  -Wall -O3  -o src/generate1drandomdata src/generate1drandomdata.cpp -Iinclude

data: USCensus1990.data.txt

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
