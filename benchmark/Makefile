# basic Makefile to build an example using the fjcore.cc|hh files.
# If you do not wish to use this Makefile, all you need to do is
#
#   g++ 01-basic.cc fjcore.cc -o 01-basic
#

CXXFLAGS=-O2

all: benchmark

benchmark: SCJet.o ParticleD.o benchmark.o fjcore.o
	$(CXX) $(LDFLAGS) -o $@ $@.o fjcore.o SCJet.o ParticleD.o

run:
	./benchmark < ../data/single-event.dat
	
clean:
	rm -vf *.o *~ benchmark

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

fjcore.o: fjcore.hh
01-basic.o: fjcore.hh Timer.h
