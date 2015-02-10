# basic Makefile to build an example using the fjcore.cc|hh files.
# If you do not wish to use this Makefile, all you need to do is
#
#   g++ 01-basic.cc fjcore.cc -o 01-basic
#

CXXFLAGS=-O2

all: copy benchmark

copy:
	cp -fv ../fastjet/fjcore.* .
	cp -fv ../scjet_cpp/KT.* .
	cp -fv ../scjet_cpp/ParticleD.* .

benchmark: KT.o ParticleD.o benchmark.o fjcore.o
	$(CXX) $(LDFLAGS) -o $@ $@.o fjcore.o KT.o ParticleD.o


run:
	./benchmark < ../data/single-event.dat
	
clean:
	rm -vf *.o *~ benchmark

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

fjcore.o: fjcore.hh
01-basic.o: fjcore.hh Timer.h
