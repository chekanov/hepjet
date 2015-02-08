# basic Makefile to build an example using the fjcore.cc|hh files.
# If you do not wish to use this Makefile, all you need to do is
#
#   g++ 01-basic.cc fjcore.cc -o 01-basic
#

CXXFLAGS=-O2

all: 01-basic

01-basic: 01-basic.o fjcore.o
	$(CXX) $(LDFLAGS) -o $@ $@.o fjcore.o

clean:
	rm -vf *.o 01-basic

distclean: clean
	rm -f fjcore.o

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

fjcore.o: fjcore.hh
01-basic.o: fjcore.hh Timer.h
