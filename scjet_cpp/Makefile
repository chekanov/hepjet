#   g++ 01-basic.cc fjcore.cc -o 01-basic
#

CXXFLAGS= -O2 -Wall -fPIC
LDFLAGS=-g

all: main

main: main.o KT.o ParticleD.o
	$(CXX) $(LDFLAGS) -o $@ $@.o KT.o ParticleD.o

clean:
	rm -vf main *.o *~ 

distclean: clean
	rm -f *.o

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

KT.o: KT.h Timer.h 
ParticleD.o: ParticleD.h 
