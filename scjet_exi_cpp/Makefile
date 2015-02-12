# Compile SCjet example
# S.Chekanov (ANL)  

CXXFLAGS= -O2 -Wall -fPIC
LDFLAGS=-g

all: main

main: main.o KT.o ParticleD.o
	$(CXX) $(LDFLAGS) -o $@ $@.o KT.o ParticleD.o

clean:
	rm -vf main *.o *~ 

run:
	./main < ../data/single-event.dat

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

KT.o: KT.h Timer.h 
ParticleD.o: ParticleD.h 
