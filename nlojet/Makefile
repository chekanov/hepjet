# S.Chekanov

# Assign or add variables
OPT2          = -O2
CXX           = g++
CXXFLAGS      = $(OPT2) -Wall -fPIC
LD            = g++
LDFLAGS       = $(OPT2) 
LIBS          = -lm 

SOURCE_FILES := $(shell ls -1 main.cc)
SOURCE_FILES += $(shell ls -1 src/*.cc)
# SOURCE_FILES += $(shell ls -1 src/bits/*.cc)


INCLUDE1= -I./src 
INCLUDE2= -I.

# build object files 
objects       = $(patsubst %.cc,%.o,$(SOURCE_FILES))


%.o: %.cc 
	$(CXX) $(OPT) $(CXXFLAGS) $(INCLUDE1) $(INCLUDE2) -o $@ -c $<

Tasks:     main

run:
	./main < ../data/single-event.dat


main: $(objects)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $(OutPutOpt)$@
clean:
	        @rm -f *.o  main *~ src/*.o src/*/*.o  src/*/*/*.o; echo "Clear.." 
