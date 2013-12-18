-include sysdeps/$(shell uname).inc

VERSION=0.1
CXXFLAGS=-Wall -I. -MMD -MP  $(CXX2011FLAGS) -ggdb # -Wno-unused-local-typedefs 
CFLAGS=-Wall -I. -O3 -MMD -MP
LDFLAGS=$(CXX2011FLAGS)  
CHEAT_ARG := $(shell ./update-git-hash-if-necessary)

PROGRAMS=hw

all: $(PROGRAMS)

-include *.d

hw: hw.o oscil.o
	g++ $^ -lCCfits -o $@

calibrate: calibrate.o oscil.o
	g++ $^ -o $@

clean:
	rm -f *~ *.o hw calibrate *.d githash githash.h

