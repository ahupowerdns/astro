-include sysdeps/$(shell uname).inc

VERSION=0.1
CXXFLAGS=-Wall -Inr_c304/code -Iodeint-v2 -MMD -MP  $(CXX2011FLAGS) -O3 # -Wno-unused-local-typedefs 
CFLAGS=-Wall -Iodeint-v2 -O3 -MMD -MP
LDFLAGS=$(CXX2011FLAGS)  
CHEAT_ARG := $(shell ./update-git-hash-if-necessary)

PROGRAMS=hw optim peaker calibrate

all: $(PROGRAMS)

-include *.d

hw: hw.o oscil.o klc.o misc.o
	g++ $^ -lCCfits -lcfitsio -lfftw3 -o $@

optim: optim.o misc.o
	g++ $^ -o $@

peaker: peaker.o misc.o
	g++ $^ -o $@



calibrate: calibrate.o oscil.o
	g++ $^ -o $@

sweep: sweep.o oscil.o
	g++ $^ -o $@


clean:
	rm -f *~ *.o $(PROGRAMS) *.d githash githash.h

