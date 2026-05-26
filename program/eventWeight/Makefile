main = $(shell basename `pwd`)
CWD=$(PWD)
CXX      = g++
OPT2       = -O2
CXXFLAGS = -O2 -Wall -fPIC

ROOTCONFIG = root-config
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs)

CXXFLAGS += $(ROOTCFLAGS)
LIBS     = $(ROOTLIBS)

CXXFLAGS += -I$(ATMNCUtil_C)/include
LIBS += -L$(ATMNCUtil_C)/lib/ -latmnc_util_c

# Action starts
# =============

all: bin/$(main)
bin/$(main): $(main).o 
	-mkdir bin
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@
$(main).o: $(main).cc
	g++ $(CXXFLAGS)  -c $(main).cc
clean:
	-rm -f bin/* *.o
.PHONY: all clean 
