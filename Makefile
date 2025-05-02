TOPDIR = $(PWD)

SRCDIR = src
SRCS = $(wildcard ${SRCDIR}/*.cc)

OBJDIR = obj
OBJS = ${patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,${SRCS}}


ROOTCONFIG = root-config
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
CXXFLAGS = -Wall -c -fPIC $(ROOTCFLAGS) -I$(TOPDIR)/include/
ROOTLIBS = $(shell  root-config --libs)
#CXXFLAGS += -I$(SATOPLOT) 
#LIBS += -L$(SATOPLOT)/lib/ -lPlotter 

LIBRARYUTIL = lib/libatmnc_util_c.so

lib: $(OBJS)
	g++ -shared $(ROOTLIBS) -o $(LIBRARYUTIL) $(OBJS)
$(OBJDIR)/%.o : $(SRCDIR)/%.cc include/%.h
	g++ $(CXXFLAGS) $< -o $@
clean:
	echo $(OBJS) 	
	rm -f *.o $(LIBRARYUTIL)
