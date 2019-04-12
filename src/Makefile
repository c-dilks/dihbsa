ARCH = linux

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS := $(shell root-config --glibs)

ifeq ($(ARCH),linux)
CXX = g++
CXXFLAGS = -g -O -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings
LD = g++
LDFLAGS = -O
SOFLAGS = -shared
INCLUDE = $(ROOTSYS)/include/
endif

CXXFLAGS += $(ROOTCFLAGS)
LIBS = $(ROOTLIBS) $(SYSLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS)

#-----------------------------------------------

CXXFILES = Observable.cxx
HEADERS = Constants.h Observable.h

#-----------------------------------------------

All: DihBsa.so

DihBsa.so: LinkDef.h $(HEADERS) $(CXXFILES)
	@echo "generating dictionary $@..."
	@rootcint -f DihBsaDict.cxx -c $(HEADERS) LinkDef.h
	g++ $(CXXFILES) DihBsaDict.cxx -shared -o DihBsa.so $(CXXFLAGS) $(GLIBS)

clean:
	@echo "removing DihBsa.so..."
	@rm DihBsa.so
