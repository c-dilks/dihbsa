ARCH = linux

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS := $(shell root-config --glibs)

CLAS12TOOL  := ../../Clas12Tool
HIPOCFLAGS  := -I$(CLAS12TOOL)/Hipo -I$(CLAS12TOOL)/Banks
HIPOLIBS    := -L$(CLAS12TOOL)/lib -lhipo -lclas12banks
LZ4LIBS     := -L$(CLAS12TOOL)/Lz4/lib -llz4
LZ4INCLUDES := -I$(CLAS12TOOL)/Lz4/lib

CXX = g++ -std=c++11

ifeq ($(ARCH),linux)
CXX = g++ -std=c++11
CXXFLAGS = -g -O -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings
LD = g++
LDFLAGS = -O
SOFLAGS = -shared
INCLUDE = $(ROOTSYS)/include/
endif

CXXFLAGS += $(ROOTCFLAGS) $(HIPOCFLAGS)
LIBS = $(ROOTLIBS) $(SYSLIBS) $(HIPOLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS) $(HIPOLIBS)

#-----------------------------------------------

CXXFILES = DIS.cxx Trajectory.cxx Dihadron.cxx EventTree.cxx Asymmetry.cxx
HEADERS = Constants.h DIS.h Trajectory.h Dihadron.h EventTree.h Asymmetry.h

#-----------------------------------------------

All: DihBsa.so

DihBsa.so: LinkDef.h $(HEADERS) $(CXXFILES)
	@echo "generating dictionary $@..."
	@rootcint -f DihBsaDict.cxx -c $(HEADERS) LinkDef.h
	$(CXX) $(CXXFILES) DihBsaDict.cxx -shared -o DihBsa.so $(CXXFLAGS) $(GLIBS)

clean:
	@echo "removing DihBsa.so..."
	@rm DihBsa.so
