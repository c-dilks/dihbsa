ARCH = linux

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS := $(shell root-config --glibs)

CLAS12TOOL  := ../../Clas12Tool
HIPOCFLAGS  := -I$(CLAS12TOOL)/Hipo3 -I$(CLAS12TOOL)/Clas12Banks3
HIPOLIBS    := -L$(CLAS12TOOL)/lib -lClas12Banks3 -lHipo3
LZ4LIBS     := -L$(CLAS12TOOL)/Lz4/lib -llz4
LZ4INCLUDES := -I$(CLAS12TOOL)/Lz4/lib

CXX = g++ -std=c++11

ifeq ($(ARCH),linux)
CXX = g++ -std=c++11
CXXFLAGS = -g -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings
LD = g++
SOFLAGS = -shared
INCLUDE = $(ROOTSYS)/include/
endif

CXXFLAGS += $(ROOTCFLAGS) $(HIPOCFLAGS)
LIBS = $(ROOTLIBS) $(SYSLIBS) $(HIPOLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS) $(HIPOLIBS)

#-----------------------------------------------

CXXFILES = DIS.cxx Trajectory.cxx Dihadron.cxx Diphoton.cxx EventTree.cxx Binning.cxx Asymmetry.cxx
HEADERS = Constants.h Tools.h DIS.h Trajectory.h Dihadron.h Diphoton.h EventTree.h Binning.h Asymmetry.h

#-----------------------------------------------

All: DihBsa.so

DihBsa.so: LinkDef.h $(HEADERS) $(CXXFILES)
	@echo "--- generating dictionary $@..."
	@rootcling -f DihBsaDict.cxx -c $(HEADERS) LinkDef.h
	@echo "--- now compiling..."
	$(CXX) $(CXXFILES) DihBsaDict.cxx -shared -o ../DihBsa.so $(CXXFLAGS) $(GLIBS)
	@echo "--- done"

clean:
	@echo "--- removing DihBsa.so..."
	@rm ../DihBsa.so
