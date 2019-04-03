CLAS12TOOL  := ../Clas12Tool

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I$(CLAS12TOOL)/Hipo -I$(CLAS12TOOL)/Banks
HIPOLIBS    := -L$(CLAS12TOOL)/lib -lhipo -lclas12banks

LZ4LIBS     := -L$(CLAS12TOOL)/Lz4/lib -llz4
LZ4INCLUDES := -I$(CLAS12TOOL)/Lz4/lib

CXX       := g++ -std=c++11
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)


all: analysis

analysis: analysis.o
	$(CXX) -o analysis.exe $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.exe

%.o: %.cc
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
