CLAS12TOOL  := ../Clas12Tool

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I$(CLAS12TOOL)/Hipo3 -I$(CLAS12TOOL)/Clas12Banks3
HIPOLIBS    := -L$(CLAS12TOOL)/lib -lClas12Banks3 -lHipo3

DIHBSACFLAGS  := -Isrc 
DIHBSALIBS    := -Lsrc ./DihBsa.so

LZ4LIBS     := -L$(CLAS12TOOL)/Lz4/lib -llz4
LZ4INCLUDES := -I$(CLAS12TOOL)/Lz4/lib

CXX       := g++ -std=c++11
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

ALLFLAGS = $(ROOTCFLAGS) $(ROOTLDFLAGS) $(DIHBSALIBS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)


all: 
	@echo "compiling src..."
	@cd src; make
	@echo ""
	@echo "now compiling main executables..."
	@echo ""
	make analysis
	@echo ""
	make diagnostics
	@echo ""
	make acceptance
	@echo ""
	make phiRtest
	@echo ""
	make asym
	@echo ""

analysis: analysis.o
	$(CXX) -o analysis.exe $< $(ALLFLAGS)
	 
diagnostics: diagnostics.o
	$(CXX) -o diagnostics.exe $< $(ALLFLAGS)

acceptance: acceptance.o
	$(CXX) -o acceptance.exe $< $(ALLFLAGS)

phiRtest: phiRtest.o
	$(CXX) -o phiRtest.exe $< $(ALLFLAGS)
	
asym: asym.o
	$(CXX) -o asym.exe $< $(ALLFLAGS)

clean:
	@cd src; make clean
	@echo 'Removing all build files'
	@rm -rf *.exe
	@rm -rf *.o

%.o: %.cpp
	$(CXX) -c $< $(ROOTCFLAGS) $(DIHBSACFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
