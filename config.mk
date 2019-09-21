# specify which particle bank to read (data vs. MC)
# 0 = REC::Particle -- for data
# 1 = MC::Lund -- for MC generated (helicity injected)
# 2 = MC::Particle -- for MC generated (helicity injected)
# 3 = REC::Particle -- for MC reconstructed (helicity set by event-match to MCgen)
PARTICLE_BANK = 3


# Clas12Tool directory, relative to this config file's `pwd`
CLAS12TOOLDIR = ../Clas12Tool

####################################

# compiler and flags
CXX = g++ -std=c++11
FLAGS = -g -Wno-deprecated -fPIC -m64 -fno-inline -Wno-write-strings
FLAGS += -DPARTICLE_BANK=$(PARTICLE_BANK)

# extra flags for valgrind
FLAGS += -O0


# ROOT
DEPS = $(shell root-config --cflags)
LIBS = $(shell root-config --glibs)
LIBS += -lMinuit -lRooFitCore -lRooFit


# Clas12Tool
THISDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
CLAS12TOOL := $(THISDIR)$(CLAS12TOOLDIR)
DEPS += -I$(CLAS12TOOL)/hipo4 -I$(CLAS12TOOL)/Clas12Banks
LIBS += -L$(CLAS12TOOL)/lib -lClas12Banks -lHipo4 -llz4


# DihBsa shared object name and source directory
DIHBSA = DihBsa
DIHBSAOBJ := lib$(DIHBSA).so

