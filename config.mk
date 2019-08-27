# HIPO version
HIPO_VERSION = 4

# Clas12Tool directory, relative to this config file's `pwd`
ifeq ($(HIPO_VERSION),3)
  CLAS12TOOLDIR = ../Clas12Tool.hipo3and4
else ifeq ($(HIPO_VERSION),4)
  CLAS12TOOLDIR = ../Clas12Tool
else
  $(error Bad HIPO version setting)
endif

####################################

# compiler and flags
CXX = g++ -std=c++11
FLAGS = -g -Wno-deprecated -fPIC -m64 -fno-inline -Wno-write-strings
FLAGS += -DHIPO_VERSION=$(HIPO_VERSION)

# extra flags for valgrind
FLAGS += -O0


# ROOT
DEPS = $(shell root-config --cflags)
LIBS = $(shell root-config --glibs)
LIBS += -lMinuit -lRooFitCore -lRooFit


# Clas12Tool
THISDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
CLAS12TOOL := $(THISDIR)$(CLAS12TOOLDIR)
ifeq ($(HIPO_VERSION),3)
  DEPS += -I$(CLAS12TOOL)/Hipo3 -I$(CLAS12TOOL)/Clas12Banks3
  LIBS += -L$(CLAS12TOOL)/lib -lClas12Banks3 -lHipo3 -llz4
else ifeq ($(HIPO_VERSION),4)
  DEPS += -I$(CLAS12TOOL)/hipo4 -I$(CLAS12TOOL)/Clas12Banks
  LIBS += -L$(CLAS12TOOL)/lib -lClas12Banks -lHipo4 -llz4
endif


# DihBsa shared object name and source directory
DIHBSA = DihBsa
DIHBSAOBJ := lib$(DIHBSA).so

