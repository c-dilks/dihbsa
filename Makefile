include config.mk


# dihbsa/src dependencies
SRCDIR = ./src
DEPS += -I$(SRCDIR)
#LIBS += -L$(SRCDIR) -l$(DIHBSA)
SOBJ := $(SRCDIR)/$(DIHBSAOBJ)


# assume each .cpp file has main and build corresponding .exe executable
SOURCES := $(basename $(wildcard *.cpp))
EXES := $(addsuffix .exe, $(SOURCES))


#--------------------------------------------


all: 
	@cd src; make
	make exe

exe: $(EXES)

%.exe: %.o
	@echo "--- make executable $@"
	$(CXX) -o $@ $< $(SOBJ) $(LIBS)

%.o: %.cpp
	@echo "----- build $@ -----"
	$(CXX) -c $^ -o $@ $(FLAGS) $(DEPS)

clean:
	@cd src; make clean
	$(RM) $(EXES)
