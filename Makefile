include makeconfig.mk


# dihbsa/src dependencies
DEPS += -Isrc
#LIBS += -L. -l$(DIHBSA)


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
	$(CXX) -o $@ $< ./$(DIHBSAOBJ) $(LIBS)

%.o: %.cpp
	@echo "----- build $@ -----"
	$(CXX) -c $^ -o $@ $(FLAGS) $(DEPS)

clean:
	@cd src; make clean
	$(RM) $(EXES) *.pcm
