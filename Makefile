include config.mk


# dihbsa/src dependencies
DIHBSA = src
DEPS += -I$(DIHBSA)
LIBS += -L$(DIHBSA)
OBJS := $(wildcard $(DIHBSA)/*.so)


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
	$(CXX) -o $@ $< $(OBJS) $(LIBS)

%.o: %.cpp
	@echo "----- build $@ -----"
	$(CXX) -c $< -o $@ $(FLAGS) $(DEPS) $(LIBS)

clean:
	@cd src; make clean
	$(RM) $(EXES)
