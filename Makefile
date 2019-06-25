include config.mk


# dihbsa/src dependencies
DEPS += -Isrc
LIBS += -Lsrc -l:DihBsa.so
FLAGS += -Wl,-rpath,$(shell pwd)/src


# assume each .cpp file has main and build corresponding .exe executable
SOURCES := $(basename $(wildcard *.cpp))
EXES := $(addsuffix .exe, $(SOURCES))


#--------------------------------------------


all: 
	@cd src; make
	make exe

exe: $(EXES)

%.exe: %.o
	$(CXX) -o $@ $< $(LIBS) $(FLAGS)

%.o: %.cpp
	@echo "--- build $< ---"
	$(CXX) -c $< -o $@ $(FLAGS) $(DEPS) $(LIBS)

clean:
	@cd src; make clean
	$(RM) $(EXES)
