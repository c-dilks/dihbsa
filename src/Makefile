include ../config.mk

# rootcling dictionary
NAME = DihBsa
DICT := $(NAME)Dict.cxx
PCM := $(NAME)Dict_rdict.pcm
LINKDEF = LinkDef.h


# source code (with $(DICT) and $(LINKDEF) at end of lists)
SOURCES := $(filter-out $(DICT), $(wildcard *.cxx)) $(DICT)
HEADERS := $(filter-out $(LINKDEF), $(wildcard *.h)) $(LINKDEF)


# shared object
OBJ := $(NAME).so

#-----------------------------------------------


$(OBJ): $(DICT) $(HEADERS) $(SOURCES)
	@echo "----- build $@ -----"
	$(CXX) $(SOURCES) -shared -o $@ $(FLAGS) $(DEPS) $(LIBS)

$(DICT): $(HEADERS)
	@echo "----- generate dictionary $@ -----"
	@rootcling -f $@ -c $^

clean:
	$(RM) $(OBJ) $(DICT) $(PCM)
