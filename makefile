#Completely stolen.
CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CETLIB_INC) \
         -I $(CETLIB_EXCEPT_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(NUSIMDATA_INC) \
         -I ${UBOONECODE_INC} \
         $(shell root-config --cflags)

CXXFLAGS=-std=c++14 -Wall -pedantic
CXX=g++
LDFLAGS=$(shell root-config --libs) \
        -L $(CANVAS_LIB) -l canvas \
        -L $(CETLIB_LIB) -l cetlib \
        -L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
        -L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase

OBJECTS = RunSelection.o TSSelection.o TSUtil.o
SOURCES = $(OBJECTS:.o=.cxx)

UNAME := $(shell uname -s)

ifeq ($(UNAME), Darwin)
  EXEC=truth_selection
else
  EXEC=truth_selection
endif

$(EXEC): $(OBJECTS)
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

%.o : %.cxx
	@$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm $(EXEC) $(OBJECTS)


