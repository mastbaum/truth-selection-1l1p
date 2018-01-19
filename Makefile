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
	-Isrc \
	$(shell root-config --cflags)

CXXFLAGS=-Wall -pedantic
CXX=g++

LDFLAGS=$(shell root-config --libs) \
	-L $(CANVAS_LIB) -l canvas \
	-L $(CETLIB_LIB) -l cetlib \
	-L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
	-L $(GALLERY_LIB) -l gallery \
	-L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
	-L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
	-L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase

SOURCES = $(wildcard src/*.cxx)
MAINS = build/RunSelection.o
OBJECTS = $(subst src,build,$(SOURCES:.cxx=.o))
LIBOBJS = $(filter-out $(MAINS), $(OBJECTS))
BINS = $(subst build,bin,$(MAINS:.o=))

vpath %.cxx src

LIB = lib/libleets.so

all: dir lib bin

$(LIB): $(LIBOBJS)

dir:
	@mkdir -p lib bin build

%.so: $(LIBOBJS)
	$(CXX) -shared -fPIC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

build/%.o: %.cxx
	$(CXX) -c -fPIC -Wall $(CPPFLAGS) -o $@ $<

lib: $(LIB)

bin/%: $(LIB) %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

bin: $(BINS)

clean:
	$(RM) $(BINS) $(OBJECTS)

