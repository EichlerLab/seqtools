CXXFLAGS = -I$(abspath include) -I/home/paudano/tools/modules/boost/1.65.1/include -I/home/paudano/tools/modules/htslib/1.6/include --std=c++14

LDFLAGS = -L/home/paudano/tools/modules/boost/1.65.1/lib -L/home/paudano/tools/modules/htslib/1.6/lib
LDLIBS = -lboost_program_options -lrt -lpthread -lbz2 -lhts -llzma

.PHONY: all
all: $(addprefix,bin/,subseqfa)

bin/subseqfa: build/subseqfa.o build/util/err.o
	mkdir -p $(dir $@)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${LDFLAGS} ${LDLIBS} -o $@ $^

build/util/%.o: src/util/%.cpp include/seqtools/util/%.h
	mkdir -p $(dir $@)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

build/%.o: src/%.cpp
	mkdir -p $(dir $@)
	${CXX} -c ${CPPFLAGS} ${CXXFLAGS} -o $@ $<
