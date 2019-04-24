#DIR_BOOST ?= /home/paudano/tools/modules/boost/1.65.1
#DIR_HTSLIB ?= /home/paudano/tools/modules/htslib/1.9

CXXFLAGS += -I$(abspath include) -I${DIR_BOOST}/include -I${DIR_HTSLIB}/include --std=c++14

override LDFLAGS += -L${DIR_BOOST}/lib -L${DIR_HTSLIB}/lib -Wl,-rpath,${DIR_BOOST}/lib,-rpath,${DIR_HTSLIB}/lib -lboost_program_options -lrt -lpthread -lbz2 -lhts -llzma

.PHONY: all
all: $(addprefix,bin/,subseqfa)

bin/subseqfa: build/subseqfa.o build/util/err.o
	mkdir -p $(dir $@)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${LDFLAGS} -o $@ $^

build/util/%.o: src/util/%.cpp include/seqtools/util/%.h
	mkdir -p $(dir $@)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

build/%.o: src/%.cpp
	mkdir -p $(dir $@)
	${CXX} -c ${CPPFLAGS} ${CXXFLAGS} -o $@ $<

clean:
	rm -rf build
