CXX ?= g++
EDEXEC= bin/ed_aligner 
EDGEXEC= bin/ed_genome_aligner

EDLIB=./edlib/edlib/include
EDLIBCPP=./edlib/edlib/src/edlib.cpp

CFLAGS = -O3 -std=c++11 -I ${EDLIB}
EDSOURCES= src/ed_aligner.cpp  ${EDLIBCPP}
EDGSOURCES= src/ed_genome_aligner.cpp  ${EDLIBCPP}
ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

.PHONY: test clean

edlib:
ifneq ($(wildcard ./edlib/.),)
	@echo "Found edlib in current directory"
else
	@echo "Did not find edlib."
	git clone https://github.com/Martinsos/edlib.git
endif


all: edlib  $(EDEXEC) $(EDGEXEC) 

clean:
	-rm -f  $(EDEXEC) $(EDGEXEC)

$(EDEXEC): ${EDSOURCES} edlib
	$(CXX) -o $@  ${EDSOURCES} $(CFLAGS)

$(EDGEXEC): ${EDGSOURCES} edlib
	$(CXX) -o $@  ${EDGSOURCES} $(CFLAGS)

