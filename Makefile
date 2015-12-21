CXX = g++

CXXFLAGS = -g -Wall -fopenmp -O2

PROGRAMS = rjFULL collapsed withinGeneCollapsedSampler

VAR1=$(PWD)

all: $(PROGRAMS)
	echo "#!/bin/bash" > cpThis.bash
	echo "#!/bin/bash" > cpThisCOLLAPSED.bash
	echo "#!/bin/bash" > cpThisCOLLAPSEDWITHIN.bash
	(echo "cp $(VAR1)/sparse.so sparse.so"; echo "cp $(VAR1)/partitionParallelCut.R partitionMerge2.R"; echo "cp $(VAR1)/getClusters.R getClusters.R" ) >> cpThis.bash
	(echo "cp $(VAR1)/sparse.so sparse.so"; echo "cp $(VAR1)/partitionParallelCutCOLLAPSED.R partitionMerge2.R"; echo "cp $(VAR1)/getClusters.R getClusters.R" ) >> cpThisCOLLAPSED.bash
	(echo "cp $(VAR1)/readProbWithinGene.R partitionMerge2.R"; echo "cp $(VAR1)/getClustersWithin.R getClustersWithin.R" ) >> cpThisCOLLAPSEDWITHIN.bash
	chmod u+x cpThis.bash
	chmod u+x cpThisCOLLAPSED.bash
	chmod u+x cpThisCOLLAPSEDWITHIN.bash
	chmod u+x rjBitSeq
	chmod u+x cjBitSeq
	chmod u+x cpThis.bash
	chmod u+x cpThisCOLLAPSED.bash
	chmod u+x cpThisCOLLAPSEDWITHIN.bash
	chmod u+x rjBitSeq
	chmod u+x cjBitSeq
	R CMD SHLIB sparse.cpp
get_k.o: get_k.h get_k.cpp infiles.h

get_k_split.o: get_k_split.h get_k_split.cpp infiles.h

get_k_within.o: get_k_within.h get_k_within.cpp infiles.h

rjFULL: get_k_split.o rjmcmcHOME.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o rjFULL get_k_split.o rjmcmcHOME.cpp

collapsed: get_k_split.o collapsed.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o collapsed get_k_split.o collapsed.cpp

# this contains geneDEprob
withinGeneCollapsedSampler: get_k_within.o withinGeneCollapsedPrior2.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o withinGeneCollapsedSampler get_k_within.o withinGeneCollapsedPrior2.cpp

# this is for the jeffreys prior
#withinGeneCollapsedSampler: get_k_within.o withinGeneCollapsed.cpp infiles.h
#	$(CXX) $(CXXFLAGS) -o withinGeneCollapsedSampler get_k_within.o withinGeneCollapsed.cpp

clean:
	rm *.o $(PROGRAMS)

