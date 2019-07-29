# The compiler to use.
CC = mpiCC

# Directories for Includes and Common clases
IDIR =include
CDIR =commons/
JIANGDIR =MatcherJiang/
MCCDIR =MCC/
PDIR =Parallel/
BINDIR =bin/

# Compiler options -Weffc++
CFLAGS= -Wall -O2 -fopenmp -I$(IDIR) -I$(JIANGDIR) -I$(MCCDIR) -DCOMPILE_USING_MPI

# Sources and Common clases sources
SOURCES= $(PDIR)genericMatching.cpp
SOURCESD= $(PDIR)DPDDFF.cpp
CSOURCES= $(CDIR)Fingerprint.cpp $(CDIR)Score.cpp $(JIANGDIR)FingerprintJiang.cpp  $(MCCDIR)MCC.cpp $(MCCDIR)Cylinder.cpp  $(CDIR)Functions.cpp $(CDIR)Minutia.cpp $(CDIR)GrahamScanConvexHull.cpp $(CDIR)Munkres.cpp $(PDIR)ParallelHandler.cpp $(PDIR)ParallelMaster.cpp $(PDIR)ParallelSlave.cpp $(PDIR)IOHandler.cpp $(CDIR)File19794.cpp


# Objects
OBJECTS=$(SOURCES:.cpp=.o)
OBJECTSD=$(SOURCESD:.cpp=.o)
COBJECTS=$(CSOURCES:.cpp=.o)

# Name of the executable
EXECUTABLE=$(BINDIR)genericMatching
EXECUTABLED=$(BINDIR)DPDDFF

all: $(EXECUTABLE) $(EXECUTABLED)

.PHONY: doc

doc:
	doxygen Doxyfile

$(EXECUTABLE): $(OBJECTS) $(COBJECTS)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(OBJECTS) $(COBJECTS) $(OBJECTFILES) -o $@ $(LDFLAGS)

$(EXECUTABLED): $(OBJECTSD) $(COBJECTS)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(OBJECTSD) $(COBJECTS) $(OBJECTFILES) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(OBJECTSD) $(COBJECTS) $(EXECUTABLE) $(EXECUTABLED)

mrproper: clean
	rm -r doc/latex doc/html
