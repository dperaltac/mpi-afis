# The compiler to use.
CC = g++

# Directories for Includes and Common clases
IDIR =../include/
CDIR =../commons/

# Compiler options
CFLAGS= -Wall -O2 -I$(IDIR) -std=c++98

# Sources and Common clases sources
SOURCES= main.cpp MCC.cpp Cylinder.cpp
CSOURCES= $(CDIR)Fingerprint.cpp $(CDIR)Functions.cpp $(CDIR)Minutia.cpp $(CDIR)GrahamScanConvexHull.cpp $(CDIR)File19794.cpp $(CDIR)Munkres.cpp

# Objects
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.cpp=.o)

# Name of the executable
EXECUTABLE=mcc

all: $(SOURCES) $(CSOURCES) $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS) $(COBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(COBJECTS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(COBJECTS) $(EXECUTABLE)
