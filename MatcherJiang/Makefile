# The compiler to use.
CC = g++

# Directories for Includes and Common clases
IDIR =../include
CDIR =../commons/


# Compiler options
CFLAGS= -Wall -I$(IDIR) -O2

# Sources and Common clases sources
SOURCES= main.cpp FingerprintJiang.cpp
CSOURCES= $(CDIR)Fingerprint.cpp $(CDIR)Functions.cpp $(CDIR)Minutia.cpp $(CDIR)File19794.cpp

# Objects
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.cpp=.o)

# Name of the executable
EXECUTABLE=jiangMatching

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(COBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(COBJECTS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(COBJECTS) $(EXECUTABLE)
