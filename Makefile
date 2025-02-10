.PHONY: all clean install uninstall

CC=g++

CFLAGS=-c -Ofast -ffast-math -std=c++17 -Wall \
	-Wcast-qual -Wpointer-arith -Wcast-align  

LDFLAGS= 

SOURCES=browanian.cpp dataio.cpp langevin.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=browanian

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -f *.o

distclean: clean 
	rm *.data