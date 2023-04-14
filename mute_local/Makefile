CC = g++
CFLAGS = -c -Wall
LDFLAGS = 

ROOTCONFIG = /usr/local/Cellar/root/6.26.06_2/bin/root-config
CFLAGS += $(shell $(ROOTCONFIG) --cflags)
LDFLAGS += $(shell $(ROOTCONFIG) --libs)

SOURCES = main.cpp molecule.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main
TILDE := $(wildcard */*~) $(wildcard *~)

all: $(SOURCES) $(EXECUTABLE)

rodri: 
	$(CC) main.cpp molecule.cpp `root-config --cflags --glibs` -o main

parallel:
	$(CC) -fopenmp main.cpp molecule.cpp `root-config --cflags --glibs` -o main

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(TILDE)
