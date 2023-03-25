CC = g++
CFLAGS = -c -Wall
LDFLAGS = 

ROOTCONFIG = /usr/local/Cellar/root/6.26.06_2/bin/root-config
CFLAGS += $(shell $(ROOTCONFIG) --cflags)
LDFLAGS += $(shell $(ROOTCONFIG) --libs)

SOURCES = main.cpp molecula.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main

all: $(SOURCES) $(EXECUTABLE)

rodri: 
	$(CC) main.cpp molecula.cpp `root-config --cflags --glibs` -o main

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)