CC=g++
CFLAGS=-g -c -Wall -I. -I/usr/X11/include
LDFLAGS=
SOURCES=main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=balloon
INCLUDE=
LIBS =  -L /usr/X11R6/lib/ -lglut -lGL -lGLU -lX11 -lXmu -lXi -lm
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $(INCLUDE) $@

clean:
	@rm $(EXECUTABLE) $(OBJECTS)
	@rm *.txt
