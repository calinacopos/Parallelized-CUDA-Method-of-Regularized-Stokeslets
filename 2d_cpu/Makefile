CC=g++
CFLAGS=-O2 -g -c -fsingle-precision-constant -Wall -I. -I/usr/X11/include
LDFLAGS=
SOURCES=main_slip.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=slip
INCLUDE=
LIBS =  -L /usr/X11R6/lib/ -lglut -lGL -lGLU -lX11 -lXmu -lXi -lm
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $(INCLUDE) $@

clean:
	@rm $(EXECUTABLE) $(OBJECTS)
cleanall: clean
	@rm *.txt
