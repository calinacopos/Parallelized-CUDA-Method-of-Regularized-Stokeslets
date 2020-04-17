CC=g++
CFLAGS= -g -c -fsingle-precision-constant -Wall -I. -I/usr/X11/include
LDFLAGS=
SOURCES=energy_slip.cpp 
#SOURCES=slip_full.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=slip
#EXECUTABLE=full
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
