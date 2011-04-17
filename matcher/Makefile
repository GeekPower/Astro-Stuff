CC	= gcc -Wall -O2 -g 
CFLAGS	= $(shell pkg-config --cflags glib-2.0)
LIBS	= $(shell pkg-config --libs glib-2.0) -lm

TARGETS = newmain groth86

all:	$(TARGETS)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $^

newmain: newmain.o
	$(CC) -o $@ $^ $(LIBS)

groth86: groth86.o
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o $(TARGETS)

