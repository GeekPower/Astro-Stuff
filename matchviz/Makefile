CC	= gcc -Wall -O2
CFLAGS	= $(shell pkg-config --cflags gtk+-2.0)
LIBS	= $(shell pkg-config --libs gtk+-2.0) -export-dynamic

TARGETS = matchviz triviz

all:	$(TARGETS)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $^

matchviz: matchviz.o
	$(CC) -o $@ $^ $(LIBS)

triviz: triviz.o
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o $(TARGETS)

