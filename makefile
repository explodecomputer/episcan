CC = gcc
CFLAGS = -c -O2 -fopenmp -Wall
LDFLAGS = -lm -Wall -O2 -fopenmp
SOURCES1 = main.c readplink.c timing.c

OBJECTS1 = $(SOURCES1:.c=.o)

EXE1 = episcan

all: $(SOURCES1) $(EXE1)


$(EXE1): $(OBJECTS1)
	$(CC) $(LDFLAGS) $(OBJECTS1) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

all:
	rm $(OBJECTS1)

clean:
	rm $(OBJECTS1) $(EXE1)



