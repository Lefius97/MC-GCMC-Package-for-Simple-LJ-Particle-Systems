
TARGET	 = GNM

OBJS     = main.o input.o output.o operate.o energy.o pair.o neigh.o

HEADERS  = GNM.h nr3.h ran.h

CC       = g++
CFLAGS   = -O3 -Wall
LDFLAGS  =
MKFILE   = Makefile

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	@-rm -f $(OBJS)

# dependencies
main.o: main.cpp GNM.h nr3.h ran.h
input.o: input.cpp GNM.h nr3.h ran.h
output.o: output.cpp GNM.h nr3.h ran.h
operate.o: operate.cpp GNM.h nr3.h ran.h
energy.o: energy.cpp GNM.h nr3.h ran.h
pair.o: pair.cpp GNM.h nr3.h ran.h
neigh.o: neigh.cpp GNM.h nr3.h ran.h
