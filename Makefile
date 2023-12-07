TARGET = BARNES
OBJS = code.o code_io.o load.o grav.o getparam.o util.o cha.o topology.o

#CC := gcc
CC := g++
#CC := clang 
#CFLAGS := -O -pthread -D_POSIX_C_SOURCE=200112 -static -integrated-as -msoft-float
#CFLAGS := -O -pthread -D_POSIX_C_SOURCE=200112 -static -msoft-float
#CFLAGS := -O -pthread -D_POSIX_C_SOURCE=200112 -integrated-as -msoft-float
#CFLAGS := -O -pthread -D_POSIX_C_SOURCE=200112 -msoft-float
CFLAGS := -O -pthread -D_POSIX_C_SOURCE=200112
#CFLAGS := -O3 -pthread -D_POSIX_C_SOURCE=200112
#CFLAGS := -g3 -pthread -D_POSIX_C_SOURCE=200112
CFLAGS := $(CFLAGS) -Wall -W -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wdisabled-optimization
CFLAGS := $(CFLAGS) -Wpadded -Winline -Wpointer-arith -Wsign-compare -Wendif-labels
LDFLAGS := -lm
#LDFLAGS := -lm -static

MACROS := ./c.m4.null
M4 := m4 -s -Ulen -Uindex

x = *

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -rf *.c *.h *.o *.cpp $(TARGET)

.SUFFIXES:
.SUFFIXES:	.o .c .C .h .H

.H.h:
	$(M4) $(MACROS) $*.H > $*.h

.C.c:
	$(M4) $(MACROS) $*.C > $*.cpp

.c.o:
	$(CC) -c $(CFLAGS) $*.cpp

.C.o:
	$(M4) $(MACROS) $*.C > $*.cpp
	$(CC) -c $(CFLAGS) $*.cpp

stdinc.h: code.h defs.h util.h vectmath.h load.h code_io.h grav.h getparam.h stdinc.H 
code.o: code.C stdinc.h cha.h
code_io.o: code_io.C stdinc.h
getparam.o: getparam.C stdinc.h
grav.o: grav.C stdinc.h
load.o: load.C stdinc.h
util.o: util.C stdinc.h
cha.o:	cha.C cha.h
topology.o:  topology.C topology.h tile.h
