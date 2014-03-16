#
# Computers with Red Hat Enterprise Linux 5 in the computer room 648, KTH Forum, Kista
#
#CC = g++44
CC = g++
CC_OPENMP = gcc-4.9 -lstdc++
MPCC =  mpicc -cc=g++44
OPENMP = -fopenmp
LIBS = -lm
CFLAGS = -O3 -std=c++11

TARGETS = serial pthreads openmp mpi

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
openmp: openmp.o common.o
	$(CC_OPENMP) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

old_serial: old_serial.o common.o
	$(CC) -o $@ $(LIBS) old_serial.o common.o
old_openmp: old_openmp.o common.o
	$(CC_OPENMP) -o $@ $(LIBS) $(OPENMP) old_openmp.o common.o
old_pthreads: old_pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread old_pthreads.o common.o

openmp.o: openmp.cpp common.h
	$(CC_OPENMP) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

old_serial.o: old_serial.cpp common.h
	$(CC) -c $(CFLAGS) old_serial.cpp
old_openmp.o: old_openmp.cpp common.h
	$(CC_OPENMP) -c $(OPENMP) $(CFLAGS) old_openmp.cpp
old_pthreads.o: old_pthreads.cpp common.h
	$(CC) -c $(CFLAGS) old_pthreads.cpp

clean:
	rm -f *.o $(TARGETS)
