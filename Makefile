PREFIX=$(HOME)
CC=g++
AR=ar
CFLAGS= -Wall -O3 -fopenmp
LDFLAGS= -L. -lfhew -lfftw3 
INCLUDE=distrib.h LWE.h FHEW.h FFT.h params.h

all: libfhew.a cmd profiler_test

cmd: cmd/gen cmd/enc cmd/nand cmd/dec

install: $(INCLUDE) libfhew.a
	install $(INCLUDE) $(PREFIX)/include
	install libfhew.a $(PREFIX)/lib

uninstall:
	rm $(PREFIX)/lib/libfhew.a
	rm $(PREFIX)/lib/{distrib,LWE,FHEW,FFT,params}.h

clean:
	rm *.o libfhew.a profiler_test cmd/gen cmd/enc cmd/dec cmd/nand || echo nothing to clean

libfhew.a: distrib.o FFT.o LWE.o FHEW.o time_profiler.o
	$(AR) -q libfhew.a distrib.o FFT.o LWE.o FHEW.o time_profiler.o

distrib.o: distrib.cpp distrib.h params.h
	$(CC) $(CFLAGS) -c distrib.cpp

FFT.o: FFT.h FFT.cpp params.h FHEW.h
	$(CC) $(CFLAGS) -c FFT.cpp

LWE.o: LWE.h LWE.cpp FFT.h params.h distrib.h
	$(CC) $(CFLAGS) -c LWE.cpp

FHEW.o: FHEW.h FHEW.cpp FFT.h LWE.h params.h
	$(CC) $(CFLAGS) -c FHEW.cpp

time_profiler.o: time_profiler.h time_profiler.cpp
	$(CC) $(CFLAGS) -c time_profiler.cpp

profiler_test: profiler_test.cpp libfhew.a
	$(CC) $(CFLAGS) -o profiler_test profiler_test.cpp $(LDFLAGS)

common.o: cmd/common.cpp cmd/common.h libfhew.a
	$(CC) $(CFLAGS) -c cmd/common.cpp 

cmd/gen: cmd/gen.cpp common.o libfhew.a
	$(CC) $(CFLAGS) -o cmd/gen cmd/gen.cpp common.o $(LDFLAGS)

cmd/enc: cmd/enc.cpp common.o libfhew.a
	$(CC) $(CFLAGS) -o cmd/enc cmd/enc.cpp common.o $(LDFLAGS)

cmd/nand: cmd/nand.cpp common.o libfhew.a
	$(CC) $(CFLAGS) -o cmd/nand cmd/nand.cpp common.o $(LDFLAGS)

cmd/dec: cmd/dec.cpp common.o libfhew.a
	$(CC) $(CFLAGS) -o cmd/dec cmd/dec.cpp common.o $(LDFLAGS)
