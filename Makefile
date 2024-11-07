CPP=CC
CFLAGS=-lm
COPTFLAGS=-O3 -ffast-math -fopenmp -flto -march=native -funroll-loops -fprefetch-loop-arrays
MPIFLAGS=-DMPI
DEBUGFLAGS=-g -pg

NVCC=nvcc
NVCCFLAGS=-DCUDA
# NVCCFLAGS=-DCUDA -O3 -dlto -use_fast_math -Xptxas -dlcm=ca -Xcompiler "-funroll-loops"

PYTHON=python3

all: mpi gpu basic_serial serial

mpi: build/mpi
gpu: build/gpu
serial: build/serial
basic_serial: build/basic_serial

build/mpi: common/main.cpp common/scenarios.cpp mpi/mpi.cpp
	$(CPP) $^ -o $@ $(MPIFLAGS) $(CFLAGS) $(OPTFLAGS)

build/gpu: common/main.cpp common/scenarios.cpp gpu/gpu.cu
	$(NVCC) $^ -o $@ $(NVCCFLAGS)

build/serial: common/main.cpp common/scenarios.cpp serial/serial.cpp
	$(CPP) $^ -o $@ $(CFLAGS) $(COPTFLAGS)

build/basic_serial: common/main.cpp common/scenarios.cpp serial/basic_serial.cpp
	$(CPP) $^ -o $@ $(CFLAGS) $(COPTFLAGS)

.PHONY: clean

clean:
	rm -f build/*.out
	rm -f build/*.o
	rm -f build/*.gif