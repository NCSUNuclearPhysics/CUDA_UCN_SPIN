
# CURRENTLY WRITING WRAPPER FUNCTION TEST SUFFIX INVOCATION IN UCN_MAIN.cpp
# Build tools

SHELL = /bin/sh

NVCC = nvcc
CXX = g++
CC = gcc

#  O0 -Xcicc -Xptxas -maxrregcount 0
NVCCFLAGS   = -lm -arch=sm_35 -lcuda -lcudadevrt -lcudart -w --relocatable-device-code true -Xcicc -Xptxas -O0 --maxrregcount 0
CXXFLAGS    =  -w
CFLAGS      =  -w
# nvcc CUDA_NEUTRON_TRANSPORT_SDM_MOD_1_0_0.cu -lm -G -g -O0 -Xcicc -Xptxas -maxrregcount 0 -arch=sm_35 -o RKIE_FAST -wnwfaeeeee 


# here are all the objects
#GPUOBJS = cuexample.o  CUDA_BFIELD_GRADIENT.o CUDA_COORD.o CUDA_RKDK_ODE_SPIN.o CUDA_RKDK_ODE_XV.o CUDA_UCN_MISC
#OBJS = cppexample.o

# make and compile


CUDA_RK45_SPIN.o: CUDA_RK45_SPIN.cu CUDA_RK45_SPIN.cuh UCN_CUDA_ALL_KERNEL.cuh UCN_ENUM.h
	$(NVCC) -c CUDA_RK45_SPIN.cu -o CUDA_RK45_SPIN.o $(NVCCFLAGS)

UCN_CUDA_ALL_KERNEL.o: UCN_CUDA_ALL_KERNEL.cu UCN_CUDA_ALL_KERNEL.cuh UCN_ENUM.h
	$(NVCC) -c UCN_CUDA_ALL_KERNEL.cu -o UCN_CUDA_ALL_KERNEL.o $(NVCCFLAGS)
	
UCN_CUDA_WRAPPER.o: UCN_CUDA_WRAPPER.cu UCN_CUDA_WRAPPER.h UCN_CUDA_ALL_KERNEL.cuh CUDA_RK45_SPIN.cuh UCN_ENUM.h
	$(NVCC) -c UCN_CUDA_WRAPPER.cu -o UCN_CUDA_WRAPPER.o $(NVCCFLAGS)

UCN_RUN.o: UCN_RUN.cpp UCN_RUN.h UCN_CUDA_WRAPPER.h UCN_ENUM.h
	$(CXX) -c UCN_RUN.cpp $(CXXFLAGS) -o UCN_RUN.o

UCN_MAIN.o: UCN_MAIN.cpp UCN_RUN.h UCN_ENUM.h
	$(CXX) -c UCN_MAIN.cpp $(CXXFLAGS) -o UCN_MAIN.o

CNT: CUDA_RK45_SPIN.o UCN_CUDA_ALL_KERNEL.o UCN_CUDA_WRAPPER.o UCN_RUN.o UCN_MAIN.o
	$(NVCC) CUDA_RK45_SPIN.o UCN_CUDA_ALL_KERNEL.o UCN_CUDA_WRAPPER.o UCN_RUN.o UCN_MAIN.o -o CNT $(NVCCFLAGS)


debug: NVCCFLAGS += -G -g 
debug: CXXFLAGS += -pg -g
debug: CFLAGS += -pg -g
debug: CNT
run: CNT
	./CNT
clean: 
	rm -rf *.o CNT
