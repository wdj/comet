all: get_info

CC = g++
NVCC = nvcc
#-ccbin
EXE = get_info.exe

CFLAGS = -I/usr/common/software/sles15_cgpu/openmpi/4.0.3/gcc/include 
#-fopenmp
#-I. -I./inc -std=c++11
#CUDA_FLAGS = -m64 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_75,code=compute_75
#CUDA_FLAGS = -m64 -maxrregcount=255 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_75,code=compute_75
#CUDA_FLAGS = -m64 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_75,code=compute_75
#CUDA_FLAGS = -m64 -gencode arch=compute_75,code=compute_75
CUDA_FLAGS = -m64 -gencode arch=compute_80,code=compute_80 -Xcompiler -fopenmp
OBJS = get_info.o 
LIBS = -L/usr/common/software/sles15_cgpu/openmpi/4.0.3/gcc/lib -lmpi 
#-L/usr/lib64 -lcublas_static -lcublasLt_static -lculibos ${MAGMA_LIB}

%.o : %.cu
	$(NVCC) ${CFLAGS} ${CUDA_FLAGS} -o $@ -c $<

%.o : %.cpp
	$(CC) ${CFLAGS} -o $@ -c $<

get_info : ${OBJS}
	${NVCC} ${CFLAGS} ${CUDA_FLAGS} ${OBJS} -o ${EXE} ${LIBS}

clean:
	rm ${EXE} *.o

