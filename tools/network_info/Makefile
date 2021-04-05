CC:=mpicc
CFLAGS:=-fopenmp

.PHONY: all
all: omp_id

omp_id: omp_id.c
	$(CC) $(CFLAGS) $< -o $@

