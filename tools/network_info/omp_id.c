#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)

#include <stdio.h>
#include <mpi.h>
#include <sched.h>
#include <omp.h>

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);
    
    fprintf(stderr, "LEGEND: 'MPI' = MPI Rank; 'CPU' = CPU Core ID, 'MCPU' = Master CPU Core ID; 'OMP' = OpenMP Thread ID, 'MOMP' = Master OpenMP Thread ID; 'F' = Fork\n");

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int thread_num_master = omp_get_thread_num();
    int cpu_num_master = sched_getcpu();
    printf("MPI\tF\tMOMP\tOMP\tMCPU\tCPU\n");
    printf("MPI PROCESS: %3d\t•\t%4d\t%3d\t%4d\t%3d\n", world_rank, thread_num_master, thread_num_master, cpu_num_master ,cpu_num_master);
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        int cpu_num = sched_getcpu();
    	printf("OMP_THREAD %3d\t└\t%4d\t%3d\t%4d\t%3d\n", world_rank, thread_num_master, thread_num, cpu_num_master ,cpu_num);
    }

    MPI_Finalize();
    return 0;
}

