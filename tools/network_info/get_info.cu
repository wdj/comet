// MPI include
#include <mpi.h>

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <sched.h>
#include <omp.h>
#include <sys/syscall.h>

//#include <cuda_runtime.h>

using std::cout;
using std::cerr;
using std::endl;

// Error handling macros
#define MPI_CHECK(call) \
    if((call) != MPI_SUCCESS) { \
        cerr << "MPI error calling \""#call"\"\n"; \
        exit(-1); \
    }

int get_cpu_id(int rank)
{
    /* Get the the current process' stat file from the proc filesystem */
    FILE* procfile = fopen("/proc/self/stat", "r");
    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof(char), to_read, procfile);
    fclose(procfile);

    printf("rank=%d stat=%s\n",rank,buffer);

    // Field with index 38 (zero-based counting) is the one we want
    char* line = strtok(buffer, " ");
    for (int i = 1; i < 38; i++)
    {
        line = strtok(NULL, " ");
    }

    line = strtok(NULL, " ");
    int cpu_id = atoi(line);
    return cpu_id;
}

inline int _ConvertSMVer2Cores(int major, int minor) {
  // Defines for GPU Architecture types (using the SM version to determine
  // the # of cores per SM
  typedef struct {
    int SM;  // 0xMm (hexidecimal notation), M = SM Major version,
    // and m = SM minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = {
      {0x30, 192},
      {0x32, 192},
      {0x35, 192},
      {0x37, 192},
      {0x50, 128},
      {0x52, 128},
      {0x53, 128},
      {0x60,  64},
      {0x61, 128},
      {0x62, 128},
      {0x70,  64},
      {0x72,  64},
      {0x75,  64},
      {0x80,  64},
      {-1, -1}};

  int index = 0;

  while (nGpuArchCoresPerSM[index].SM != -1) {
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchCoresPerSM[index].Cores;
    }

    index++;
  }

  // If we don't find the values, we default use the previous one
  // to run properly
  printf(
      "MapSMtoCores for SM %d.%d is undefined."
      "  Default to use %d Cores/SM\n",
      major, minor, nGpuArchCoresPerSM[index - 1].Cores);
  return nGpuArchCoresPerSM[index - 1].Cores;
}

void printGPUDetails(int dev, cudaDeviceProp deviceProp) {
  int driverVersion = 0, runtimeVersion = 0;

  printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

  // Console log
  cudaDriverGetVersion(&driverVersion);
  cudaRuntimeGetVersion(&runtimeVersion);
  printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
         driverVersion / 1000, (driverVersion % 100) / 10,
         runtimeVersion / 1000, (runtimeVersion % 100) / 10);
  printf("  CUDA Capability Major/Minor version number:    %d.%d\n",
         deviceProp.major, deviceProp.minor);

  printf("  Total amount of global memory:                 %.0f MBytes "
         "(%llu bytes)\n",
         static_cast<float>(deviceProp.totalGlobalMem / 1048576.0f),
         (unsigned long long)deviceProp.totalGlobalMem);

  printf("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
         deviceProp.multiProcessorCount,
         _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
         _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
         deviceProp.multiProcessorCount);
  printf("  GPU Max Clock rate:                            %.0f MHz (%0.2f "
         "GHz)\n",
         deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);

  printf("  Memory Clock rate:                             %.0f Mhz\n",
         deviceProp.memoryClockRate * 1e-3f);
  printf("  Memory Bus Width:                              %d-bit\n",
         deviceProp.memoryBusWidth);

  if (deviceProp.l2CacheSize) {
    printf("  L2 Cache Size:                                 %d bytes\n",
           deviceProp.l2CacheSize);
  }

  printf("  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, "
         "%d), 3D=(%d, %d, %d)\n",
         deviceProp.maxTexture1D, deviceProp.maxTexture2D[0],
         deviceProp.maxTexture2D[1], deviceProp.maxTexture3D[0],
         deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
  printf("  Maximum Layered 1D Texture Size, (num) layers  1D=(%d), %d layers\n",
         deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1]);
  printf("  Maximum Layered 2D Texture Size, (num) layers  2D=(%d, %d), %d "
         "layers\n",
         deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1],
         deviceProp.maxTexture2DLayered[2]);

  printf("  Total amount of constant memory:               %zu bytes\n",
         deviceProp.totalConstMem);
  printf("  Total amount of shared memory per block:       %zu bytes\n",
         deviceProp.sharedMemPerBlock);
  printf("  Total number of registers available per block: %d\n",
         deviceProp.regsPerBlock);
  printf("  Warp size:                                     %d\n",
         deviceProp.warpSize);
  printf("  Maximum number of threads per multiprocessor:  %d\n",
         deviceProp.maxThreadsPerMultiProcessor);
  printf("  Maximum number of threads per block:           %d\n",
         deviceProp.maxThreadsPerBlock);
  printf("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
         deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
         deviceProp.maxThreadsDim[2]);
  printf("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
         deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
         deviceProp.maxGridSize[2]);
  printf("  Maximum memory pitch:                          %zu bytes\n",
         deviceProp.memPitch);
  printf("  Texture alignment:                             %zu bytes\n",
         deviceProp.textureAlignment);
  printf("  Concurrent copy and kernel execution:          %s with %d copy "
         "engine(s)\n",
         (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
  printf("  Run time limit on kernels:                     %s\n",
         deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
  printf("  Integrated GPU sharing Host Memory:            %s\n",
         deviceProp.integrated ? "Yes" : "No");
  printf("  Support host page-locked memory mapping:       %s\n",
         deviceProp.canMapHostMemory ? "Yes" : "No");
  printf("  Alignment requirement for Surfaces:            %s\n",
         deviceProp.surfaceAlignment ? "Yes" : "No");
  printf("  Device has ECC support:                        %s\n",
         deviceProp.ECCEnabled ? "Enabled" : "Disabled");

  printf("  Device supports Unified Addressing (UVA):      %s\n",
         deviceProp.unifiedAddressing ? "Yes" : "No");
  printf("  Device supports Managed Memory:                %s\n",
         deviceProp.managedMemory ? "Yes" : "No");
  printf("  Device supports Compute Preemption:            %s\n",
         deviceProp.computePreemptionSupported ? "Yes" : "No");
  printf("  Supports Cooperative Kernel Launch:            %s\n",
         deviceProp.cooperativeLaunch ? "Yes" : "No");
  printf("  Supports MultiDevice Co-op Kernel Launch:      %s\n",
         deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
  printf("  Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n",
         deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);

  const char *sComputeMode[] = {
      "Default (multiple host threads can use ::cudaSetDevice() with device "
      "simultaneously)",
      "Exclusive (only one host thread in one process is able to use "
      "::cudaSetDevice() with this device)",
      "Prohibited (no host thread can use ::cudaSetDevice() with this "
      "device)",
      "Exclusive Process (many threads in one process is able to use "
      "::cudaSetDevice() with this device)",
      "Unknown",
      NULL};
  printf("  Compute Mode:\n");
  printf("     < %s >\n", sComputeMode[deviceProp.computeMode]);
}

/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d,", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  *ptr = 0;
  return(str);
}

// Host code
int main(int argc, char *argv[])
{
  // Initialize MPI state
  MPI_CHECK(MPI_Init(&argc, &argv));

  // Get our MPI node number and node count
  int commSize, commRank;
  MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &commSize));
  MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &commRank));

  // Get core number
  int core_num = get_cpu_id(commRank);
  MPI_Barrier(MPI_COMM_WORLD);

  // Get core information
  int thread;
  cpu_set_t coremask;
  char clbuf[7 * CPU_SETSIZE], hnbuf[64];
  memset(clbuf, 0, sizeof(clbuf));
  memset(hnbuf, 0, sizeof(hnbuf));
  (void)gethostname(hnbuf, sizeof(hnbuf));
  #pragma omp parallel private(thread, coremask, clbuf)
  {
    thread = omp_get_thread_num();
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);
    #pragma omp barrier
    printf("Hello from rank %d, thread %d, on %s. (core affinity = %s)\n",
            commRank, thread, hnbuf, clbuf);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int processor_name_len;
  MPI_Get_processor_name(processor_name, &processor_name_len);

  printf("rank=%d/%d core=%d proc=%s\n",commRank,commSize,core_num,processor_name);
  MPI_Barrier(MPI_COMM_WORLD);

  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);

  // This function call returns 0 if there are no CUDA capable devices.
  printf("rank=%d core=%d Detected %d CUDA Capable device(s)\n",
    commRank,core_num,deviceCount);

  int print_all = false;
  int dev;
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    printf("rank=%d core=%d Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n",
      commRank,core_num,deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);
    if(print_all) printGPUDetails(dev,deviceProp);
  }

  if(print_all) {
    // *****************************
    // exe and CUDA driver name
    printf("\n");
    std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
    char cTemp[16];

    // driver version
    int driverVersion = 0;
    cudaDriverGetVersion(&driverVersion);
    sProfileString += ", CUDA Driver Version = ";
    snprintf(cTemp, sizeof(cTemp), "%d.%d", driverVersion / 1000,
             (driverVersion % 100) / 10);
    sProfileString += cTemp;

    // Runtime version
    int runtimeVersion = 0;
    cudaRuntimeGetVersion(&runtimeVersion);
    sProfileString += ", CUDA Runtime Version = ";
    snprintf(cTemp, sizeof(cTemp), "%d.%d", runtimeVersion / 1000,
             (runtimeVersion % 100) / 10);
    sProfileString += cTemp;

    // Device count
    sProfileString += ", NumDevs = ";
    snprintf(cTemp, sizeof(cTemp), "%d", deviceCount);
    sProfileString += cTemp;
    sProfileString += "\n";
    printf("%s", sProfileString.c_str());
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_CHECK(MPI_Finalize());
  return 0;
}

