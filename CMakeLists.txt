#==============================================================================
#
# CMakeLists.txt
#
#==============================================================================

cmake_minimum_required(VERSION 2.8)

project("genomics_gpu")

SET(INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src
  )

SET(LIB_SOURCES
  src/env.c
  src/vectors.c
  src/vector_sums.c
  src/metrics.c
  src/compute_metrics.c
  src/compute_metrics_utils.c
  src/compute_metrics_ccc.c
  src/compute_metrics_czekanowski_2way.c
  src/compute_metrics_czekanowski_3way.c
  src/compute_metrics_sorenson.c
  )

SET(LIB_INCLUDES
  src/env.h
  src/vectors.h
  src/vector_sums.h
  src/metrics.h
  src/compute_metrics.h
  src/compute_metrics_utils.h
  src/compute_metrics_ccc.h
  src/compute_metrics_czekanowski_2way.h
  src/compute_metrics_czekanowski_3way.h
  src/compute_metrics_sorenson.h
  )

find_package(MPI REQUIRED)
find_package(CUDA REQUIRED)

ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES})
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

ADD_EXECUTABLE(genomics_metric driver/genomics_metric.c)
TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})

SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)

install(TARGETS genomics_metric DESTINATION bin)
install(TARGETS genomics_metrics DESTINATION lib)
install(FILES ${LIB_INCLUDES} DESTINATION include)

MESSAGE("Compiler ID: ${CMAKE_C_COMPILER_ID}")

#==============================================================================
