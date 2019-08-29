#==============================================================================
#
# CMakeLists.txt
#
#==============================================================================

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0074 NEW)

project("genomics_gpu")

MESSAGE("--------------------")
MESSAGE("Compiler ID: ${CMAKE_C_COMPILER_ID}")

#--------------------
# Dirs, Files.

SET(INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/driver
  )

SET(LIB_SOURCES
  src/env.cc
  src/assertions.cc
  src/decomp_mgr.cc
  src/mirrored_buf.cc
  src/linalg.cc
  src/linalg_accel.cc
  src/metrics.cc
  src/checksum.cc
  src/vectors.cc
  src/vector_sums.cc
  src/comm_xfer_utils.cc
  src/compute_metrics_2way_block_nums.cc
  src/compute_metrics_2way_block_combine.cc
  src/compute_metrics_3way_block_nongpu.cc
  src/compute_metrics_3way_block_gpu.cc
  src/compute_metrics_3way_block.cc
  src/compute_metrics_2way.cc
  src/compute_metrics_3way.cc
  src/compute_metrics.cc
  )

SET(LIB_INCLUDES
  src/env.hh
  src/assertions.hh
  src/types.hh
  src/decomp_mgr.hh
  src/mirrored_buf.hh
  src/linalg.hh
  src/linalg_accel.hh
  src/metrics.hh
  src/metrics_2way_indexing.hh
  src/metrics_2way_accessors.hh
  src/metrics_3way_indexing.hh
  src/metrics_3way_accessors.hh
  src/checksum.hh
  src/vectors.hh
  src/vector_sums.hh
  src/comm_xfer_utils.hh
  src/compute_metrics_2way_block_nums.hh
  src/compute_metrics_2way_block_combine.hh
  src/compute_metrics_3way_block_nongpu.hh
  src/compute_metrics_3way_block_gpu.hh
  src/compute_metrics_3way_block.hh
  src/compute_metrics_2way.hh
  src/compute_metrics_3way.hh
  src/compute_metrics.cc
  )

SET(DRIVER_SOURCES
  driver/driver.cc
  driver/test_problems.cc
  driver/input_output.cc
  )

#--------------------
# Preliminaries.

if(NOT NOMPI)
  find_package(MPI REQUIRED)
endif()

if(NOT NOCUDA)
  find_package(CUDA REQUIRED)
  CUDA_INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  # see https://stackoverflow.com/questions/36021233/
  #       cmake-cuda-compile-cpp-files-in-cuda-mode-x-cu
  # Interpret cc file as cu file.
  # NOTE: order of next two lines matters.
  set_source_files_properties(src/linalg_accel.cc
                              PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
  CUDA_ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES} ${CUDA_SOURCES})
else()
  ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES} ${CUDA_SOURCES})
endif()

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)

#--------------------
# Main executable, lib.

ADD_EXECUTABLE(genomics_metric driver/genomics_metric.cc ${DRIVER_SOURCES})
TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})

if(TESTING) 
  TARGET_LINK_LIBRARIES(genomics_metric ${GTEST_LD_FLAGS})
endif() 

install(FILES ${LIB_INCLUDES} DESTINATION include)
install(TARGETS genomics_metrics DESTINATION lib)
install(TARGETS genomics_metric DESTINATION bin)

#--------------------
# Tester.

if(TESTING) 
  enable_testing()

  foreach(TEST_TARGET driver_test metrics_test)
    ADD_EXECUTABLE(${TEST_TARGET} testing/${TEST_TARGET}.cc ${DRIVER_SOURCES})
    TARGET_LINK_LIBRARIES(${TEST_TARGET} genomics_metrics)
    TARGET_LINK_LIBRARIES(${TEST_TARGET} ${CMAKE_EXE_LINKER_FLAGS})
    TARGET_LINK_LIBRARIES(${TEST_TARGET} ${GTEST_LD_FLAGS})

    add_test(${TEST_TARGET} bash -c
             "cd ${CMAKE_INSTALL_PREFIX}; \
              ${TEST_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET}")
  endforeach()
endif() 

MESSAGE("--------------------")

#==============================================================================
