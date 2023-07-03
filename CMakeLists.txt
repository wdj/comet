#==============================================================================
#
# CMakeLists.txt
#
#==============================================================================

cmake_minimum_required(VERSION 2.8.12)
#cmake_policy(SET CMP0004 OLD)
#cmake_policy(SET CMP0074 NEW)
string(STRIP "${CMAKE_EXE_LINKER_FLAGS}" CMAKE_EXE_LINKER_FLAGS)

if(USE_HIPINTEL)
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
endif()

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
  src/assertions.cc
  src/utils.cc
  src/env.cc
  src/magma_wrapper.cc
  src/mirrored_buf.cc
  src/compressed_buf.cc
  src/histograms.cc
  src/tc.cc
  src/decomp_mgr.cc
  src/linalg.cc
  src/vectors.cc
  src/metrics.cc
  src/vector_sums.cc
  src/checksum.cc
  src/comm_xfer_utils.cc
  src/compute_metrics_2way_block_nums.cc
  src/compute_metrics_2way_block_finalize.cc
  src/compute_metrics_3way_block_nonlinalg.cc
  src/compute_metrics_3way_block_linalg.cc
  src/compute_metrics_3way_block.cc
  src/compute_metrics_2way.cc
  src/compute_metrics_3way.cc
  src/compute_metrics.cc
)

SET(ACCEL_SOURCES
  src/tc.cc
  src/compressed_buf.cc
  src/vector_sums.cc
)

SET(LIB_INCLUDES
  src/assertions.hh
  src/assertions.i.hh
  src/utils.hh
  src/types.hh
  src/env.hh
  src/formulas.hh
  src/magma_wrapper.hh
  src/mirrored_buf.hh
  src/compressed_buf.hh
  src/compressed_buf.i.hh
  src/histograms.hh
  src/tc.hh
  src/tc_helpers.i.hh
  src/tc_in.i.hh
  src/tc_solve.i.hh
  src/tc_out.i.hh
  src/decomp_mgr.hh
  src/linalg.hh
  src/vectors.hh
  src/metrics.hh
  src/metrics_2way_indexing.i.hh
  src/metrics_2way_accessors.i.hh
  src/metrics_3way_indexing.i.hh
  src/metrics_3way_accessors.i.hh
  src/vector_sums.hh
  src/checksum.hh
  src/comm_xfer_utils.hh
  src/compute_metrics_2way_block.hh
  src/compute_metrics_3way_block.hh
  src/compute_metrics_2way.hh
  src/compute_metrics_3way.hh
  src/compute_metrics.cc
)

SET(DRIVER_SOURCES
  driver/driver.cc
  driver/vectors_io.cc
  driver/metrics_io.cc
  driver/test_problems.cc
)

SET(DRIVER_INCLUDES
  driver/driver.hh
  driver/vectors_io.hh
  driver/metrics_io.hh
  driver/test_problems.hh
)

SET(TEST_TARGETS
  env_test
  metrics_test
  driver_test
  perf_test
)

#--------------------
# Preliminaries.

if(USE_MPI)
  find_package(MPI REQUIRED)
endif()

if(USE_CUDA)
  find_package(CUDA REQUIRED)
  CUDA_INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  # see https://stackoverflow.com/questions/36021233/
  #       cmake-cuda-compile-cpp-files-in-cuda-mode-x-cu
  # Interpret cc file as cu file.
  # NOTE: order of next lines matters.
#  set_source_files_properties(src/tc.cc
#                              PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
#  set_source_files_properties(src/assertions.cc
#                              PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
  foreach(ACCEL_SOURCE ${ACCEL_SOURCES})
    set_source_files_properties(${ACCEL_SOURCE}
                                PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
  endforeach()
  CUDA_ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES} ${CUDA_SOURCES})
elseif(USE_HIPINTEL) #=========================================
  set(CMAKE_MODULE_PATH "$ENV{HIPCL_DIR}/cmake" ${CMAKE_MODULE_PATH})
  #set(CMAKE_MODULE_PATH "$ENV{HOME}/cmake" ${CMAKE_MODULE_PATH})
  #find_package(HIP)
  #set(HIPLOCAL_FOUND TRUE)
  #find_package(HIPLOCAL)
  #set(HIP_DIR "${HIPCL_DIR}/cmake")
  #set(CMAKE_HIP_COMPILER "${HIPCL_CLANG_DIR}/bin/clang++")
  find_package(HIP)
  if(HIP_FOUND)
     message(STATUS "Found HIP: " ${HIP_VERSION})
  else()
     message(FATAL_ERROR "Could not find HIP. Ensure that the variable HIP_DIR is set to point to the right location.")
  endif()
  #set(CMAKE_LINKER "clang++-link --gcc-toolchain=/usr/bin --target=x86_64-linux-gnu")
  #set(CMAKE_LINKER /soft/modulefiles/hipcl/20210506_wrapper//clang++-link)
  set_source_files_properties(${LIB_SOURCES} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
  #hip_add_library(genomics_metrics ${LIB_SOURCES} HIPCC_OPTIONS ${CMAKE_CXX_FLAGS})
  #hip_add_library(genomics_metrics STATIC ${LIB_SOURCES} HIPCC_OPTIONS ${CMAKE_CXX_FLAGS})
  set_source_files_properties(driver/genomics_metric.cc ${DRIVER_SOURCES} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
  #XXXset(HIP_HOST_FLAGS ${CMAKE_CXX_FLAGS})
elseif(USE_HIP) #=========================================
  set(HIPLOCAL_DIR "${CMAKE_SOURCE_DIR}/scripts")
  if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
      message(FATAL_ERROR "Could not find HIP..")
    else()
      set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to HIP install")
    endif()
  endif()
  set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
  find_package(HIPLOCAL QUIET)
  if(HIPLOCAL_FOUND)
    message(STATUS "Found HIP: " ${HIP_VERSION})
  else()
    message(FATAL_ERROR "Could not find HIP.")
  endif()
  set(CMAKE_LINKER "${ROCM_PATH}/bin/hipcc")
  set_source_files_properties(${LIB_SOURCES} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
  hip_add_library(genomics_metrics ${LIB_SOURCES} HIPCC_OPTIONS ${CMAKE_CXX_FLAGS})
  set_source_files_properties(driver/genomics_metric.cc ${DRIVER_SOURCES} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
  #set_property(TARGET tgt PROPERTY HIP_ARCHITECTURES gfx801 gfx900)
else()
  ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES} ${CUDA_SOURCES})
endif()


#--------------------
# Main executable, lib.

if(USE_HIPINTEL)
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS} ${CMAKE_EXTRA_INCLUDE_DIRS})
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
  ADD_EXECUTABLE(genomics_metric driver/genomics_metric.cc ${DRIVER_SOURCES} ${LIB_SOURCES})
  set_property(TARGET genomics_metric PROPERTY HIP_ARCHITECTURES gfx90a)
  TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})
elseif(USE_HIP)
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
  hip_add_executable(genomics_metric driver/genomics_metric.cc ${DRIVER_SOURCES} HIPCC_OPTIONS ${CMAKE_CXX_FLAGS})
  TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
  TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})
#elseif(USE_HIP)
#  hip_add_executable(genomics_metric driver/genomics_metric.cc ${DRIVER_SOURCES} ${LIB_SOURCES} HIPCC_OPTIONS ${CMAKE_CXX_FLAGS})
#  #hip_add_executable(genomics_metric driver/genomics_metric.cc ${DRIVER_SOURCES} HIPCC_OPTIONS ${CMAKE_CXX_FLAGS})
#  #set_property(TARGET genomics_metric PROPERTY HIP_ARCHITECTURES ${COMET_HIP_ARCHITECTURES})
#  TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
#  TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})
else()
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
  ADD_EXECUTABLE(genomics_metric driver/genomics_metric.cc ${DRIVER_SOURCES})
  TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
  TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})
endif()

install(FILES ${LIB_INCLUDES} DESTINATION include)
if(USE_HIPINTEL)
  install(FILES ${DRIVER_INCLUDES} DESTINATION include)
else()
  install(TARGETS genomics_metrics DESTINATION lib)
endif()
install(TARGETS genomics_metric DESTINATION bin)

#--------------------
# Tester.

if(TESTING) 
  enable_testing()

  foreach(TEST_TARGET ${TEST_TARGETS})
    if(USE_HIPINTEL)
      ADD_EXECUTABLE(${TEST_TARGET} testing/${TEST_TARGET}.cc ${DRIVER_SOURCES} ${LIB_SOURCES} ${ACCEL_SOURCES})
    else()
      ADD_EXECUTABLE(${TEST_TARGET} testing/${TEST_TARGET}.cc ${DRIVER_SOURCES})
      TARGET_LINK_LIBRARIES(${TEST_TARGET} genomics_metrics)
    endif()
    TARGET_LINK_LIBRARIES(${TEST_TARGET} ${CMAKE_EXE_LINKER_FLAGS})

    if("${TEST_TARGET}" STREQUAL "perf_test")
      if(NOT "${TEST_COMMAND_PERF}" STREQUAL "")
      if("${FP_PRECISION}" STREQUAL "DOUBLE")
        add_test(${TEST_TARGET} bash -c
          "date && cd ${CMAKE_INSTALL_PREFIX} && ${TEST_COMMAND_PERF} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET} && date")
      endif()
      endif()
    else()
      add_test(${TEST_TARGET} bash -c
        "date && cd ${CMAKE_INSTALL_PREFIX} && ${TEST_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET} && date")
      set(UNUSED_NO_WARN_1, ${FP_PRECISION})
    endif()

  endforeach()

else()
  # turn off unused var warnings
  set(UNUSED_NO_WARN_1, ${TEST_COMMAND})
  set(UNUSED_NO_WARN_1, ${TEST_COMMAND_PERF})
  set(UNUSED_NO_WARN_1, ${FP_PRECISION})
endif() 

MESSAGE("--------------------")

#==============================================================================
