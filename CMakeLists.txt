#==============================================================================
#
# CMakeLists.txt
#
#==============================================================================

cmake_minimum_required(VERSION 2.8)

project("genomics_gpu")

#--------------------

SET(INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/driver
  )

SET(LIB_SOURCES
  src/env.cc
  src/vectors.cc
  src/vector_sums.cc
  src/metrics.cc
  src/compute_metrics.cc
  src/compute_metrics_utils.cc
  src/compute_metrics_ccc.cc
  src/compute_metrics_czekanowski_2way.cc
  src/compute_metrics_czekanowski_3way.cc
  src/compute_metrics_sorenson.cc
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

#--------------------

if(TESTING) 
  set(C_AND_CXX_FLAGS "${C_AND_CXX_FLAGS} -isystem ${GTEST_DIR}/include -pthread -DTESTING")
endif(TESTING) 

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_AND_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_AND_CXX_FLAGS}")

find_package(MPI REQUIRED)
find_package(CUDA REQUIRED)

ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES})
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

ADD_EXECUTABLE(genomics_metric driver/genomics_metric.cc)
TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})
if(TESTING) 
  TARGET_LINK_LIBRARIES(genomics_metric "-L${GTEST_DIR}/lib -lgtest")
endif(TESTING) 

SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)

install(TARGETS genomics_metric DESTINATION bin)
install(TARGETS genomics_metrics DESTINATION lib)
install(FILES ${LIB_INCLUDES} DESTINATION include)

MESSAGE("Compiler ID: ${CMAKE_C_COMPILER_ID}")

#--------------------

if(TESTING) 

  enable_testing()

  foreach(test_target system_test)
    ADD_EXECUTABLE(${test_target} testing/${test_target}.cc)
    TARGET_LINK_LIBRARIES(${test_target} genomics_metrics)
    TARGET_LINK_LIBRARIES(${test_target} ${CMAKE_EXE_LINKER_FLAGS})
    TARGET_LINK_LIBRARIES(${test_target} "-L${GTEST_DIR}/lib -lgtest")

    add_test(${test_target} aprun -n6 -N1 ${CMAKE_CURRENT_BINARY_DIR}/${test_target})
  endforeach(test_target)
endif(TESTING) 

#==============================================================================
