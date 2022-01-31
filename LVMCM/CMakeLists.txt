project(LVMCM)

# Set cmake and compiler options
cmake_minimum_required(VERSION 3.5.1)

# Pick the correct C++ standard option
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  set(CMAKE_CXX_STANDARD 11)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# if(NOT CMAKE_BUILD_TYPE) set(CMAKE_BUILD_TYPE Release) endif()

# Find libraries Find and include Armadillo
find_package(Armadillo REQUIRED)
if(Armadillo_FOUND)
  include_directories(${Armadillo_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Could not find Armadillo!")
endif()

# Find and include Boost
#set(BOOST_ROOT /usr/lib/x86_64-linux-gnu)
find_package(Boost COMPONENTS system filesystem REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  add_definitions("-DHAS_BOOST")
else()
  message(FATAL_ERROR "Could not find Boost!")
endif()

# Find and include Sundials nvecserial, kinsol, cvode, cvodes,
find_library(SUNDIALS_NVECSERIAL_LIB sundials_nvecserial)
if(NOT SUNDIALS_NVECSERIAL_LIB)
  message(FATAL_ERROR "Could not find sundials_nvecserial!")
endif()

find_library(SUNDIALS_KINSOL_LIB sundials_kinsol)
if(NOT SUNDIALS_KINSOL_LIB)
  message(FATAL_ERROR "Could not find sundials_kinsol!")
endif()

find_library(SUNDIALS_IDA_LIB sundials_ida)
if(NOT SUNDIALS_IDA_LIB)
  message(FATAL_ERROR "Could not find sundials_ida!")
endif()

find_library(SUNDIALS_CVODE_LIB sundials_cvode)
if(NOT SUNDIALS_CVODE_LIB)
  message(FATAL_ERROR "Could not find sundials_cvode!")
endif()

find_library(SUNDIALS_CVODES_LIB sundials_cvodes)
if(NOT SUNDIALS_CVODES_LIB)
  message(FATAL_ERROR "Could not find sundials_cvodes!")
endif()

find_library(SUNDIALS_SUNLINSOLSPGMR_LIB sundials_sunlinsolspgmr)
if(NOT SUNDIALS_SUNLINSOLSPGMR_LIB)
  message(FATAL_ERROR "Could not find sundials_sunlinsolspgmr!")
else()
  message(STATUS "Found sundials_sunlinsolspgmr: ${SUNDIALS_SUNLINSOLSPGMR_LIB}")
endif()

if (NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  # find_library(OPENBLAS_LIB openblasp) # might be better sometimes
  find_library(OPENBLAS_LIB openblas)
  if(NOT OPENBLAS_LIB)
    message(FATAL_ERROR "Could not find openblas!")
  endif()
endif()

set(SUNDIALS_LIBRARIES ${SUNDIALS_NVECSERIAL_LIB} ${SUNDIALS_KINSOL_LIB}
        ${SUNDIALS_CVODE_LIB} ${SUNDIALS_SUNLINSOLSPGMR_LIB}
        ${SUNDIALS_IDA_LIB})
include_directories(${SUNDIALS_INCLUDE_DIRS})

# find OpenMPI
#find_package(MPI REQUIRED)
#if(MPI_FOUND)
#  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
#else(MPI_FOUND)
#  message(SEND_ERROR "This application cannot compile without MPI")
#endif(MPI_FOUND)

# Link to source files
set(SOURCE_FILES main.cpp)
add_executable(
        LVMCM
        main.cpp
        ODE.cc
        Topography.cpp
        Species.cpp
        Metacommunity.cpp
        CommunityDynamics.cpp
        CommunityDynamics.h
        LVMCM_rng.cpp
        LVMCM_rng.h)
target_link_libraries(LVMCM ${ARMADILLO_LIBRARIES})
target_link_libraries(LVMCM ${SUNDIALS_LIBRARIES})
target_link_libraries(LVMCM ${Boost_LIBRARIES})
if (NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  target_link_libraries(LVMCM ${OPENBLAS_LIB})
endif()
#target_link_libraries(LVMCM ${MPI_LIBRARIES})

# For testing only:
add_executable(testODE testODE.cpp ODE.cc)
target_link_libraries(testODE ${SUNDIALS_LIBRARIES})
target_link_libraries(testODE ${ARMADILLO_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(LVMCM PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(LVMCM PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# Uncomment if running HPCtoolkit
#SET(CMAKE_CXX_FLAGS "-g")

# For native compiling opimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native")

# For additional debugging
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -check=conversions,stack,uninit")
