## Set cmake and compiler options
cmake_minimum_required(VERSION 3.5.1)
set(CMAKE_CXX_STANDARD 11)
project(LVMCM)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

## Find libraries
#Find and include Armadillo
find_package(Armadillo REQUIRED)
IF (Armadillo_FOUND)
    INCLUDE_DIRECTORIES(${Armadillo_INCLUDE_DIR})
ELSE ()
    MESSAGE(FATAL_ERROR "Could not find Armadillo!")
ENDIF()

#Find and include Boost
set(BOOST_ROOT /usr/lib/x86_64-linux-gnu)
find_package(Boost COMPONENTS system filesystem REQUIRED)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ELSE ()
    MESSAGE(FATAL_ERROR "Could not find Boost!")
ENDIF()

#Find and include Sundials nvecserial, kinsol, cvode, cvodes,
find_library(SUNDIALS_NVECSERIAL_LIB sundials_nvecserial)
IF(NOT SUNDIALS_NVECSERIAL_LIB)
    MESSAGE(FATAL_ERROR "Could not find sundials_nvecserial!")
ENDIF()

find_library(SUNDIALS_KINSOL_LIB sundials_kinsol)
IF(NOT SUNDIALS_KINSOL_LIB)
    MESSAGE(FATAL_ERROR "Could not find sundials_kinsol!")
ENDIF()

find_library(SUNDIALS_CVODE_LIB sundials_cvode)
IF(NOT SUNDIALS_CVODE_LIB)
    MESSAGE(FATAL_ERROR "Could not find sundials_cvode!")
ENDIF()

find_library(SUNDIALS_CVODES_LIB sundials_cvodes)
IF(NOT SUNDIALS_KINSOL_LIB)
    MESSAGE(FATAL_ERROR "Could not find sundials_cvodes!")
ENDIF()

SET(SUNDIALS_LIBRARIES ${SUNDIALS_NVECSERIAL_LIB} ${SUNDIALS_KINSOL_LIB} ${SUNDIALS_CVODE_LIB})
include_directories(${SUNDIALS_INCLUDE_DIRS})

# find_library(OPENBLAS_LIB openblasp) # might be better sometimes
find_library(OPENBLAS_LIB openblas)
IF(NOT OPENBLAS_LIB)
    MESSAGE(FATAL_ERROR "Could not find openblas!")
ENDIF()

# find OpenMPI
find_package(MPI REQUIRED)
IF(MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
ELSE(MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
ENDIF(MPI_FOUND)

## Link to source files
set(SOURCE_FILES main.cpp)
add_executable(LVMCM main.cpp ODE.cc Topography.cpp Species.cpp Metacommunity.cpp CommunityDynamics.cpp CommunityDynamics.h LVMCM_rng.cpp LVMCM_rng.h)
target_link_libraries(LVMCM ${ARMADILLO_LIBRARIES})
target_link_libraries(LVMCM ${SUNDIALS_LIBRARIES})
target_link_libraries(LVMCM ${Boost_LIBRARIES})
target_link_libraries(LVMCM ${OPENBLAS_LIB})
target_link_libraries(LVMCM ${MPI_LIBRARIES})

IF(MPI_COMPILE_FLAGS)
    set_target_properties(LVMCM PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
ENDIF()

IF(MPI_LINK_FLAGS)
    set_target_properties(LVMCM PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
ENDIF()

