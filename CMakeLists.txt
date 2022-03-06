#
# Main CMakeLists.txt file to build library "linear_algebra"
#

cmake_minimum_required( VERSION 3.10 )

#
# compiler settings
#

# set compiler type
if ( NOT DEFINED CMAKE_Fortran_COMPILER)
  if ( DEFINED My_Fortran_Compiler)
    SET( CMAKE_Fortran_COMPILER ${My_Fortran_Compiler})
  else ()
    SET( CMAKE_Fortran_COMPILER /usr/bin/gfortran )
  endif()
  MESSAGE( STATUS "Compiling with ${CMAKE_Fortran_COMPILER}")
endif()

message("**")
message("** PROJECT: linear_algebra")

# set default 
if ( NOT DEFINED BUILD_TYPE)
  set (BUILD_TYPE "RELEASE")
endif()
set (CMAKE_Fortran_FLAGS_BASE "-fPIC -ffree-line-length-none")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_BASE} -O5")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_BASE} -g -C -Wall -fcheck=all -O")
set (CMAKE_Fortran_FLAGS_RelWithIEEE "${CMAKE_Fortran_FLAGS_REALESE} -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")

# set build type (DEBUG/RELEASE)
if ( ${BUILD_TYPE} MATCHES "DEBUG" )
   message("** Compiling in debug mode")
   set ( CMAKE_BUILD_TYPE Debug )
elseif ( ${BUILD_TYPE} MATCHES "RELEASE"   )
   message("** Compiling in release mode")
   set ( CMAKE_BUILD_TYPE RELEASE )
elseif ( ${BUILD_TYPE} MATCHES "IEEE"   )
   message("** Compiling in release with IEEE mode")
   set ( CMAKE_BUILD_TYPE RelWithIEEE )
 else ()
   set ( CMAKE_BUILD_TYPE user )
   set ( CMAKE_Fortran_FLAGS_user  "${CMAKE_Fortran_FLAGS_BASE} ${BUILD_TYPE}" )
   message("** Compiling in user mode:  ${CMAKE_Fortran_FLAGS_user}")
 endif()

# set the project name
project( linear_algebra LANGUAGES Fortran VERSION 1.0 )


#
# set default type of library to look for
#
# set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)
# LOOK FOR BLAS LIBRARY
# SETTING VARIABLES:
# BLAS_LIB, BLAS_DIR
if ( NOT DEFINED BLAS_DIR)
  message(STATUS "LOOKING FOR BLAS LIBRARIES")
  find_package(BLAS REQUIRED)
  if (BLAS_FOUND)
    message(STATUS "BLAS LIBRARIES: ${BLAS_LIBRARIES}")
    list(GET BLAS_LIBRARIES 0 BLAS_LIB)
    message(STATUS "BLAS LIB SELECTED: ${BLAS_LIB}")
    get_filename_component(BLAS_DIR "${BLAS_LIB}" DIRECTORY)
  else()
    # print error message
    MESSAGE( SEND_ERROR  "+ BLAS NOT FOUND" )
    MESSAGE( SEND_ERROR  "+ Pass the full path of the directory containg it with" )
    MESSAGE( SEND_ERROR  "+ -DBLAS_DIR=absolute path")
    MESSAGE( SEND_ERROR  "+ Check README for instructions" )
    MESSAGE( FATAL_ERROR "++++++++++++++++++++++++++++++++++++++++++++++" )
  endif()
else()
  message(STATUS "DIRECTORY PATH  FOR BLAS ${BLAS_DIR}")
endif()

# # LOOK FOR LAPACK LIBRARY
# # SETTING VARIABLES:
# # LAPACK_LIB, LAPACK_DIR
# if ( NOT DEFINED LAPACK_DIR)
#   message(STATUS "LOOKING FOR LAPACK LIBRARIES")
#   find_package(LAPACK REQUIRED)
#   if (LAPACK_FOUND)
#     message(VERBOSE "LAPACK LIBRARIES: ${LAPACK_LIBRARIES}")
#     list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
#     message(VERBOSE "LAPACK LIB SELECTED: ${LAPACK_LIB}")
#     get_filename_component(LAPACK_DIR "${LAPACK_LIB}" DIRECTORY)
#   else()
#     # print error message
#     MESSAGE( SEND_ERROR  "+ LAPACK LIBRARY NOT FOUND: Pass the full path of the directory containg it with" )
#     MESSAGE( SEND_ERROR  "+ -DLAPACK_DIR= absolute path for lapack ")
#     MESSAGE( FATAL_ERROR "++++++++++++++++++++++++++++++++++++++++++++++" )
#   endif()
# else()
#   message(STATUS "DIRECTORY PATH FOR LAPACK: ${LAPACK_DIR}")
# endif()


###########
SET( LIBS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
SET( LIBS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} )
SET( MODS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/mods )
SET( CMAKE_Fortran_MODULE_DIRECTORY ${MODS_BINARY_DIR} )



# get all sources
FILE( GLOB SOURCES src_abs_vect/*.f90 )
SET( linalg_SOURCES ${SOURCES} )
# create objects
MESSAGE("srs files ${linalg_SOURCES} ")
ADD_LIBRARY( linalg STATIC ${linalg_SOURCES} )


add_executable ( test test/main_dvector.f90)
target_link_libraries( test  linalg)

message("**")



