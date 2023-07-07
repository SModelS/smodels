# CMake file for Resummino.
#
# Copyright 2008-2010 Jonathan Debove.
# Copyright 2011-2014 David R. Lamprea.
# Copyright 2011-2014 Marcel Rothering.
#
# Licensed under the terms of the EUPL version 1.1 or later.
# See the LICENCE file for more information.

cmake_minimum_required(VERSION 2.9)
project(resummino C CXX Fortran)


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()





set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/external/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
include_directories(${CMAKE_SOURCE_DIR}/external/include)
include_directories(${CMAKE_SOURCE_DIR}/external/include/Fastor)

include(ExternalProject)

if(CI)
	set(BUILD_LOOPTOOLS TRUE)
	set(BUILD_LHAPDF TRUE)
else()
	set(BUILD_LOOPTOOLS TRUE)
	#set(BUILD_LHAPDF FALSE)
endif()

## compiles looptools 213 (configure script changed so that the folder is always calles lib)
if(BUILD_LOOPTOOLS)
	ExternalProject_Add(lt
		SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/LoopTools-2.13
		CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FFLAGS=-std=legacy CCFLAGS=-fPIC CPPFLAGS=-fPIC ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/LoopTools-2.13/configure --prefix=<INSTALL_DIR>
		BUILD_COMMAND ${MAKE})
	ExternalProject_Get_Property(lt install_dir)
	set( LT_DIR ${install_dir})
	include_directories(${LT_DIR}/include)
	set(LOOPTOOLS_LIBRARIES ${LT_DIR}/lib/libooptools.a)
else()
	find_package(LoopTools REQUIRED)
	include_directories(${LOOPTOOLS_INCLUDE_DIRS})
endif()

# LHAPDF
if(BUILD_LHAPDF)
	# This is needed since otherwise LHAPDF autoconf files get new timestamp and
	# this screws the compilation up.
	if(NOT EXISTS ${CMAKE_SOURCE_DIR}/external/lib/LHAPDF-6.2.3)
		execute_process(COMMAND tar xvfz
		${CMAKE_SOURCE_DIR}/external/lib/LHAPDF-6.2.3.tar.gz -C
		${CMAKE_SOURCE_DIR}/external/lib/)
	endif()
	ExternalProject_Add(lhapdf
		SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/LHAPDF-6.2.3
		CONFIGURE_COMMAND
		${CMAKE_CURRENT_SOURCE_DIR}/external/lib/LHAPDF-6.2.3/configure
		--prefix=<INSTALL_DIR> --disable-python
		BUILD_COMMAND ${MAKE}) 
	ExternalProject_Get_Property(lhapdf install_dir)
	set( LHAPDF_DIR ${install_dir})
	include_directories(${LHAPDF_DIR}/include)
	set(LHAPDF_LIBRARIES ${LHAPDF_DIR}/lib/libLHAPDF.a)
else()
	find_package(LHAPDF REQUIRED)
	include_directories(${LHAPDF_INCLUDE_DIRS})
endif()



# Boost headers
find_package(Boost 1.67.0 REQUIRED)


# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


add_subdirectory(src)
