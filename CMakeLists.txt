#
# CMakeLists.txt
#
# ‘created’ 3 Jan 2018 GCRich, though i take no credit for its authorship
# based on a CMakeList made (at least in part) by J. Newby of ORNL
#  
#
# no effort put into ensuring this is robust and widely usable
#

cmake_minimum_required(VERSION 2.6)
project(waveformTools)

if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" "$ENV{ROOTSYS}/lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" "$ENV{ROOTSYS}/lib")
    endif()
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/cmake/modules)
# include(FindROOT)
# include(ROOT_GENERATE_DICTIONARY)

include(CheckIncludeFile)

find_package(ROOT REQUIRED)
message(${ROOT_INCLUDE_DIRS})
message(${ROOT_USE_FILE})
include(${ROOT_USE_FILE})

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


set(WAVEBASE $ENV{HOME}/Dev/COHERENT/CASA/TWaveform-casa)

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR}
${CMAKE_SOURCE_DIR}
${WAVEBASE}/WaveBase
$ENV{PWD}
)

include_directories( ${INCLUDE_DIRECTORIES})
set( SOURCES
waveformTools.cc
)

set( HEADERS
waveformTools.hh
)

set(TWAVEFORMLIB ${WAVEBASE}/lib/libWaveWaveBase.dylib)

set(LINKDEF LinkDef.h)
set(DICTIONARY ${PROJECT_NAME}_dict)
set(PCMFILE ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm)
set(RMAPFILE ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}.rootmap)
ROOT_GENERATE_DICTIONARY("${DICTIONARY}" "${HEADERS}" LINKDEF "${LINKDEF}")
SET(SOURCES ${HEADERS} ${SOURCES} ${DICTIONARY})

add_library(${PROJECT_NAME}Lib SHARED ${SOURCES})
target_link_libraries(${PROJECT_NAME}Lib ${ROOT_LIBRARIES} ${TWAVEFORMLIB})
set_target_properties(${PROJECT_NAME}Lib PROPERTIES SUFFIX ".so")