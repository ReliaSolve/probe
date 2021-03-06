cmake_minimum_required(VERSION 3.10.0)
project(probe)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

# Find the packages we require and can use
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CCTBX REQUIRED)
find_package(Boost REQUIRED)

set(probe_SOURCES
	DotSpheres.cpp
)

set(probe_HEADERS
  DotSpheres.h
)

add_library(probeLib ${probe_SOURCES} ${probe_HEADERS})
target_include_directories(probeLib PUBLIC SYSTEM 
  ${CCTBX_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR})
if (NOT WIN32)
  target_link_libraries(probeLib PUBLIC m)
endif()

add_executable(tst_probe tst_probe.cpp)
target_link_libraries(tst_probe probeLib)
if (NOT WIN32)
  target_link_libraries(tst_probe m)
endif()
