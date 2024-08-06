# Common boilerplate
cmake_minimum_required(VERSION 3.5)

message(${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prevent in-source builds
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory and run CMake from there")
endif()

# Determine if this is built as a subproject (using add_subdirectory)
# or if it is the master project.
set(MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MASTER_PROJECT ON)
  message(STATUS "Feaspump is MASTER PROJECT")
  message(STATUS "CMake version: ${CMAKE_VERSION}")

  project(dominiqs-feaspump VERSION 1.0)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  # # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE KP_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE KP_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

# Find CPLEX library
find_package(CPLEX)

# Find XPRESS library
find_package(XPRESS)

# Find SCIP library
find_package(SCIP)

# Find ORTOOLS library
find_package(ORTOOLS)

# Find Threads
find_package(Threads)

# Define libkp
add_library(libkp STATIC ../src/feaspump.cpp ../src/transformers.cpp ../src/ranking.cpp ../src/solution.cpp ../src/kernelpump.cpp)
target_link_libraries(libkp PUBLIC Utils::Lib fmt::fmt Prop::Lib)
add_library(Kp::Lib ALIAS libkp)

target_include_directories(libkp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
  $<INSTALL_INTERFACE:include>
)

# Define results executable
add_executable(results src/csv_reader.hpp src/main.cpp)


list(APPEND includePath "${CMAKE_CURRENT_SOURCE_DIR}/../include" "${CMAKE_CURRENT_SOURCE_DIR}/..")

target_include_directories(results PUBLIC
  "$<BUILD_INTERFACE:${includePath}>"
  $<INSTALL_INTERFACE:include>
)

if (APPLE)
	target_link_libraries(results -Wl,-force_load Prop::Lib -Wl,-force_load Kp::Lib Utils::Lib fmt::fmt)
else()
	target_link_libraries(results  -Wl,--whole-archive Prop::Lib Kp::Lib -Wl,--no-whole-archive Utils::Lib fmt::fmt)
endif()


# Deal with optional dependencies
if (CPLEX_FOUND)
  target_compile_definitions(libkp PUBLIC HAS_CPLEX=1)
  target_compile_definitions(results PUBLIC HAS_CPLEX=1)
  target_sources(libkp PRIVATE ../src/cpxmodel.cpp)
  target_link_libraries(libkp PUBLIC Cplex::Cplex)
endif()

if (XPRESS_FOUND)
  target_compile_definitions(libkp PUBLIC HAS_XPRESS=1)
  target_compile_definitions(results PUBLIC HAS_XPRESS=1)
  target_sources(libkp PRIVATE ../src/xprsmodel.cpp)
  target_link_libraries(libkp PUBLIC Xpress::Xpress)
endif()

if (SCIP_FOUND)
  target_compile_definitions(libkp PUBLIC HAS_SCIP=1)
  target_compile_definitions(results PUBLIC HAS_SCIP=1)
  target_sources(libkp PRIVATE ../src/scipmodel.cpp)
  target_link_libraries(libkp PUBLIC Scip::Scip)
endif()

if (ORTOOLS_FOUND)
  target_compile_definitions(libkp PUBLIC HAS_ORTOOLS=1)
  target_compile_definitions(results PUBLIC HAS_ORTOOLS=1)
  target_sources(libkp PRIVATE ../src/pdlpmodel.cpp)
  target_link_libraries(libkp PUBLIC Pdlp::Pdlp)
endif()

# Generate version.h
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/kernelpump/version.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/kernelpump/version.h
)

include_directories ("/usr/include/eigen3")

# Add subprojects
if (GIT_FOUND AND MASTER_PROJECT)
  add_subdirectory(../extern/fmt extern/fmt)
  add_subdirectory(../extern/utils extern/utils)
  add_subdirectory(../extern/propagator extern/propagator)
endif()
