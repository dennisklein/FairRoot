# CMake build for the imported Boost subset. Upstream builds Boost with b2.
# This file installs the headers and builds Boost.ProgramOptions as a
# shared library, so that the only host tools are CMake and the compiler.
cmake_minimum_required(VERSION 3.24)
project(Boost VERSION 1.90.0 LANGUAGES CXX)

set(BOOST_SOURCE_DIR "" CACHE PATH "Boost source tree")
include(GNUInstallDirs)

add_library(headers INTERFACE)
target_include_directories(headers INTERFACE
  $<BUILD_INTERFACE:${BOOST_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

file(GLOB po_sources "${BOOST_SOURCE_DIR}/libs/program_options/src/*.cpp")
list(FILTER po_sources EXCLUDE REGEX "winmain\\.cpp$")
add_library(program_options SHARED ${po_sources})
set_target_properties(program_options PROPERTIES OUTPUT_NAME boost_program_options)
target_link_libraries(program_options PUBLIC headers)
target_compile_definitions(program_options PUBLIC BOOST_PROGRAM_OPTIONS_DYN_LINK PRIVATE BOOST_PROGRAM_OPTIONS_SOURCE)

install(DIRECTORY "${BOOST_SOURCE_DIR}/boost" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS headers program_options EXPORT BoostTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT BoostTargets NAMESPACE Boost::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FlatsimBoost)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/FlatsimBoostConfig.cmake"
  "include(\"\${CMAKE_CURRENT_LIST_DIR}/BoostTargets.cmake\")\n")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/FlatsimBoostConfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FlatsimBoost)
