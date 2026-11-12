# CMake build for Pythia 8. Upstream builds with its own configure script and
# GNU make. This file builds the same library from the same sources with the
# flags of the upstream configure script, so that the only host tools are
# CMake and the compiler.
cmake_minimum_required(VERSION 3.24)
project(Pythia8 LANGUAGES CXX)

set(PYTHIA8_SOURCE_DIR "" CACHE PATH "Pythia 8 source tree")
file(GLOB sources CONFIGURE_DEPENDS "${PYTHIA8_SOURCE_DIR}/src/*.cc")

add_library(pythia8 SHARED ${sources})
target_include_directories(pythia8 PUBLIC
  $<BUILD_INTERFACE:${PYTHIA8_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_options(pythia8 PRIVATE -pedantic -W -Wall -Wshadow)
target_link_libraries(pythia8 PUBLIC ${CMAKE_DL_LIBS})
# Per-file definitions of the upstream Makefile.
include(GNUInstallDirs)
set_source_files_properties("${PYTHIA8_SOURCE_DIR}/src/Pythia.cc" PROPERTIES
  COMPILE_DEFINITIONS "XMLDIR=\"${CMAKE_INSTALL_FULL_DATADIR}/Pythia8/xmldoc\"")
set_source_files_properties("${PYTHIA8_SOURCE_DIR}/src/FJcore.cc" PROPERTIES
  COMPILE_DEFINITIONS "FJCORE_HAVE_LIMITED_THREAD_SAFETY")

include(GNUInstallDirs)
install(TARGETS pythia8 EXPORT Pythia8Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY "${PYTHIA8_SOURCE_DIR}/include/"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY "${PYTHIA8_SOURCE_DIR}/share/Pythia8"
  DESTINATION ${CMAKE_INSTALL_DATADIR})
install(EXPORT Pythia8Targets NAMESPACE Pythia8::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Pythia8)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/Pythia8Config.cmake"
  "include(\"\${CMAKE_CURRENT_LIST_DIR}/Pythia8Targets.cmake\")\n"
  "set(Pythia8_FOUND TRUE)\n")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Pythia8Config.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Pythia8)
