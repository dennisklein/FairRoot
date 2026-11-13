cmake_minimum_required(VERSION 3.5...3.28 FATAL_ERROR)
project(TGeo2VecGeom VERSION 0.1)

if(NOT CMAKE_BUILD_TYPE)
    message(WARNING "CMAKE_BUILD_TYPE not set : will use RelWithDebInfo")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

find_package(ROOT REQUIRED)
find_package(VecGeom REQUIRED)

#-------------------------------------------------q
# Locate sources for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

# define the main target : The TGeo2VecGeom library
add_library(tgeo2vecgeom SHARED ${sources})
add_library(TGeo2VecGeom::tgeo2vecgeom ALIAS tgeo2vecgeom)

# The dependency include directories are deliberately confined to BUILD_INTERFACE.
# Anything listed unguarded here is written verbatim, as an absolute path, into the
# installed TGeo2VecGeomTargets.cmake -- which makes the install non-relocatable: a
# binary package moved to another prefix (or unpacked by a build-cache consumer on a
# different machine) then points at a ROOT/VecGeom include directory that does not
# exist, and every consumer fails at configure time with
#   Imported target "TGeo2VecGeom::TGeo2VecGeom" includes non-existent path ...
# Consumers get the ROOT and VecGeom headers through the imported targets below,
# which each resolve their own prefix at find_package() time.
target_include_directories(
        tgeo2vecgeom
        PUBLIC
        "$<BUILD_INTERFACE:${ROOT_INCLUDE_DIRS}>"
        "$<BUILD_INTERFACE:${VECGEOM_INCLUDE_DIR}>"
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# we depend on VecGeom / ROOT. These are imported target names (VECGEOM_LIBRARIES expands
# to VecGeom::vecgeom, plus VecGeom::vecgeomcuda for a CUDA-enabled VecGeom), so the export
# records the names and the consumer's own find_dependency resolves their usage requirements
# -- include directories included -- against the ROOT and VecGeom it actually has.
target_link_libraries(tgeo2vecgeom
        PUBLIC
        ROOT::Geom ${VECGEOM_LIBRARIES}
        )

# installation specific
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/TGeo2VecGeom)

install(TARGETS tgeo2vecgeom
        EXPORT tgeo2vecgeom-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

# This is required so that the exported target
# has the name TGeo2VecGeom and not tgeo2vecgeom
set_target_properties(tgeo2vecgeom PROPERTIES EXPORT_NAME TGeo2VecGeom)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/TGeo2VecGeom)

#Export the targets to a script
install(EXPORT tgeo2vecgeom-targets
        FILE
        TGeo2VecGeomTargets.cmake
        NAMESPACE
        TGeo2VecGeom::
        DESTINATION
        ${INSTALL_CONFIGDIR}
        )

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/TGeo2VecGeomConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/TGeo2VecGeomConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
        )

# Install the config, configversion and custom find modules
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/TGeo2VecGeomConfig.cmake
        # ${CMAKE_CURRENT_BINARY_DIR}/TGeo2VecGeomConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

##############################################
## Exporting from the build tree

export(EXPORT tgeo2vecgeom-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/TGeo2VecGeomTargets.cmake NAMESPACE TGeo2VecGeom::)

# NOT registered in the user package registry (~/.cmake/packages): an entry there makes a
# stale build tree win over the installed package in a later, unrelated find_package(), which
# is exactly the kind of surprise a build-cache/CI setup cannot afford. Consumers point at the
# install with CMAKE_PREFIX_PATH or TGeo2VecGeom_DIR.

##############################################
## Tests

option(TGEO2VECGEOM_BUILD_TESTS "Build tgeo2vecgeom test/smoke executables" ON)
if(TGEO2VECGEOM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
