cmake_minimum_required(VERSION 3.24)
project(flatsim LANGUAGES C CXX)

# Baseline superbuild: each upstream package keeps its own build system and
# its own shared libraries, as the organizations ship them. The packages are
# built in dependency order into one staging prefix.

include(ExternalProject)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

set(FLATSIM_STAGE "${CMAKE_BINARY_DIR}/stage" CACHE PATH
  "Install prefix of all packages")
set(FLATSIM_G4DATA_DIR "$ENV{HOME}/.cache/flatsim/geant4-data" CACHE PATH
  "Geant4 data sets. Outside of the build tree, so that a clean build reuses them")
set(FLATSIM_JOBS 4 CACHE STRING "Parallel jobs for each package build")

set(ext "${CMAKE_SOURCE_DIR}/extern")
set(common_args
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_INSTALL_PREFIX=${FLATSIM_STAGE}
  -DCMAKE_PREFIX_PATH=${FLATSIM_STAGE}
  -DCMAKE_INSTALL_LIBDIR=lib
  -DCMAKE_CXX_STANDARD=20
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
)
set(build_cmd ${CMAKE_COMMAND} --build <BINARY_DIR> --parallel ${FLATSIM_JOBS})

function(flatsim_package name)
  cmake_parse_arguments(P "" "SOURCE_DIR" "DEPENDS;ARGS" ${ARGN})
  if(NOT P_SOURCE_DIR)
    set(P_SOURCE_DIR "${ext}/${name}")
  endif()
  ExternalProject_Add(${name}
    SOURCE_DIR "${P_SOURCE_DIR}"
    BINARY_DIR "${CMAKE_BINARY_DIR}/pkg/${name}"
    PREFIX "${CMAKE_BINARY_DIR}/ep/${name}"
    INSTALL_DIR "${FLATSIM_STAGE}"
    CMAKE_ARGS ${common_args} ${P_ARGS}
    BUILD_COMMAND ${build_cmd}
    DEPENDS ${P_DEPENDS}
    LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON LOG_OUTPUT_ON_FAILURE ON
  )
endfunction()

flatsim_package(fmt ARGS
  -DBUILD_SHARED_LIBS=ON -DFMT_TEST=OFF -DFMT_DOC=OFF)

flatsim_package(fairlogger DEPENDS fmt ARGS
  -DUSE_EXTERNAL_FMT=ON -DBUILD_TESTING=OFF)

flatsim_package(ms_gsl ARGS -DGSL_TEST=OFF)

# ROOT: only the parts a simulation needs. All third-party code comes from
# ROOT's own builtins, so that the only host tools are CMake and GCC.
flatsim_package(root ARGS
  -Dminimal=ON
  -Dgeom=ON
  -Dbuiltin_zlib=ON -Dbuiltin_lzma=ON -Dbuiltin_lz4=ON -Dbuiltin_zstd=ON
  -Dbuiltin_xxhash=ON -Dbuiltin_pcre=ON -Dbuiltin_freetype=ON
  -Dbuiltin_nlohmannjson=ON -Dbuiltin_openui5=OFF
  -Dbuiltin_llvm=ON -Dbuiltin_clang=ON -Dbuiltin_cling=ON
  -Dfail-on-missing=ON
  -Dsoversion=ON
  -Drpath=ON
  -Dtesting=OFF)

flatsim_package(vmc DEPENDS root)

flatsim_package(geant4 ARGS
  -DBUILD_SHARED_LIBS=ON
  -DGEANT4_BUILD_MULTITHREADED=ON
  -DGEANT4_BUILD_TLS_MODEL=initial-exec
  -DGEANT4_ENABLE_TESTING=OFF
  -DGEANT4_INSTALL_EXAMPLES=OFF
  -DGEANT4_USE_SYSTEM_EXPAT=OFF
  -DGEANT4_USE_SYSTEM_ZLIB=OFF
  -DGEANT4_USE_SYSTEM_CLHEP=OFF
  -DGEANT4_USE_G3TOG4=ON
  -DGEANT4_USE_GDML=OFF
  -DGEANT4_INSTALL_DATA=ON
  -DGEANT4_INSTALL_DATADIR=${FLATSIM_G4DATA_DIR}
  -DGEANT4_INSTALL_DATA_TIMEOUT=2000)

flatsim_package(geant4_vmc DEPENDS geant4 vmc root ARGS
  -DGeant4VMC_USE_VGM=OFF
  -DGeant4VMC_BUILD_EXAMPLES=OFF
  -DGeant4VMC_INSTALL_EXAMPLES=OFF)

flatsim_package(pythia8 SOURCE_DIR "${CMAKE_SOURCE_DIR}/cmake/pythia8" ARGS
  -DPYTHIA8_SOURCE_DIR=${ext}/pythia8)

flatsim_package(boost SOURCE_DIR "${CMAKE_SOURCE_DIR}/cmake/boost" ARGS
  -DBOOST_SOURCE_DIR=${ext}/boost)

flatsim_package(tbb ARGS
  -DTBB_TEST=OFF -DTBB_EXAMPLES=OFF -DTBB_STRICT=OFF)

flatsim_package(o2 SOURCE_DIR "${CMAKE_SOURCE_DIR}/cmake/o2"
  DEPENDS root vmc geant4 geant4_vmc pythia8 fairlogger fmt tbb boost ms_gsl
  ARGS
  -DO2_SOURCE_DIR=${ext}/o2
  -DFAIRMQ_SOURCE_DIR=${ext}/fairmq
  -DRAPIDJSON_SOURCE_DIR=${ext}/rapidjson
  -DMSGSL_SOURCE_DIR=${ext}/ms_gsl
  -DFLATSIM_APP_DIR=${CMAKE_SOURCE_DIR}/app
  -DFLATSIM_HARNESS_DIR=${CMAKE_SOURCE_DIR}/harness)
