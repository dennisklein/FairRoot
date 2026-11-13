# Determine the alignment required by the configured host vector backend.
#
# AlignedBase is used by raw-copied helper types such as Vector3D.  Host and
# CUDA declarations must therefore use the same stride for these types.

function(vecgeom_determine_alignedbase_alignment out_var)
  set(_override_var "${out_var}_OVERRIDE")
  set(${_override_var} "" CACHE STRING
    "Override AlignedBase alignment for cross builds where the host vector-backend probe cannot run")
  mark_as_advanced(${_override_var})

  set(_alignment 1)
  string(TOLOWER "${VECGEOM_BACKEND}" _backend_lower)

  if(NOT "${${_override_var}}" STREQUAL "")
    set(_alignment "${${_override_var}}")
  elseif(NOT _backend_lower STREQUAL "scalar")
    if(CMAKE_CROSSCOMPILING)
      message(FATAL_ERROR
        "Cannot run the host vector-backend alignment probe while cross-compiling. "
        "Set -D${_override_var}=<alignment> explicitly.")
    endif()

    # Generate a temporary Config.h so Global.h can select the configured host
    # backend before the final configured header is generated.
    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/VecGeom/base")
    set(${out_var} 1)
    configure_file("${PROJECT_SOURCE_DIR}/VecGeom/base/Config.h.in"
      "${PROJECT_BINARY_DIR}/VecGeom/base/Config.h" @ONLY)

    set(_probe_dir "${PROJECT_BINARY_DIR}/CMakeFiles/VecGeomAlignmentProbe")
    file(MAKE_DIRECTORY "${_probe_dir}")
    set(_probe_src "${_probe_dir}/VecGeomAlignmentProbe.cpp")
    file(WRITE "${_probe_src}"
"#include \"VecGeom/base/Global.h\"
#include <iostream>

int main()
{
  std::cout << alignof(vecgeom::VectorBackend::Real_v) << '\\n';
  return 0;
}
")
    separate_arguments(_probe_compile_flags NATIVE_COMMAND
      "${CMAKE_CXX_FLAGS} ${VECGEOM_CXX_FLAGS} ${VECGEOM_COMPILATION_FLAGS}")

    try_run(_probe_run_result _probe_compile_result
      "${_probe_dir}"
      "${_probe_src}"
      CMAKE_FLAGS
        "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}"
        "-DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}"
        "-DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}"
        "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} ${VECGEOM_CXX_FLAGS} ${VECGEOM_COMPILATION_FLAGS}"
      COMPILE_DEFINITIONS
        "-I${PROJECT_SOURCE_DIR}"
        ${_probe_compile_flags}
        "-I${PROJECT_BINARY_DIR}"
      LINK_LIBRARIES VecCore::VecCore
      RUN_OUTPUT_VARIABLE _probe_output
      COMPILE_OUTPUT_VARIABLE _probe_compile_output)

    if(NOT _probe_compile_result OR NOT _probe_run_result EQUAL 0)
      message(FATAL_ERROR "Failed to determine host vector-backend alignment.\n"
        "Compile output:\n${_probe_compile_output}\n"
        "Run output:\n${_probe_output}")
    endif()
    string(STRIP "${_probe_output}" _alignment)
  endif()

  if(NOT "${_alignment}" MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "Invalid AlignedBase alignment value: '${_alignment}'")
  endif()
  set(${out_var} "${_alignment}" PARENT_SCOPE)
  message(STATUS "AlignedBase host vector alignment: ${_alignment}")
endfunction()
