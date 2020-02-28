################################################################################
# Copyright (C) 2020 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH       #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

################################################################################
# Exchange file extention of file in list from ext1 to ext2,
# only if the resulting file exists in CMAKE_CURRENT_SOURCE_DIR,
# and assign the newly created list to 'output'.
# The input list is not changed at all
# Ex: fair_change_file_extension_ifexists(*.cxx *.h "${TRD_SRCS}" TRD_HEADERS)
################################################################################
function(fairroot_change_extensions_if_exists ext1 ext2)
  cmake_parse_arguments(ARGS "" "OUTVAR" "FILES" ${ARGN})

  set(required_args "FILES;OUTVAR")
  foreach(required_arg IN LISTS required_args)
    if(NOT ARGS_${required_arg})
      message(FATAL_ERROR "fair_change_extensions_if_exists is missing a required argument: ${required_arg}")
    endif()
  endforeach()

  if(${ext1} MATCHES "^[*][.]+.*$")
    string(REGEX REPLACE "^[*]+([.].*)$" "\\1" ext1new ${ext1})
  else()
    set(ext1new ${ext1})
  endif()
  if(${ext2} MATCHES "^[*][.]+.*$")
    string(REGEX REPLACE "^[*]+([.].*)$" "\\1" ext2new ${ext2})
  else()
    set(ext2new ${ext2})
  endif()
  foreach(file ${ARGS_FILES})
    set(newFile "")
    string(REGEX REPLACE "^(.*)${ext1new}$" "\\1${ext2new}" newFile ${file})
    if(NOT ${file} STREQUAL ${newFile} AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${newFile})
      list(APPEND result ${newFile})
    endif()
  endforeach()
  set(${ARGS_OUTVAR} ${result} PARENT_SCOPE)
endfunction()

macro(fairroot_remove_system_paths path_var)
  set(__excluded_system_prefixes)
  list(APPEND __excluded_system_prefixes
    ${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}
    /bin
    /usr/bin)
  list(REMOVE_DUPLICATES __excluded_system_prefixes)

  list(JOIN __excluded_system_prefixes "|" __excluded_system_prefixes)
  set(__excluded_system_prefixes "^(${__excluded_system_prefixes})")

  set(__tmp ${${path_var}})
  set(${path_var})
  foreach(dir IN LISTS __tmp)
    if(NOT dir MATCHES ${__excluded_system_prefixes})
      list(APPEND ${path_var} ${dir})
    endif()
  endforeach()
  unset(__tmp)
  unset(__excluded_system_prefixes)
endmacro()

# extendable map as global properties ROOT_CONFIG_LIBRARY_TARGETS_<pkg> = list of targets
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_ROOT ROOT::Core)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_Geant3 geant321)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_Geant4 G4run)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_Geant4VMC geant4vmc)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_VMC VMCLibrary)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_Pythia6 Pythia6)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_Pythia8 Pythia8)
set_property(GLOBAL APPEND PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_VGM VGM::BaseVGM)

function(fairroot_get_root_config_library_targets pkg)
  cmake_parse_arguments(ARGS "" "OUTVAR" "" ${ARGN})
  if(NOT ARGS_OUTVAR)
    set(ARGS_OUTVAR ROOT_CONFIG_LIBRARY_TARGETS)
  endif()

  set(${ARGS_OUTVAR})
  get_property(tgts GLOBAL PROPERTY ROOT_CONFIG_LIBRARY_TARGETS_${pkg})
  foreach(tgt IN LISTS tgts)
    if(TARGET ${tgt})
      list(APPEND ${ARGS_OUTVAR} ${tgt})
    endif()
  endforeach()

  set(${ARGS_OUTVAR} ${${ARGS_OUTVAR}} PARENT_SCOPE)
endfunction()

function(fairroot_generate_config_sh file)
  cmake_parse_arguments(ARGS "BUILD;INSTALL" "VMCWORKDIR" "DEPENDENCIES;PATH" ${ARGN})

  file(WRITE ${file} "#!/bin/bash\n# Generated by ${PROJECT_NAME} ${PROJECT_GIT_VERSION}\n")

  if(ARGS_BUILD AND ARGS_INSTALL)
    message(FATAL_ERROR "BUILD and INSTALL are mutually exclusive")
  endif()

  list(APPEND ARGS_DEPENDENCIES ${PROJECT_NAME})
  foreach(dep IN LISTS ARGS_DEPENDENCIES)
    fairroot_get_root_config_library_targets(${dep} OUTVAR targets)
    if(targets)
      list(APPEND target_deps ${targets})
    endif()
  endforeach()

  if(ARGS_BUILD)
    set(lib config_sh_build_interface)
  else()
    set(lib config_sh_install_interface)
  endif()
  add_library(${lib} INTERFACE)
  target_link_libraries(${lib} INTERFACE ${target_deps})

  # ld library path
  foreach(target_dep IN LISTS target_deps)
    get_property(imported TARGET ${target_dep} PROPERTY IMPORTED)
    if(imported)
      get_property(ld_path TARGET ${target_dep} PROPERTY LOCATION)
      get_filename_component(ld_path ${ld_path} DIRECTORY)
      list(APPEND ld_library_path ${ld_path})
    endif()
  endforeach()

  if(ARGS_BUILD)
    list(APPEND ld_library_path "${CMAKE_BINARY_DIR}/${PROJECT_INSTALL_LIBDIR}")
  else()
    list(APPEND ld_library_path "${CMAKE_INSTALL_PREFIX}/${PROJECT_INSTALL_LIBDIR}")
  endif()
  fairroot_remove_system_paths(ld_library_path)
  if(ld_library_path)
    list(REMOVE_DUPLICATES ld_library_path)
    file(TO_NATIVE_PATH "${ld_library_path}" ld_library_path)
    string(REPLACE ";" ":" ld_library_path "${ld_library_path}")
    if(APPLE)
      file(APPEND ${file} "export DYLD_LIBRARY_PATH=\"${ld_library_path}:\$DYLD_LIBRARY_PATH\"\n")
    else()
      file(APPEND ${file} "export LD_LIBRARY_PATH=\"${ld_library_path}:\$LD_LIBRARY_PATH\"\n")
    endif()
  endif()

  # path
  if(ARGS_PATH)
    fairroot_remove_system_paths(ARGS_PATH)
    file(TO_NATIVE_PATH "${ARGS_PATH}" path)
    string(REPLACE ";" ":" path "${path}")
    file(APPEND ${file} "export PATH=\"${path}:\$PATH\"\n")
  endif()

  # data dirs
  foreach(dep IN LISTS ARGS_DEPENDENCIES)
    if(dep STREQUAL Geant4)
      geant4_generate_data_config(OUTVAR Geant4_DATA_CONFIG_SH)
      file(APPEND ${file} "${Geant4_DATA_CONFIG_SH}")
    elseif(dep STREQUAL Geant3)
      file(APPEND ${file} "export G3SYS=\"${Geant3_SYSTEM_DIR}\"\n")
    elseif(dep STREQUAL Pythia8)
      file(APPEND ${file} "export PYTHIA8DATA=\"${Pythia8_DATA}\"\n")
    elseif(dep STREQUAL ROOT)
      get_filename_component(ROOT_PREFIX ${ROOT_BINARY_DIR}/.. ABSOLUTE)
      file(APPEND ${file} "export ROOTSYS=\"${ROOT_PREFIX}\"\n")
    endif()
  endforeach()
  if(ARGS_VMCWORKDIR)
    file(APPEND ${file} "export VMCWORKDIR=\"${ARGS_VMCWORKDIR}\"\n")
  endif()

  # root include dirs
  if(ARGS_BUILD)
    set(rgenfile ${file}.gen_root_include_path)
    file(GENERATE OUTPUT ${rgenfile} CONTENT
      "export ROOT_INCLUDE_PATH=\"$<JOIN:$<REMOVE_DUPLICATES:$<TARGET_PROPERTY:${lib},INTERFACE_INCLUDE_DIRECTORIES>>,:>:$ROOT_INCLUDE_PATH\"\n")
  else()
  endif()
  if(ARGS_BUILD)
    set(name update_root_include_path_build_interface)
  else()
    set(name update_root_include_path_install_interface)
  endif()
  add_custom_target(${name} ALL
    sed -i '/^export ROOT_INCLUDE_PATH=/d' ${file}
    COMMAND cat ${rgenfile} >> ${file}
    DEPENDS ${file} ${rgenfile})
endfunction()

function(fairroot_copy)
  cmake_parse_arguments(ARGS "" "TO" "HEADERS" ${ARGN})

  if(NOT ARGS_TO)
    set(ARGS_TO "${CMAKE_BINARY_DIR}/${PROJECT_INSTALL_INCDIR}")
  endif()
  file(MAKE_DIRECTORY ${ARGS_TO})

  foreach(header IN LISTS ARGS_HEADERS)
    message(STATUS "COPY ${header} TO ${ARGS_TO}")
    configure_file(${header} ${ARGS_TO} COPYONLY)
  endforeach()
endfunction()
