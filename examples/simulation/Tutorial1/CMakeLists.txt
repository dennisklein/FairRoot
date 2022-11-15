################################################################################
# Copyright (C) 2014-2022 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

#
# add_ex_tut1_test(MCENGINE <mce> VARIANT <variant> N <n>
#                 [BASE_SCRIPT] [PASS_REGEX <regex>]
#                 [ARGS <arg1> [<arg2> ...]])
#
# Defines a new test case with name
#
#    ex_sim_tutorial1_<variant>_<mce>
#
# that runs the command
#
#    macros/run_tutorial1[_<variant>].C <n> "<mce>" [<arg1> [<arg2> ...]]
#                        ^~~~~~~~~~~~
#                        omitted if BASE_SCRIPT option set
#
function(add_ex_tut1_test)
  cmake_parse_arguments(PARSE_ARGV 0 ARG
    "BASE_SCRIPT" "VARIANT;MCENGINE;N;PASS_REGEX" "ARGS")
  set(test_prefix "ex_sim_tutorial1")
  if(NOT ARG_BASE_SCRIPT)
    set(script "_${ARG_VARIANT}")
  endif()
  set(test "${test_prefix}_${ARG_VARIANT}_${ARG_MCENGINE}")
  add_test(NAME ${test}
           COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_tutorial1${script}.sh
                   ${ARG_N} \"${ARG_MCENGINE}\" ${ARG_ARGS}
  )
  if(ARG_PASS_REGEX)
    set(regex ${ARG_PASS_REGEX})
  else()
    set(regex "Macro finished successfully")
  endif()
  set_tests_properties(${test} PROPERTIES
                       TIMEOUT 60
                       PASS_REGULAR_EXPRESSION "${regex}"
  )
endfunction()

add_subdirectory(src)
add_subdirectory(macros)
