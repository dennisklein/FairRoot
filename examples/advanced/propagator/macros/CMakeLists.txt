 ################################################################################
 # Copyright (C) 2019-2022 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  #
 #                                                                              #
 #              This software is distributed under the terms of the             #
 #              GNU Lesser General Public Licence (LGPL) version 3,             #
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################

 function(add_ex_prop_test)
  cmake_parse_arguments(PARSE_ARGV 0 ARG "" "STAGE;DEPENDS" "")
  set(test_prefix "ex_propagator")
  set(test "${test_prefix}_${ARG_STAGE}")
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run${ARG_STAGE}.sh)
  set_tests_properties(${test} PROPERTIES
    FIXTURES_SETUP ${test}
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "Macro finished successfully"
  )
  if(ARG_DEPENDS)
    set_tests_properties(${test} PROPERTIES
      FIXTURES_REQUIRED "${test_prefix}_${ARG_DEPENDS}"
    )
  endif()
endfunction()

GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/runMC.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/runProp.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/runPull.C)

add_ex_prop_test(STAGE "MC")
add_ex_prop_test(STAGE "Prop" DEPENDS "MC")
add_ex_prop_test(STAGE "Pull" DEPENDS "Prop")
