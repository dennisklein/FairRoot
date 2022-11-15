################################################################################
# Copyright (C) 2014-2022 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

function(add_ex_tut3_test)
  cmake_parse_arguments(PARSE_ARGV 0 ARG
    "LONGER" "STAGE;MCENGINE;DEPENDS;N" "")
  set(test_prefix "ex_tutorial3")
  set(test "${test_prefix}_${stage}_${ARG_MCENGINE}")
  add_test(NAME ${test}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_${stage}.sh
            ${ARG_N} \"${ARG_MCENGINE}\"
  )
  set_tests_properties(${test} PROPERTIES
    FIXTURES_SETUP ${test}
    PASS_REGULAR_EXPRESSION "Macro finished successfully"
  )
  set(timeout 60)
  if(ARG_LONGER_TIMEOUT)
    math(EXPR timeout 3*${timeout})
  endif()
  set_tests_properties(${test} PROPERTIES TIMEOUT ${timeout})
  if(ARG_DEPENDS)
    set_tests_properties(${test} PROPERTIES
      FIXTURES_REQUIRED "${test_prefix}_${ARG_DEPENDS}_${ARG_MCENGINE}"
    )
  endif()
endfunction()

set(stages "sim" "digi" "reco" "digi_timebased" "reco_timebased")
set(           sim_args N 100 LONGER)
set(          digi_args DEPENDS "sim")
set(          reco_args DEPENDS "digi")
set(digi_timebased_args DEPENDS "sim")
set(reco_timebased_args DEPENDS "digi_timebased")
if(BUILD_PROOF_SUPPORT)
  list(APPEND stages "digi_reco_proof")
  set(digi_reco_proof_args N 1 DEPENDS "sim")
endif()

set(tut3_macros eventDisplay.C)
foreach(stage IN LISTS stages)
  GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/run_${stage}.C)
  list(APPEND tut3_macros "run_${stage}.C")
  foreach(mc_engine IN LISTS mcEngine_list)
    add_ex_tut3_test(STAGE ${stage} MCENGINE ${mc_engine} ${${stage}_args})
  endforeach()
endforeach()

set(tut3_install_dir "${PROJECT_INSTALL_DATADIR}/examples/advanced/Tutorial3")
install(FILES ${tut3_macros} DESTINATION ${tut3_install_dir})
install(FILES data/dummy.txt DESTINATION ${tut3_install_dir}/macro/data)
