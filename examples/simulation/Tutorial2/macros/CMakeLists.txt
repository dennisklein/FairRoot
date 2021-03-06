################################################################################
# Copyright (C) 2014-2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/create_digis_mixed.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/create_digis.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/read_digis.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/run_background.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/run_signal.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/run_tutorial2.C)

set(maxTestTime 60)

add_test(ex_sim_tutorial2 ${CMAKE_CURRENT_BINARY_DIR}/run_tutorial2.sh)
set_tests_properties(ex_sim_tutorial2 PROPERTIES
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

add_test(ex_sim_tutorial2_create_digis ${CMAKE_CURRENT_BINARY_DIR}/create_digis.sh)
set_tests_properties(ex_sim_tutorial2_create_digis PROPERTIES
  DEPENDS ex_sim_tutorial2
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

add_test(ex_sim_tutorial2_read_digis ${CMAKE_CURRENT_BINARY_DIR}/read_digis.sh)
set_tests_properties(ex_sim_tutorial2_read_digis PROPERTIES
  DEPENDS ex_sim_tutorial2_create_digis
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

add_test(ex_sim_tutorial2_run_background ${CMAKE_CURRENT_BINARY_DIR}/run_background.sh)
set_tests_properties(ex_sim_tutorial2_run_background PROPERTIES
  DEPENDS ex_sim_tutorial2_read_digis
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

add_test(ex_sim_tutorial2_run_signal1 ${CMAKE_CURRENT_BINARY_DIR}/run_signal.sh 1 10)
set_tests_properties(ex_sim_tutorial2_run_signal1 PROPERTIES
  DEPENDS ex_sim_tutorial2_run_background
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

add_test(ex_sim_tutorial2_run_signal2 ${CMAKE_CURRENT_BINARY_DIR}/run_signal.sh 2 20)
set_tests_properties(ex_sim_tutorial2_run_signal2 PROPERTIES
  DEPENDS ex_sim_tutorial2_run_signal1
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

add_test(ex_sim_tutorial2_create_digis_mixed ${CMAKE_CURRENT_BINARY_DIR}/create_digis_mixed.sh)
set_tests_properties(ex_sim_tutorial2_create_digis_mixed PROPERTIES
  DEPENDS ex_sim_tutorial2_run_signal2
  TIMEOUT ${maxTestTime}
  PASS_REGULAR_EXPRESSION "Macro finished successfully"
)

install(FILES
  create_digis_mixed.C
  create_digis.C
  read_digis.C
  run_background.C
  run_signal.C
  run_tutorial2.C
  DESTINATION ${PROJECT_INSTALL_DATADIR}/examples/simulation/Tutorial2
)
