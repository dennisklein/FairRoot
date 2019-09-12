################################################################################
#    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/runMC.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/runGeane.C)
GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/runPull.C)

set(MaxTestTime 60)
set(PassExpr "Macro finished successfully")

add_test(runGeaneMC ${CMAKE_CURRENT_BINARY_DIR}/runMC.sh)
set_tests_properties(runGeaneMC PROPERTIES
  TIMEOUT ${MaxTestTime}
  PASS_REGULAR_EXPRESSION ${PassExpr})

add_test(runGeane ${CMAKE_CURRENT_BINARY_DIR}/runGeane.sh)
set_tests_properties(runGeane PROPERTIES
  DEPENDS runGeaneMC
  TIMEOUT ${MaxTestTime}
  PASS_REGULAR_EXPRESSION ${PassExpr})

add_test(runGeanePull ${CMAKE_CURRENT_BINARY_DIR}/runPull.sh)
set_tests_properties(runGeanePull PROPERTIES
  DEPENDS runGeane
  TIMEOUT ${MaxTestTime}
  PASS_REGULAR_EXPRESSION ${PassExpr})
