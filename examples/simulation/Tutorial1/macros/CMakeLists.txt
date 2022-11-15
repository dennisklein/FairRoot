################################################################################
# Copyright (C) 2014-2022 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

foreach(macro
        "run_tutorial1_mesh.C"
        "run_tutorial1_pythia6.C"
        "run_tutorial1_pythia8.C"
        "run_tutorial1_urqmd.C"
        "run_tutorial1.C")
  GENERATE_ROOT_TEST_SCRIPT(${CMAKE_CURRENT_SOURCE_DIR}/${macro})
  install(FILES ${macro}
          DESTINATION ${PROJECT_INSTALL_DATADIR}/examples/simulation/Tutorial1)
endforeach()

foreach(mce IN LISTS mcEngine_list)
  add_ex_tut1_test(VARIANT "mesh"    MCENGINE ${mce} N 10)
  add_ex_tut1_test(VARIANT "pythia6" MCENGINE ${mce} N 10)
  add_ex_tut1_test(VARIANT "pythia8" MCENGINE ${mce} N 10)
  add_ex_tut1_test(VARIANT "urqmd"   MCENGINE ${mce} N 2)
  add_ex_tut1_test(VARIANT "plain"   MCENGINE ${mce} N 10 BASE_SCRIPT)
endforeach()
add_ex_tut1_test(VARIANT "PostInitConfig" MCENGINE "TGeant4" N 1 BASE_SCRIPT
                 ARGS false true PASS_REGEX "Loading Geant4 PostInit Config.")
