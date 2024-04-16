################################################################################
#   Copyright (C) 2024 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH     #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

# Test suite for the `FairRoot::GeoBase` library

set(suite GeoBase)
# set(exe "${test_prefix}_${suite}")
# add_executable(${exe}
#   test_.cxx
# )
# target_link_libraries(${exe} PRIVATE FairRoot::GeoBase Catch2::Catch2WithMain)
# catch_discover_tests(${exe} TEST_PREFIX "${suite}::")

set(checks_lib "${suite}_checks")
add_library(${checks_lib} INTERFACE)
target_link_libraries(${checks_lib} INTERFACE FairRoot::GeoBase Catch2::Catch2 fmt::fmt)
target_include_directories(${checks_lib} INTERFACE ${CMAKE_SOURCE_DIR})
fairroot_library_settings(${checks_lib})
