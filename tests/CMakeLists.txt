################################################################################
# Copyright (C) 2014-2024 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH  #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

# Shall be prefixed to all test executables
set(test_prefix test)

add_subdirectory(base)
add_subdirectory(tgeo)
add_subdirectory(examples/simulation/Tutorial1)
add_subdirectory(examples/simulation/Tutorial2)
add_subdirectory(examples/simulation/Tutorial4)
add_subdirectory(examples/advanced/Tutorial3)
add_subdirectory(examples/advanced/propagator)
add_subdirectory(examples/MQ/pixelDetector)
add_subdirectory(geobase)
