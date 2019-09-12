################################################################################
#    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
#                                                                              #
#              This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (LGPL) version 3,             #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/examples/advanced/geane/src
  ${CMAKE_SOURCE_DIR}/examples/common/mcstack
  ${CMAKE_SOURCE_DIR}/geane
  ${CMAKE_SOURCE_DIR}/trackbase
  ${FairLogger_INCDIR}
)

set(SYSTEM_INCLUDE_DIRECTORIES
  ${SYSTEM_INCLUDE_DIRECTORIES}
  ${Geant3_INCLUDE_DIRS}
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
  FairTutGeaneContFact.cxx
  FairTutGeaneDet.cxx
  FairTutGeaneField.cxx
  FairTutGeaneGeo.cxx
  FairTutGeaneGeoPar.cxx
  FairTutGeanePoint.cxx
  FairTutGeaneTr.cxx
)

set(HEADERS
  FairTutGeaneContFact.h
  FairTutGeaneDet.h
  FairTutGeaneField.h
  FairTutGeaneGeo.h
  FairTutGeaneGeoPar.h
  FairTutGeanePoint.h
  FairTutGeaneTr.h
)

set(LINKDEF TutGeaneLinkDef.h)
set(LIBRARY_NAME GeaneExample)
set(DEPENDENCIES
    Base
    GeoBase
    ParBase
    Geane
    MCStack
    Gen
    TrkBase
)

GENERATE_LIBRARY()
