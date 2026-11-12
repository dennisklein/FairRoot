// Reduced from Generators/src/GeneratorsLinkDef.h by tools/gen-o2-cmake.py.
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
// -------------------------------------------------------------------------
// -----                  M. Al-Turany   June 2014                     -----
// -------------------------------------------------------------------------

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::eventgen::Generator + ;
#pragma link C++ class o2::eventgen::GeneratorGeantinos + ;
#ifdef GENERATORS_WITH_TPCLOOPERS
#pragma link C++ class o2::eventgen::GenTPCLoopers + ;
#endif
#ifdef GENERATORS_WITH_HEPMC3
#endif
#ifdef GENERATORS_WITH_PYTHIA6
#endif
#ifdef GENERATORS_WITH_PYTHIA8
#pragma link C++ class o2::eventgen::DecayerPythia8 + ;
#pragma link C++ class o2::eventgen::DecayerPythia8Param + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::eventgen::DecayerPythia8Param> + ;
#pragma link C++ class o2::eventgen::EventDecayerParam + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::eventgen::EventDecayerParam> + ;

#endif
#if defined(GENERATORS_WITH_PYTHIA8) && defined(GENERATORS_WITH_HEPMC3)
#endif
#pragma link C++ class o2::eventgen::PrimaryGenerator + ;
#pragma link C++ class o2::eventgen::PrimaryGeneratorParam + ;

#pragma link C++ class o2::eventgen::BoxGunParam + ;
#pragma link C++ class o2::eventgen::BoxGenConfig + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::eventgen::BoxGunParam> + ;

#pragma link C++ class o2::eventgen::BoxGenerator + ;

#endif
