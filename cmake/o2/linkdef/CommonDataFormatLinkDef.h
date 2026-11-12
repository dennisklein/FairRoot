// Reduced from DataFormats/common/src/CommonDataFormatLinkDef.h by tools/gen-o2-cmake.py.
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

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::dataformats::TimeStamp < float> + ;
#pragma link C++ class o2::dataformats::TimeStamp < double> + ;
#pragma link C++ class o2::dataformats::TimeStamp < int> + ;
#pragma link C++ class o2::dataformats::TimeStamp < uint32_t> + ;
#pragma link C++ class o2::dataformats::TimeStamp < Float16_t> + ;
#pragma link C++ class o2::dataformats::TimeStampWithError < float, float> + ;
#pragma link C++ class o2::dataformats::TimeStampWithError < double, double> + ;
#pragma link C++ class o2::dataformats::TimeStampWithError < int, int> + ;
#pragma link C++ class o2::dataformats::TimeStampWithError < uint32_t, uint32_t> + ;






#pragma link C++ class o2::InteractionRecord + ;
#pragma link C++ class o2::InteractionTimeRecord + ;
#pragma link C++ class o2::BunchFilling + ;
#pragma link C++ class std::pair < long, o2::BunchFilling> + ;







#endif
