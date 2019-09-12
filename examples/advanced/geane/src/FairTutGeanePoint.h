/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#ifndef FAIRTUTGEANEPOINT_H
#define FAIRTUTGEANEPOINT_H 1

class TVector3;

#include "FairMCPoint.h"
#include <Rtypes.h>

class FairTutGeanePoint : public FairMCPoint
{
  public:
    // FIXME FairMCPoint should delete its copy ctor,
    // then this whole section can be dropped
    FairTutGeanePoint() = default;
    FairTutGeanePoint(const FairTutGeanePoint&) = delete;
    FairTutGeanePoint& operator=(const FairTutGeanePoint&) = delete;
    FairTutGeanePoint(FairTutGeanePoint&&) = delete;
    FairTutGeanePoint& operator=(FairTutGeanePoint&&) = delete;
    ~FairTutGeanePoint() = default;

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param pos      Ccoordinates at entrance to active volume [cm]
     *@param mom      Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     **/
    FairTutGeanePoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom,
                      Double_t tof, Double_t length, Double_t eLoss);

    /** Output to screen **/
    void Print(const Option_t* opt) const override;

    ClassDefOverride(FairTutGeanePoint,1)
};

#endif
