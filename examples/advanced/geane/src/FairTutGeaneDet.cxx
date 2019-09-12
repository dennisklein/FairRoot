/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "FairDetectorList.h" // for DetectorId::kTutDet
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoVolume.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairStack.h"
#include "FairTutGeaneDet.h"
#include "FairTutGeanePoint.h"
#include "FairTutGeaneGeo.h"
#include "FairTutGeaneGeoPar.h"
#include "FairVolume.h"

#include <iostream>
#include <memory>
#include <TClonesArray.h>
#include <TList.h> // for TListIter, TList (ptr only)
#include <TObjArray.h>
#include <TString.h>
#include <TVirtualMC.h>
#include <TVirtualMCStack.h>
#include <TVector3.h>
#include <TGeoPhysicalNode.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>

using std::cout;
using std::endl;

FairTutGeaneDet::FairTutGeaneDet()
  : FairTutGeaneDet("FairTutGeaneDet", kTRUE)
{
}

FairTutGeaneDet::FairTutGeaneDet(const char* name, Bool_t active)
  : FairDetector(name, active, kTutGeane),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fMisalignDetector(kFALSE),
    fFairTutGeanePointCollection(new TClonesArray("FairTutGeanePoint"))
{
}

FairTutGeaneDet::~FairTutGeaneDet()
{
    if (fFairTutGeanePointCollection) {
        fFairTutGeanePointCollection->Delete();
        delete fFairTutGeanePointCollection;
    }
}

Bool_t FairTutGeaneDet::ProcessHits(FairVolume* vol)
{
    /** This method is called from the MC stepping */

    // Set parameters at entrance of volume. Reset ELoss.
    if (TVirtualMC::GetMC()->IsTrackEntering() ) {
        fELoss  = 0.;
        fTime   = TVirtualMC::GetMC()->TrackTime() * 1.0e09;
        fLength = TVirtualMC::GetMC()->TrackLength();
        TVirtualMC::GetMC()->TrackPosition(fPos);
        TVirtualMC::GetMC()->TrackMomentum(fMom);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += TVirtualMC::GetMC()->Edep();

    // Create FairTutGeanePoint at exit of active volume
    if (TVirtualMC::GetMC()->IsTrackExiting()    ||
        TVirtualMC::GetMC()->IsTrackStop()       ||
        TVirtualMC::GetMC()->IsTrackDisappeared())
    {
        fTrackID  = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();

        if (fELoss == 0.) { return kFALSE; }

        // FIXME Is this still needed?
        // Taking stationNr and sectorNr from string is almost effortless.
        // Simulation of 100k events with 5 pions without magnetic field takes:
        // - Real time 142.366 s, CPU time 140.32s WITH USING VolPath TO GET fVolumeID
        // - Real time 142.407 s, CPU time 140.64s WITHOUT THE FOLLOWING TString OPERATIONS
        // {
            // TString detPath = TVirtualMC::GetMC()->CurrentVolPath();
            // detPath.Remove (0,detPath.Last('/')+1);
            // detPath.Remove (0,detPath.First("Pixel")+5);
            // Int_t stationNr = detPath.Atoi();
            // detPath.Remove (0,detPath.First("_")+1);
            // Int_t sectorNr  = detPath.Atoi();
            // fVolumeID = stationNr*256+sectorNr;
        // }

        AddHit(fTrackID, fVolumeID, TVector3(fPos.X(), fPos.Y(), fPos.Z()),
               TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,
               fELoss);

        // Increment number of FairTutGeaneDet points in TParticle
        FairStack* stack = static_cast<FairStack*>(TVirtualMC::GetMC()->GetStack());
        stack->AddPoint(kTutGeane);
    }

    return kTRUE;
}

void FairTutGeaneDet::EndOfEvent()
{
    fFairTutGeanePointCollection->Clear();
}

void FairTutGeaneDet::Register()
{
  /** This will create a branch in the output tree called
      FairTutGeanePoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("FairTutGeanePoint", "FairTutGeaneDet",
                                        fFairTutGeanePointCollection, kTRUE);
}

TClonesArray* FairTutGeaneDet::GetCollection(Int_t iColl) const
{
    if (iColl == 0) {
        return fFairTutGeanePointCollection;
    }
    return nullptr;
}

void FairTutGeaneDet::Reset()
{
    fFairTutGeanePointCollection->Clear();
}

void FairTutGeaneDet::ConstructGeometry()
{
    /** If you are using the standard ASCII input for the geometry
        just copy this and use it for your detector, otherwise you can
        implement here you own way of constructing the geometry. */

    ConstructASCIIGeometry<FairTutGeaneGeo, FairTutGeaneGeoPar>(nullptr, "FairTutGeaneGeoPar");
}

FairTutGeanePoint *FairTutGeaneDet::AddHit(Int_t trackID, Int_t detID,
                                           TVector3 pos, TVector3 mom,
                                           Double_t time, Double_t length,
                                           Double_t eLoss)
{
    Int_t size = fFairTutGeanePointCollection->GetEntriesFast();
    return new ((*fFairTutGeanePointCollection)[size])
        FairTutGeanePoint(trackID, detID, pos, mom, time, length, eLoss);
}

ClassImp(FairTutGeaneDet)
