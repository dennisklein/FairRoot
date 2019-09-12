/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#ifndef FAIRTUTGEANEDET_H
#define FAIRTUTGEANEDET_H

class FairTutGeanePoint;
class FairVolume;
class TClonesArray;
class TVector3;

#include "FairDetector.h"

#include <Rtypes.h> // for Int_t, Double32_t, Double_t, etc
#include <TLorentzVector.h>

class FairTutGeaneDet : public FairDetector
{
  public:
    /** @arg Name Detector Name
     *  @arg Active kTRUE for active detectors (ProcessHits() will be called)
     *              kFALSE for inactive detectors
    */
    FairTutGeaneDet(const char* Name, Bool_t Active);

    FairTutGeaneDet();
    FairTutGeaneDet(const FairTutGeaneDet&) = delete;
    FairTutGeaneDet& operator=(const FairTutGeaneDet&) = delete;
    FairTutGeaneDet(FairTutGeaneDet&&) = delete;
    FairTutGeaneDet& operator=(FairTutGeaneDet&&) = delete;
    ~FairTutGeaneDet();

    /** Initialization of the detector is done here */
    // void Initialize() override;

    /** this method is called for each step during simulation
     *  (see FairMCApplication::Stepping())
    */
    Bool_t ProcessHits(FairVolume* v = nullptr) override;

    /** Registers the produced collections in FairRootManager. */
    void Register() override;

    /** Gets the produced collections */
    TClonesArray* GetCollection(Int_t iColl) const override;

    /** has to be called after each event to reset the containers */
    void Reset() override;

    /** Create the detector geometry */
    void ConstructGeometry() override;

    /** This method is an example of how to add your own point
     *  of type FairTutGeaneDetPoint to the clones array
    */
    FairTutGeanePoint* AddHit(Int_t trackID, Int_t detID,
                              TVector3 pos, TVector3 mom,
                              Double_t time, Double_t length,
                              Double_t eLoss);

    /** The following methods can be implemented if you need to make
     *  any optional action in your detector during the transport.
    */
    virtual void EndOfEvent() override;

    // void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset) override {}
    // void SetSpecialPhysicsCuts() override {}
    // void FinishPrimary() override {}
    // void FinishRun() override {}
    // void BeginPrimary() override {}
    // void PostTrack() override {}
    // void PreTrack() override {}
    // void BeginEvent() override {}

  private:
    /** Track information to be stored until the track leaves the
    active volume.
    */
    Int_t          fTrackID;           //!  track index
    Int_t          fVolumeID;          //!  volume id
    TLorentzVector fPos;               //!  position at entrance
    TLorentzVector fMom;               //!  momentum at entrance
    Double32_t     fTime;              //!  time
    Double32_t     fLength;            //!  length
    Double32_t     fELoss;             //!  energy loss

    /** container for data points */
    Bool_t         fMisalignDetector;
    TClonesArray*  fFairTutGeanePointCollection;

    ClassDefOverride(FairTutGeaneDet, 1)
};

#endif //FAIRTUTGEANEDET_H
