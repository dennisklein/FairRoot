/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
/**
 * FairTutGeaneTr.h
 *
 * @since 2019-09-10
 * @author R. Karabowicz
 */

#ifndef FAIRTUTGEANETR_H
#define FAIRTUTGEANETR_H 1

class FairGeanePro;
class FairTutGeanePoint;

#include "FairTask.h"

#include <memory>
#include <Rtypes.h>
#include <vector>

class FairTutGeaneTr : public FairTask
{
 public:
  FairTutGeaneTr();
  FairTutGeaneTr(const FairTutGeaneTr&) = delete;
  FairTutGeaneTr& operator=(const FairTutGeaneTr&) = delete;
  FairTutGeaneTr(FairTutGeaneTr&&) = delete;
  FairTutGeaneTr& operator=(FairTutGeaneTr&&) = delete;
  ~FairTutGeaneTr();

  InitStatus Init() override;

  void Exec(Option_t* opt) override;
 private:
  void Reset();

  /** Finish at the end of run **/
  void Finish() override;

  /** Input array of Points **/
  std::vector<FairTutGeanePoint> fPointArray;

  /** Output array of Hits **/
  std::vector<FairTrackParP> fTrackParIni;   // initial MC track parameters
  std::vector<FairTrackParP> fTrackParFinal; // final MC track parameters
  std::vector<FairTrackParP> fTrackParGeane; // calculated track parameters
  std::vector<FairTrackParP> fTrackParWrong; // demonstrate effect of wrong charge

  Int_t fEvent;
  std::unique_ptr<FairGeanePro> fPro;

  ClassDefOverride(FairTutGeaneTr, 1);
};

#endif //FAIRTUTGEANETR_H
