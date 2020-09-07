/********************************************************************************
 *    Copyright (C) 2020 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "FairEveTransparencyControl.h"

#include <TGNumberEntry.h>

FairEveTransparencyControl::FairEveTransparencyControl(TGFrame const *parent, char const *label)
    : TGHorizontalFrame(parent)
    , fCheck(new TGCheckButton(this, label))
    , fNumber(new TGNumberEntry(this,
                                80.,   // initial number
                                6,     // digitwidth
                                -1,
                                TGNumberFormat::kNESInteger,
                                TGNumberFormat::kNEANonNegative,
                                TGNumberFormat::kNELLimitMinMax,
                                0,      // min
                                100))   // max
{
    SetCleanup(kDeepCleanup);

    // display
    AddFrame(fCheck);    // takes ownership
    AddFrame(fNumber);   // takes ownership

    // wire up observers
    fCheck->Connect("Toggled(Bool_t)", this->ClassName(), this, "Toggled(Bool_t)");
    fNumber->Connect("ValueSet(Long_t)", this->ClassName(), this, "ValueSet()");
}

void FairEveTransparencyControl::Toggled(Bool_t on) { SetSignal(on ? fNumber->GetIntNumber() : 0); }

void FairEveTransparencyControl::ValueSet()
{
    if (fCheck->IsOn()) {
        SetSignal(fNumber->GetIntNumber());
    }
}

void FairEveTransparencyControl::SetSignal(Int_t transparency) { Emit("SetSignal(Int_t)", transparency); }
