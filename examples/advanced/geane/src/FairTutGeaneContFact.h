/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#ifndef FAIRTUTGEANECONTFACT_H_
#define FAIRTUTGEANECONTFACT_H_

#include "FairContFact.h"

#include <Rtypes.h> // for ClassDef

// Factory for all FairTutGeane parameter containers
class FairTutGeaneContFact : public FairContFact
{
  public:
    FairTutGeaneContFact();

    /// @return Transfers ownership of the returned object to the caller
    FairParSet* createContainer(FairContainer*) override;

  private:
    void setAllContainers();

    ClassDefOverride(FairTutGeaneContFact, 0)
};

#endif /* FAIRTUTGEANECONTFACT_H_ */
