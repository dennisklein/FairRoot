/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "FairRuntimeDb.h"
#include "FairTutGeaneContFact.h"
#include "FairTutGeaneGeoPar.h"

#include <memory>
#include <string.h> // for strcmp

static FairTutGeaneContFact gFairTutGeaneContFact;

/** Constructor (called when the library is loaded) */
FairTutGeaneContFact::FairTutGeaneContFact()
    : FairContFact()
{
    fName = "FairTutGeaneContFact";
    fTitle = "Factory for parameter containers in libFairTutGeane";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

/** Creates the Container objects with all accepted
    contexts and adds them to
    the list of containers for the FairTutGeane library.
*/
void FairTutGeaneContFact::setAllContainers()
{
    std::unique_ptr<FairContainer> p(new FairContainer(
        "FairTutGeaneGeoPar",
        "FairTutGeane Geometry Parameters",
        "TestDefaultContext"));
    p->addContext("TestNonDefaultContext");
    containers->Add(p.release());
}

/** Calls the constructor of the corresponding parameter container.
    For an actual context, which is not an empty string and not
    the default context of this container, the name is concatinated
    with the context.
*/
FairParSet* FairTutGeaneContFact::createContainer(FairContainer* c)
{
    const char* name = c->GetName();
    if (strcmp(name, "FairTutGeaneGeoPar") == 0)
    {
        return new FairTutGeaneGeoPar(c->getConcatName().Data(), c->GetTitle(),
                                      c->getContext());
    }
    return nullptr;
}

ClassImp(FairTutGeaneContFact)
