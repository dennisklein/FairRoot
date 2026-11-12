// LHEFHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Christian T. Preuss, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Christian T. Preuss.

// This class implements an interface to write LHE files that can be loaded
// via the plugin structure. It can be run with PythiaParallel.

#ifndef Pythia8_LHEFHooks_H
#define Pythia8_LHEFHooks_H

// Pythia includes.
#include "Pythia8/Pythia.h"
#include "Pythia8/Plugins.h"

namespace Pythia8 {

//==========================================================================

// UserHook to write LHE files.

class LHEFHooks : public UserHooks {

public:

  // Constructors and destructor.
  LHEFHooks() {}
  LHEFHooks(Pythia* pythiaPtrIn, Settings*, Logger*) :
    pythiaPtr(pythiaPtrIn) {}
  ~LHEFHooks() {if (lhefWriterPtr != nullptr) delete lhefWriterPtr;}

  //--------------------------------------------------------------------------

  // Print event to LHE file.
  void onEndEvent(Status) {

    // Create the LHEF converter.
    if (lhefWriterPtr == nullptr) {
      // Create a LHEF writer.
      lhefWriterPtr = new LHEF3FromPythia8(&pythiaPtr->event,
        &pythiaPtr->info);

      // If run in parallel, every thread writes its own file.
      string filename = word("LHEF:filename");
      int idx = mode("Parallelism:index");
      if (idx >= 0 && mode("Parallelism:numThreads") > 1) {
        size_t iSuffix = filename.find(".lhe");
        if (iSuffix != string::npos)
          filename = filename.substr(0, iSuffix);
        filename = filename + "_" + to_string(idx) + ".lhe";
      }
      lhefWriterPtr->openLHEF(filename);
      lhefWriterPtr->setInit();
    }

    // Write event.
    lhefWriterPtr->setEvent();

  }

  //--------------------------------------------------------------------------

  // Finalise.
  void onStat() {
    if (lhefWriterPtr != nullptr) lhefWriterPtr->closeLHEF(true);
  }

  //--------------------------------------------------------------------------

  // Close files.
  void onStat(vector<PhysicsBase*> hookPtrs, Pythia* pythiaPtrIn) {
    // TODO: merge the LHE files into one.
    for (int iPtr = 0; iPtr < (int)hookPtrs.size(); ++iPtr) {
      LHEFHooks* hookNow = dynamic_cast<LHEFHooks*>(hookPtrs[iPtr]);
      if (hookNow == nullptr) {
        loggerPtr->ERROR_MSG("could not retrieve LHEFHooks for thread ",
          toString(iPtr));
        return;
      }
      hookNow->onStat();
    }

    onStat();
  }

  //--------------------------------------------------------------------------

 private:

  Pythia* pythiaPtr{};
  LHEF3FromPythia8* lhefWriterPtr{nullptr};

};

//--------------------------------------------------------------------------

// Register LHEF settings.

void lhefSettings(Settings *settingsPtr) {
  settingsPtr->addWord("LHEF:fileName", "events.lhe");
}

//--------------------------------------------------------------------------

// Declare the plugin.

PYTHIA8_PLUGIN_CLASS(UserHooks, LHEFHooks, true, false, false)
PYTHIA8_PLUGIN_SETTINGS(lhefSettings)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_LHEFHooks_H
