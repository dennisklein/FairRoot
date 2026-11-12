// ResonanceDecayFilterHooks.h is part of the PYTHIA event generator.
// Copyright (C) 2026 Stephen Mrenna, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Stephen Mrenna, January 2021.

// This class is used to filter out resonance decay products to obtain
// a desired parton-level final state.  It copies an implementation
// already in the CMS software originally authored by Josh Bendavid,
// with some minor modifications.

// The following settings are available:
// flag ResonanceDecayFilter:filter = false
//   Activate filter if true.
// flag ResonanceDecayFilter:exclusive = false
//   Demand exactly the particles requested by the filter if true.
//   If false, additional particles can be present in the resonance decays.
// flag ResonanceDecayFilter:eMuAsEquivalent = false
//   If true, treat the particle id for electron and muon as the same.
// flag ResonanceDecayFilter:eMuTauAsEquivalent = false
//  As above, but include tau leptons.
// flag ResonanceDecayFilter:allNuAsEquivalent = false
// flag ResonanceDecayFilter:udscAsEquivalent = false
//   If true, treat the particle id for all up, down, strange, and charm
//   quarks as the same.</flag>
// flag ResonanceDecayFilter:udscbAsEquivalent = false
//   As above, but include bottom quarks.
// flag ResonanceDecayFilter:wzAsEquivalent = false
//   If true, treat the particle id for W and Z bosons as the same.
// mvec ResonanceDecayFilter:mothers = {}
//   If provided, the set of mothers that will be filtered.
//   Otherwise, all resonance decay products are considered.
// mvec ResonanceDecayFilter:daughters = {}
//   A list of particle ids for resonance decay products that
//   will pass the filter.   If empty, all particles will pass the
//   filter. Beware: an infinite loop is possible if the daughter never
//   appears in resonance decays.

#ifndef Pythia8_ResonanceDecayFilterHooks_H
#define Pythia8_ResonanceDecayFilterHooks_H

// Includes.
#include "Pythia8/Pythia.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/Event.h"
#include "Pythia8/Plugins.h"

namespace Pythia8 {

class ResonanceDecayFilterHooks : public UserHooks {
public:

  // Constructor.
  ResonanceDecayFilterHooks(Pythia*, Settings*, Logger*) : counter(0) {;}

  // Override base class methods.
  bool canVetoResonanceDecays() override {return true;}
  bool doVetoResonanceDecays(Event& process) override;
  bool initAfterBeams() override;

private:

  // Return the particle ID category.
  int idCat(int id);

  // Data members.
  bool filter, exclusive, eMuAsEquivalent, eMuTauAsEquivalent,
    allNuAsEquivalent, udscAsEquivalent, udscbAsEquivalent, wzAsEquivalent;
  unsigned long int counter;
  set<int> mothers;
  vector<int> daughters;
  unordered_map<int, int> requestedDaughters, observedDaughters;

};

//--------------------------------------------------------------------------

// Return a particle ID given equivalence user settings.

int ResonanceDecayFilterHooks::idCat(int id) {
  id = abs(id);
  if (id == 13 && (eMuAsEquivalent || eMuTauAsEquivalent)) id = 11;
  else if (id == 15 && eMuTauAsEquivalent) id = 11;
  else if ((id == 14 || id == 16) && allNuAsEquivalent) id = 12;
  else if ((id == 2 || id == 3 || id == 4) && udscAsEquivalent) id = 1;
  else if ((id == 2 || id == 3 || id == 4 || id == 5) &&
           udscbAsEquivalent) id = 1;
  else if ((id == 23 || id == 24) && wzAsEquivalent) id = 23;
  return id;
}

//--------------------------------------------------------------------------

// Intialize the user hook after the beams.

bool ResonanceDecayFilterHooks::initAfterBeams() {
  filter                = flag("ResonanceDecayFilter:filter");
  exclusive             = flag("ResonanceDecayFilter:exclusive");
  eMuAsEquivalent       = flag("ResonanceDecayFilter:eMuAsEquivalent");
  eMuTauAsEquivalent    = flag("ResonanceDecayFilter:eMuTauAsEquivalent");
  allNuAsEquivalent     = flag("ResonanceDecayFilter:allNuAsEquivalent");
  udscAsEquivalent      = flag("ResonanceDecayFilter:udscAsEquivalent");
  udscbAsEquivalent     = flag("ResonanceDecayFilter:udscbAsEquivalent");
  wzAsEquivalent        = flag("ResonanceDecayFilter:wzAsEquivalent");
  vector<int> mothersIn = mvec("ResonanceDecayFilter:mothers");
  mothers.clear();
  mothers.insert(mothersIn.begin(), mothersIn.end());
  daughters = mvec("ResonanceDecayFilter:daughters");
  requestedDaughters.clear();

  // Loop over the daughters.
  for (int id : daughters) ++requestedDaughters[idCat(id)];
  return true;

}

//--------------------------------------------------------------------------

// Return true of the resonance decays are vetoed.

bool ResonanceDecayFilterHooks::doVetoResonanceDecays(Event &process) {
  if (!filter) return false;

  // Count the number of times hook is called.
  counter++;
  settingsPtr->mode("ResonanceDecayFilter:counter", counter);
  observedDaughters.clear();

  // Loop over particles and determine equivalent types.
  for (int i = 0; i < process.size(); ++i) {
    const Particle &p = process[i];
    int mid = p.mother1() > 0 ? abs(process[p.mother1()].id()) : 0;

    // If no list of mothers is provided, then all particles
    // in hard process and resonance decays are counted together
    if (mothers.empty() || mothers.count(mid) || mothers.count(-mid))
      ++observedDaughters[idCat(p.id())];
  }

  // Check if criteria is satisfied.
  // inclusive mode: at least as many decay products as requested
  // exclusive mode: exactly as many decay products as requested
  // (but additional particle types not appearing in the list of requested
  // daughter id's are ignored)
  for (const auto &reqDau : requestedDaughters) {
    int reqId = reqDau.first;
    int reqCount = reqDau.second;
    auto obsItr = observedDaughters.find(reqId);
    int obsCount = ( obsItr != observedDaughters.end() ) ? obsItr->second : 0;

    // Inclusive criteria not satisfied, veto event
    if (obsCount < reqCount) return true;

    // Exclusive criteria not satisfied, veto event
    if (exclusive && obsCount > reqCount) return true;
  }

  // All criteria satisfied, don't veto
  return false;

}

//--------------------------------------------------------------------------

// Register settings.

void resonanceSettings(Settings *settingsPtr) {
  settingsPtr->addFlag("ResonanceDecayFilter:filter", false);
  settingsPtr->addFlag("ResonanceDecayFilter:exclusive", false);
  settingsPtr->addFlag("ResonanceDecayFilter:eMuAsEquivalent", false);
  settingsPtr->addFlag("ResonanceDecayFilter:eMuTauAsEquivalent", false);
  settingsPtr->addFlag("ResonanceDecayFilter:allNuAsEquivalent", false);
  settingsPtr->addFlag("ResonanceDecayFilter:udscAsEquivalent", false);
  settingsPtr->addFlag("ResonanceDecayFilter:udscbAsEquivalent", false);
  settingsPtr->addFlag("ResonanceDecayFilter:wzAsEquivalent", false);
  settingsPtr->addMVec("ResonanceDecayFilter:mothers", {}, false, false, 0, 0);
  settingsPtr->addMVec("ResonanceDecayFilter:daughters", {},
    false, false, 0, 0);
  settingsPtr->addMode("ResonanceDecayFilter:counter", 0, false, false, 0, 0);
}

//--------------------------------------------------------------------------

// Declare the plugin.

PYTHIA8_PLUGIN_CLASS(UserHooks, ResonanceDecayFilterHooks, false, false, false)
PYTHIA8_PLUGIN_SETTINGS(resonanceSettings)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_ResonanceDecayFilterHooks_H
