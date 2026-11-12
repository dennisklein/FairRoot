// HepMC3.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
//
// Author: HepMC 3 Collaboration, hepmc-dev@.cern.ch,
//         Philip Ilten, March 2026.
//
// Based on the HepMC2 interface by Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch.
// Header file and function definitions for the Pythia8ToHepMC class,
// which converts a PYTHIA event record to the standard HepMC format.
//
// The Pythia8ToFilteredHepMC3 class which converts a PYTHIA event
// record to the standard HepMC format, but does not use
// GenEvent::add_tree and allows for filtering while keeping the
// correct vertex structure.

#ifndef Pythia8_HepMC3_H
#define Pythia8_HepMC3_H

#ifdef Pythia8_HepMC2_H
#error Cannot include HepMC3.h if HepMC2.h has already been included.
#endif

#include <vector>
#include "Pythia8/Pythia.h"
#include "Pythia8/HIInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/GenHeavyIon.h"
#include "HepMC3/GenPdfInfo.h"

namespace HepMC3 {

//==========================================================================

// A class to convert a Pythia 8 event into a HepMC event. This class
// keeps the entire event structure and is defined in the HepMC3 namespace.
// Typically, users do not directly interact with this class but rather
// the Pythia8::Pythia8ToHepMC3 class instead. This class is based on the
// HepMC2 interface by Mikhail Kirsanov.

class Pythia8ToHepMC3 {

public:

  // Constructor and destructor.
  Pythia8ToHepMC3(): m_internal_event_number(0), m_print_inconsistency(true),
    m_free_parton_warnings(true), m_crash_on_problem(false),
    m_convert_gluon_to_0(false), m_store_pdf(true), m_store_proc(true),
    m_store_xsec(true), m_store_weights(true) {}
  virtual ~Pythia8ToHepMC3() {}

  // The recommended method to convert Pythia events into HepMC3 ones.
  bool fill_next_event( Pythia8::Pythia& pythia, GenEvent* evt,
    int ievnum = -1 ) { return fill_next_event( pythia.event, evt,
    ievnum, &pythia.info, &pythia.settings); }
  bool fill_next_event( Pythia8::Pythia& pythia, GenEvent& evt) {
    return fill_next_event( pythia, &evt); }

  // Alternative method to convert Pythia events into HepMC3 ones.
  bool fill_next_event( Pythia8::Event& pyev, GenEvent&evt, int ievnum = -1,
    const Pythia8::Info* pyinfo = 0, Pythia8::Settings* pyset = 0) {
    return fill_next_event(pyev, &evt, ievnum, pyinfo, pyset); }
  bool fill_next_event( Pythia8::Event& pyev, GenEvent* evt, int ievnum = -1,
    const Pythia8::Info* pyinfo = 0, Pythia8::Settings* pyset = 0) {

    // 1. Error if no event passed.
    if (evt == nullptr) return warning(pyinfo,
      "Pythia8ToHepMC::fill_next_event", "passed null event");

    // Event number counter.
    if ( ievnum >= 0 ) {
      evt->set_event_number(ievnum);
      m_internal_event_number = ievnum;
    }
    else {
      evt->set_event_number(m_internal_event_number);
      ++m_internal_event_number;
    }

    // Set units to be GeV and mm, to agree with Pythia ones.
    evt->set_units(Units::GEV,Units::MM);

    // 1a. If there is a HIInfo object fill info from that.
    if ( pyinfo && pyinfo->hiInfo ) {
      auto ion = make_shared<HepMC3::GenHeavyIon>();
      ion->Ncoll_hard = pyinfo->hiInfo->nCollND();
      ion->Ncoll = pyinfo->hiInfo->nCollTot();
      ion->Npart_proj = pyinfo->hiInfo->nAbsProj() +
                        pyinfo->hiInfo->nDiffProj();
      ion->Npart_targ = pyinfo->hiInfo->nAbsTarg() +
                        pyinfo->hiInfo->nDiffTarg();
      ion->impact_parameter = pyinfo->hiInfo->b();
      evt->set_heavy_ion(ion);
    }

    // 2. Fill particle information.
    vector<GenParticlePtr> hepevt_particles;
    hepevt_particles.reserve( pyev.size() );
    for(int i = 0; i < pyev.size(); ++i) {
      hepevt_particles.push_back( make_shared<GenParticle>(
        FourVector( pyev[i].px(), pyev[i].py(), pyev[i].pz(), pyev[i].e() ),
        pyev[i].id(), pyev[i].statusHepMC() ) );
      hepevt_particles[i]->set_generated_mass( pyev[i].m() );
    }

    // 3. Fill vertex information.
    vector<GenVertexPtr> vertex_cache;
    vector<GenParticlePtr> beam_particles;
    for (int i = 1; i < pyev.size(); ++i) {
      vector<int> mothers = pyev[i].motherList();
      sort(mothers.begin(),mothers.end());
      for (;;) {
        if (!mothers.empty() && mothers.front() == 0)
          mothers.erase(mothers.begin());
        else break;
      }
      if (mothers.size()) {
        GenVertexPtr prod_vtx = hepevt_particles[mothers[0]]->end_vertex();
        if (!prod_vtx) {
          prod_vtx = make_shared<GenVertex>();
          vertex_cache.push_back(prod_vtx);
          for (unsigned int j = 0; j < mothers.size(); ++j)
            prod_vtx->add_particle_in( hepevt_particles[mothers[j]] );
        }
        FourVector prod_pos( pyev[i].xProd(), pyev[i].yProd(),pyev[i].zProd(),
          pyev[i].tProd() );

        // Update vertex position if necessary.
        if (!prod_pos.is_zero() && prod_vtx->position().is_zero())
          prod_vtx->set_position( prod_pos );
        prod_vtx->add_particle_out( hepevt_particles[i] );
      } else beam_particles.push_back(hepevt_particles[i]);
    }

    // Reserve memory for the event.
    evt->reserve( hepevt_particles.size(), vertex_cache.size() );

    // Add particles and vertices in topological order.
    evt->add_tree( beam_particles );

    // Attributes should be set after adding the particles to event.
    for (int i = 0; i < pyev.size(); ++i) {
      /* TODO: Set polarization */
      // Colour flow uses index 1 and 2.
      int colType = pyev[i].colType();
      if (colType ==  -1 ||colType ==  1 || colType == 2) {
        int flow1 = 0, flow2 = 0;
        if (colType ==  1 || colType == 2) flow1 = pyev[i].col();
        if (colType == -1 || colType == 2) flow2 = pyev[i].acol();
        hepevt_particles[i]->add_attribute("flow1",
          make_shared<IntAttribute>(flow1));
        hepevt_particles[i]->add_attribute("flow2",
          make_shared<IntAttribute>(flow2));
      }
    }

    // If hadronization switched on then no final coloured particles.
    bool doHadr = (pyset == 0) ? m_free_parton_warnings
      : pyset->flag("HadronLevel:all") && pyset->flag("HadronLevel:Hadronize");

    // 4. Check for particles which come from nowhere, i.e. are without
    // mothers or daughters. These need to be attached to a vertex, or else
    // they will never become part of the event.
    for (int i = 1; i < pyev.size(); ++i) {

      // Check for particles not added to the event.
      // NOTE: We have to check if this step makes any sense in
      // the HepMC event standard.
      if ( hepevt_particles[i] == nullptr ||
        !hepevt_particles[i]->in_event()) {
        warning(pyinfo, "Pythia8ToHepMC::fill_next_event",
          "found orphan particle", "i = " + Pythia8::toString(i));
        GenVertexPtr prod_vtx = make_shared<GenVertex>();
        prod_vtx->add_particle_out( hepevt_particles[i] );
        evt->add_vertex(prod_vtx);
      }

      // Also check for free partons (= gluons and quarks; not diquarks?).
      if ( doHadr && m_free_parton_warnings ) {
        if ( hepevt_particles[i]->pid() == 21
           && hepevt_particles[i]->end_vertex() == nullptr ) {
          warning(pyinfo, "Pythia8ToHepMC::fill_next_event",
            "found gluon without end vertex", "i = " + Pythia8::toString(i));
          if ( m_crash_on_problem ) exit(1);
        }
        if ( abs(hepevt_particles[i]->pid()) <= 6
          && hepevt_particles[i]->end_vertex() == nullptr ) {
          warning(pyinfo, "Pythia8ToHepMC::fill_next_event",
            "found quark without end vertex", "i = " + Pythia8::toString(i));
          if ( m_crash_on_problem ) exit(1);
        }
      }
    }

    // 5. Store PDF, weight, cross section and other event information.
    // Flavours of incoming partons.
    if (m_store_pdf && pyinfo != 0) {
      int id1pdf = pyinfo->id1pdf();
      int id2pdf = pyinfo->id2pdf();
      if ( m_convert_gluon_to_0 ) {
        if (id1pdf == 21) id1pdf = 0;
        if (id2pdf == 21) id2pdf = 0;
      }

      // Store PDF information.
      GenPdfInfoPtr pdfinfo = make_shared<GenPdfInfo>();
      pdfinfo->set(id1pdf, id2pdf, pyinfo->x1pdf(), pyinfo->x2pdf(),
        pyinfo->QFac(), pyinfo->pdf1(), pyinfo->pdf2() );
      evt->set_pdf_info( pdfinfo );
    }

    // Store process code, scale, alpha_em, alpha_s.
    if (m_store_proc && pyinfo != 0) {
      evt->add_attribute("signal_process_id",
        make_shared<IntAttribute>( pyinfo->code()));
      evt->add_attribute("mpi",
        make_shared<IntAttribute>( pyinfo->nMPI()));
      evt->add_attribute("event_scale",
        make_shared<DoubleAttribute>(pyinfo->QRen()));
      evt->add_attribute("alphaQCD",
        make_shared<DoubleAttribute>(pyinfo->alphaS()));
      evt->add_attribute("alphaQED",
        make_shared<DoubleAttribute>(pyinfo->alphaEM()));
    }

    // Store event weights.
    if (m_store_weights && pyinfo != 0) {
      evt->weights().clear();
      for (int iWeight = 0; iWeight < pyinfo->numberOfWeights(); ++iWeight)
        evt->weights().push_back(pyinfo->weightValueByIndex(iWeight));
    }

    // Store cross-section information in pb.
    if (m_store_xsec && pyinfo != 0) {
      // First set attribute to event, such that
      // GenCrossSection::set_cross_section knows how many weights the
      // event has and sets the number of cross sections accordingly.
      GenCrossSectionPtr xsec = make_shared<GenCrossSection>();
      evt->set_cross_section(xsec);
      xsec->set_cross_section( pyinfo->sigmaGen() * 1e9,
        pyinfo->sigmaErr() * 1e9);
      // If multiweights with possibly different xsec, overwrite central value
      vector<double> xsecVec = pyinfo->weightContainerPtr->getTotalXsec();
      if (xsecVec.size() > 0) {
        for (unsigned int iXsec = 0; iXsec < xsecVec.size(); ++iXsec) {
          xsec->set_xsec(iXsec, xsecVec[iXsec]*1e9);
        }
      }
    }

    // Done.
    return true;
  }

  // Read out values for some switches.
  bool print_inconsistency()  const { return m_print_inconsistency; }
  bool free_parton_warnings() const { return m_free_parton_warnings; }
  bool crash_on_problem()     const { return m_crash_on_problem; }
  bool convert_gluon_to_0()   const { return m_convert_gluon_to_0; }
  bool store_pdf()            const { return m_store_pdf; }
  bool store_proc()           const { return m_store_proc; }
  bool store_xsec()           const { return m_store_xsec; }
  bool store_weights()        const { return m_store_weights; }

  // Set values for some switches.
  void set_print_inconsistency(bool b = true)  { m_print_inconsistency  = b; }
  void set_free_parton_warnings(bool b = true) { m_free_parton_warnings = b; }
  void set_crash_on_problem(bool b = false)    { m_crash_on_problem     = b; }
  void set_convert_gluon_to_0(bool b = false)  { m_convert_gluon_to_0   = b; }
  void set_store_pdf(bool b = true)            { m_store_pdf            = b; }
  void set_store_proc(bool b = true)           { m_store_proc           = b; }
  void set_store_xsec(bool b = true)           { m_store_xsec           = b; }
  void set_store_weights(bool b = true)        { m_store_weights        = b; }

private:

    // Try to send warning message to the logger if present, otherwise
  // send it to cout if print_inconsistency().
  bool warning(const Pythia8::Info * pyinfo, string loc,
               string message, string extraInfo = "") {
    if ( pyinfo )
      pyinfo->loggerPtr->warningMsg(loc, message, extraInfo);
    else if ( print_inconsistency() )
      cout << "Warning in " << loc << ": " << message << extraInfo << endl;
    return false;
  }

  // Data members.
  int  m_internal_event_number;
  bool m_print_inconsistency, m_free_parton_warnings, m_crash_on_problem,
       m_convert_gluon_to_0, m_store_pdf, m_store_proc, m_store_xsec,
       m_store_weights;

};

//==========================================================================

// This class converts the Pythia event records into the HepMC3 event
// format, with the possibility of filtering the event.
//
// Traversing the Pythia event record in the "forward" direction by
// following daughters does not necsssarily produce the same graph
// structure as traversing in the "backward" direction by following
// mothers. From the Pythia HTML manual: "the mother-daughter relation
// normally is reciprocal, but not always. An example is hadron beams
// (indices 1 and 2), where each beam remnant and the initiator of
// each multiparton interaction has the respective beam as mother, but
// the beam itself only has the initiator of the hardest interaction
// as daughter."
//
// This class takes the "backward" approach.

class Pythia8ToFilteredHepMC3 {

public:

  // Constructor.
  Pythia8ToFilteredHepMC3(): m_internal_event_number(0), m_vertex_status(0),
    m_store_color(true), m_store_pdf(true), m_store_proc(true),
    m_store_xsec(true), m_store_weights(true), m_selector(nullptr) {}

  // The recommended method to convert Pythia events into HepMC3 ones.
  bool fill_next_event( Pythia8::Pythia& pythia, GenEvent* evt,
    int ievnum = -1 ) { return fill_next_event( pythia.event, evt,
    ievnum, &pythia.info); }
  bool fill_next_event( Pythia8::Pythia& pythia, GenEvent& evt) {
    return fill_next_event( pythia, &evt); }

  // Alternative method to convert Pythia events into HepMC3 ones.
  bool fill_next_event( Pythia8::Event& pyev, GenEvent&evt, int ievnum = -1,
    const Pythia8::Info* pyinfo = 0) {
    return fill_next_event(pyev, &evt, ievnum, pyinfo); }
  bool fill_next_event( Pythia8::Event& pyev, GenEvent* evt, int ievnum = -1,
    const Pythia8::Info* pyinfo = 0) {

    // Error if no event passed.
    if (evt == nullptr) return warning(pyinfo,
      "Pythia8ToFilteredHepMC::fill_next_event", "passed null event");

    // Event number counter.
    if ( ievnum >= 0 ) {
      evt->set_event_number(ievnum);
      m_internal_event_number = ievnum;
    }
    else {
      evt->set_event_number(m_internal_event_number);
      ++m_internal_event_number;
    }

    // Set units to be GeV and mm, to agree with Pythia ones.
    evt->set_units(Units::GEV, Units::MM);

    // If there is a HIInfo object fill info from that.
    if ( pyinfo && pyinfo->hiInfo ) {
      auto ion = make_shared<HepMC3::GenHeavyIon>();
      ion->Ncoll_hard = pyinfo->hiInfo->nCollND();
      ion->Ncoll = pyinfo->hiInfo->nCollTot();
      ion->Npart_proj = pyinfo->hiInfo->nAbsProj() +
                        pyinfo->hiInfo->nDiffProj();
      ion->Npart_targ = pyinfo->hiInfo->nAbsTarg() +
                        pyinfo->hiInfo->nDiffTarg();
      ion->impact_parameter = pyinfo->hiInfo->b();
      evt->set_heavy_ion(ion);
    }

    // Build all the vertices from the particle record.
    // An ordered map is used for reproducibility of the record.
    map<pair<int, int>, GenVertexPtr> vrts;
    vector<GenParticlePtr> prts(pyev.size(), nullptr);
    for (int iPrt = pyev.size() - 1; iPrt > 0; --iPrt)
      add(pyev, iPrt, vrts, prts);

    // Reconnect intermediate particles without an end vertex.
    for (int iPrt = 0; iPrt < pyev.size(); ++iPrt) {
      if (pyev[iPrt].isFinal()) continue;
      GenParticlePtr prt = prts[iPrt];
      if (prt == nullptr) continue;
      if (prt->end_vertex() != nullptr) continue;
      // Find the first daughter that has a production vertex.
      for (const int& iDtr : pyev[iPrt].daughterListRecursive()) {
        GenParticlePtr dtr = prts[iDtr];
        if (dtr != nullptr && dtr->production_vertex() != nullptr) {
          dtr->production_vertex()->add_particle_in(prt);
        }
      }
    }

    // Reserve memory for the event.
    evt->reserve(prts.size(), vrts.size());

    // Add the vertices.
    for (auto& vrt : vrts)
      evt->add_vertex(vrt.second);

    // Store color flow.
    // This has to be performed after the particles have been added to
    // the event, otherwise the attribute will not be persisted.
    if (m_store_color) {
      for (int iPrt = 0; iPrt < pyev.size(); ++iPrt) {
        if (prts[iPrt] == nullptr) continue;
        int colType = pyev[iPrt].colType();
        if (colType ==  -1 || colType ==  1 || colType == 2) {
          int flow1 = 0, flow2 = 0;
          if (colType ==  1 || colType == 2) flow1 = pyev[iPrt].col();
          if (colType == -1 || colType == 2) flow2 = pyev[iPrt].acol();
          prts[iPrt]->add_attribute("flow1", make_shared<IntAttribute>(flow1));
          prts[iPrt]->add_attribute("flow2", make_shared<IntAttribute>(flow2));
        }
      }
    }

    // Store PDF, weight, cross section and other event information.
    // Flavours of incoming partons.
    if (m_store_pdf && pyinfo != 0) {
      int id1pdf = pyinfo->id1pdf();
      int id2pdf = pyinfo->id2pdf();

      // Store PDF information.
      GenPdfInfoPtr pdfinfo = make_shared<GenPdfInfo>();
      pdfinfo->set(id1pdf, id2pdf, pyinfo->x1pdf(), pyinfo->x2pdf(),
        pyinfo->QFac(), pyinfo->pdf1(), pyinfo->pdf2() );
      evt->set_pdf_info( pdfinfo );
    }

    // Store process code, scale, alpha_em, alpha_s.
    if (m_store_proc && pyinfo != 0) {
      evt->add_attribute("signal_process_id",
        make_shared<IntAttribute>( pyinfo->code()));
      evt->add_attribute("mpi",
        make_shared<IntAttribute>( pyinfo->nMPI()));
      evt->add_attribute("event_scale",
        make_shared<DoubleAttribute>(pyinfo->QRen()));
      evt->add_attribute("alphaQCD",
        make_shared<DoubleAttribute>(pyinfo->alphaS()));
      evt->add_attribute("alphaQED",
        make_shared<DoubleAttribute>(pyinfo->alphaEM()));
    }

    // Store event weights.
    if (m_store_weights && pyinfo != 0) {
      evt->weights().clear();
      for (int iWeight = 0; iWeight < pyinfo->numberOfWeights(); ++iWeight)
        evt->weights().push_back(pyinfo->weightValueByIndex(iWeight));
    }

    // Store cross-section information in pb.
    if (m_store_xsec && pyinfo != 0) {
      // First set attribute to event, such that
      // GenCrossSection::set_cross_section knows how many weights the
      // event has and sets the number of cross sections accordingly.
      GenCrossSectionPtr xsec = make_shared<GenCrossSection>();
      evt->set_cross_section(xsec);
      xsec->set_cross_section( pyinfo->sigmaGen() * 1e9,
        pyinfo->sigmaErr() * 1e9);
      // If multiweights with possibly different xsec, overwrite central value
      vector<double> xsecVec = pyinfo->weightContainerPtr->getTotalXsec();
      if (xsecVec.size() > 0) {
        for (unsigned int iXsec = 0; iXsec < xsecVec.size(); ++iXsec) {
          xsec->set_xsec(iXsec, xsecVec[iXsec]*1e9);
        }
      }
    }

    // Done.
    return true;
  }

  // Read out values for some switches.
  bool store_color()          const { return m_store_color; }
  bool store_pdf()            const { return m_store_pdf; }
  bool store_proc()           const { return m_store_proc; }
  bool store_xsec()           const { return m_store_xsec; }
  bool store_weights()        const { return m_store_weights; }

  // Set values for some switches.
  void set_store_color(bool b = true)          { m_store_color   = b; }
  void set_store_pdf(bool b = true)            { m_store_pdf     = b; }
  void set_store_proc(bool b = true)           { m_store_proc    = b; }
  void set_store_xsec(bool b = true)           { m_store_xsec    = b; }
  void set_store_weights(bool b = true)        { m_store_weights = b; }

  // Get and set the particle selector.
  const function<bool(const Pythia8::Event&, size_t)>& selector() const {
    return m_selector;}
  void set_selector(function<bool(const Pythia8::Event&, size_t)>
    selector) {m_selector = std::move(selector);}

  // Get and set the default vertex status.
  int vertex_status() const { return m_vertex_status; }
  void set_vertex_status(int i = 0) { m_vertex_status = i; }

private:

  // Returns a HepMC particle for a given Pythia particle index.
  GenParticlePtr particle(const Pythia8::Event& pyev, size_t idx,
    vector<GenParticlePtr>& prts) {
    GenParticlePtr prt = prts[idx];
    if (prt == nullptr) {
      prt = make_shared<GenParticle>(
        FourVector(pyev[idx].px(), pyev[idx].py(), pyev[idx].pz(),
          pyev[idx].e()), pyev[idx].id(), pyev[idx].statusHepMC());
      prt->set_generated_mass(pyev[idx].m());
      prts[idx] = prt;
    }
    return prt;
  }

  // Add the vertex for a particle.
  GenVertexPtr add(const Pythia8::Event& pyev, size_t iPrt,
    map<pair<int, int>, GenVertexPtr>& vrts,
    vector<GenParticlePtr>& prts, bool check = true) {
    // Create the particle and production vertex.
    GenVertexPtr vrt = nullptr;
    GenParticlePtr prt = nullptr;
    // Filter and create the particle if it does not already exist.
    if (check) {
      if (m_selector == nullptr || m_selector(pyev, iPrt))
        prt = particle(pyev, iPrt, prts);
      else
        return vrt;
    }
    // Check if the vertex exists.
    pair<int, int> key = make_pair(pyev[iPrt].mother1(), pyev[iPrt].mother2());
    auto result = vrts.emplace(key, nullptr);
    auto itr = result.first;
    bool added = result.second;
    // Set the vertex if it exists.
    if (!added) vrt = itr->second;
    // Set the vertex if it does not exist.
    else {
      // Check if there are linking mothers.
      vector<int> iMoms = pyev[iPrt].motherList();
      bool link = false;
      for (int& iMom : iMoms) {
        if (m_selector == nullptr || m_selector(pyev, iMom)) link = true;
        else iMom = -iMom;
      }
      // Create the vertex if linking mothers exist.
      if (link) {
        vrt = make_shared<GenVertex>(
          FourVector(pyev[iPrt].xProd(), pyev[iPrt].yProd(),
            pyev[iPrt].zProd(), pyev[iPrt].tProd()));
        vrt->set_status(m_vertex_status);
        for (const int& iMom : iMoms) {
          if (iMom > 0) vrt->add_particle_in(particle(pyev, iMom, prts));
        }
      // Get the upstream vertex if no linking mothers exist.
      } else {
        for (const int& iMom : iMoms) {
          if (iMom != 0) {
            vrt = add(pyev, abs(iMom), vrts, prts, check = false);
            break;
          }
        }
      }
      // Set the vertex.
      itr->second = vrt;
    }
    // Add the outgoing particle.
    if (prt != nullptr && vrt != nullptr) vrt->add_particle_out(prt);
    // Return the vertex.
    return vrt;
  }

  // Try to send warning message to the logger if present, otherwise
  // send it to cout.
  bool warning(const Pythia8::Info * pyinfo, string loc,
               string message, string extraInfo = "") {
    if (pyinfo)
      pyinfo->loggerPtr->warningMsg(loc, message, extraInfo);
    else
      cout << "Warning in " << loc << ": " << message << extraInfo << endl;
    return false;
  }

  // Data members.
  int  m_internal_event_number, m_vertex_status;
  bool m_store_color, m_store_pdf, m_store_proc, m_store_xsec, m_store_weights;
  function<bool(const Pythia8::Event&, size_t)> m_selector;

};

//==========================================================================

} // end namespace HepMC3

namespace Pythia8 {

//==========================================================================

// This is a wrapper around HepMC3::Pythia8ToHepMC in the Pythia8
// namespace that simplifies the most common use cases. It stores the
// current GenEvent and output stream internally to avoid cluttering
// of user code. This class is also defined in HepMC2.h with the same
// signatures, and the user can therefore switch between HepMC version
// 2 and 3, by simply changing the include file.

class Pythia8ToHepMC : public HepMC3::Pythia8ToHepMC3 {

public:

  // We can either have standard ascii output version 2 or three or
  // none at all.
  enum OutputType { none, ascii2, ascii3 };

  // Typedef for the version 3 specific classes used.
  typedef HepMC3::GenEvent GenEvent;
  typedef shared_ptr<GenEvent> EventPtr;
  typedef HepMC3::Writer Writer;
  typedef shared_ptr<Writer> WriterPtr;

  // The empty constructor does not creat an aoutput stream.
  Pythia8ToHepMC() : runinfo(make_shared<HepMC3::GenRunInfo>()) {}

  // Construct an object with an internal output stream.
  Pythia8ToHepMC(string filename, OutputType ft = ascii3)
    : runinfo(make_shared<HepMC3::GenRunInfo>()) {
    setNewFile(filename, ft);
  }

  // Open a new external output stream.
  bool setNewFile(string filename, OutputType ft = ascii3) {
    switch ( ft ) {
    case ascii3:
      writerPtr = make_shared<HepMC3::WriterAscii>(filename);
      break;
    case ascii2:
      writerPtr = make_shared<HepMC3::WriterAsciiHepMC2>(filename);
      break;
    case none:
      break;
    }
    return writerPtr != nullptr;
  }

  // Create a new GenEvent object and fill it with information from
  // the given Pythia object.
  bool fillNextEvent(Pythia & pythia) {
    geneve = make_shared<HepMC3::GenEvent>(runinfo);
    if (runinfo->weight_names().size() == 0)
      setWeightNames(pythia.info.weightNameVector());
    return fill_next_event(pythia, *geneve);
  }

  // Write out the current GenEvent to the internal stream.
  void writeEvent() {
    writerPtr->write_event(*geneve);
  }

  // Create a new GenEvent object and fill it with information from
  // the given Pythia object and write it out directly to the
  // internal stream.
  bool writeNextEvent(Pythia & pythia) {
    if ( !fillNextEvent(pythia) ) return false;
    writeEvent();
    return !writerPtr->failed();
  }

  // Get a reference to the current GenEvent.
  GenEvent & event() {
    return *geneve;
  }

  // Get a pointer to the current GenEvent.
   EventPtr getEventPtr() {
    return geneve;
  }

  // Get a reference to the internal stream.
  Writer & output() {
    return *writerPtr;
  }

  // Get a pointer to the internal stream.
  WriterPtr outputPtr() {
    return writerPtr;
  }

  // Set cross section information in the current GenEvent.
  void setXSec(double xsec, double xsecerr) {
    auto xsecptr = geneve->cross_section();
    if ( !xsecptr ) {
      xsecptr = make_shared<HepMC3::GenCrossSection>();
      geneve->set_cross_section(xsecptr);
    }
    xsecptr->set_cross_section(xsec, xsecerr);
  }

  // Update all weights in the current GenEvent.
  void setWeights(const vector<double> & wv) {
    geneve->weights() = wv;
  }

  // Set all weight names in the current run.
  void setWeightNames(const vector<string> &wnv) {
    runinfo->set_weight_names(wnv);
  }

  // Update the PDF information in the current GenEvent
  void setPdfInfo(int id1, int id2, double x1, double x2,
                  double scale, double xf1, double xf2,
                  int pdf1 = 0, int pdf2 = 0) {
    auto pdf = make_shared<HepMC3::GenPdfInfo>();
    pdf->set(id1, id2, x1, x2, scale, xf1, xf2, pdf1, pdf2);
    geneve->set_pdf_info(pdf);
  }

  // Add an additional attribute derived from HepMC3::Attribute
  // to the current event.
  template<class T>
  void addAttribute(const string& name, T& attribute) {
    shared_ptr<HepMC3::Attribute> att = make_shared<T>(attribute);
    geneve->add_attribute(name, att);
  }

  // Add an attribute of double type.
  template<class T=double>
  void addAttribute(const string& name, double& attribute) {
    auto dAtt = HepMC3::DoubleAttribute(attribute);
    shared_ptr<HepMC3::Attribute> att =
      make_shared<HepMC3::DoubleAttribute>(dAtt);
    geneve->add_attribute(name, att);
  }

  // Add an attribute of integer type.
  template<class T=int>
  void addAttribute(const string& name, int& attribute) {
    auto iAtt = HepMC3::IntAttribute(attribute);
    shared_ptr<HepMC3::Attribute> att =
      make_shared<HepMC3::IntAttribute>(iAtt);
    geneve->add_attribute(name, att);
  }

  // Remove an attribute from the current event.
  void removeAttribute(const string& name) {
    geneve->remove_attribute(name);
  }

private:

  // The current GenEvent
  EventPtr geneve = nullptr;

  // The output stream.
  WriterPtr writerPtr = nullptr;

  // The current run info.
  shared_ptr<HepMC3::GenRunInfo> runinfo;

};

}

#endif // end Pythia8_HepMC3_H
