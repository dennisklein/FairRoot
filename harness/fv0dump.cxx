// fv0dump: write the physics observables of an fv0sim output file as CSV.
//
// Usage: fv0dump <o2sim.root> <prefix>
// Output:
//   <prefix>.hits.csv    one line per FV0 hit
//   <prefix>.events.csv  one line per event (sums and counts)
//
// The harness compares these files between a reference and a candidate
// build. It reads the files only through ROOT I/O and the O2 dictionaries,
// so it does not depend on the simulation code it tests.

#include <DataFormatsFV0/Hit.h>
#include <SimulationDataFormat/MCTrack.h>
#include <SimulationDataFormat/TrackReference.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::fprintf(stderr, "usage: fv0dump <o2sim.root> <prefix>\n");
    return 2;
  }
  std::unique_ptr<TFile> file(TFile::Open(argv[1]));
  if (!file || file->IsZombie()) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 1;
  }
  // o2-sim names the tree "o2sim"; the serial FairRoot sink keeps the FairRoot default "cbmsim".
  auto tree = file->Get<TTree>("o2sim");
  if (!tree) {
    tree = file->Get<TTree>("cbmsim");
  }
  if (!tree) {
    std::fprintf(stderr, "no tree o2sim or cbmsim in %s\n", argv[1]);
    return 1;
  }
  std::vector<o2::fv0::Hit>* hits = nullptr;
  std::vector<o2::MCTrack>* tracks = nullptr;
  std::vector<o2::TrackReference>* refs = nullptr;
  tree->SetBranchAddress("FV0Hit", &hits);
  tree->SetBranchAddress("MCTrack", &tracks);
  if (tree->GetBranch("TrackRefs")) {
    tree->SetBranchAddress("TrackRefs", &refs);
  }

  std::string prefix = argv[2];
  auto fh = std::fopen((prefix + ".hits.csv").c_str(), "w");
  auto fe = std::fopen((prefix + ".events.csv").c_str(), "w");
  // dphi, dr: hit entry point relative to the straight line of its track from the track vertex.
  // The field bends the line (dphi), and the material scatters it (dphi, dr).
  std::fprintf(fh, "event,cell,pdg,x,y,z,xs,ys,zs,t,edep,e,primary,dphi,dr\n");
  std::fprintf(fe, "event,nhits,edep,ntracks,nrefs\n");
  for (Long64_t ev = 0; ev < tree->GetEntries(); ++ev) {
    tree->GetEntry(ev);
    double edep = 0;
    for (auto const& h : *hits) {
      edep += h.GetEnergyLoss();
      int primary = -1;
      double dphi = 0, dr = 0;
      auto id = h.GetTrackID();
      if (id >= 0 && id < static_cast<int>(tracks->size())) {
        auto const& t = (*tracks)[id];
        primary = t.isPrimary() ? 1 : 0;
        double pz = t.GetStartVertexMomentumZ();
        if (std::abs(pz) > 1e-9) {
          double s = (h.GetStartZ() - t.GetStartVertexCoordinatesZ()) / pz;
          double lx = t.GetStartVertexCoordinatesX() + s * t.GetStartVertexMomentumX();
          double ly = t.GetStartVertexCoordinatesY() + s * t.GetStartVertexMomentumY();
          dphi = std::remainder(std::atan2(h.GetStartY(), h.GetStartX()) - std::atan2(ly, lx), 2 * M_PI);
          dr = std::hypot(h.GetStartX(), h.GetStartY()) - std::hypot(lx, ly);
        }
      }
      std::fprintf(fh, "%lld,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%.9g,%.9g\n", ev,
                   h.GetDetectorID(), h.GetParticlePdg(), h.GetX(), h.GetY(), h.GetZ(),
                   h.GetStartX(), h.GetStartY(), h.GetStartZ(), h.GetTime(),
                   h.GetEnergyLoss(), h.GetE(), primary, dphi, dr);
    }
    std::fprintf(fe, "%lld,%zu,%.9g,%zu,%zu\n", ev, hits->size(), edep, tracks->size(),
                 refs ? refs->size() : 0);
  }
  std::fclose(fh);
  std::fclose(fe);
  return 0;
}
