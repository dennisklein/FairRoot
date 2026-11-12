// fv0sim: single-process O2 simulation of the FV0 detector.
//
// This program is a reduced copy of the O2 program o2-sim-serial
// (run/o2sim.cxx, macro/o2sim.C and macro/build_geometry.C of AliceO2).
// It keeps the simulation path of o2-sim-serial unchanged and removes only
// the parts that do not affect the transport:
//  - CCDB access, GRP objects, alignment and the hit file migration
//    (the magnetic field comes from the command line, as with --field <kG>),
//  - all modules other than CAVE, PIPE and FV0,
//  - all event generators other than the box generator.
//
// The command line options are the o2-sim options, for example:
//   fv0sim -m CAVE PIPE FV0 -e TGeant4 -g boxgen -n 100 --seed 1 \
//          --configKeyValues "BoxGun.pdg=211;BoxGun.eta[0]=2.2;BoxGun.eta[1]=5.1"

#include <CommonUtils/ConfigurableParam.h>
#include <CommonUtils/RngHelper.h>
#include <DetectorsPassive/Cave.h>
#include <DetectorsPassive/Pipe.h>
#include <FV0Simulation/Detector.h>
#include <Field/MagneticField.h>
#include <Generators/BoxGenerator.h>
#include <Generators/BoxGunParam.h>
#include <Generators/PrimaryGenerator.h>
#include <SimConfig/SimConfig.h>
#include <SimSetup/SimSetup.h>
#include <SimSetup/StaticTLSReserve.h>

#ifdef FLATSIM_STATIC
// T2b: the Geant4 setup is linked into the program, so it is called directly,
// not loaded from libO2G4Setup.so with dlopen.
namespace o2
{
void SetCuts();
namespace g4config
{
void G4Config();
void G4Terminate();
} // namespace g4config
} // namespace o2
#endif
#include <SimulationDataFormat/MCEventHeader.h>
#include <Steer/O2RunSim.h>

#include <FairRootFileSink.h>
#include <FairSystemInfo.h>
#include <fairlogger/Logger.h>

#include <TGeoManager.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include <algorithm>
#include <vector>
#include <string>

namespace
{

bool isActivated(std::string const& name)
{
  auto& modules = o2::conf::SimConfig::Instance().getActiveModules();
  return std::find(modules.begin(), modules.end(), name) != modules.end();
}

bool isReadout(std::string const& name)
{
  auto& readout = o2::conf::SimConfig::Instance().getReadoutDetectors();
  return std::find(readout.begin(), readout.end(), name) != readout.end();
}

// From macro/build_geometry.C, reduced to CAVE, PIPE and FV0.
void buildGeometry(FairRunSim* run)
{
  TString dir = getenv("VMCWORKDIR");
  gSystem->Setenv("GEOMPATH", (dir + "/Detectors/Geometry/").Data());
  gSystem->Setenv("CONFIG_DIR", (dir + "/Detectors/gconfig").Data());

  run->SetMaterials("media.geo");
  // From Detectors/Base/src/SimFieldUtils.cxx, without the CCDB branch.
  auto& conf = o2::conf::SimConfig::Instance().getConfigData();
  run->SetField(o2::field::MagneticField::createNominalField(
    conf.mField, conf.mFieldMode == o2::conf::SimFieldMode::kUniform));

  auto cave = new o2::passive::Cave("CAVE");
  cave->includeZDC(false);
  cave->SetGeometryFileName("cave.geo");
  run->AddModule(cave);

  if (isActivated("PIPE")) {
    run->AddModule(new o2::passive::Pipe("PIPE", "Beam pipe"));
  }
  // Dense index of the active detectors, for the hit bits of the MC tracks (as in build_geometry.C).
  std::vector<int> detId2RunningId(o2::detectors::DetID::nDetectors, -1);
  if (isActivated("FV0")) {
    auto* fv0 = new o2::fv0::Detector(isReadout("FV0"));
    run->AddModule(fv0);
    if (fv0->IsActive()) {
      detId2RunningId[fv0->GetDetId()] = 0;
    }
  }
  o2::base::Detector::setDetId2HitBitIndex(detId2RunningId);
}

} // namespace

int main(int argc, char* argv[])
{
#ifndef FLATSIM_STATIC
  // Geant4 is dlopened with initial-exec TLS; glibc takes the static TLS reserve only at program start.
  o2::ensureStaticTLSReserve(argv);
#endif
  auto& conf = o2::conf::SimConfig::Instance();
  if (!conf.resetFromArguments(argc, argv)) {
    return 1;
  }
  FairLogger::GetLogger()->SetLogScreenLevel(conf.getLogSeverity().c_str());
  FairLogger::GetLogger()->SetLogVerbosityLevel(conf.getLogVerbosity().c_str());

  TStopwatch timer;
  timer.Start();

  unsigned int nEvents = conf.getNEvents();
  o2::eventgen::Generator::setTotalNEvents(nEvents);
  o2::conf::ConfigurableParam::updateFromFile(conf.getConfigFile());
  o2::conf::ConfigurableParam::updateFromString(conf.getKeyValueString());
  auto seed = o2::utils::RngHelper::setGRandomSeed(conf.getStartSeed());
  LOG(info) << "RNG INITIAL SEED " << seed;

  auto run = new o2::steer::O2RunSim(false, false);
  run->SetImportTGeoToVMC(false);
#ifdef FLATSIM_STATIC
  run->SetSimSetup([]() {
    o2::g4config::G4Config();
    o2::SetCuts();
  });
#else
  run->SetSimSetup([&conf]() { o2::SimSetup::setup(conf.getMCEngine().c_str()); });
#endif
  run->SetRunId(conf.getTimestamp());
  run->SetSink(new FairRootFileSink((conf.getOutPrefix() + ".root").c_str()));
  run->SetName(conf.getMCEngine().c_str());
  run->SetIsMT(conf.getIsMT());
  run->SetMCEventHeader(new o2::dataformats::MCEventHeader());

  auto lock = TGeoManager::LockDefaultUnits(false);
  TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
  TGeoManager::LockDefaultUnits(lock);

  buildGeometry(run);

  // From Generators/src/GeneratorFactory.cxx, generator "boxgen".
  auto primGen = new o2::eventgen::PrimaryGenerator();
  auto& box = o2::eventgen::BoxGunParam::Instance();
  primGen->AddGenerator(new o2::eventgen::BoxGenerator(box.pdg, box.number, box.eta[0], box.eta[1],
                                                       box.prange[0], box.prange[1],
                                                       box.phirange[0], box.phirange[1]));
  run->SetGenerator(primGen);

  run->Init();
  timer.Stop();
  FairSystemInfo sysinfo;
  LOG(info) << "Init: Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s";
  LOG(info) << "Init: Memory used " << sysinfo.GetMaxMemory() << " MB";

  timer.Start();
  run->Run(conf.getNEvents());
  timer.Stop();
  LOG(info) << "Run: Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s";
  LOG(info) << "Run: Memory used " << sysinfo.GetMaxMemory() << " MB";

#ifdef FLATSIM_STATIC
  o2::g4config::G4Terminate();
#else
  o2::SimSetup::shutdown();
#endif
  delete run;
  return 0;
}
