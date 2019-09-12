/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

enum class pdgCode : Int_t {
  eplus = 11,
  pion = 211,
  proton = 2212
};

int runMC(Int_t nEvents = 1000, TString mcEngine = "TGeant3", Bool_t isMT = true)
{
  // -----   Parameters   ---------------------------------------------------
  TString dir = getenv("VMCWORKDIR");

  TString tutGeomDir = dir + "/common/geometry";
  gSystem->Setenv("GEOMPATH", tutGeomDir.Data());

  TString tutConfigDir = dir + "/common/gconfig";
  gSystem->Setenv("CONFIG_DIR",tutConfigDir.Data());

  Double_t momentum = 2.;

  Double_t theta    = 10.;

  TString outDir = "./";

  TString outFile = Form("%s/geane.mc.root", outDir.Data()); // Output file name

  TString parFile = Form("%s/geane.par.root", outDir.Data()); // Parameter file name

  TString geoFile = "geofile_" + mcEngine + "_full.root";

  gDebug = 0;

  // -----   Create simulation run   ----------------------------------------
  FairRunSim* run = new FairRunSim(); // singleton! FIXME
  // FairRunSim* run = FairRunSim::Instance(); // should be this
  run->SetName(mcEngine);              // Transport engine
  // run->SetSimulationConfig(new FairVMCConfig());
  run->SetIsMT(isMT);                  // Multi-threading mode (Geant4 only)
  run->SetSink(new FairRootFileSink(outFile));

  // -----   Create media   -------------------------------------------------
  run->SetMaterials("media.geo");

  // -----   Create geometry   ----------------------------------------------
  unique_ptr<FairModule> cave(new FairCave("CAVE"));
  cave->SetGeometryFileName("cave_vacuum.geo");
  run->AddModule(cave.release());

  unique_ptr<FairTutGeaneDet> det(new FairTutGeaneDet("TutGeaneDetector", kTRUE));
  det->SetGeometryFileName("pixel.geo");
  run->AddModule(det.release());

  // -----   Create PrimaryGenerator   --------------------------------------
  unique_ptr<FairPrimaryGenerator> primGen(new FairPrimaryGenerator());

  unique_ptr<FairBoxGenerator> boxGen(new FairBoxGenerator(static_cast<Int_t>(pdgCode::pion), 1));
  boxGen->SetThetaRange (   theta,    theta+0.01);
  boxGen->SetPRange     (momentum, momentum+0.01);
  boxGen->SetPhiRange   (0., 360.);
  boxGen->SetDebug      (kTRUE);

  primGen->AddGenerator(boxGen.release());

  run->SetGenerator(primGen.release());

  // -----   Create magnetic field   ----------------------------------------
  unique_ptr<FairConstField> magField(new FairConstField());
  magField->SetField      (0., 0., 20.); // values are in kG
  magField->SetFieldRegion(-150, 150, -150, 150, -250, 250);// values are in cm (xmin,xmax,ymin,ymax,zmin,zmax)

  run->SetField(magField.release());

  // -----   Initialize simulation run   ------------------------------------
  run->Init();

  // -----   Runtime database   ---------------------------------------------
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  Bool_t kParameterMerged = kTRUE;
  unique_ptr<FairParRootFileIo> parOut(new FairParRootFileIo(kParameterMerged));
  parOut->open(parFile.Data());
  rtdb->setOutput(parOut.release());
  rtdb->saveOutput();
  rtdb->print();

  // -----   Run simulation   -----------------------------------------------
  gRandom->SetSeed(123456);
  // gRandom.SetSeed(); // Use this to get a random seed

  TStopwatch timer;
  timer.Start();
  run->Run(nEvents);
  timer.Stop();

  // -----   Summary   -------------------------------------------------------
  run->CreateGeometryFile(geoFile);

  cout << endl << endl;

  // Extract the maximal used memory an add is as Dart measurement
  // This line is filtered by CTest and the value send to CDash
  FairSystemInfo sysInfo;
  Float_t maxMemory=sysInfo.GetMaxMemory();
  cout << "<DartMeasurement name=\"MaxMemory\" type=\"numeric/double\">";
  cout << maxMemory;
  cout << "</DartMeasurement>" << endl;

  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  Float_t cpuUsage=ctime/rtime;
  cout << "<DartMeasurement name=\"CpuLoad\" type=\"numeric/double\">";
  cout << cpuUsage;
  cout << "</DartMeasurement>" << endl;

  cout << endl << endl;
  cout << "Output file is "    << outFile << endl;
  cout << "Parameter file is " << parFile << endl;
  cout << "Real time " << rtime << "s, CPU time "
                       << ctime << "s" << endl << endl;
  cout << "Macro finished successfully." << endl;

  return 0;
}
