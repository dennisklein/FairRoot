/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <memory>

void printSummary(TStopwatch& timer, const TString& parFile, const TString& outFile)
{
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

  Float_t cpuUsage = ctime/rtime;
  cout << "<DartMeasurement name=\"CpuLoad\" type=\"numeric/double\">";
  cout << cpuUsage;
  cout << "</DartMeasurement>" << endl;

  cout << endl << endl;
  cout << "Output file is "    << outFile << endl;
  cout << "Parameter file is " << parFile << endl;
  cout << "Real time " << rtime << "s, CPU time "
                       << ctime << "s" << endl << endl;

  cout << "Macro finished successfully." << endl;
}

int runGeane()
{
  // -----   Parameters   ---------------------------------------------------
  TString dir = getenv("VMCWORKDIR");

  TString tutGeomDir = dir + "/common/geometry";
  gSystem->Setenv("GEOMPATH", tutGeomDir.Data());

  TString tutConfigDir = dir + "/common/gconfig";
  gSystem->Setenv("CONFIG_DIR", tutConfigDir.Data());

  TString inFile = "geane.mc.root"; // Input file (MC events)
  TString parFile = "geane.par.root"; // Parameter file
  TString outFile = "geane.cal.root"; // Output file

  // -----   Create run   ---------------------------------------------------
  FairRunAna* run = new FairRunAna(); // singleton! FIXME
  // FairRunAna* run = FairRunAna::Instance(); // should be this
  run->SetSource(new FairFileSource  (inFile));
  run->SetSink  (new FairRootFileSink(outFile));

  FairRuntimeDb* rtdb = run->GetRuntimeDb();

  unique_ptr<FairParRootFileIo> parInput1(new FairParRootFileIo());
  parInput1->open(parFile.Data());
  rtdb->setFirstInput(parInput1.release());

  unique_ptr<FairConstField> magField(new FairConstField());
  magField->SetField      (0., 0., 20.); // values are in kG
  magField->SetFieldRegion(-150, 150, -150, 150, -250, 250); // values are in cm (xmin,xmax,ymin,ymax,zmin,zmax)
  run->SetField(magField.release());

  // -----   TorinoDetector hit producers   ---------------------------------
  run->AddTask(new FairGeane());
  run->AddTask(new FairTutGeaneTr());

  // -----   Initialize reconstruction run   --------------------------------
  run->Init();

  // -----   Run   ----------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  run->Run();
  timer.Stop();

  // -----   Summary  -------------------------------------------------------
  printSummary(timer, parFile, outFile);

  return 0;
}
