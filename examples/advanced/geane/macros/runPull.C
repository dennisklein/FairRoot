/********************************************************************************
 *    Copyright (C) 2019 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include <vector>

int runPull(bool drawHist = false)
{
  gROOT->Reset();
  gStyle->SetOptFit(1);

  TFile f("geane.cal.root");
  TTree* simtree = static_cast<TTree*>(f.Get("cbmsim"));

  std::vector<FairTrackParP>  trackParGeane;
  std::vector<FairTrackParP>* trackParGeanePtr; // required by TTree::SetBranchAddress
  std::vector<FairTrackParP>  trackParFinal;
  std::vector<FairTrackParP>* trackParFinalPtr; // required by TTree::SetBranchAddress
  simtree->SetBranchAddress("GeaneTrackPar",   &trackParGeanePtr);
  simtree->SetBranchAddress("GeaneTrackFinal", &trackParFinalPtr);

  TH1F hQP("hQP","charge over momentum",100,-10.,10.);
  TH1F hX ("hX", "position X",          100,-10,10);
  TH1F hY ("hY", "position Y",          100,-10,10);
  TH1F hPx("hPx","momentum X",          100,-10,10);
  TH1F hPy("hPy","momentum Y",          100,-10,10);
  TH1F hPz("hPz","momentum Z",          100,-10,10);

  Int_t nEvents = simtree->GetEntriesFast();
  cout << nEvents << endl;

  for (Int_t i = 0; i < nEvents; ++i) {
    simtree->GetEntry(i);
    for (Int_t k = 0; k < trackParGeane.size(); ++k)	{
      try {
        auto trkG(trackParGeane.at(k));
        auto trkF(trackParFinal.at(k));
        if(trkG.GetDQp()) hQP.Fill( ( trkF.GetQp() - trkG.GetQp() ) / trkG.GetDQp() );
        if(trkG.GetDX() ) hX .Fill( ( trkF.GetX()  - trkG.GetX()  ) / trkG.GetDX()  );
        if(trkG.GetDY() ) hY .Fill( ( trkF.GetY()  - trkG.GetY()  ) / trkG.GetDY()  );
        if(trkG.GetDPx()) hPx.Fill( ( trkF.GetPx() - trkG.GetPx() ) / trkG.GetDPx() );
        if(trkG.GetDPy()) hPy.Fill( ( trkF.GetPy() - trkG.GetPy() ) / trkG.GetDPy() );
        if(trkG.GetDPz()) hPz.Fill( ( trkF.GetPz() - trkG.GetPz() ) / trkG.GetDPz() );
      } catch (const std::out_of_range& ex) {
        break;
      }
    }
  }

  if (drawHist) {
    TCanvas c("c", "c", 900, 600);
    c.Divide(3, 2);
    c.cd(1);
    hQP.Draw();
    hQP.Fit("gaus");
    c.cd(2);
    hX.Draw();
    hX. Fit("gaus");
    c.cd(3);
    hY.Draw();
    hY. Fit("gaus");
    c.cd(4);
    hPx.Draw();
    hPx. Fit("gaus");
    c.cd(5);
    hPy.Draw();
    hPy. Fit("gaus");
    c.cd(6);
    hPz.Draw();
    hPz. Fit("gaus");
    c.cd();
  }

  TF1 fitX("fitX", "gaus", -5., 5.);
  hX.Fit("fitX", "QN");

  if ( fitX.GetParameter(1) > -0.03 && fitX.GetParameter(1) < 0.03 &&
       fitX.GetParameter(2) >  0.   && fitX.GetParameter(2) < 1.5 )
      cout << "Macro finished successfully. Mean (" << fitX.GetParameter(1) << ") and sigma (" << fitX.GetParameter(2) << ") inside limits." << endl;
  else
      cout << "Macro failed. Mean (" << fitX.GetParameter(1) << ") or sigma (" << fitX.GetParameter(2) << ") too far off." << endl;

  return 0;
}
