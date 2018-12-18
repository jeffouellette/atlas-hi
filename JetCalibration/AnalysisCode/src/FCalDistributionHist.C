#include "FCalDistributionHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>

#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasUtils.h>
#include <AtlasStyle.h>

namespace JetCalibration {

void FCalDistributionHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("FCalDistribution/", "JetCalibration/");

  TH1D* fCal_p_et[2];
  fCal_p_et[0] = new TH1D ("fCal_p_et_data", "", 350, -50, 300);
  fCal_p_et[0]->Sumw2 ();
  fCal_p_et[1] = new TH1D ("fCal_p_et_mc", "", 350, -50, 300);
  fCal_p_et[1]->Sumw2 ();

  TH1D* fCal_Pb_et[2];
  fCal_Pb_et[0] = new TH1D ("fCal_Pb_et_data", "", 350, -50, 300);
  fCal_Pb_et[0]->Sumw2 ();
  fCal_Pb_et[1] = new TH1D ("fCal_Pb_et_mc", "", 350, -50, 300);
  fCal_Pb_et[1]->Sumw2 ();

  //////////////////////////////////////////////////////////////////////////////
  // Load analyzed TTrees
  //////////////////////////////////////////////////////////////////////////////
  double evtWeight = 0, fcal_p_et = 0, fcal_Pb_et = 0;
  bool isMC = false, isPeriodA = false;

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* inTree = (TTree*)inFile->Get ("FCalTree");

  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("fCal_p_et", &fcal_p_et);
  inTree->SetBranchAddress ("fCal_Pb_et", &fcal_Pb_et);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);

  //////////////////////////////////////////////////////////////////////////////
  // Fill desired histograms
  //////////////////////////////////////////////////////////////////////////////
  const long nEvents = inTree->GetEntries ();
  for (long iEvent = 0; iEvent < nEvents; iEvent++) {
    inTree->GetEntry (iEvent);

    //const short iPer = isPeriodA ? 0 : 1;
    const short iMC = isMC ? 1 : 0;

    fCal_p_et[iMC]->Fill (fcal_p_et, evtWeight);
    fCal_Pb_et[iMC]->Fill (fcal_Pb_et, evtWeight);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Save histograms for interactive access
  //////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");
  fCal_p_et[0]->Write ();
  fCal_p_et[1]->Write ();
  fCal_Pb_et[0]->Write ();
  fCal_Pb_et[1]->Write ();
  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }

  //////////////////////////////////////////////////////////////////////////////
  // Canvas definitions
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  canvas->Draw ();

  Style_t markers[2] = {kFullCircle, kOpenCircle};
  for (int iData = 0; iData < 2; iData++) {
   fCal_p_et[iData]->SetLineColor (kBlue);
   fCal_p_et[iData]->SetMarkerColor (kBlue);
   fCal_Pb_et[iData]->SetLineColor (kBlack);
   fCal_Pb_et[iData]->SetMarkerColor (kBlack);

   fCal_p_et[iData]->SetMarkerStyle (markers[iData]);
   fCal_Pb_et[iData]->SetMarkerStyle (markers[iData]);

   fCal_p_et[iData]->GetXaxis ()->SetTitle ("Forward #SigmaE_{T} #left[GeV#right]");
   fCal_Pb_et[iData]->GetXaxis ()->SetTitle ("Forward #SigmaE_{T} #left[GeV#right]");
   fCal_p_et[iData]->GetYaxis ()->SetTitle ("Counts / Total");
   fCal_Pb_et[iData]->GetYaxis ()->SetTitle ("Counts / Total");

   fCal_p_et[iData]->Scale (1./fCal_p_et[iData]->Integral ());
   fCal_Pb_et[iData]->Scale (1./fCal_Pb_et[iData]->Integral ());

   if (iData == 0)
    fCal_p_et[iData]->Draw ("e1");
   else
    fCal_p_et[iData]->Draw ("same e1");
   fCal_Pb_et[iData]->Draw ("same e1");
  }

  myMarkerText (0.55, 0.85, kBlack, kFullCircle, "Pb-going direction", 1.25, 0.05);
  myMarkerText (0.55, 0.78, kBlue, kFullCircle, "p-going direction", 1.25, 0.05);

  canvas->SaveAs (Form ("%s/fcal_et.pdf", plotPath.Data ()));

  return;
}

} // end namespace
