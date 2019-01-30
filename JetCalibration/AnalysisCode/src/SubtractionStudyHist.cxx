#include "SubtractionStudyHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <TTree.h>
#include <TVirtualFitter.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace atlashi;

namespace JetCalibration {

void SubtractionStudyHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("SubtractionStudy/", "JetCalibration/");

  //////////////////////////////////////////////////////////////////////////////
  // Create desired histograms
  //////////////////////////////////////////////////////////////////////////////
  TH2D*** jetSubHist = Get2DArray <TH2D*> (3, 2); // iPer, iData

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* per = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
      const char* data = (iData == 0 ? "data" : "mc");

      //jetSubHist[iPer][iData] = new TH2D (Form ("jetSubHist_%s_%s", per, data), "", numphibins, phibins, 60, -5, 25);
      jetSubHist[iPer][iData] = new TH2D (Form ("jetSubHist_%s_%s", per, data), "", 45, -4.5, 4.5, 60, -5, 25);
      jetSubHist[iPer][iData]->Sumw2 ();
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Load analyzed TTrees
  //////////////////////////////////////////////////////////////////////////////
  float subjpt = 0, subjeta = 0, subjphi = 0, subje = 0, unsubjpt = 0, unsubjeta = 0, unsubjphi = 0, unsubje = 0, fcalScaleFactor = 0, ppt = 0, peta = 0, pphi = 0; 
  double evtWeight = 0;
  bool isMC = false, isPeriodA = false;

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* inTree = (TTree*)inFile->Get ("SubtractionTree");

  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("sub_jet_pt", &subjpt);
  inTree->SetBranchAddress ("sub_jet_eta", &subjeta);
  inTree->SetBranchAddress ("sub_jet_phi", &subjphi);
  inTree->SetBranchAddress ("sub_jet_e", &subje);
  inTree->SetBranchAddress ("unsub_jet_pt", &unsubjpt);
  inTree->SetBranchAddress ("unsub_jet_eta", &unsubjeta);
  inTree->SetBranchAddress ("unsub_jet_phi", &unsubjphi);
  inTree->SetBranchAddress ("unsub_jet_e", &unsubje);
  inTree->SetBranchAddress ("photon_pt", &ppt);
  inTree->SetBranchAddress ("photon_eta", &peta);
  inTree->SetBranchAddress ("photon_phi", &pphi);
  inTree->SetBranchAddress ("fcalScaleFactor", &fcalScaleFactor);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);

  //////////////////////////////////////////////////////////////////////////////
  // Fill desired histograms
  //////////////////////////////////////////////////////////////////////////////
  const int nJets = inTree->GetEntries ();
  for (int iJet = 0; iJet < nJets; iJet++) {
    inTree->GetEntry (iJet);

    if (ppt < 40)
      continue;

    if (subjpt < 20)
      continue;

    if (-pi < subjphi && subjphi < -pi/2)
      continue;

    const short iPer = isPeriodA ? 0 : 1;
    const short iMC = isMC ? 1 : 0;

    //jetSubHist[iPer][iMC]->Fill (subjphi, unsubjpt-subjpt);
    //jetSubHist[2][iMC]->Fill (subjphi, unsubjpt-subjpt);
    jetSubHist[iPer][iMC]->Fill (subjeta, unsubjpt-subjpt, evtWeight*fcalScaleFactor);
    jetSubHist[2][iMC]->Fill (subjeta, unsubjpt-subjpt, evtWeight*fcalScaleFactor);
  }
  inTree = NULL;
  inFile->Close ();
  if (inFile) { delete inFile; inFile = NULL; }

  //////////////////////////////////////////////////////////////////////////////
  // Save histograms for interactive access
  //////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");
  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iData = 0; iData < 2; iData++) {
      jetSubHist[iPer][iData]->Write ();
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Canvas definitions
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  FormatTH2Canvas (canvas, true);
  gStyle->SetPalette (kRainBow);

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
    const char* per = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    //jetSubHist[iPer][0]->GetXaxis ()->SetTitle ("Jet #phi");
    jetSubHist[iPer][0]->GetXaxis ()->SetTitle ("Jet #eta");
    jetSubHist[iPer][0]->GetYaxis ()->SetTitle ("#it{p}_{T}^{Constituent} - #it{p}_{T}^{EM} #left[GeV#right]");
    jetSubHist[iPer][0]->GetZaxis ()->SetTitle ("Counts");

    jetSubHist[iPer][0]->GetXaxis ()->SetTitleOffset (1);
    jetSubHist[iPer][0]->GetYaxis ()->SetTitleOffset (1);

    TProfile* _px = jetSubHist[iPer][0]->ProfileX ();
    TH1D* jetSubHist_px = TProfile2TH1D (Form ("jetSubHist_px_%s_data", per), _px, 45, linspace (-4.5, 4.5, 45));
    if (_px) { delete _px; _px = NULL; }

    jetSubHist_px->SetLineColor (kWhite);

    TH1D* jetSubHist_px_flipped = (TH1D*) jetSubHist_px->Clone (Form ("jetSubHist_px_flipped_%s_data", per));
    GetReflectionX (jetSubHist_px_flipped, jetSubHist_px_flipped->GetNbinsX () / 2);

    jetSubHist[iPer][0]->Draw ("colz");
    jetSubHist_px->Draw ("same hist");
    jetSubHist_px_flipped->Draw ("same hist");

    jetSubHist_px->Write ();
    jetSubHist_px_flipped->Write ();

    canvas->SaveAs (Form ("%s/Period%s/jetSubHist_data.pdf", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"))); 
   
    //jetSubHist[iPer][1]->GetXaxis ()->SetTitle ("Jet #phi");
    jetSubHist[iPer][1]->GetXaxis ()->SetTitle ("Jet #eta");
    jetSubHist[iPer][1]->GetYaxis ()->SetTitle ("#it{p}_{T}^{Constituent} - #it{p}_{T}^{EM} #left[GeV#right]");
    jetSubHist[iPer][1]->GetZaxis ()->SetTitle ("Counts");

    jetSubHist[iPer][1]->GetXaxis ()->SetTitleOffset (1);
    jetSubHist[iPer][1]->GetYaxis ()->SetTitleOffset (1);

     _px = jetSubHist[iPer][1]->ProfileX ();
    jetSubHist_px = TProfile2TH1D (Form ("jetSubHist_px_%s_mc", per), _px, 45, linspace (-4.5, 4.5, 45));
    if (_px) { delete _px; _px = NULL; }

    jetSubHist_px->SetLineColor (kWhite);

    jetSubHist_px_flipped = (TH1D*) jetSubHist_px->Clone (Form ("jetSubHist_px_flipped_%s_mc", per));
    GetReflectionX (jetSubHist_px_flipped, jetSubHist_px_flipped->GetNbinsX () / 2);

    jetSubHist[iPer][1]->Draw ("colz");
    jetSubHist_px->Draw ("same hist");
    jetSubHist_px_flipped->Draw ("same hist");

    jetSubHist_px->Write ();
    jetSubHist_px_flipped->Write ();

    canvas->SaveAs (Form ("%s/Period%s/jetSubHist_mc.pdf", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"))); 
  } // end loop over periods

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }

  return;
}

} // end namespace
