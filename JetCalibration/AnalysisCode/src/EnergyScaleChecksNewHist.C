#include "EnergyScaleChecksNewHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace JetCalibration {

//const Color_t colors[15] = {kAzure, kViolet, kMagenta, kPink+10, kPink, kOrange+10, kOrange, kSpring+10, kSpring, kTeal+10, kTeal, kAzure+10, kGray, kBlack, kGray+2};
const Color_t colors[12] = {kAzure, kOrange+7, kSpring+10, kTeal+10, kAzure+10, kRed, kPink+6, kViolet, kSpring, kTeal, kBlack, kGray};

void EnergyScaleChecksNewHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("EnergyScaleChecks/", "JetCalibration/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* JetTree = (TTree*)inFile->Get ("jeffsjets");
  //TTree* PhotonTree = (TTree*)inFile->Get ("jeffsphotons");

  float evtWeight;
  float jpt, jeta, jphi, je, jes, pjes, jpts, pjpts;
  float ppt, peta, pphi, pes;

  JetTree->SetBranchAddress ("evt_weight", &evtWeight);
  JetTree->SetBranchAddress ("jet_pt", &jpt);
  JetTree->SetBranchAddress ("jet_eta", &jeta);
  JetTree->SetBranchAddress ("jet_phi", &jphi);
  JetTree->SetBranchAddress ("jet_e", &je);
  JetTree->SetBranchAddress ("jet_energy_scale", &jes);
  JetTree->SetBranchAddress ("precalib_jet_energy_scale", &pjes);
  JetTree->SetBranchAddress ("jet_pt_scale", &jpts);
  JetTree->SetBranchAddress ("precalib_jet_pt_scale", &pjpts);

  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");

  TH2D* jetEnergyRespDist_phi = new TH2D ("jetEnergyRespDist_phi", "", numphibins, phibins, 200, 0, 2);
  TH3D* jetEnergyRespDist_pt_eta = new TH3D ("jetEnergyRespDist_pt_eta", "", numpbins, pbins, numetabins, etabins, 200, linspace (0, 2, 200));

  const int njets = JetTree->GetEntries ();

  for (int iJ = 0; iJ < njets; iJ++) {
   JetTree->GetEntry (iJ);
   if (jeta < 1.5 || 3.2 < jeta)
    continue;
   //if (jpt < 60)
   // continue;

   jetEnergyRespDist_phi->Fill (jphi, jes, evtWeight);
   jetEnergyRespDist_pt_eta->Fill (jpt, jeta, jes, evtWeight);
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plot energy scale & resolution as a function of phi
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D* jetEnergyScale_phi = new TH1D ("jetEnergyScale_phi", "", numphibins, phibins);
  TH1D* jetEnergyRes_phi = new TH1D ("jetEnergyRes_phi", "", numphibins, phibins);
  for (int iPhi = 1; iPhi <= numphibins; iPhi++) {
   TH1D* projy = jetEnergyRespDist_phi->ProjectionY ("_py", iPhi, iPhi);

   TF1* fit = new TF1 ("fit", "gaus(0)", projy->GetMean ()-2*projy->GetStdDev (), projy->GetMean ()+2*projy->GetStdDev ());
   projy->Fit (fit, "R0Q");

   float m = fit->GetParameter (1);
   float me = fit->GetParError (1);
   float s = fit->GetParameter (2);
   float se = fit->GetParError (2);

   if (projy) { delete projy; projy = NULL; }
   if (fit) { delete fit; fit = NULL; }

   jetEnergyScale_phi->SetBinContent (iPhi, m);
   jetEnergyScale_phi->SetBinError (iPhi, me);

   jetEnergyRes_phi->SetBinContent (iPhi, s);
   jetEnergyRes_phi->SetBinError (iPhi, se);
  }

  TCanvas* jetEnergyScale_phi_canvas = new TCanvas ("jetEnergyScale_phi_canvas", "", 800, 600);
  jetEnergyScale_phi_canvas->cd ();

  jetEnergyScale_phi->GetXaxis ()->SetTitle ("#phi");
  jetEnergyScale_phi->GetYaxis ()->SetTitle ("#mu, JES");
  jetEnergyScale_phi->GetYaxis ()->SetRangeUser (0.8, 1.2);

  jetEnergyScale_phi->Draw ("e1 x0");
  
  jetEnergyScale_phi_canvas->SaveAs (Form ("%s/jetEnergyScale_phi.pdf", plotPath.Data ()));


  TCanvas* jetEnergyRes_phi_canvas = new TCanvas ("jetEnergyRes_phi_canvas", "", 800, 600);
  jetEnergyRes_phi_canvas->cd ();

  jetEnergyRes_phi->GetXaxis ()->SetTitle ("#phi");
  jetEnergyRes_phi->GetYaxis ()->SetTitle ("#sigma, JER");

  jetEnergyRes_phi->Draw ("e1 x0");
  
  jetEnergyRes_phi_canvas->SaveAs (Form ("%s/jetEnergyRes_phi.pdf", plotPath.Data ()));

  jetEnergyRespDist_phi->Write ();
  jetEnergyScale_phi->Write ();
  jetEnergyRes_phi->Write ();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plot energy scale & resolution as a function of pT
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* jetEnergyScale_pt_canvas = new TCanvas ("jetEnergyScale_pt_canvas", "", 800, 600);
  TCanvas* jetEnergyRes_pt_canvas = new TCanvas ("jetEnergyRes_pt_canvas", "", 800, 600);
  for (short iEta = 1; iEta <= numetabins; iEta++) {
   TH2D* proj2d = Project2D ("proj2d", jetEnergyRespDist_pt_eta, "x", "z", iEta, iEta, false);

   TH1D* jetEnergyScale_pt = new TH1D (Form ("jetEnergyScale_pt_iEta%i", iEta), "", numpbins, pbins);
   TH1D* jetEnergyRes_pt = new TH1D (Form ("jetEnergyRes_pt_iEta%i", iEta), "", numpbins, pbins);

   for (short iP = 1; iP <= numpbins; iP++) {
    TH1D* projy = proj2d->ProjectionY ("_py", iP, iP);

    TF1* fit = new TF1 ("fit", "gaus(0)", projy->GetMean ()-2*projy->GetStdDev (), projy->GetMean ()+2*projy->GetStdDev ());
    projy->Fit (fit, "R0Q");

    float m = fit->GetParameter (1);
    float me = fit->GetParError (1);
    float s = fit->GetParameter (2);
    float se = fit->GetParError (2);

    if (projy) { delete projy; projy = NULL; }
    if (fit) { delete fit; fit = NULL; }

    jetEnergyScale_pt->SetBinContent (iP, m);
    jetEnergyScale_pt->SetBinError (iP, me);

    jetEnergyRes_pt->SetBinContent (iP, s);
    jetEnergyRes_pt->SetBinError (iP, se);
   }
   if (proj2d) { delete proj2d; proj2d = NULL; }


   jetEnergyScale_pt_canvas->cd ();
   gPad->SetLogx ();

   TGraphAsymmErrors* jetEnergyScale_pt_graph = make_graph (jetEnergyScale_pt);
   deltaize (jetEnergyScale_pt_graph, (numetabins/2-iEta), false);

   jetEnergyScale_pt_graph->SetLineColor (colors[iEta]);
   jetEnergyScale_pt_graph->SetMarkerColor (colors[iEta]);
   jetEnergyScale_pt_graph->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet}");
   jetEnergyScale_pt_graph->GetYaxis ()->SetTitle ("#mu, JES");
   jetEnergyScale_pt_graph->GetYaxis ()->SetRangeUser (0.9, 1.2);

   ( (TGraphAsymmErrors*)jetEnergyScale_pt_graph->Clone ())->Draw (iEta==1?"ap":"p");


   jetEnergyRes_pt_canvas->cd ();
   gPad->SetLogx ();

   TGraphAsymmErrors* jetEnergyRes_pt_graph = make_graph (jetEnergyRes_pt);
   deltaize (jetEnergyScale_pt_graph, (numetabins/2-iEta), false);

   jetEnergyRes_pt_graph->SetLineColor (colors[iEta]);
   jetEnergyRes_pt_graph->SetMarkerColor (colors[iEta]);
   jetEnergyRes_pt_graph->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet}");
   jetEnergyRes_pt_graph->GetYaxis ()->SetTitle ("#sigma, JER");
   ( (TGraphAsymmErrors*)jetEnergyRes_pt_graph->Clone ())->Draw (iEta==1?"ap":"p");

   jetEnergyScale_pt->Write ();
   jetEnergyRes_pt->Write ();

   if (jetEnergyScale_pt) { delete jetEnergyScale_pt; jetEnergyScale_pt = NULL; }
   if (jetEnergyRes_pt) { delete jetEnergyRes_pt; jetEnergyRes_pt = NULL; }
   if (jetEnergyScale_pt_graph) { delete jetEnergyScale_pt_graph; jetEnergyScale_pt_graph = NULL; }
   if (jetEnergyRes_pt_graph) { delete jetEnergyRes_pt_graph; jetEnergyRes_pt_graph = NULL; }
  }

  jetEnergyScale_pt_canvas->SaveAs (Form ("%s/jetEnergyScale_pt.pdf", plotPath.Data ()));
  jetEnergyRes_pt_canvas->SaveAs (Form ("%s/jetEnergyRes_pt.pdf", plotPath.Data ()));

  jetEnergyRespDist_pt_eta->Write ();

  outFile->Close ();
  

}

} // end namespace
