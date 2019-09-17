#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematic.h"

const bool doSys = true;
// nominal analyses
FullAnalysis* data18 = nullptr;
FullAnalysis* data15 = nullptr;
MinbiasAnalysis* bkg18 = nullptr;
MinbiasAnalysis* bkg15 = nullptr;
MCAnalysis* mc = nullptr;
TruthAnalysis* truth = nullptr;

// master systematics objects
Systematic* combSys = nullptr;
Systematic* bkgSys = nullptr;
Systematic* trkSys = nullptr;
Systematic* trkEffSys = nullptr;
Systematic* trkPurSys = nullptr;
Systematic* leptonRejSys = nullptr;
Systematic* electronPtSys = nullptr;
Systematic* muonPtSys = nullptr;

// variations for systematics
PhysicsAnalysis* data_trackHItight = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr;
PhysicsAnalysis* data_trackEff = nullptr;
MinbiasAnalysis* bkg_trackEff = nullptr;
PhysicsAnalysis* data_trackPurity = nullptr;
MinbiasAnalysis* bkg_trackPurity = nullptr;
PhysicsAnalysis* data_electronPtUp = nullptr, *data_electronPtDown = nullptr;
MinbiasAnalysis* bkg_electronPtUp = nullptr, *bkg_electronPtDown = nullptr;
PhysicsAnalysis* data_muonPtUp = nullptr, *data_muonPtDown = nullptr;
MinbiasAnalysis* bkg_muonPtUp = nullptr, *bkg_muonPtDown = nullptr;
PhysicsAnalysis* data_leptonRejVar = nullptr;

MinbiasAnalysis* bkg_statUpVar = nullptr, *bkg_statDownVar = nullptr;
PhysicsAnalysis* data_bkgStatUpVar = nullptr, *data_bkgStatDownVar = nullptr;

void Run () {

  data18    = new FullAnalysis ("data18");
  data15  = new FullAnalysis ("data15");
  data15->is2015Conds = true;
  data15->useHijingEffs = true;
  bkg18     = new MinbiasAnalysis ();
  bkg15     = new MinbiasAnalysis ();
  bkg15->is2015Conds = true;
  bkg15->useHijingEffs = true;
  mc      = new MCAnalysis ();
  truth   = new TruthAnalysis ();

  if (doSys) {
    data_trackHItight       = new PhysicsAnalysis ("data_trackHITightVar");
    data_trackHItight->useHITight = true;
    bkg_trackHItight        = new MinbiasAnalysis ("bkg_trackHITightVar");
    bkg_trackHItight->useHITight = true;

    data_trackEff        = new PhysicsAnalysis ("data_trackEffVar");
    data_trackEff->doTrackEffVar = true;
    bkg_trackEff         = new MinbiasAnalysis ("bkg_trackEffVar");
    bkg_trackEff->doTrackEffVar = true;

    data_trackPurity        = new PhysicsAnalysis ("data_trackPurityVar");
    data_trackPurity->doTrackPurVar = true;
    bkg_trackPurity         = new MinbiasAnalysis ("bkg_trackPurityVar");
    bkg_trackPurity->doTrackPurVar = true;

    //data_leptonRejVar       = new PhysicsAnalysis ("data_leptonRejVar");
    //data_leptonRejVar->doLeptonRejVar = true;

    //data_electronPtUp       = new PhysicsAnalysis ("data_electronPtUpVar");
    //data_electronPtDown     = new PhysicsAnalysis ("data_electronPtDownVar");
    //data_muonPtUp           = new PhysicsAnalysis ("data_muonPtUpVar");
    //data_muonPtDown         = new PhysicsAnalysis ("data_muonPtDownVar");
    //bkg_electronPtUp       = new MinbiasAnalysis ("bkg");//_electronPtUpVar");
    //bkg_electronPtDown     = new MinbiasAnalysis ("bkg");//_electronPtDownVar");
    //bkg_muonPtUp           = new MinbiasAnalysis ("bkg");//_muonPtUpVar");
    //bkg_muonPtDown         = new MinbiasAnalysis ("bkg");//_muonPtDownVar");

    //data_bkgStatUpVar       = new PhysicsAnalysis ("data_bkgStatUpVar");
    //data_bkgStatDownVar     = new PhysicsAnalysis ("data_bkgStatDownVar");
    //bkg_statUpVar           = new MinbiasAnalysis ("bkg_statUpVar");
    //bkg_statDownVar         = new MinbiasAnalysis ("bkg_statDownVar");
  }

  //data18->Execute   ("DataAnalysis/Nominal/data18hi.root",   "DataAnalysis/Nominal/data18hi_hists.root");
  //data15->Execute ("DataAnalysis/Nominal/data15hi.root",  "DataAnalysis/Nominal/data15hi_hists.root");

  if (doSys) {
    //data_electronPtUp->Execute      ("DataAnalysis/Variations/ElectronPtUpVariation/data18hi.root",        "DataAnalysis/Variations/ElectronPtUpVariation/data18hi_hists.root");
    //data_electronPtDown->Execute    ("DataAnalysis/Variations/ElectronPtDownVariation/data18hi.root",      "DataAnalysis/Variations/ElectronPtDownVariation/data18hi_hists.root");
    //data_muonPtUp->Execute          ("DataAnalysis/Variations/MuonPtUpVariation/data18hi.root",            "DataAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    //data_muonPtDown->Execute        ("DataAnalysis/Variations/MuonPtDownVariation/data18hi.root",          "DataAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");
    //data_leptonRejVar->Execute      ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");
    //data_trackHItight->Execute      ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi.root",      "DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    //data_trackEff->Execute          ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    //data_trackPurity->Execute       ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
  }


  data18->LoadHists   ("DataAnalysis/Nominal/data18hi_hists.root");
  data15->LoadHists   ("DataAnalysis/Nominal/data15hi_hists.root");
  bkg18->LoadHists      ("MinbiasAnalysis/Nominal/data18hi_hists.root");
  bkg15->LoadHists      ("MinbiasAnalysis/Nominal/data15hi_hists.root");
  mc->LoadHists       ("MCAnalysis/Nominal/savedHists.root");
  truth->LoadHists  ("TruthAnalysis/Nominal/savedHists.root");

  //if (doSys) {
  //  data_bkgStatUpVar->CopyAnalysis (data18);
  //  data_bkgStatDownVar->CopyAnalysis (data18);
  //  bkg_statUpVar->CopyAnalysis (bkg18);
  //  bkg_statUpVar->ConvertToStatVariation (true, 1);
  //  bkg_statDownVar->CopyAnalysis (bkg18);
  //  bkg_statDownVar->ConvertToStatVariation (false, 1);
  //  data_bkgStatUpVar->SubtractBackground (bkg_statUpVar);
  //  data_bkgStatDownVar->SubtractBackground (bkg_statDownVar);
  //  data_bkgStatUpVar->CalculateIAA ();
  //  data_bkgStatUpVar->CalculateICP ();
  //  data_bkgStatDownVar->CalculateIAA ();
  //  data_bkgStatDownVar->CalculateICP ();
  //}

  data18->SubtractBackground (bkg18);
  data18->CalculateIAA ();
  data18->CalculateICP ();

  data15->SubtractBackground (bkg15);
  data15->CalculateIAA ();
  data15->CalculateICP ();

  data18->CalculateZPtDistRatio (truth);
  data18->CalculateZEtaDistRatio ();
  data18->CalculateZYDistRatio (truth);
  data18->CalculateZMassSpectraRatio (mc);
  data15->CalculateZMassSpectraRatio (mc);

  mc->SubtractBackground (bkg18);
  mc->CalculateIAA ();
  mc->CalculateICP ();
  truth->SubtractBackground ();

  if (doSys) {
    data_trackHItight->LoadHists      ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    bkg_trackHItight->LoadHists       ("MinbiasAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    data_trackEff->LoadHists          ("DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    bkg_trackEff->LoadHists           ("MinbiasAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    data_trackPurity->LoadHists       ("DataAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    bkg_trackPurity->LoadHists        ("MinbiasAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    //data_leptonRejVar->LoadHists      ("DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");
    //data_electronPtUp->LoadHists      ("DataAnalysis/Variations/ElectronPtUpVariation/data18hi_hists.root");
    //bkg_electronPtUp->LoadHists       ("MinbiasAnalysis/Variations/ElectronPtUpVariation/data18hi_hists.root");
    //data_electronPtDown->LoadHists    ("DataAnalysis/Variations/ElectronPtDownVariation/data18hi_hists.root");
    //bkg_electronPtDown->LoadHists     ("MinbiasAnalysis/Variations/ElectronPtDownVariation/data18hi_hists.root");
    //data_muonPtUp->LoadHists          ("DataAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    //bkg_muonPtUp->LoadHists           ("MinbiasAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    //data_muonPtDown->LoadHists        ("DataAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");
    //bkg_muonPtDown->LoadHists         ("MinbiasAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");

    data_trackHItight->SubtractBackground (bkg_trackHItight);
    data_trackEff->SubtractBackground (bkg_trackEff);
    data_trackPurity->SubtractBackground (bkg_trackPurity);
    //data_leptonRejVar->SubtractBackground (bkg18);
    //data_electronPtUp->SubtractBackground (bkg18);
    //data_electronPtDown->SubtractBackground (bkg18);
    //data_muonPtUp->SubtractBackground (bkg18);
    //data_muonPtDown->SubtractBackground (bkg18);

    //data_electronPtUp->SubtractBackground (bkg_electronPtUp);
    //data_electronPtDown->SubtractBackground (bkg_electronPtDown);
    //data_muonPtUp->SubtractBackground (bkg_muonPtUp);
    //data_muonPtDown->SubtractBackground (bkg_muonPtDown);
/*
    bkgSys = new Systematic (data18, "bkgSys", "Background");
    bkgSys->AddVariation (data_bkgStatUpVar, false);
    bkgSys->AddVariation (data_bkgStatDownVar, false);
    bkgSys->AddVariations ();

    trkSys = new Systematic (data18, "trkSys", "Track ID");
    trkSys->AddVariation (data_trackHItight);
    trkSys->AddVariations ();

    trkEffSys = new Systematic (data18, "trkEffSys", "Tracking Efficiency");
    trkEffSys->AddVariation (data_trackEff);
    trkEffSys->AddVariations ();

    trkPurSys = new Systematic (data18, "trkPurSys", "Tracking Purity");
    trkPurSys->AddVariation (data_trackPurity);
    trkPurSys->AddVariations ();

    leptonRejSys = new Systematic (data18, "leptonRejSys", "Lepton Rejection");
    leptonRejSys->AddVariation (data_leptonRejVar, -1);
    leptonRejSys->AddVariations ();

    electronPtSys = new Systematic (data18, "electronPtSys", "Electron ES");
    electronPtSys->AddVariation (data_electronPtUp);
    electronPtSys->AddVariation (data_electronPtDown);
    electronPtSys->AddVariations ();

    muonPtSys = new Systematic (data18, "muonPtSys", "Muon ES");
    muonPtSys->AddVariation (data_muonPtUp);
    muonPtSys->AddVariation (data_muonPtDown);
    muonPtSys->AddVariations ();

    combSys = new Systematic (data18, "combSys", "Total");
    combSys->AddSystematic (trkSys);
    combSys->AddSystematic (trkEffSys);
    combSys->AddSystematic (trkPurSys);
    combSys->AddSystematic (bkgSys);
    combSys->AddSystematic (leptonRejSys);
    combSys->AddSystematic (electronPtSys);
    combSys->AddSystematic (muonPtSys);
    combSys->AddSystematics ();
*/
  }


  SetupDirectories ("", "ZTrackAnalysis/");

}


void MakePhysicsPlots () {

  data18->PlotTrkYield (1, 0, 2, 4);
  combSys->PlotTrkYield (1, 1, 2, 4);
  bkg18->PlotTrkYield (1, 0, 2, 4);
  data18->PlotTrkYield (1, 0, 2, 4);

  max_iaa=4.4;
  combSys->PlotIAAdPhi (1, 1, 2, 4);
  data18->PlotIAAdPhi (1, 0, 2, 4);
  max_icp=3.;
  combSys->PlotICPdPhi (1, 1, 2, 4);
  data18->PlotICPdPhi (1, 0, 2, 4);
  

  combSys->PlotTrkYieldZPt (1, 1, 2);
  data18->PlotTrkYieldZPt (1, 0, 2);

  combSys->PlotTrkYieldZPt (0, 1, 2);
  data18->PlotTrkYieldZPt (0, 0, 2);

  max_iaa = 10;
  combSys->PlotSingleIAAdPtZ (1, 1);
  data18->PlotSingleIAAdPtZ (1, 0);

  combSys->PlotSingleIAAdPtZ (0, 1);
  data18->PlotSingleIAAdPtZ (0, 0);

  max_iaa = 3.2;
  combSys->PlotIAAdPtZ (1, 1);
  data18->PlotIAAdPtZ (1, 0);

  max_iaa = 3.2;
  combSys->PlotIAAdPtZ (0, 1);
  data18->PlotIAAdPtZ (0, 0);

  max_icp = 3.;
  combSys->PlotICPdPtZ (1, 1);
  data18->PlotICPdPtZ (1, 0);
  
}




void ComparePbPbSubYields (const short iSpc = 2, const short iPtZ = nPtZBins-1) {
  TCanvas* c = new TCanvas ("c", "", 1600, 800);
  const double dPadY = 0.4;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h1 = nullptr, *h2 = nullptr, *hrat = nullptr;
  TGraphAsymmErrors* g = nullptr;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* uPadName = Form ("uPad_%i", iCent);
    const char* dPadName = Form ("dPad_%i", iCent);

    TPad* uPad = new TPad (uPadName, "", (1./(numCentBins))*(iCent), dPadY, (1./(numCentBins))*(iCent+1), 1);
    TPad* dPad = new TPad (dPadName, "", (1./(numCentBins))*(iCent), 0, (1./(numCentBins))*(iCent+1), dPadY);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (0.17);
    uPad->SetRightMargin (0.06);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    dPad->SetLeftMargin (0.17);
    dPad->SetRightMargin (0.06);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    h1 = bkg18->h_z_trk_zpt[iSpc][iPtZ][iCent];
    h2 = bkg_trackHItight->h_z_trk_zpt[iSpc][iPtZ][iCent];
    const float min = fmin (h1->GetMinimum (0), h2->GetMinimum (0));
    const float max = fmax (h1->GetMaximum (),  h2->GetMaximum ());

    g = data18->GetTGAE (h1);
    if (iCent == 3)
      g->Print ("ALL");

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);
    g->SetFillColorAlpha (fillColors[1], 0.3);

    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("dY / #Delta#phi d#it{p}_{T} [GeV^{-1}]");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->Draw ("AP");

    g = data_trackHItight->GetTGAE (h2);
    if (iCent == 3)
      g->Print ("ALL");

    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);
    g->SetFillColorAlpha (fillColors[2], 0.3);

    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("dY / #Delta#phi d#it{p}_{T} [GeV^{-1}]");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->Draw ("P");
    if (iCent == 0)
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.04/uPadY);
    else
      myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/uPadY);

    dPad->cd ();
    dPad->SetLogx ();

    hrat = (TH1D*) h2->Clone ();
    hrat->Divide (h1);

    g = make_graph (hrat);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);
    g->SetFillColorAlpha (fillColors[2], 0.3);

    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (0.6, 1.4);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("HITight / HILoose"); 

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->GetYaxis ()->CenterTitle ();

    g->Draw ("AP");

    TF1* fit = new TF1 ("fit", "[0]", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    fit->SetParameter (0, 1);
    g->Fit (fit, "RQN0");

    fit->SetLineColor (kPink-8);
    fit->Draw ("same");

    myText (0.22, 0.3, kBlack, Form ("Avg. = %g#pm%g", fit->GetParameter (0), fit->GetParError (0)), 0.02/dPadY);
  }

}

#endif
