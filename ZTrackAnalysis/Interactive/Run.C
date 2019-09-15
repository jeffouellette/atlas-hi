#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematic.h"

const bool doSys = false;

// nominal analyses
FullAnalysis* data18 = nullptr;
FullAnalysis* data18_1_0_0 = nullptr;
FullAnalysis* data15 = nullptr;
MCAnalysis* mc = nullptr;
MinbiasAnalysis* bkg = nullptr;
MinbiasAnalysis* bkg_1_0_0 = nullptr;
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
Systematic* electronLHMedSys = nullptr;
Systematic* muonTightSys = nullptr;

// variations for systematics
FullAnalysis* data_trigEff = nullptr;
PhysicsAnalysis* data_trackHItight = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr;
PhysicsAnalysis* data_trackEff = nullptr;
MinbiasAnalysis* bkg_trackEff = nullptr;
PhysicsAnalysis* data_trackPurity = nullptr;
MinbiasAnalysis* bkg_trackPurity = nullptr;
PhysicsAnalysis* data_electronPtUp = nullptr, *data_electronPtDown = nullptr;
PhysicsAnalysis* data_muonPtUp = nullptr, *data_muonPtDown = nullptr;
PhysicsAnalysis* data_electronLHMedium = nullptr;
PhysicsAnalysis* data_muonTight = nullptr;
PhysicsAnalysis* data_leptonRejVar = nullptr;

MinbiasAnalysis* bkg_statUpVar = nullptr, *bkg_statDownVar = nullptr;
PhysicsAnalysis* data_bkgStatUpVar = nullptr, *data_bkgStatDownVar = nullptr;

void Run () {

  data18    = new FullAnalysis ("data18");
  data15  = new FullAnalysis ("data15");
  data15->is2015Conds = true;
  data15->useHijingEffs = true;
  mc      = new MCAnalysis ();
  bkg     = new MinbiasAnalysis ();
  //truth   = new TruthAnalysis ();

  if (doSys) {
    data_trigEff            = new FullAnalysis ("data_trigEff");

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

    data_leptonRejVar       = new PhysicsAnalysis ("data_leptonRejVar");
    data_leptonRejVar->doLeptonRejVar = true;

    data_electronPtUp       = new PhysicsAnalysis ("data_electronPtUpVar");
    data_electronPtDown     = new PhysicsAnalysis ("data_electronPtDownVar");
    data_muonPtUp           = new PhysicsAnalysis ("data_muonPtUpVar");
    data_muonPtDown         = new PhysicsAnalysis ("data_muonPtDownVar");

    data_electronLHMedium   = new PhysicsAnalysis ("data_electronLHMediumVar");
    data_muonTight          = new PhysicsAnalysis ("data_muonTightVar");

    data_bkgStatUpVar       = new PhysicsAnalysis ("data_bkgStatUpVar");
    data_bkgStatDownVar     = new PhysicsAnalysis ("data_bkgStatDownVar");
    bkg_statUpVar           = new MinbiasAnalysis ("bkg_statUpVar");
    bkg_statDownVar         = new MinbiasAnalysis ("bkg_statDownVar");
  }

  data18->Execute   ("DataAnalysis/Nominal/outFile.root",   "DataAnalysis/Nominal/savedHists.root");
  //data15->Execute ("DataAnalysis/Nominal/data15hi.root",  "DataAnalysis/Nominal/data15hi_hists.root");

  if (doSys) {
    //data_trackHItight->Execute      ("DataAnalysis/Variations/TrackHITightWPVariation/outFile.root",      "DataAnalysis/Variations/TrackHITightWPVariation/savedHists.root");
    //data_trackEff->Execute          ("DataAnalysis/Nominal/outFile.root",                                 "DataAnalysis/Variations/TrackEffPionsVariation/savedHists.root");
    //data_trackPurity->Execute       ("DataAnalysis/Nominal/outFile.root",                                 "DataAnalysis/Variations/TrackPurityVariation/savedHists.root");
    //data_leptonRejVar->Execute      ("DataAnalysis/Nominal/outFile.root",                                 "DataAnalysis/Variations/LeptonRejVariation/savedHists.root");
    //data_electronPtUp->Execute      ("DataAnalysis/Variations/ElectronPtUpVariation/outFile.root",        "DataAnalysis/Variations/ElectronPtUpVariation/savedHists.root");
    //data_electronPtDown->Execute    ("DataAnalysis/Variations/ElectronPtDownVariation/outFile.root",      "DataAnalysis/Variations/ElectronPtDownVariation/savedHists.root");
    //data_muonPtUp->Execute          ("DataAnalysis/Variations/MuonPtUpVariation/outFile.root",            "DataAnalysis/Variations/MuonPtUpVariation/savedHists.root");
    //data_muonPtDown->Execute        ("DataAnalysis/Variations/MuonPtDownVariation/outFile.root",          "DataAnalysis/Variations/MuonPtDownVariation/savedHists.root");
    //data_electronLHMedium->Execute  ("DataAnalysis/Variations/ElectronLHMediumWPVariation/outFile.root",  "DataAnalysis/Variations/ElectronLHMediumWPVariation/savedHists.root");
    //data_muonTight->Execute         ("DataAnalysis/Variations/MuonTightWPVariation/outFile.root",         "DataAnalysis/Variations/MuonTightWPVariation/savedHists.root");
  }


/*
  data18->LoadHists   ("DataAnalysis/Nominal/savedHists.root");
  data15->LoadHists   ("DataAnalysis/Nominal/data15hi_hists.root");
  mc->LoadHists       ("MCAnalysis/Nominal/savedHists.root");
  bkg->LoadHists      ("MCAnalysis/Nominal/savedHists.root");
  //truth->LoadHists  ("TruthAnalysis/Nominal/savedHists.root");

  if (doSys) {
    data_bkgStatUpVar->CopyAnalysis (data18);
    data_bkgStatDownVar->CopyAnalysis (data18);
    bkg_statUpVar->CopyAnalysis (bkg);
    bkg_statUpVar->ConvertToStatVariation (true, 1);
    bkg_statDownVar->CopyAnalysis (bkg);
    bkg_statDownVar->ConvertToStatVariation (false, 1);
    data_bkgStatUpVar->SubtractBackground (bkg_statUpVar);
    data_bkgStatDownVar->SubtractBackground (bkg_statDownVar);
    data_bkgStatUpVar->CalculateIAA ();
    data_bkgStatUpVar->CalculateICP ();
    data_bkgStatDownVar->CalculateIAA ();
    data_bkgStatDownVar->CalculateICP ();
  }

  data18->SubtractBackground (bkg);
  data18->CalculateIAA ();
  data18->CalculateICP ();

  //data18_1_0_0->SubtractBackground (bkg_1_0_0);
  //data18_1_0_0->CalculateIAA ();
  //data18_1_0_0->CalculateICP ();

  //data15->SubtractBackground (bkg);
  //data15->CalculateIAA ();
  //data15->CalculateICP ();

  data18->CalculateZPtDistRatio (mc);
  data18->CalculateZEtaDistRatio ();
  data18->CalculateZYDistRatio (truth);
  data18->CalculateZMassSpectraRatio (mc);

  mc->SubtractBackground (bkg);
  //mc->CalculateIAA ();
  //mc->CalculateICP ();
  //truth->SubtractBackground ();

  if (doSys) {
    data_trigEff->LoadHists           ("DataAnalysis/Variations/TriggerEfficiencyCorrected/savedHists.root");
    data_trackHItight->LoadHists      ("DataAnalysis/Variations/TrackHITightWPVariation/savedHists.root");
    bkg_trackHItight->LoadHists       ("MinbiasAnalysis/Variations/TrackHITightWPVariation/savedHists.root");
    data_trackEff->LoadHists          ("DataAnalysis/Variations/TrackEffPionsVariation/savedHists.root");
    bkg_trackEff->LoadHists           ("MinbiasAnalysis/Variations/TrackEffPionsVariation/savedHists.root");
    data_trackPurity->LoadHists       ("DataAnalysis/Variations/TrackPurityVariation/savedHists.root");
    bkg_trackPurity->LoadHists        ("MinbiasAnalysis/Variations/TrackPurityVariation/savedHists.root");
    data_leptonRejVar->LoadHists      ("DataAnalysis/Variations/LeptonRejVariation/savedHists.root");
    data_electronPtUp->LoadHists      ("DataAnalysis/Variations/ElectronPtUpVariation/savedHists.root");
    data_electronPtDown->LoadHists    ("DataAnalysis/Variations/ElectronPtDownVariation/savedHists.root");
    data_muonPtUp->LoadHists          ("DataAnalysis/Variations/MuonPtUpVariation/savedHists.root");
    data_muonPtDown->LoadHists        ("DataAnalysis/Variations/MuonPtDownVariation/savedHists.root");
    data_electronLHMedium->LoadHists  ("DataAnalysis/Variations/ElectronLHMediumWPVariation/savedHists.root");
    data_muonTight->LoadHists         ("DataAnalysis/Variations/MuonTightWPVariation/savedHists.root");

    data_trigEff->SubtractBackground (bkg);
    data_trackHItight->SubtractBackground (bkg_trackHItight);
    data_trackEff->SubtractBackground (bkg_trackEff);
    data_trackPurity->SubtractBackground (bkg_trackPurity);
    data_leptonRejVar->SubtractBackground (bkg);
    data_electronPtUp->SubtractBackground (bkg);
    data_electronPtDown->SubtractBackground (bkg);
    data_muonPtUp->SubtractBackground (bkg);
    data_muonPtDown->SubtractBackground (bkg);
    data_electronLHMedium->SubtractBackground (bkg);
    data_muonTight->SubtractBackground (bkg);

    data_trigEff->CalculateZYDistRatio (truth);

    bkgSys = new Systematic (data18, "bkgSys", "Background");
    bkgSys->AddVariation (data_bkgStatUpVar, -1);
    bkgSys->AddVariation (data_bkgStatDownVar, 1);
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

    electronLHMedSys = new Systematic (data18, "electronLHMedium", "Electron ID");
    electronLHMedSys->AddVariation (data_electronLHMedium);
    electronLHMedSys->AddVariations ();

    muonTightSys = new Systematic (data18, "muonTight", "Muon ID");
    muonTightSys->AddVariation (data_muonTight);
    muonTightSys->AddVariations ();

    combSys = new Systematic (data18, "combSys", "Total");
    combSys->AddSystematic (trkSys);
    combSys->AddSystematic (trkEffSys);
    combSys->AddSystematic (trkPurSys);
    combSys->AddSystematic (bkgSys);
    combSys->AddSystematic (leptonRejSys);
    combSys->AddSystematic (electronLHMedSys);
    combSys->AddSystematic (muonTightSys);
    combSys->AddSystematic (electronPtSys);
    combSys->AddSystematic (muonPtSys);
    combSys->AddSystematics ();
  }
*/


  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}


void MakePhysicsPlots () {

  combSys->PlotTrkYieldZPt (1, 1, 2);
  data18->PlotTrkYieldZPt (1, 0, 2);

  max_iaa = 10;
  combSys->PlotSingleIAAdPtZ (1, 1);
  data18->PlotSingleIAAdPtZ (1, 0);

  combSys->PlotSingleIAAdPtZ (0, 1);
  data18->PlotSingleIAAdPtZ (0, 0);

  max_iaa = 3.2;
  combSys->PlotIAAdPtZ (1, 1);
  data18->PlotIAAdPtZ (1, 0);

  max_icp = 3.;
  combSys->PlotICPdPtZ (1, 1);
  data18->PlotICPdPtZ (1, 0);
  
}


void ComparePbPbSubYields (const short iSpc = 2, const short iPtZ = nPtZBins-1) {
  TCanvas* c = new TCanvas ("c", "", 1600, 800);
  const double dPadY = 0.4;
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

    h1 = data18->h_z_trk_zpt_sub[iSpc][iPtZ][iCent];
    h2 = data15->h_z_trk_zpt_sub[iSpc][iPtZ][iCent];
    const float min = fmin (h1->GetMinimum (0), h2->GetMinimum (0));
    const float max = fmax (h1->GetMaximum (),  h2->GetMaximum ());

    g = data18->GetTGAE (h1);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);
    g->SetFillColorAlpha (fillColors[1], 0.3);

    g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[nPtZBins-1]]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]");

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

    g = data15->GetTGAE (h2);

    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);
    g->SetFillColorAlpha (fillColors[2], 0.3);

    g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[nPtZBins-1]]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]");

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

    g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[nPtZBins-1]]);
    g->GetYaxis ()->SetRangeUser (0, 4);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("2015 / 2018");

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
  }

}

#endif
