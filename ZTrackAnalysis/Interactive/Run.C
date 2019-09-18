#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematic.h"
#include "TrackIDSystematic.h"
#include "ReweightingSystematic.h"

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
TrackIDSystematic* trkSys = nullptr;
ReweightingSystematic* trkEffSys = nullptr;
ReweightingSystematic* trkPurSys = nullptr;
Systematic* leptonRejSys = nullptr;
Systematic* electronPtSys = nullptr;
Systematic* muonPtSys = nullptr;

// variations for systematics
PhysicsAnalysis* data_trackHItight = nullptr, *data_trkIDUpVar = nullptr, *data_trkIDDownVar = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr, *bkg_trkIDUpVar = nullptr, *bkg_trkIDDownVar = nullptr;
PhysicsAnalysis* data_partComp = nullptr, *data_partCompUpVar = nullptr, *data_partCompDownVar = nullptr;
MinbiasAnalysis* bkg_partComp = nullptr, *bkg_partCompUpVar = nullptr, *bkg_partCompDownVar = nullptr;
PhysicsAnalysis* data_trackPurity = nullptr, *data_noPurUpVar = nullptr, *data_noPurDownVar = nullptr;
MinbiasAnalysis* bkg_trackPurity = nullptr, *bkg_noPurUpVar = nullptr, *bkg_noPurDownVar = nullptr;
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
    //bkg_trackHItight        = new MinbiasAnalysis ("bkg_trackHITightVar");
    //bkg_trackHItight->useHITight = true;
    data_trkIDUpVar         = new PhysicsAnalysis ("data_trkIDUpVar");
    data_trkIDUpVar->useHITight = true;
    bkg_trkIDUpVar          = new MinbiasAnalysis ("bkg_trkIDUpVar");
    bkg_trkIDUpVar->useHITight = true;
    data_trkIDDownVar       = new PhysicsAnalysis ("data_trkIDDownVar");
    data_trkIDDownVar->useHITight = true;
    bkg_trkIDDownVar        = new MinbiasAnalysis ("bkg_trkIDDownVar");
    bkg_trkIDDownVar->useHITight = true;
    

    data_partComp        = new PhysicsAnalysis ("data_trackEffVar");
    data_partComp->doTrackEffVar = true;
    //bkg_partComp         = new MinbiasAnalysis ("bkg_trackEffVar");
    //bkg_partComp->doTrackEffVar = true;
    data_partCompUpVar   = new PhysicsAnalysis ("data_partCompUpVar");
    data_partCompUpVar->doTrackEffVar = true;
    bkg_partCompUpVar    = new MinbiasAnalysis ("bkg_partCompUpVar");
    bkg_partCompUpVar->doTrackEffVar = true;
    data_partCompDownVar   = new PhysicsAnalysis ("data_partCompDownVar");
    data_partCompDownVar->doTrackEffVar = true;
    bkg_partCompDownVar    = new MinbiasAnalysis ("bkg_partCompDownVar");
    bkg_partCompDownVar->doTrackEffVar = true;


    data_trackPurity        = new PhysicsAnalysis ("data_trackPurityVar");
    data_trackPurity->doTrackPurVar = true;
    //bkg_trackPurity         = new MinbiasAnalysis ("bkg_trackPurityVar");
    //bkg_trackPurity->doTrackPurVar = true;
    data_noPurUpVar   = new PhysicsAnalysis ("data_noPurUpVar");
    data_noPurUpVar->doTrackPurVar = true;
    bkg_noPurUpVar    = new MinbiasAnalysis ("bkg_noPurUpVar");
    bkg_noPurUpVar->doTrackPurVar = true;
    data_noPurDownVar   = new PhysicsAnalysis ("data_noPurDownVar");
    data_noPurDownVar->doTrackPurVar = true;
    bkg_noPurDownVar    = new MinbiasAnalysis ("bkg_noPurDownVar");
    bkg_noPurDownVar->doTrackPurVar = true;

    data_leptonRejVar       = new PhysicsAnalysis ("data_leptonRejVar");
    data_leptonRejVar->doLeptonRejVar = true;

    data_electronPtUp       = new PhysicsAnalysis ("data_electronPtUpVar");
    data_electronPtDown     = new PhysicsAnalysis ("data_electronPtDownVar");
    bkg_electronPtUp       = new MinbiasAnalysis ("bkg");//_electronPtUpVar");
    bkg_electronPtDown     = new MinbiasAnalysis ("bkg");//_electronPtDownVar");

    data_muonPtUp           = new PhysicsAnalysis ("data_muonPtUpVar");
    data_muonPtDown         = new PhysicsAnalysis ("data_muonPtDownVar");
    bkg_muonPtUp           = new MinbiasAnalysis ("bkg");//_muonPtUpVar");
    bkg_muonPtDown         = new MinbiasAnalysis ("bkg");//_muonPtDownVar");

    data_bkgStatUpVar       = new PhysicsAnalysis ("data_bkgStatUpVar");
    data_bkgStatDownVar     = new PhysicsAnalysis ("data_bkgStatDownVar");
    bkg_statUpVar           = new MinbiasAnalysis ("bkg_statUpVar");
    bkg_statDownVar         = new MinbiasAnalysis ("bkg_statDownVar");
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
    //data_partComp->Execute          ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    //data_trackPurity->Execute       ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
  }


  data18->LoadHists   ("DataAnalysis/Nominal/data18hi_hists.root");
  data15->LoadHists   ("DataAnalysis/Nominal/data15hi_hists.root");
  bkg18->LoadHists      ("MinbiasAnalysis/Nominal/data18hi_hists.root");
  bkg15->LoadHists      ("MinbiasAnalysis/Nominal/data15hi_hists.root");
  mc->LoadHists       ("MCAnalysis/Nominal/savedHists.root");
  truth->LoadHists  ("TruthAnalysis/Nominal/savedHists.root");

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
    cout << "Initializing systematic objects. " << endl;
    bkgSys = new Systematic (data18, "bkgSys", "Bkg. Stat.");
    trkSys = new TrackIDSystematic (data18, "trkSys", "Track ID Cuts");
    trkEffSys = new ReweightingSystematic (data18, "trkEffSys", "Particle Composition");
    trkPurSys = new ReweightingSystematic (data18, "trkPurSys", "Purity Corr.");
    leptonRejSys = new Systematic (data18, "leptonRejSys", "Lepton Rejection");
    electronPtSys = new Systematic (data18, "electronPtSys", "Electron ES");
    muonPtSys = new Systematic (data18, "muonPtSys", "Muon ES");
    combSys = new Systematic (data18, "combSys", "Total");

    cout << "Calculating bkg. stat. systematic errors." << endl;
    data_bkgStatUpVar->CopyAnalysis (data18, false);
    data_bkgStatDownVar->CopyAnalysis (data18, false);
    bkg_statUpVar->CopyAnalysis (bkg18, false);
    bkg_statUpVar->ConvertToStatVariation (true, 1);
    bkg_statDownVar->CopyAnalysis (bkg18, false);
    bkg_statDownVar->ConvertToStatVariation (false, 1);
    data_bkgStatUpVar->SubtractBackground (bkg_statUpVar);
    data_bkgStatDownVar->SubtractBackground (bkg_statDownVar);
    data_bkgStatUpVar->CalculateIAA ();
    data_bkgStatUpVar->CalculateICP ();
    data_bkgStatDownVar->CalculateIAA ();
    data_bkgStatDownVar->CalculateICP ();
    bkgSys->AddVariation (data_bkgStatUpVar, false);
    bkgSys->AddVariation (data_bkgStatDownVar, false);
    bkgSys->AddVariations ();

    cout << "Calculating track ID relative systematic errors." << endl;
    data_trackHItight->LoadHists      ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    trkSys->GetRelativeVariation (data18, data_trackHItight);
    //bkg_trackHItight->LoadHists       ("MinbiasAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    //data_trackHItight->SubtractBackground (bkg_trackHItight);
    cout << "Applying track ID systematic errors." << endl;
    data_trkIDUpVar->CopyAnalysis (data18, false);
    data_trkIDDownVar->CopyAnalysis (data18, false);
    bkg_trkIDUpVar->CopyAnalysis (bkg18, false);
    bkg_trkIDDownVar->CopyAnalysis (bkg18, false);
    data_trkIDUpVar->ApplyRelativeVariation (trkSys->relVar, true); // track yields go up
    data_trkIDDownVar->ApplyRelativeVariation (trkSys->relVar, false); // track yields go down
    bkg_trkIDUpVar->ApplyRelativeVariation (trkSys->relVar, true); // track yields go up
    bkg_trkIDDownVar->ApplyRelativeVariation (trkSys->relVar, false); // track yields go down
    data_trkIDUpVar->SubtractBackground (bkg_trkIDUpVar);
    data_trkIDDownVar->SubtractBackground (bkg_trkIDDownVar);
    trkSys->AddVariation (data_trkIDUpVar, true);
    trkSys->AddVariation (data_trkIDDownVar, true);
    trkSys->AddVariations ();

    cout << "Calculating particle composition systematic errors." << endl;
    data_partComp->LoadHists          ("DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    //bkg_partComp->LoadHists           ("MinbiasAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    //data_partComp->SubtractBackground (bkg_partComp);
    trkEffSys->GetRelativeVariations (data18, data_partComp);
    data_partCompUpVar->CopyAnalysis (data18, false);
    data_partCompDownVar->CopyAnalysis (data18, false);
    bkg_partCompUpVar->CopyAnalysis (bkg18, false);
    bkg_partCompDownVar->CopyAnalysis (bkg18, false);
    trkEffSys->ApplyRelativeVariations (data_partCompUpVar, true);
    trkEffSys->ApplyRelativeVariations (bkg_partCompUpVar, true);
    trkEffSys->ApplyRelativeVariations (data_partCompDownVar, false);
    trkEffSys->ApplyRelativeVariations (bkg_partCompDownVar, false);
    data_partCompUpVar->SubtractBackground (bkg_partCompUpVar);
    data_partCompDownVar->SubtractBackground (bkg_partCompDownVar);
    trkEffSys->AddVariation (data_partCompUpVar);
    trkEffSys->AddVariation (data_partCompDownVar);
    trkEffSys->AddVariations ();


    cout << "Calculating track purity systematic errors." << endl;
    data_trackPurity->LoadHists       ("DataAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    //bkg_trackPurity->LoadHists        ("MinbiasAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    //data_trackPurity->SubtractBackground (bkg_trackPurity);
    trkPurSys->GetRelativeVariations (data18, data_trackPurity);
    data_noPurUpVar->CopyAnalysis (data18, false);
    data_noPurDownVar->CopyAnalysis (data18, false);
    bkg_noPurUpVar->CopyAnalysis (bkg18, false);
    bkg_noPurDownVar->CopyAnalysis (bkg18, false);
    trkPurSys->ApplyRelativeVariations (data_noPurUpVar, true);
    trkPurSys->ApplyRelativeVariations (bkg_noPurUpVar, true);
    trkPurSys->ApplyRelativeVariations (data_noPurDownVar, false);
    trkPurSys->ApplyRelativeVariations (bkg_noPurDownVar, false);
    data_noPurUpVar->SubtractBackground (bkg_noPurUpVar);
    data_noPurDownVar->SubtractBackground (bkg_noPurDownVar);
    trkPurSys->AddVariation (data_noPurUpVar);
    trkPurSys->AddVariation (data_noPurDownVar);
    trkPurSys->AddVariations ();


    cout << "Calculating lepton rejection systematic errors." << endl;
    data_leptonRejVar->LoadHists      ("DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");
    data_leptonRejVar->SubtractBackground (bkg18);
    leptonRejSys->AddVariation (data_leptonRejVar, -1);
    leptonRejSys->AddVariations ();

    cout << "Calculating electron ES systematic errors." << endl;
    data_electronPtUp->LoadHists      ("DataAnalysis/Variations/ElectronPtUpVariation/data18hi_hists.root");
    bkg_electronPtUp->LoadHists       ("MinbiasAnalysis/Variations/ElectronPtUpVariation/data18hi_hists.root");
    data_electronPtUp->SubtractBackground (bkg_electronPtUp);
    data_electronPtDown->LoadHists    ("DataAnalysis/Variations/ElectronPtDownVariation/data18hi_hists.root");
    bkg_electronPtDown->LoadHists     ("MinbiasAnalysis/Variations/ElectronPtDownVariation/data18hi_hists.root");
    data_electronPtDown->SubtractBackground (bkg_electronPtDown);
    electronPtSys->AddVariation (data_electronPtUp);
    electronPtSys->AddVariation (data_electronPtDown);
    electronPtSys->AddVariations ();

    cout << "Calculating muon ES systematic errors." << endl;
    data_muonPtUp->LoadHists          ("DataAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    bkg_muonPtUp->LoadHists           ("MinbiasAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    data_muonPtUp->SubtractBackground (bkg_muonPtUp);
    data_muonPtDown->LoadHists        ("DataAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");
    bkg_muonPtDown->LoadHists         ("MinbiasAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");
    data_muonPtDown->SubtractBackground (bkg_muonPtDown);
    muonPtSys->AddVariation (data_muonPtUp);
    muonPtSys->AddVariation (data_muonPtDown);
    muonPtSys->AddVariations ();

    cout << "Adding errors in quadrature." << endl;
    combSys->AddSystematic (trkSys);
    combSys->AddSystematic (trkEffSys);
    combSys->AddSystematic (trkPurSys);
    combSys->AddSystematic (bkgSys);
    combSys->AddSystematic (leptonRejSys);
    combSys->AddSystematic (electronPtSys);
    combSys->AddSystematic (muonPtSys);
    combSys->AddSystematics ();
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

  TH1D* h1 = nullptr, *h2 = nullptr, *hrat = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr, *purrat = nullptr;
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

    PhysicsAnalysis* a1, *a2 = nullptr;

    a1 = data18;
    a2 = data_trackPurity;

    //h1 = a1->h_z_trk_zpt[iSpc][iPtZ][iCent];
    //h2 = a2->h_z_trk_zpt[iSpc][iPtZ][iCent];
    h1 = (TH1D*) a1->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h1");
    h1->Add (a1->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    h2 = (TH1D*) a2->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h2");
    h2->Add (a2->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    const float min = fmin (h1->GetMinimum (0), h2->GetMinimum (0));
    const float max = fmax (h1->GetMaximum (),  h2->GetMaximum ());

    g = a1->GetTGAE (h1);
    //if (iCent == 3)
    //  g->Print ("ALL");

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
    g->GetYaxis ()->SetTitle ("N_{ch}^{total}");

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


    g = a2->GetTGAE (h2);
    //if (iCent == 3)
    //  g->Print ("ALL");

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
    g->GetYaxis ()->SetTitle ("N_{ch}^{total}");

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

    if (iCent == 3) {
      myMarkerTextNoLine (0.50, 0.9, colors[1], kFullCircle, "HILoose tracks", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.82, colors[2], kOpenCircle, "HITight tracks", 1.4, 0.04/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    hrat = (TH1D*) h2->Clone ("hrat");
    hrat->Divide (h1);
    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / h1->GetBinContent (ix)));
    //}

    a1->LoadTrackingEfficiencies (true);
    a2->LoadTrackingEfficiencies (true);
    a1->LoadTrackingPurities (true);
    a2->LoadTrackingPurities (true);

    //eff1 = (TH1D*) a1->h2_num_trk_effs[iCent]->ProjectionY ("1");
    //effrat = (TH1D*) a1->h2_den_trk_effs[iCent]->ProjectionY ("2");
    //eff1->Divide (effrat);
    //delete effrat;
    //eff2 = (TH1D*) a2->h2_num_trk_effs[iCent]->ProjectionY ("3");
    //effrat = (TH1D*) a2->h2_den_trk_effs[iCent]->ProjectionY ("4");
    //eff2->Divide (effrat);
    //delete effrat;

    //effrat = (TH1D*) eff2->Clone ("effrat");
    //effrat->Divide (eff1);
    //delete eff1, eff2;

    eff1 = (TH1D*) a1->h2_num_trk_purs[iCent]->ProjectionY ("1");
    eff1->Divide ((TH1D*) a1->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    eff2 = (TH1D*) a2->h2_num_trk_purs[iCent]->ProjectionY ("3");
    eff2->Divide ((TH1D*) a2->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    purrat = (TH1D*) eff2->Clone ("purrat");
    purrat->Divide (eff1);
    delete eff1, eff2;

    const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
    const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
    for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
      hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1));// / effrat->GetBinContent (iy+bin1-1));

    g = make_graph (hrat);
    delete hrat, effrat;

    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);
    g->SetFillColorAlpha (fillColors[2], 0.3);

    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (0.8, 1.2);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("Pions / Inclusive");

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

    TF1* inv_fit = new TF1 ("inv_fit", "[0]", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    inv_fit->SetParameter (0, 1./fit->GetParameter (0));

    inv_fit->SetLineColor (kPink-8);
    inv_fit->SetLineStyle (2);
    inv_fit->Draw ("same");

    myText (0.22, 0.3, kBlack, Form ("Avg. = %s", FormatMeasurement (fit->GetParameter (0), fit->GetParError (0), 2)), 0.02/dPadY);

    //delete h1, h2;
  }

}

#endif
