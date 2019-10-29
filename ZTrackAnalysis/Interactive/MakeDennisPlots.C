#ifndef __MakeDennisPlots_C__
#define __MakeDennisPlots_C__

#include "Params.h"
#include "PlotHybridModel.h"

#include "ArrayTemplates.h"

#include <AtlasUtils.h>

#include <TFile.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TColor.h>

#include <iostream>
#include <string>

using namespace atlashi;
using namespace std;

typedef TGraphErrors TGE;
typedef TGraphAsymmErrors TGAE;

TColor* kirchCol = new TColor ();
Color_t starColor = kirchCol->GetColor(240,197,40);//kOrange-7;//kMagenta-6;//kirchCol->GetColor(240,197,40);
Color_t cmsColor = kirchCol->GetColor (153,55,55);//kirchCol->GetColor(71,97,130);//kGreen-2;//kAzure+7;//kirchCol->GetColor(145, 103, 65);
Color_t phenixColor = kMagenta-6;//kOrange-7;
Color_t atlasColor = kAzure+2;//kRed+1;//kMagenta+1;
Color_t atlasFillColor = kAzure+8;//kRed-9;

void MakeDennisPlots () {

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Load our histograms and systematics graphs
  ////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* resultsFile = new TFile ("../rootFiles/Results/finalResults.root", "read");
  TFile* sysFile = new TFile ("../rootFiles/Systematics/CombinedSys.root", "read");

  TH1D*** h_ztrk_pt_iaa_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_ztrk_xhz_iaa_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_ztrk_pt_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_ztrk_xhz_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);

  for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      h_ztrk_pt_iaa_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_ztrk_pt_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_ztrk_xhz_iaa_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_ztrk_xhz_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    }
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      h_ztrk_pt_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_ztrk_pt_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_ztrk_xhz_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_ztrk_xhz_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    }
  }

  TGAE*** g_ztrk_pt_iaa_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE*** g_ztrk_xhz_iaa_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE*** g_ztrk_pt_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE*** g_ztrk_xhz_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);

  for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      g_ztrk_pt_iaa_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_z_trk_zpt_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_ztrk_xhz_iaa_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_z_trk_zxzh_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    }
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      g_ztrk_pt_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_z_trk_zpt_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_ztrk_xhz_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_z_trk_zxzh_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    }
  }


  TGAE** g_hybridModel_xhz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_hybridModel_pt  = Get1DArray <TGAE*> (nPtZBins);
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    g_hybridModel_xhz[iPtZ] = new TGAE ();
    g_hybridModel_pt[iPtZ] = new TGAE ();
  }
  {
    const bool useTrkPt = false;
    //const string modelFileName = (useTrkPt ? "../HybridModel/IAAs/010_IAA_pt_wake_0_ignore_neg_0.dat" : "../HybridModel/IAAs/010_IAA_z_wake_0_ignore_neg_0.dat"); // means no medium response.data"
    string modelFileName = "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_0.dat"; // means medium response including only the positive contribution from the wake
    //const string modelFileName = (useTrkPt ? "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_1.dat" : "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_1.dat"); // means full medium response, including also the negative contribution from the wake
    ifstream f;
    f.open (modelFileName.c_str ());
    float dummy = 0, x = 0, y1 = 0, y2 = 0, y = 0, yerr = 0;
    short ix = 0;
    while (f) {
      f >> x;
      for (int iPtZ = 2; iPtZ <= 4; iPtZ++) {
        f >> dummy >> dummy >> dummy >> dummy >> y1 >> y2;
        if (x*zPtBins[iPtZ+1] < 1)
          continue;
        y = 0.5*(y1+y2);
        yerr = 0.5*(y1-y2); 
        g_hybridModel_xhz[iPtZ]->SetPoint (ix, x, y);
        g_hybridModel_xhz[iPtZ]->SetPointEYhigh (ix, yerr);
        g_hybridModel_xhz[iPtZ]->SetPointEYlow (ix, yerr);
      }
      ix++;
    }
    f.close ();

    modelFileName = "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_0.dat";
    f.open (modelFileName.c_str ());
    ix = 0;
    while (f) {
      f >> x;
      for (int iPtZ = 2; iPtZ <= 4; iPtZ++) {
        f >> dummy >> dummy >> dummy >> dummy >> y1 >> y2;
        if (x > zPtBins[iPtZ])
          continue;
        y = 0.5*(y1+y2);
        yerr = 0.5*(y1-y2); 
        g_hybridModel_pt[iPtZ]->SetPoint (ix, x, y);
        g_hybridModel_pt[iPtZ]->SetPointEYhigh (ix, yerr);
        g_hybridModel_pt[iPtZ]->SetPointEYlow (ix, yerr);
      }
      ix++;
    }
    f.close ();
  }


  TGAE** g_vitevModel = Get1DArray <TGAE*> (nPtZBins);
  for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++)
    g_vitevModel[iPtZ] = new TGAE ();

  {
    ifstream f;
    double dummy = 0, x = 0, y = 0;

    TGAE* g = g_vitevModel[3];
    vector<double> xarr (0);
    vector<double> yarr_g1_8 (0), yarr_g2_0 (0), yarr_g2_2 (0);

    string modelFileName = "../VitevModel/R_SigZ0Had5020_010.g1.8LOZ30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../VitevModel/R_SigZ0Had5020_010.g2.0LOZ30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../VitevModel/R_SigZ0Had5020_010.g2.2LOZ30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();

    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }


    g = g_vitevModel[4];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();
    
    modelFileName = "../VitevModel/R_SigZ0Had5020_010.g1.8LOZ60-300";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../VitevModel/R_SigZ0Had5020_010.g2.0LOZ60-300";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../VitevModel/R_SigZ0Had5020_010.g2.2LOZ60-300";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();
    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }

  }




  TFile* dataCompFile = new TFile ("../DataComparisons/Zhadron_IAA_data_comparisons.root", "read");
  TGE* tg_PHENIX_IAA_stat = (TGE*) dataCompFile->Get ("tg_PHENIX_IAA_stat");
  TGE* tg_PHENIX_IAA_syst = (TGE*) dataCompFile->Get ("tg_PHENIX_IAA_syst");
  TGAE* tg_STAR_IAA_stat = (TGAE*) dataCompFile->Get ("tg_STAR_IAA_stat");
  TGE* tg_STAR_IAA_syst = (TGE*) dataCompFile->Get ("tg_STAR_IAA_syst");
  TGE* tg_CMS_IAA_stat = (TGE*) dataCompFile->Get ("tg_CMS_IAA_stat");
  TGE* tg_CMS_IAA_syst = (TGE*) dataCompFile->Get ("tg_CMS_IAA_syst");




  {
    TCanvas* c1 = new TCanvas ("c1", "", 1000, 1000);

    double llMargin = 0.17;
    double lrMargin = 0.024;
    double rlMargin = 0.035;
    double rrMargin = 0.02;

    double a = (double) 1./(2. + (llMargin+lrMargin)/(1.-llMargin-lrMargin) + (rlMargin+rrMargin)/(1.-rlMargin-rrMargin));
    double xPadMiddle = a * (1 + (llMargin+lrMargin)/(1.-llMargin-lrMargin));

    double yPadMiddle = 0.5;

    TPad* luPad = new TPad ("luPad", "", 0, 0.5, xPadMiddle, 1);
    TPad* ldPad = new TPad ("ldPad", "", 0, 0, xPadMiddle, 0.5);
    TPad* ruPad = new TPad ("ruPad", "", xPadMiddle, 0.5, 1, 1);
    TPad* rdPad = new TPad ("rdPad", "", xPadMiddle, 0, 1, 0.5);

    luPad->SetLeftMargin (llMargin);
    luPad->SetRightMargin (lrMargin);
    ruPad->SetLeftMargin (rlMargin);
    ruPad->SetRightMargin (rrMargin);
    ldPad->SetLeftMargin (llMargin);
    ldPad->SetRightMargin (lrMargin);
    rdPad->SetLeftMargin (rlMargin);
    rdPad->SetRightMargin (rrMargin);

    luPad->Draw ();
    ruPad->Draw ();
    ldPad->Draw ();
    rdPad->Draw ();

    short iPtZ = 3;
    short iCent = 1;
    luPad->cd ();
    luPad->SetLogx ();

    TH1D* h = new TH1D ("", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
    h->GetXaxis ()->SetRangeUser (ptTrkBins[iPtZ][0], 1.);
    h->GetYaxis ()->SetRangeUser (0, 3.1);

    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    h->GetYaxis ()->SetTitle ("#it{I}_{AA}");

    h->GetXaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelFont (43);

    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetXaxis ()->SetTitleOffset (2.5);
    h->GetYaxis ()->SetTitleOffset (2.5);
    h->GetXaxis ()->SetLabelSize (27);
    h->GetYaxis ()->SetLabelSize (27);

    h->Draw ("");

    //Color_t _colors[4] = {kBlack, kAzure-3, kViolet-3, kRed+1};
    //Color_t _fillColors[4] = {kGray, kAzure-9, kViolet-9, kRed-9};
    Color_t _colors[4] = {kBlack, kAzure-1, kGreen+2, kRed+1};
    Color_t _fillColors[4] = {kGray, kAzure-9, kGreen-7, kRed-9};

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    //for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_syst = (TGAE*) g_ztrk_pt_iaa_syst[iPtZ][iCent]->Clone ();

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetMarkerColor (_colors[iCent]);
      g_syst->SetLineColor (_colors[iCent]);
      g_syst->SetFillColorAlpha (_fillColors[iCent], 0.3);

      g_syst->Draw ( "5P");

    }

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_stat = make_graph (h_ztrk_pt_iaa_stat[iPtZ][iCent]);

      const Style_t markerStyle = markerStyles[iCent-1];
      const float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.5);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 1 + 0.03*((iCent-1) - 0.5*(numCentBins-1)), true);
      ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (2);
      g_stat->SetMarkerColor (_colors[iCent]);
      g_stat->SetLineColor (_colors[iCent]);

      g_stat->Draw ("P");
    }

    TLine* l = new TLine (1, 1, 30, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");

    iCent = numCentBins-1;
    iPtZ = 4;
    ldPad->cd ();
    ldPad->SetLogx ();

    h = new TH1D ("", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
    h->GetXaxis ()->SetRangeUser (ptTrkBins[iPtZ][0], 1.);
    h->GetYaxis ()->SetRangeUser (0, 3.8);

    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    h->GetYaxis ()->SetTitle ("#it{I}_{AA}");

    h->GetXaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelFont (43);

    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetXaxis ()->SetTitleOffset (2.5);
    h->GetYaxis ()->SetTitleOffset (2.5);
    h->GetXaxis ()->SetLabelSize (27);
    h->GetYaxis ()->SetLabelSize (27);

    h->Draw ("");

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_syst = (TGAE*) g_ztrk_pt_iaa_syst[iPtZ][iCent]->Clone ();

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetMarkerColor (_colors[iCent]);
      g_syst->SetLineColor (_colors[iCent]);
      g_syst->SetFillColorAlpha (_fillColors[iCent], 0.3);

      g_syst->Draw ( "5P");

    }

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_stat = make_graph (h_ztrk_pt_iaa_stat[iPtZ][iCent]);

      const Style_t markerStyle = markerStyles[iCent-1];
      const float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.5);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 1 + 0.03*((iCent-1) - 0.5*(numCentBins-1)), true);
      ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (2);
      g_stat->SetMarkerColor (_colors[iCent]);
      g_stat->SetLineColor (_colors[iCent]);

      g_stat->Draw ("P");
    }

    l = new TLine (1, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");

    iCent = numCentBins-1;
    iPtZ = 3;
    ruPad->cd ();
    ruPad->SetLogx ();

    h = new TH1D ("", "", nXHZBins[iPtZ], xHZBins[iPtZ]);
    h->GetXaxis ()->SetRangeUser (xHZBins[iPtZ][0], 1.);
    h->GetYaxis ()->SetRangeUser (0, 3.1);
    h->GetYaxis ()->SetLabelOffset (h->GetYaxis ()->GetLabelOffset () * 5.);

    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
    //h->GetYaxis ()->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");

    h->GetXaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelFont (43);

    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetXaxis ()->SetTitleOffset (2.5);
    h->GetYaxis ()->SetTitleOffset (2.5);
    h->GetXaxis ()->SetLabelSize (27);
    h->GetYaxis ()->SetLabelSize (0);

    h->Draw ("");

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_syst = (TGAE*) g_ztrk_xhz_iaa_syst[iPtZ][iCent]->Clone ();

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetMarkerColor (_colors[iCent]);
      g_syst->SetLineColor (_colors[iCent]);
      g_syst->SetFillColorAlpha (_fillColors[iCent], 0.3);

      g_syst->Draw ("5P");
    }

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_stat = make_graph (h_ztrk_xhz_iaa_stat[iPtZ][iCent]);

      const Style_t markerStyle = markerStyles[iCent-1];
      const float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.5);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 1 + 0.03*(iCent-1 - 0.5*(numCentBins-1)), true);
      ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (2);
      g_stat->SetMarkerColor (_colors[iCent]);
      g_stat->SetLineColor (_colors[iCent]);

      g_stat->Draw ("P");
    }

    l = new TLine (1./30., 1, 1, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");


    iCent = numCentBins-1;
    iPtZ = 4;
    rdPad->cd ();
    rdPad->SetLogx ();

    h = new TH1D ("", "", nXHZBins[iPtZ], xHZBins[iPtZ]);
    h->GetXaxis ()->SetRangeUser (xHZBins[iPtZ][0], 1.);
    h->GetYaxis ()->SetRangeUser (0, 3.8);
    h->GetYaxis ()->SetLabelOffset (h->GetYaxis ()->GetLabelOffset () * 5.);

    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
    //h->GetYaxis ()->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");

    h->GetXaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelFont (43);

    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetXaxis ()->SetTitleOffset (2.5);
    h->GetYaxis ()->SetTitleOffset (2.5);
    h->GetXaxis ()->SetLabelSize (27);
    h->GetYaxis ()->SetLabelSize (0);

    h->Draw ("");

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_syst = (TGAE*) g_ztrk_xhz_iaa_syst[iPtZ][iCent]->Clone ();

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetMarkerColor (_colors[iCent]);
      g_syst->SetLineColor (_colors[iCent]);
      g_syst->SetFillColorAlpha (_fillColors[iCent], 0.3);

      g_syst->Draw ("5P");
    }

    //for (iCent = numCentBins-1; iCent >= 1; iCent--) {
    for (iCent = 2; iCent >= 1; iCent--) {
    //for (iCent = 1; iCent >= 1; iCent--) {
      TGAE* g_stat = make_graph (h_ztrk_xhz_iaa_stat[iPtZ][iCent]);

      const Style_t markerStyle = markerStyles[iCent-1];
      const float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.5);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 1 + 0.03*(iCent-1 - 0.5*(numCentBins-1)), true);
      ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (2);
      g_stat->SetMarkerColor (_colors[iCent]);
      g_stat->SetLineColor (_colors[iCent]);

      g_stat->Draw ("P");
    }

    l = new TLine (1./60., 1, 1, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");

    luPad->cd ();
    myText (0.35, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.35, 0.78, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.050);
    myText (0.35, 0.71, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02, 1.4-1.7 nb^{-1}", 0.050);

    ldPad->cd ();
    myText (0.35, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.35, 0.78, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.050);
    myText (0.35, 0.71, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02, 1.4-1.7 nb^{-1}", 0.050);

    ruPad->cd ();
    myText (0.10, 0.85, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.06);
    rdPad->cd ();
    myText (0.10, 0.85, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.06);

    ruPad->cd ();
    myText             (0.70+0.010, 0.88-0.020, kBlack, "30-80\%", 0.062);
    myText             (0.70+0.010, 0.82-0.020, kBlack, "10-30\%", 0.062);
    myText             (0.70+0.010, 0.76-0.020, kBlack, "0-10\%", 0.062);
    myOnlyBoxText      (0.70, 0.88, 1.2, _fillColors[1], _colors[1], 0, "", 0.080, 1001, 0.40);
    myOnlyBoxText      (0.70, 0.82, 1.2, _fillColors[2], _colors[2], 0, "", 0.080, 1001, 0.40);
    myOnlyBoxText      (0.70, 0.76, 1.2, _fillColors[3], _colors[3], 0, "", 0.080, 1001, 0.40);
    myMarkerTextNoLine (0.70-0.0125, 0.88, _colors[1], markerStyles[0], "", 1.5, 0.07);
    myMarkerTextNoLine (0.70-0.0125, 0.82, _colors[2], markerStyles[1], "", 1.5, 0.07);
    myMarkerTextNoLine (0.70-0.0125, 0.76, _colors[3], markerStyles[2], "", 2.2, 0.07);

    rdPad->cd ();
    myText             (0.70+0.010, 0.88-0.020, kBlack, "30-80\%", 0.062);
    myText             (0.70+0.010, 0.82-0.020, kBlack, "10-30\%", 0.062);
    myText             (0.70+0.010, 0.76-0.020, kBlack, "0-10\%", 0.062);
    myOnlyBoxText      (0.70, 0.88, 1.2, _fillColors[1], _colors[1], 0, "", 0.080, 1001, 0.40);
    myOnlyBoxText      (0.70, 0.82, 1.2, _fillColors[2], _colors[2], 0, "", 0.080, 1001, 0.40);
    myOnlyBoxText      (0.70, 0.76, 1.2, _fillColors[3], _colors[3], 0, "", 0.080, 1001, 0.40);
    myMarkerTextNoLine (0.70-0.0125, 0.88, _colors[1], markerStyles[0], "", 1.5, 0.07);
    myMarkerTextNoLine (0.70-0.0125, 0.82, _colors[2], markerStyles[1], "", 1.5, 0.07);
    myMarkerTextNoLine (0.70-0.0125, 0.76, _colors[3], markerStyles[2], "", 2.2, 0.07);
  }
  



  //{
  //  TCanvas* c3 = new TCanvas ("c3", "", 1600, 800);
  //  c3->SetLogx ();
  //  c3->Divide (2, 1);
  //  short iPtZ = 3;
  //  const short iCent = 3;



  //  for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
  //    c3->cd (iPtZ - 3 + 1);
  //    TGAE* g_syst = (TGAE*) g_ztrk_xhz_iaa_syst[iPtZ][iCent]->Clone ();

  //    g_syst->GetXaxis ()->SetRangeUser (1, 60);
  //    g_syst->GetYaxis ()->SetRangeUser (0, 3.8);

  //    g_syst->GetXaxis ()->SetMoreLogLabels ();

  //    g_syst->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
  //    g_syst->GetYaxis ()->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ ch})");

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetMarkerColor (colors[iPtZ-2]);
  //    g_syst->SetLineColor (colors[iPtZ-2]);
  //    g_syst->SetFillColorAlpha (fillColors[iPtZ-2], 0.3);

  //    g_syst->Draw ("A5P" : "5P");

  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
  //    TGAE* g_stat = make_graph (h_ztrk_xhz_iaa_stat[iPtZ][iCent]);

  //    const Style_t markerStyle = markerStyles[iPtZ-2];
  //    const float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.5);

  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    deltaize (g_stat, 1 + 0.012*(iPtZ - 3.5), true);

  //    g_stat->GetXaxis ()->SetRangeUser (1, 60);
  //    g_stat->GetYaxis ()->SetRangeUser (0, 3.8);

  //    g_stat->GetXaxis ()->SetMoreLogLabels ();

  //    g_stat->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
  //    g_stat->GetYaxis ()->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ ch})");

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (2);
  //    g_stat->SetMarkerColor (colors[iPtZ-2]);
  //    g_stat->SetLineColor (colors[iPtZ-2]);

  //    g_stat->Draw ("P");

  //    if (iPtZ == nPtZBins-1)
  //      myMarkerTextNoLine (0.650, 0.75-0.06*(iPtZ-3), colors[iPtZ-2], markerStyle, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 1.1 * markerSize, 0.038);
  //    else
  //      myMarkerTextNoLine (0.650, 0.75-0.06*(iPtZ-3), colors[iPtZ-2], markerStyle, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 1.1 * markerSize, 0.038);
  //  }

  //  myText (0.56, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
  //  myText (0.56, 0.82, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.038);

  //  PlotHybridModel (true);
  //}




  {
    TCanvas* c4 = new TCanvas ("c4", "", 1500, 800);
    double llMargin = 0.17;
    double lrMargin = 0.025;
    double rlMargin = 0.028;
    double rrMargin = 0.02;

    double a = (double) 1./(2. + (llMargin+lrMargin)/(1.-llMargin-lrMargin) + (rlMargin+rrMargin)/(1.-rlMargin-rrMargin));
    double xPadMiddle = a * (1 + (llMargin+lrMargin)/(1.-llMargin-lrMargin));

    TPad* lPad = new TPad ("lPad", "", 0, 0, xPadMiddle, 1);
    TPad* rPad = new TPad ("rPad", "", xPadMiddle, 0, 1, 1);

    lPad->SetLeftMargin (llMargin);
    lPad->SetRightMargin (lrMargin);
    rPad->SetLeftMargin (rlMargin);
    rPad->SetRightMargin (rrMargin);

    lPad->Draw ();
    rPad->Draw ();

    TPad* pads[2] = {lPad, rPad};

    short iPtZ = 3;
    const short iCent = 3;

    const Color_t hybridColor = kRed+2;//kViolet-6;
    const Color_t vitevColor = kCyan-3;

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TH1D* h = new TH1D ("", "", nXHZBins[iPtZ], xHZBins[iPtZ]);
      h->GetXaxis ()->SetRangeUser (xHZBins[iPtZ][0], 1.);
      h->GetYaxis ()->SetRangeUser (0, 2.5);

      h->GetXaxis ()->SetMoreLogLabels ();

      h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
      h->GetYaxis ()->SetTitle ("#it{I}_{AA}");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelFont (43);

      h->GetXaxis ()->SetTitleSize (36);
      h->GetYaxis ()->SetTitleSize (36);
      h->GetXaxis ()->SetTitleOffset (1.5);
      h->GetYaxis ()->SetTitleOffset (1.5);
      h->GetXaxis ()->SetLabelSize (36);

      if (gPad == lPad)
        h->GetYaxis ()->SetLabelSize (32);
      else
        h->GetYaxis ()->SetLabelSize (0);

      h->Draw ("");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      gPad->SetLogx ();

      TGAE* g = g_hybridModel_xhz[iPtZ];

      g->SetFillColorAlpha (hybridColor, 0.8);
      g->Draw ("3");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = g_vitevModel[iPtZ];
  
      g->SetFillColorAlpha (vitevColor, 0.5);
      g->Draw ("3");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_syst = (TGAE*) g_ztrk_xhz_iaa_syst[iPtZ][iCent]->Clone ();

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (kBlack);
      g_syst->SetFillColorAlpha (kGray, 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_stat = make_graph (h_ztrk_xhz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 1 + 0.012*(iPtZ - 3.5), true);

      g_stat->SetMarkerStyle (kOpenCircle);
      g_stat->SetMarkerSize (1.5);
      g_stat->SetLineWidth (2);
      g_stat->SetMarkerColor (kBlack);
      g_stat->SetLineColor (kBlack);

      g_stat->Draw ("P");
    }

    lPad->cd ();
    myText (0.35, 0.876, kBlack, "#bf{#it{ATLAS}} Internal", 0.060);
    myText (0.35, 0.815, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.040);
    myText (0.35, 0.755, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}", 0.040);
    myText (0.22, 0.21, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.046 * 2. * (1.-xPadMiddle));
  
    TLine* l = new TLine (1./30., 1, 1, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");

    rPad->cd ();
    myText             (-0.05+0.56,  0.10+0.800-0.0130, kBlack, "ATLAS 0-10\% Pb+Pb", 0.048);
    myText             (-0.05+0.56,  0.10+0.745-0.0130, kBlack, "Hybrid Model", 0.048);
    myText             (-0.05+0.56,  0.10+0.690-0.0130, kBlack, "Li & Vitev", 0.048);
    myOnlyBoxText      (-0.05+0.550, 0.10+0.800, 1.2, kGray, kBlack, 0, "", 0.050, 1001, 0.30);
    myMarkerTextNoLine (-0.05+0.54,  0.10+0.800+0.0001, kBlack, kOpenCircle, "", 1.1 * 1.5, 0.040);
    myOnlyBoxText      (-0.05+0.550, 0.10+0.745, 1.2, hybridColor, hybridColor, 0, "", 0.050, 1001, 0.80);
    myOnlyBoxText      (-0.05+0.550, 0.10+0.690, 1.2, vitevColor, vitevColor,  0, "", 0.050, 1001, 0.50);
    myText (0.12, 0.21, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.044 * 2.*xPadMiddle);

    l = new TLine (1./60., 1, 1, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");
  }




  {
    TCanvas* c5 = new TCanvas ("c5", "", 1100,  1000);
    c5->SetLogx ();
    c5->SetLogy ();
    short iPtZ = nPtZBins-1;
    const short iCent = 3;

    TH1D* h = new TH1D ("", "", nXHZBins[iPtZ], xHZBins[iPtZ]);
    h->GetXaxis ()->SetRangeUser (xHZBins[iPtZ][0], 1.);
    h->GetYaxis ()->SetRangeUser (0.1, 5.8);
    //h->GetYaxis ()->SetRangeUser (0.0, 2.5);

    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitle ("#it{x}_{h,#gamma/Z} = #it{p}_{T}^{ ch} / #it{p}_{T}^{#gamma/Z}");
    h->GetYaxis ()->SetTitle ("#it{I}_{AA}");
    h->Draw ("");

    TGAE* g_syst = (TGAE*) g_ztrk_xhz_iaa_syst[iPtZ][iCent]->Clone ();

    //Color_t starColor = kOrange-7;//kMagenta-6;//kirchCol->GetColor(240,197,40);
    //Color_t cmsColor = kirchCol->GetColor(71,97,130);//kGreen-2;//kAzure+7;//kirchCol->GetColor(145, 103, 65);
    //Color_t phenixColor = kMagenta-6;//kOrange-7;
    //Color_t atlasColor = kAzure+2;//kRed+1;//kMagenta+1;
    //Color_t atlasFillColor = kAzure+8;//kRed-9;

    TGE* g = nullptr;

    /*
    g = tg_PHENIX_IAA_syst;
    g->SetMarkerSize (0);
    g->SetLineColor (phenixColor);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (phenixColor, 0.3);
    g->Draw ("5P");

    g = tg_PHENIX_IAA_stat;
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1.5);
    g->SetMarkerColor (phenixColor);
    g->SetLineColor (phenixColor);
    g->SetLineWidth (2);
    g->Draw ("P");

    g = tg_STAR_IAA_syst;
    g->SetMarkerSize (0);
    g->SetLineColor (starColor);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (starColor, 0.3);
    g->Draw ("5P");

    g = (TGE*) tg_STAR_IAA_stat;
    g->SetMarkerStyle (kOpenCrossX);
    g->SetMarkerSize (1.8);
    g->SetMarkerColor (starColor);
    g->SetLineColor (starColor);
    g->SetLineWidth (2);
    g->Draw ("P");
    */

    g = tg_CMS_IAA_syst;
    g->SetMarkerSize (0);
    g->SetLineColor (cmsColor);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (cmsColor, 0.3);
    g->Draw ("5P");

    g = tg_CMS_IAA_stat;
    g->SetMarkerStyle (kOpenDiamond);
    g->SetMarkerSize (2.0);
    g->SetMarkerColor (cmsColor);
    g->SetLineColor (cmsColor);
    g->SetLineWidth (2);
    g->Draw ("P");

    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (atlasColor);
    g_syst->SetLineColor (atlasColor);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (atlasFillColor, 0.3);

    g_syst->Draw ("5P");

    TGAE* g_stat = make_graph (h_ztrk_xhz_iaa_stat[iPtZ][iCent]);

    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
    //deltaize (g_stat, 1 + 0.01*(iCent - (numCentBins-1)), true);

    g_stat->SetMarkerStyle (kFullCircle);
    g_stat->SetMarkerSize (1.5);
    g_stat->SetLineWidth (2);
    g_stat->SetMarkerColor (atlasColor);
    g_stat->SetLineColor (atlasColor);

    g_stat->Draw ("P");

    myText             (-0.30+0.56,  0.10+0.80-0.012, kBlack, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 0-10% Pb+Pb", 0.04);
    myText             (-0.30+0.56,  0.10+0.75-0.012, kBlack, "CMS #it{p}_{T}^{#gamma} > 60 GeV, #it{p}_{T}^{jet} > 30, 0-10\% Pb+Pb", 0.04);
    //myText             (-0.30+0.56,  0.10+0.75-0.012, kBlack, "PHENIX 5 < #it{p}_{T}^{#gamma} < 9 GeV, 0-40\% Au+Au", 0.04);
    //myText             (-0.30+0.56,  0.10+0.70-0.012, kBlack, "STAR 12 < #it{p}_{T}^{#gamma} < 20 GeV, 0-12\% Au+Au", 0.04);
    myOnlyBoxText      (-0.30+0.550, 0.10+0.80, 1.2, atlasFillColor, atlasColor, 0, "", 0.050, 1001, 0.3);
    myOnlyBoxText      (-0.30+0.550, 0.10+0.75, 1.2, cmsColor, cmsColor,  0, "", 0.050, 1001, 0.3);
    //myOnlyBoxText      (-0.30+0.550, 0.10+0.75, 1.2, phenixColor, phenixColor, 0, "", 0.050, 1001, 0.3);
    //myOnlyBoxText      (-0.30+0.550, 0.10+0.70, 1.2, starColor, starColor, 0, "", 0.050, 1001, 0.3);
    myMarkerTextNoLine (-0.30+0.54,  0.10+0.80+0.0001, atlasColor, kFullCircle, "", 1.1 * 1.5, 0.038);
    myMarkerTextNoLine (-0.30+0.54,  0.10+0.75+0.0001, cmsColor, kOpenDiamond, "", 1.1 * 2.2, 0.038);
    //myMarkerTextNoLine (-0.30+0.54,  0.10+0.75+0.0001, phenixColor, kOpenSquare, "", 1.1 * 1.5, 0.038);
    //myMarkerTextNoLine (-0.30+0.54,  0.10+0.70+0.0001, starColor, kOpenCrossX, "", 1.1 * 1.8, 0.038);

    TLine* l = new TLine (1./60., 1, 1, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");

    myText (0.20, 0.320, kBlack, "#bf{#it{ATLAS}} Internal", 0.052);
    myText (0.20, 0.265, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.044);
    myText (0.20, 0.210, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}", 0.044);
  }




  {
    TCanvas* c6 = new TCanvas ("c6", "", 1000, 1000);

    double llMargin = 0.18;
    double lrMargin = 0.023;
    double rlMargin = 0.032;
    double rrMargin = 0.02;

    double a = (double) 1./(2. + (llMargin+lrMargin)/(1.-llMargin-lrMargin) + (rlMargin+rrMargin)/(1.-rlMargin-rrMargin));
    double xPadMiddle = a * (1 + (llMargin+lrMargin)/(1.-llMargin-lrMargin));

    double yPadMiddle = 0.5;

    TPad* luPad = new TPad ("luPad", "", 0, 0.5, xPadMiddle, 1);
    TPad* ldPad = new TPad ("ldPad", "", 0, 0, xPadMiddle, 0.5);
    TPad* ruPad = new TPad ("ruPad", "", xPadMiddle, 0.5, 1, 1);
    TPad* rdPad = new TPad ("rdPad", "", xPadMiddle, 0, 1, 0.5);

    luPad->SetLeftMargin (llMargin);
    luPad->SetRightMargin (lrMargin);
    ruPad->SetLeftMargin (rlMargin);
    ruPad->SetRightMargin (rrMargin);
    ldPad->SetLeftMargin (llMargin);
    ldPad->SetRightMargin (lrMargin);
    rdPad->SetLeftMargin (rlMargin);
    rdPad->SetRightMargin (rrMargin);

    luPad->Draw ();
    ruPad->Draw ();
    ldPad->Draw ();
    rdPad->Draw ();

    TPad* pads[4] = {luPad, ruPad, ldPad, rdPad};

    Color_t _colors[4] = {kBlack, kAzure-1, kGreen+2, kRed+1};
    Color_t _fillColors[4] = {kGray, kAzure-9, kGreen-7, kRed-9};

    short iPtZ = 3;
    short iCent = 3;

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
      h->GetXaxis ()->SetRangeUser (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]);
      h->GetYaxis ()->SetRangeUser (6e-4, 20);

      //if (gPad == ruPad)
      //  h->GetYaxis ()->SetLabelOffset (h->GetYaxis ()->GetLabelOffset () * 500.);

      h->GetXaxis ()->SetMoreLogLabels ();

      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
      h->GetYaxis ()->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelFont (43);

      h->GetXaxis ()->SetTitleSize (30);
      h->GetYaxis ()->SetTitleSize (30);
      h->GetXaxis ()->SetTitleOffset (2.5);
      h->GetYaxis ()->SetTitleOffset (2.5);
      h->GetXaxis ()->SetLabelSize (27);
      if (gPad == ruPad || gPad == rdPad)
        h->GetYaxis ()->SetLabelSize (0);
      else
        h->GetYaxis ()->SetLabelSize (27);

      h->Draw ("");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();
      //for (iCent = 0; iCent < 1; iCent++) {
      //for (iCent = 0; iCent < 2; iCent++) {
      //for (iCent = 0; iCent < 3; iCent++) {
      for (iCent = 0; iCent < numCentBins; iCent++) {
        TGAE* g_syst = (TGAE*) g_ztrk_pt_sub_syst[iPtZ][iCent]->Clone ();

        RecenterGraph (g_syst);
        for (int ix = 0; ix < g_syst->GetN (); ix++) {
          double xlo = 0, xhi = 0, x = 0, y = 0;

          g_syst->GetPoint (ix, x, y);
          xlo = fabs (x - g_syst->GetErrorXlow (ix));
          xhi = fabs (x + g_syst->GetErrorXhigh (ix));

          g_syst->SetPoint (ix, x*(1 + 0.040*(iCent-1.5)), y);
          g_syst->SetPointEXlow (ix, fabs (x*(1 + 0.040*(iCent-1.5)) - xlo));
          g_syst->SetPointEXhigh (ix, fabs (x*(1 + 0.040*(iCent-1.5)) - xhi));
        }

        g_syst->SetMarkerSize (0);
        //g_syst->SetMarkerStyle (kDot);
        g_syst->SetLineWidth (1);
        g_syst->SetLineColor (_colors[iCent]);
        g_syst->SetFillColorAlpha (_fillColors[iCent], 0.3);

        g_syst->Draw ("5P");
      }
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();
      //for (iCent = 0; iCent < 1; iCent++) {
      //for (iCent = 0; iCent < 2; iCent++) {
      //for (iCent = 0; iCent < 3; iCent++) {
      for (iCent = 0; iCent < numCentBins; iCent++) {
        TGAE* g_stat = make_graph (h_ztrk_pt_sub_stat[iPtZ][iCent]);

        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 1 + 0.040*(iCent-1.5), true);
        ResetXErrors (g_stat);

        g_stat->SetMarkerStyle (markerStyles[iCent]);
        g_stat->SetMarkerSize (markerStyles[iCent] == kOpenDiamond ? 1.2 : 1);
        g_stat->SetLineWidth (2);
        g_stat->SetMarkerColor (_colors[iCent]);
        g_stat->SetLineColor (_colors[iCent]);

        g_stat->Draw ("P");
      }
    }


    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3+2]->cd ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXHZBins[iPtZ], xHZBins[iPtZ]);
      h->GetXaxis ()->SetRangeUser (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[iPtZ]]);
      h->GetYaxis ()->SetRangeUser (5e-3, 1e3);

      if (gPad == rdPad)
        h->GetYaxis ()->SetLabelOffset (h->GetYaxis ()->GetLabelOffset () * 500.);

      h->GetXaxis ()->SetMoreLogLabels ();

      h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
      h->GetYaxis ()->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{x} d#Delta#phi)");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelFont (43);

      h->GetXaxis ()->SetTitleSize (30);
      h->GetYaxis ()->SetTitleSize (30);
      h->GetXaxis ()->SetTitleOffset (2.5);
      h->GetYaxis ()->SetTitleOffset (2.5);
      h->GetXaxis ()->SetLabelSize (27);
      if (gPad == ruPad || gPad == rdPad)
        h->GetYaxis ()->SetLabelSize (0);
      else
        h->GetYaxis ()->SetLabelSize (27);

      h->Draw ("");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3+2]->cd ();
      //for (iCent = 0; iCent < 1; iCent++) {
      //for (iCent = 0; iCent < 2; iCent++) {
      //for (iCent = 0; iCent < 3; iCent++) {
      for (iCent = 0; iCent < numCentBins; iCent++) {
        TGAE* g_syst = (TGAE*) g_ztrk_xhz_sub_syst[iPtZ][iCent]->Clone ();

        RecenterGraph (g_syst);
        for (int ix = 0; ix < g_syst->GetN (); ix++) {
          double xlo = 0, xhi = 0, x = 0, y = 0;

          g_syst->GetPoint (ix, x, y);
          xlo = fabs (x - g_syst->GetErrorXlow (ix));
          xhi = fabs (x + g_syst->GetErrorXhigh (ix));

          g_syst->SetPoint (ix, x*(1 + 0.040*(iCent-1.5)), y);
          g_syst->SetPointEXlow (ix, fabs (x*(1 + 0.040*(iCent-1.5)) - xlo));
          g_syst->SetPointEXhigh (ix, fabs (x*(1 + 0.040*(iCent-1.5)) - xhi));
        }

        //RecenterGraph (g_syst);
        //deltaize (g_syst, 1 + 0.040*(iCent-0.5*(numCentBins+1)), true);

        g_syst->SetMarkerSize (0);
        //g_syst->SetMarkerStyle (kDot);
        g_syst->SetLineWidth (1);
        g_syst->SetLineColor (_colors[iCent]);
        g_syst->SetFillColorAlpha (_fillColors[iCent], 0.3);

        g_syst->Draw ("5P");
      }
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3+2]->cd ();
      //for (iCent = 0; iCent < 1; iCent++) {
      //for (iCent = 0; iCent < 2; iCent++) {
      //for (iCent = 0; iCent < 3; iCent++) {
      for (iCent = 0; iCent < numCentBins; iCent++) {
        TGAE* g_stat = make_graph (h_ztrk_xhz_sub_stat[iPtZ][iCent]);

        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 1 + 0.040*(iCent-1.5), true);
        ResetXErrors (g_stat);

        g_stat->SetMarkerStyle (markerStyles[iCent]);
        g_stat->SetMarkerSize (1);
        g_stat->SetLineWidth (2);
        g_stat->SetMarkerColor (_colors[iCent]);
        g_stat->SetLineColor (_colors[iCent]);

        g_stat->Draw ("P");
      }
    }

    const double textSF = 0.91;

    luPad->cd ();
    myText (0.44, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);

    //ruPad->cd ();
    //myText (0.53, 0.870, kBlack, "#it{pp}, 260 pb^{-1}", 0.054);
    //myText (0.53, 0.805, kBlack, "#sqrt{s} = 5.02 TeV", 0.054);
    //myText (0.53, 0.730, kBlack, "Pb+Pb, 1.4-1.7 nb^{-1}", 0.054);
    //myText (0.53, 0.665, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.054);
    ruPad->cd ();
    myText (0.09, 0.28, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", textSF*0.050);
    myText (0.09, 0.22, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}", textSF*0.050);

    luPad->cd ();
    myText (0.60, 0.765, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.8*0.065);
    ruPad->cd ();
    myText (0.63, 0.765, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.8*0.060 * (xPadMiddle)/(1.-xPadMiddle));

    luPad->cd ();
    myText             (-0.37+0.66,         -0.080+0.500-0.018, kBlack, "#it{pp}", 0.052);
    myText             (-0.37+0.66,         -0.080+0.440-0.018, kBlack, "Pb+Pb 30-80\%", 0.052);
    myText             (-0.37+0.66,         -0.080+0.380-0.018, kBlack, "Pb+Pb 10-30\%", 0.052);
    myText             (-0.37+0.66,         -0.080+0.320-0.018, kBlack, "Pb+Pb 0-10\%", 0.052);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.500, 1.2, _fillColors[0], _colors[0], 0, "", 0.068, 1001, 0.30);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.440, 1.2, _fillColors[1], _colors[1], 0, "", 0.068, 1001, 0.30);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.380, 1.2, _fillColors[2], _colors[2], 0, "", 0.068, 1001, 0.30);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.320, 1.2, _fillColors[3], _colors[3], 0, "", 0.068, 1001, 0.30);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.500+0.0000, _colors[0], markerStyles[0], "", 1.2 * 1.0, 0.040);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.440+0.0000, _colors[1], markerStyles[1], "", 1.2 * 1.0, 0.040);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.380-0.0005, _colors[2], markerStyles[2], "", 1.2 * 1.2, 0.040);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.320+0.0000, _colors[3], markerStyles[3], "", 1.2 * 1.0, 0.040);
    //myText             (-0.000+0.60,         0.390+0.500-0.018, kBlack, "#it{pp}", 0.052);
    //myText             (-0.000+0.60,         0.390+0.440-0.018, kBlack, "Pb+Pb 30-80\%", 0.052);
    //myText             (-0.000+0.60,         0.390+0.380-0.018, kBlack, "Pb+Pb 10-30\%", 0.052);
    //myText             (-0.000+0.60,         0.390+0.320-0.018, kBlack, "Pb+Pb 0-10\%", 0.052);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.500, 1.2, _fillColors[0], _colors[0], 0, "", 0.068, 1001, 0.30);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.440, 1.2, _fillColors[1], _colors[1], 0, "", 0.068, 1001, 0.30);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.380, 1.2, _fillColors[2], _colors[2], 0, "", 0.068, 1001, 0.30);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.320, 1.2, _fillColors[3], _colors[3], 0, "", 0.068, 1001, 0.30);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.500+0.0000, _colors[0], markerStyles[0], "", 1.2 * 1.0, 0.040);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.440+0.0000, _colors[1], markerStyles[1], "", 1.2 * 1.0, 0.040);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.380-0.0005, _colors[2], markerStyles[2], "", 1.2 * 1.2, 0.040);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.320+0.0000, _colors[3], markerStyles[3], "", 1.2 * 1.0, 0.040);


    ldPad->cd ();
    myText (0.44, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);

    //rdPad->cd ();
    //myText (0.53, 0.870, kBlack, "#it{pp}, 260 pb^{-1}", 0.054);
    //myText (0.53, 0.805, kBlack, "#sqrt{s} = 5.02 TeV", 0.054);
    //myText (0.53, 0.730, kBlack, "Pb+Pb, 1.4-1.7 nb^{-1}", 0.054);
    //myText (0.53, 0.665, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.054);
    rdPad->cd ();
    myText (0.09, 0.28, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", textSF*0.050);
    myText (0.09, 0.22, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}", textSF*0.050);

    ldPad->cd ();
    myText (0.60, 0.765, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.8*0.065);
    rdPad->cd ();
    myText (0.63, 0.765, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.8*0.060 * (xPadMiddle)/(1.-xPadMiddle));

    ldPad->cd ();
    myText             (-0.37+0.66,         -0.080+0.500-0.018, kBlack, "#it{pp}", 0.052);
    myText             (-0.37+0.66,         -0.080+0.440-0.018, kBlack, "Pb+Pb 30-80\%", 0.052);
    myText             (-0.37+0.66,         -0.080+0.380-0.018, kBlack, "Pb+Pb 10-30\%", 0.052);
    myText             (-0.37+0.66,         -0.080+0.320-0.018, kBlack, "Pb+Pb 0-10\%", 0.052);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.500, 1.2, _fillColors[0], _colors[0], 0, "", 0.068, 1001, 0.30);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.440, 1.2, _fillColors[1], _colors[1], 0, "", 0.068, 1001, 0.30);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.380, 1.2, _fillColors[2], _colors[2], 0, "", 0.068, 1001, 0.30);
    myOnlyBoxText      (-0.375+0.650,        -0.080+0.320, 1.2, _fillColors[3], _colors[3], 0, "", 0.068, 1001, 0.30);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.500+0.0000, _colors[0], markerStyles[0], "", 1.2 * 1.0, 0.040);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.440+0.0000, _colors[1], markerStyles[1], "", 1.2 * 1.0, 0.040);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.380-0.0005, _colors[2], markerStyles[2], "", 1.2 * 1.2, 0.040);
    myMarkerTextNoLine (-0.378+0.64-0.0050,  -0.080+0.320+0.0000, _colors[3], markerStyles[3], "", 1.2 * 1.0, 0.040);
    //myText             (-0.000+0.60,         0.390+0.500-0.018, kBlack, "#it{pp}", 0.052);
    //myText             (-0.000+0.60,         0.390+0.440-0.018, kBlack, "Pb+Pb 30-80\%", 0.052);
    //myText             (-0.000+0.60,         0.390+0.380-0.018, kBlack, "Pb+Pb 10-30\%", 0.052);
    //myText             (-0.000+0.60,         0.390+0.320-0.018, kBlack, "Pb+Pb 0-10\%", 0.052);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.500, 1.2, _fillColors[0], _colors[0], 0, "", 0.068, 1001, 0.30);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.440, 1.2, _fillColors[1], _colors[1], 0, "", 0.068, 1001, 0.30);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.380, 1.2, _fillColors[2], _colors[2], 0, "", 0.068, 1001, 0.30);
    //myOnlyBoxText      (-0.005+0.590,        0.390+0.320, 1.2, _fillColors[3], _colors[3], 0, "", 0.068, 1001, 0.30);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.500+0.0000, _colors[0], markerStyles[0], "", 1.2 * 1.0, 0.040);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.440+0.0000, _colors[1], markerStyles[1], "", 1.2 * 1.0, 0.040);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.380-0.0005, _colors[2], markerStyles[2], "", 1.2 * 1.2, 0.040);
    //myMarkerTextNoLine (-0.007+0.58-0.0050,  0.390+0.320+0.0000, _colors[3], markerStyles[3], "", 1.2 * 1.0, 0.040);
  }





  {
    TCanvas* c7 = new TCanvas ("c7", "", 1200, 1200);
    c7->SetLogx ();
    //c7->Divide (2, 1);
    short iPtZ = 4;
    const short iCent = 3;

    //const Color_t vitevColor = kCyan+3;

    TH1D* h = new TH1D ("", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
    h->GetXaxis ()->SetRangeUser (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]);
    h->GetYaxis ()->SetRangeUser (0, 3.8);

    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetYaxis ()->SetTitle ("#it{I}_{AA}");
    h->Draw ("");

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      c7->cd (iPtZ-3+1);
      gPad->SetLogx ();

      TGAE* g = g_hybridModel_pt[iPtZ];

      g->SetFillColorAlpha (modelFillColors[iPtZ-2], 0.8);
      g->Draw ("3");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      //c7->cd (iPtZ-3+1);
      TGAE* g_syst = (TGAE*) g_ztrk_pt_iaa_syst[iPtZ][iCent]->Clone ();

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (colors[iPtZ-2]);
      g_syst->SetFillColorAlpha (fillColors[iPtZ-2], 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      //c7->cd (iPtZ-3+1);
      TGAE* g_stat = make_graph (h_ztrk_pt_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 1 + 0.030*(iPtZ - 3.5), true);

      g_stat->SetMarkerStyle (markerStyles[iPtZ-3]);
      g_stat->SetMarkerSize (1.5);
      g_stat->SetLineWidth (2);
      g_stat->SetMarkerColor (colors[iPtZ-2]);
      g_stat->SetLineColor (colors[iPtZ-2]);

      g_stat->Draw ("P");
    }

    c7->cd (1);
    myText (0.32, 0.880, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
    myText (0.32, 0.830, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.040);
    myText (0.32, 0.780, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}", 0.040);

    c7->cd (2);
    myText             (0.05-0.150+0.50,   -0.08+0.80-0.01, kBlack, "30-60", 0.036);
    myText             (0.05-0.050+0.50,   -0.08+0.80-0.01, kBlack, "60+", 0.036);
    myText             (0.05-0.050+0.58,   -0.08+0.80-0.01, kBlack, "#it{p}_{T}^{Z} [GeV]", 0.036);
    myText             (0.05-0.050+0.58,   -0.08+0.74-0.01, kBlack, "ATLAS 0-10\% Pb+Pb", 0.036);
    myText             (0.05-0.050+0.58,   -0.08+0.68-0.01, kBlack, "Hybrid Model", 0.036);
    myOnlyBoxText      (0.05-0.050+0.5566, -0.08+0.74, 1.2, fillColors[2], colors[2], 0, "", 0.060, 1001, 0.36);
    myOnlyBoxText      (0.05-0.133+0.5566, -0.08+0.74, 1.2, fillColors[1], colors[1], 0, "", 0.060, 1001, 0.46);
    myOnlyBoxText      (0.05-0.050+0.5566, -0.08+0.68, 1.2, modelFillColors[2], modelFillColors[2], 0, "", 0.060, 1001, 0.80);
    myOnlyBoxText      (0.05-0.133+0.5566, -0.08+0.68, 1.2, modelFillColors[1], modelFillColors[1], 0, "", 0.060, 1001, 0.80);
    myMarkerTextNoLine (0.05-0.050+0.54,   -0.08+0.74+0.0001, colors[2], markerStyles[1], "", 1.3 * 1.5, 0.036);
    myMarkerTextNoLine (0.05-0.133+0.54,   -0.08+0.74+0.0001, colors[1], markerStyles[0], "", 1.3 * 1.5, 0.036);

    TLine* l = new TLine (ptTrkBins[nPtZBins-1][0], 1, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]], 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    //l->SetLineColor (kPink-8);
    l->Draw ("same");
  }


}

#endif

