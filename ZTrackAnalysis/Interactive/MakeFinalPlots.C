#ifndef __MakeCONFPlots_C__
#define __MakeCONFPlots_C__

#include "Params.h"
//#include "PlotHybridModel.h"

#include "ArrayTemplates.h"

#include <AtlasUtils.h>

#include <TFile.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TColor.h>
#include <TLatex.h>

#include <iostream>
#include <string>

using namespace atlashi;
using namespace std;

typedef TGraphErrors TGE;
typedef TGraphAsymmErrors TGAE;

const Color_t starColor     = (Color_t) tcolor->GetColor (50, 175, 50); //(Color_t) tcolor->GetColor (240, 197,  40);
const Style_t starMarker    = kFullCrossX;
//const Color_t cmsColor      = (Color_t) tcolor->GetColor (153,  55,  55);
const Color_t cmsColor      = (Color_t) tcolor->GetColor ( 71,  97, 130);
const Style_t cmsMarker     = kFullDiamond;
const Color_t phenixColor   = kMagenta-6;
const Style_t phenixMarker  = kFullSquare;
const Color_t atlasColor    = finalColors[3];
const Color_t atlasFillColor = finalFillColors[3];

const Color_t hybridColor = (Color_t) tcolor->GetColor (82, 82, 173);
const double  hybridAlpha = 0.7;
const Color_t vitevColor  = (Color_t) tcolor->GetColor (39, 180, 66);
const double  vitevAlpha  = 0.6;
const Color_t jewelColor  = tcolor->GetColor (255, 170, 50);
const double  jewelAlpha  = 0.5;

const bool plotXhZ = true;
const double minModelUnc = 0.08;


void MakeTheoryBox (const double x, const double y, const Color_t color, const double colorAlpha, const double boxMultiplier = 1.) {
  const double ytsize = 0.07;
  const double xtsize = 0.18;
  const double y1 = y - 0.25*ytsize;
  const double y2 = y + 0.25*ytsize;
  const double x2 = x - 0.15*xtsize;
  const double x1 = x - 0.55*xtsize*boxMultiplier;
  TPave *mbox = new TPave (x1, y1, x2, y2, 0, "NDC");
  mbox->SetFillColorAlpha (color, colorAlpha);
  mbox->SetFillStyle (1001);
  mbox->Draw ();

  TLine mline;
  mline.SetLineWidth (1);
  mline.SetLineColor (color);
  //mline.SetLineStyle (lstyle);
  mline.SetLineStyle (0);
  Double_t y_new = (y1+y2)/2.;
  //mline.DrawLineNDC (x1, y_new, x2, y_new);
  mline.DrawLineNDC (x1, y1, x2, y1);
  mline.DrawLineNDC (x1, y2, x2, y2);
  mline.DrawLineNDC (x1, y1, x1, y2);
  mline.DrawLineNDC (x2, y1, x2, y2);
  return;
}


void MakeDataBox (const double x, const double y, const Color_t color, const double colorAlpha, const Style_t mstyle, const double msize) {
  MakeTheoryBox (x, y, color, colorAlpha);

  const double ytsize = 0.07;
  const double xtsize = 0.18;

  const double y1 = y - 0.25*ytsize;
  const double y2 = y + 0.25*ytsize;
  const double x2 = x - 0.15*xtsize;
  const double x1 = x - 0.55*xtsize;

  TLine* ml = new TLine ();
  ml->SetNDC();
  ml->SetLineColor (color);
  ml->SetLineStyle (1);
  ml->SetLineWidth (2);

  ml->DrawLineNDC (0.9*x1+0.1*x2, 0.5*(y1+y2), 0.1*x1+0.9*x2, 0.5*(y1+y2));
  ml->DrawLineNDC (0.5*(x1+x2), 0.9*y1+0.1*y2, 0.5*(x1+x2), 0.1*y1+0.9*y2);

  TMarker* marker = new TMarker (x-0.35*0.18, y, 0);
  marker->SetNDC();
  marker->SetMarkerColor (color);
  marker->SetMarkerStyle (mstyle);
  marker->SetMarkerSize (msize);
  marker->Draw ();

  if (IsFullMarker (mstyle)) {
    TMarker* marker2 = new TMarker (x-0.35*0.18, y, 0);
    marker2->SetNDC();
    marker2->SetMarkerColor (kBlack);
    marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
    marker2->SetMarkerSize (msize);
    marker2->Draw();
  }
  return;
}

void MakeFinalPlots () {

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Load our histograms and systematics graphs
  ////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* resultsFile = new TFile ("../rootFiles/Results/finalResults.root", "read");
  TFile* sysFile = new TFile ("../rootFiles/Systematics/CombinedSys.root", "read");

  TH1D*** h_trk_pt_ptz_iaa_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_trk_pt_ptz_iaa_stat_2015proj = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_trk_xhz_ptz_iaa_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_trk_xhz_ptz_iaa_stat_2015proj = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_trk_pt_ptz_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
  TH1D*** h_trk_xhz_ptz_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);

  TGAE** g_trk_avg_pt_ptz_stat = Get1DArray <TGAE*> (numCentBins);
  TGAE** g_trk_avg_xhz_ptz_stat = Get1DArray <TGAE*> (numCentBins);

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      h_trk_pt_ptz_iaa_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_pt_ptz_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_pt_ptz_iaa_2015proj_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_xhz_ptz_iaa_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_xhz_ptz_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_xhz_ptz_iaa_stat_2015proj[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_xhz_ptz_iaa_2015proj_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    }
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      h_trk_pt_ptz_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_xhz_ptz_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    }
  }
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    g_trk_avg_pt_ptz_stat[iCent] = (TGAE*) resultsFile->Get (Form ("g_trk_avg_pt_ptz_comb_iCent%i_data18", iCent));
    g_trk_avg_xhz_ptz_stat[iCent] = (TGAE*) resultsFile->Get (Form ("g_trk_avg_xhz_ptz_comb_iCent%i_data18", iCent));
  }

  TGAE*** g_trk_pt_ptz_iaa_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE*** g_trk_xhz_ptz_iaa_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE*** g_trk_pt_ptz_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE*** g_trk_xhz_ptz_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
  TGAE** g_trk_avg_pt_ptz_syst = Get1DArray <TGAE*> (numCentBins);
  TGAE** g_trk_avg_xhz_ptz_syst = Get1DArray <TGAE*> (numCentBins);

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      g_trk_pt_ptz_iaa_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_pt_ptz_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_trk_xhz_ptz_iaa_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_xhz_ptz_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    }
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      g_trk_pt_ptz_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_trk_xhz_ptz_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    }
  }
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    g_trk_avg_pt_ptz_syst[iCent] = (TGAE*) sysFile->Get (Form ("g_trk_avg_pt_ptz_comb_iCent%i_combSys", iCent));
    g_trk_avg_xhz_ptz_syst[iCent] = (TGAE*) sysFile->Get (Form ("g_trk_avg_xhz_ptz_comb_iCent%i_combSys", iCent));
  }


  TGAE** g_pythia_pth_ptz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_pythia_xhz_ptz = Get1DArray <TGAE*> (nPtZBins);
  {
    TFile* pythiaFile = new TFile ("../rootFiles/TruthAnalysis/Nominal/pythiaCompare.root", "read");
    TH1D* h_z_counts = (TH1D*) pythiaFile->Get ("h_z_counts");
    TH1D* h = nullptr, *hSig = nullptr;
    TH2D* h2 = nullptr;
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      const float nTotalEvents = h_z_counts->GetBinContent (iPtZ+1);

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_pth_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_pth_cov_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (4./pi);

      hSig = h;

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_pth_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_pth_cov_bkg_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (8./pi);

      hSig->Add (h, -1);
      g_pythia_pth_ptz[iPtZ] = make_graph (hSig);

      TGAE* g = g_pythia_pth_ptz[iPtZ];
      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (y > 0 && g->GetErrorYhigh (i) / y < 0.10) g->SetPointEYhigh (i, 0.10*y);
        if (y > 0 && g->GetErrorYlow (i) / y < 0.10) g->SetPointEYlow (i, 0.10*y);
      }

      
      h = (TH1D*) pythiaFile->Get (Form ("h_trk_xhz_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_xhz_cov_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (4./pi);

      hSig = h;

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_xhz_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_xhz_cov_bkg_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (8./pi);

      hSig->Add (h, -1);
      g_pythia_xhz_ptz[iPtZ] = make_graph (hSig);

      g = g_pythia_xhz_ptz[iPtZ];
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (y > 0 && g->GetErrorYhigh (i) / y < 0.10) g->SetPointEYhigh (i, 0.10*y);
        if (y > 0 && g->GetErrorYlow (i) / y < 0.10) g->SetPointEYlow (i, 0.10*y);
      }
    }
  }


  TGAE** g_hybrid_xhz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_hybrid_pth  = Get1DArray <TGAE*> (nPtZBins);
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    g_hybrid_xhz[iPtZ] = new TGAE ();
    g_hybrid_pth[iPtZ] = new TGAE ();
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
        g_hybrid_xhz[iPtZ]->SetPoint (ix, x, y);
        g_hybrid_xhz[iPtZ]->SetPointEYhigh (ix, yerr);
        g_hybrid_xhz[iPtZ]->SetPointEYlow (ix, yerr);
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
        g_hybrid_pth[iPtZ]->SetPoint (ix, x, y);
        g_hybrid_pth[iPtZ]->SetPointEYhigh (ix, yerr);
        g_hybrid_pth[iPtZ]->SetPointEYlow (ix, yerr);
      }
      ix++;
    }
    f.close ();
  }

  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    TGAE* g = g_hybrid_pth[iPtZ];
    double x, y;
    for (int i = 0; i < g->GetN (); i++) {
      g->GetPoint (i, x, y);
      if (y != 0 && g->GetErrorYhigh (i) < minModelUnc) g->SetPointEYhigh (i, minModelUnc);
      if (y != 0 && g->GetErrorYlow (i) < minModelUnc) g->SetPointEYlow (i, minModelUnc);
    }
    g = g_hybrid_xhz[iPtZ];
    for (int i = 0; i < g->GetN (); i++) {
      g->GetPoint (i, x, y);
      if (y != 0 && g->GetErrorYhigh (i) < minModelUnc) g->SetPointEYhigh (i, minModelUnc);
      if (y != 0 && g->GetErrorYlow (i) < minModelUnc) g->SetPointEYlow (i, minModelUnc);
    }
  } // end loop over iPtZ


  TGAE** g_vitev_xhz = Get1DArray <TGAE*> (nPtZBins);
  for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++)
    g_vitev_xhz[iPtZ] = new TGAE ();

  {
    ifstream f;
    double dummy = 0, x = 0, y = 0;

    TGAE* g = g_vitev_xhz[3];
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


    g = g_vitev_xhz[4];
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


  TGAE** g_jewel_xhz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_jewel_pth = Get1DArray <TGAE*> (nPtZBins);
  {
    TFile* jewelFile = new TFile ("../rootFiles/Jewel/hists.root", "read");
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_pth_medium_iSpc2_iPtZ%i", iPtZ));
      TH1D* h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_pth_vacuum_iSpc2_iPtZ%i", iPtZ));

      h_num->Divide (h_den);
      g_jewel_pth[iPtZ] = make_graph (h_num);

      TGAE* g = g_jewel_pth[iPtZ];
      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (y != 0 && g->GetErrorYhigh (i) < minModelUnc) g->SetPointEYhigh (i, minModelUnc);
        if (y != 0 && g->GetErrorYlow (i) < minModelUnc) g->SetPointEYlow (i, minModelUnc);
      }

      h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_medium_iSpc2_iPtZ%i", iPtZ));
      h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_vacuum_iSpc2_iPtZ%i", iPtZ));

      h_num->Divide (h_den);
      g_jewel_xhz[iPtZ] = make_graph (h_num);
      g = g_jewel_xhz[iPtZ];
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (y != 0 && g->GetErrorYhigh (i) < minModelUnc) g->SetPointEYhigh (i, minModelUnc);
        if (y != 0 && g->GetErrorYlow (i) < minModelUnc) g->SetPointEYlow (i, minModelUnc);
      }
    }
    jewelFile->Close ();
  }




  TFile* dataCompFile = new TFile ("../DataComparisons/Zhadron_IAA_data_comparisons.root", "read");
  TGE* tg_PHENIX_IAA_stat = (TGE*) dataCompFile->Get ("tg_PHENIX_IAA_stat");
  TGE* tg_PHENIX_IAA_syst = (TGE*) dataCompFile->Get ("tg_PHENIX_IAA_syst");
  TGAE* tg_STAR_IAA_stat = (TGAE*) dataCompFile->Get ("tg_STAR_IAA_stat");
  TGE* tg_STAR_IAA_syst = (TGE*) dataCompFile->Get ("tg_STAR_IAA_syst");
  TGE* tg_CMS_IAA_stat = (TGE*) dataCompFile->Get ("tg_CMS_IAA_stat");
  TGE* tg_CMS_IAA_syst = (TGE*) dataCompFile->Get ("tg_CMS_IAA_syst");

  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();

  {
    TCanvas* c1 = new TCanvas ("c1", "", 1600, plotXhZ ? 1200 : 600);

    const double llMargin = 0.17;
    const double lrMargin = 0.032;
    const double clMargin = 0.032;
    const double crMargin = 0.032;
    const double rlMargin = 0.032;
    const double rrMargin = 0.040;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    const double deltaL = (1. - llMargin - lrMargin);
    const double deltaC = (1. - clMargin - crMargin);
    const double deltaR = (1. - rlMargin - rrMargin);

    const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
    const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

    const double xPadLCMiddle = a;
    const double xPadCRMiddle = a+b;

    const double yPadMiddle = plotXhZ ? 0.5 : 0.0;

    TPad* luPad = nullptr;
    TPad* cuPad = nullptr;
    TPad* ruPad = nullptr;
    TPad* ldPad = nullptr;
    TPad* cdPad = nullptr;
    TPad* rdPad = nullptr;

    luPad = new TPad ("luPad", "", 0, yPadMiddle, xPadLCMiddle, 1);
    cuPad = new TPad ("cuPad", "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
    ruPad = new TPad ("ruPad", "", xPadCRMiddle, yPadMiddle, 1, 1);
    if (plotXhZ) {
      ldPad = new TPad ("ldPad", "", 0, 0, xPadLCMiddle, yPadMiddle);
      cdPad = new TPad ("cdPad", "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
      rdPad = new TPad ("rdPad", "", xPadCRMiddle, 0, 1, yPadMiddle);
    }

    luPad->SetLeftMargin (llMargin);
    luPad->SetRightMargin (lrMargin);
    cuPad->SetLeftMargin (clMargin);
    cuPad->SetRightMargin (crMargin);
    ruPad->SetLeftMargin (rlMargin);
    ruPad->SetRightMargin (rrMargin);
    luPad->SetBottomMargin (bMargin);
    luPad->SetTopMargin (tMargin);
    cuPad->SetBottomMargin (bMargin);
    cuPad->SetTopMargin (tMargin);
    ruPad->SetBottomMargin (bMargin);
    ruPad->SetTopMargin (tMargin);
    if (plotXhZ) {
      ldPad->SetLeftMargin (llMargin);
      ldPad->SetRightMargin (lrMargin);
      cdPad->SetLeftMargin (clMargin);
      cdPad->SetRightMargin (crMargin);
      rdPad->SetLeftMargin (rlMargin);
      rdPad->SetRightMargin (rrMargin);
      ldPad->SetBottomMargin (bMargin);
      ldPad->SetTopMargin (tMargin);
      cdPad->SetBottomMargin (bMargin);
      cdPad->SetTopMargin (tMargin);
      rdPad->SetBottomMargin (bMargin);
      rdPad->SetTopMargin (tMargin);
    }

    TPad* uPads[3] = {luPad, cuPad, ruPad};
    TPad* dPads[3] = {ldPad, cdPad, rdPad};

    for (int i = 0; i < 3; i++) {
      uPads[i]->Draw ();
      if (plotXhZ)
        dPads[i]->Draw ();
    }
    for (int i = 0; i < 3; i++) {
      uPads[i]->SetLogx ();
      uPads[i]->SetLogy ();
      if (plotXhZ) {
        dPads[i]->SetLogx ();
        dPads[i]->SetLogy ();
      }
    }


    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      uPads[iPtZ-2]->cd ();

      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      xax->SetTitleFont (43);
      xax->SetTitleSize (30);
      xax->SetTitleOffset (plotXhZ ? 2.5 : 1.25);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      const double ymin = 0.10;
      const double ymax = 10;
      yax->SetRangeUser (ymin, ymax);
      yax->SetTitleFont (43);
      yax->SetTitleSize (30);
      yax->SetTitleOffset (plotXhZ ? 2.6 : 1.30);
      yax->SetLabelFont (43);
      if (iPtZ == 2) yax->SetLabelSize (28);
      else yax->SetLabelSize (0);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (28);
      tl->SetTextAlign (21);
      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
    } // end loop over iPtZ

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      uPads[iPtZ-2]->cd ();

      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.040, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        g_syst->Draw ( "5P");

      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        ResetXErrors (g_stat);

        Style_t markerStyle = markerStyles[iCent-1];
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ


    if (plotXhZ) {
      for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
        dPads[iPtZ-2]->cd ();

        TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

        TAxis* xax = h->GetXaxis ();
        TAxis* yax = h->GetYaxis ();

        xax->SetTitle ("#it{x}_{hZ}");
        xax->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
        xax->SetTitleFont (43);
        xax->SetTitleSize (30);
        xax->SetTitleOffset (plotXhZ ? 2.5 : 1.25);
        xax->SetLabelSize (0);

        yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
        const double ymin = 0.05;
        const double ymax = 10;
        yax->SetRangeUser (ymin, ymax);
        yax->SetTitleFont (43);
        yax->SetTitleSize (30);
        yax->SetTitleOffset (plotXhZ ? 2.6 : 1.30);
        yax->SetLabelFont (43);
        if (iPtZ == 2) yax->SetLabelSize (28);
        else yax->SetLabelSize (0);

        h->SetLineWidth (0);

        h->DrawCopy ("");
        SaferDelete (&h);

        tl->SetTextFont (43);
        tl->SetTextSize (28);
        tl->SetTextAlign (21);
        const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
        if (iPtZ > 2) {
          if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
          else tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
        }
        tl->DrawLatex (1e-1,  yoff, "10^{-1}");
        tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
        tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
        tl->DrawLatex (1,     yoff, "1");

        l->SetLineStyle (2);
        l->SetLineWidth (2);
        //l->SetLineColor (kPink-8);
        l->DrawLine (xhZBins[iPtZ][0], 1, xhZBins[iPtZ][nXhZBins[iPtZ]], 1);
      } // end loop over iPtZ

      for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
        dPads[iPtZ-2]->cd ();

        for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
          TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

          RecenterGraph (g_syst);
          ResetXErrors (g_syst);
          deltaize (g_syst, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
          ResetXErrors (g_syst);
          SetConstantXErrors (g_syst, 0.040, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

          g_syst->SetMarkerSize (0);
          g_syst->SetLineWidth (1);
          g_syst->SetMarkerColor (finalColors[iCent]);
          g_syst->SetLineColor (finalColors[iCent]);
          g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

          g_syst->Draw ("5P");
        } // end loop over iCent

        for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
          TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

          Style_t markerStyle = markerStyles[iCent-1];
          float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

          RecenterGraph (g_stat);
          ResetXErrors (g_stat);
          deltaize (g_stat, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
          ResetXErrors (g_stat);

          g_stat->SetMarkerStyle (markerStyle);
          g_stat->SetMarkerSize (markerSize);
          g_stat->SetLineWidth (3);
          g_stat->SetMarkerColor (finalColors[iCent]);
          g_stat->SetLineColor (finalColors[iCent]);

          ((TGAE*) g_stat->Clone ())->Draw ("P");

          markerStyle = FullToOpenMarker (markerStyle);
          markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
          
          g_stat->SetMarkerStyle (markerStyle);
          g_stat->SetMarkerSize (markerSize);
          g_stat->SetLineWidth (0);
          g_stat->SetMarkerColor (kBlack);

          ((TGAE*) g_stat->Clone ())->Draw ("P");

          SaferDelete (&g_stat);
        } // end loop over iCent
      } // end loop over iPtZ
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    luPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    if (plotXhZ) {
      ldPad->cd ();
      tl->SetTextSize (32);
      tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
      tl->DrawLatexNDC (0.22, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
    }

    tl->SetTextSize (30);
    luPad->cd ();
    tl->DrawLatexNDC (0.62, 0.86, "#it{p}_{T}^{Z} = 15-30 GeV");
    cuPad->cd ();
    tl->DrawLatexNDC (0.56, 0.86, "#it{p}_{T}^{Z} = 30-60 GeV");
    ruPad->cd ();
    tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} = 60+ GeV");

    if (plotXhZ) {
      ldPad->cd ();
      tl->DrawLatexNDC (0.62, 0.86, "#it{p}_{T}^{Z} = 15-30 GeV");
      cdPad->cd ();
      tl->DrawLatexNDC (0.56, 0.86, "#it{p}_{T}^{Z} = 30-60 GeV");
      rdPad->cd ();
      tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} = 60+ GeV");
    }

    ruPad->cd ();
    myMarkerAndBoxAndLineText (0.60, 0.780, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "30-80\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
    myMarkerAndBoxAndLineText (0.60, 0.715, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "10-30\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
    myMarkerAndBoxAndLineText (0.60, 0.650, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));

    if (plotXhZ) {
      rdPad->cd ();
      myMarkerAndBoxAndLineText (0.60, 0.780, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "30-80\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
      myMarkerAndBoxAndLineText (0.60, 0.715, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "10-30\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
      myMarkerAndBoxAndLineText (0.60, 0.650, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
    }

    c1->SaveAs (plotXhZ ? "../Plots/FinalPlots/iaa_allptz.pdf" : "../Plots/FinalPlots/iaa_allptz_pTchOnly.pdf");
  }




  {
    TCanvas* c2 = new TCanvas ("c2", "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c2->SetLeftMargin (lMargin);
    c2->SetRightMargin (rMargin);
    c2->SetBottomMargin (0.15);
    c2->SetTopMargin (0.04);

    c2->SetLogx ();
    c2->SetLogy ();

    TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
    //xax->SetTitleFont (43);
    //xax->SetTitleSize (30);
    //xax->SetTitleOffset (1.25);
    xax->SetLabelSize (0);

    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
    const double ymin = 0.10;
    const double ymax = 1000;
    yax->SetRangeUser (ymin, ymax);
    //yax->SetTitleFont (43);
    //yax->SetTitleSize (30);
    //yax->SetTitleOffset (1.30);
    //yax->SetLabelFont (43);
    //yax->SetLabelSize (28);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    tl->SetTextFont (43);
    tl->SetTextSize (36);
    tl->SetTextAlign (21);
    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));

    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (pTchBins[nPtZBins-1][0], 1, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1);
    l->DrawLine (pTchBins[nPtZBins-1][0], 10, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 10);
    l->DrawLine (pTchBins[nPtZBins-1][0], 100, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 100);

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

        OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.040, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        g_syst->Draw ( "5P");

      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

        OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_stat);

        Style_t markerStyle = markerStyles[iCent-1];
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.19, 0.19, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->DrawLatex (13, 1.25, "#it{p}_{T}^{Z} = 15-30 GeV (#times 1)");
    tl->DrawLatex (13, 12.5, "#it{p}_{T}^{Z} = 30-60 GeV (#times 10)");
    tl->DrawLatex (13, 125, "#it{p}_{T}^{Z} = 60+ GeV (#times 10^{2})");

    myMarkerAndBoxAndLineText (0.77, 0.280, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\% #/#it{pp}", 0.036);
    myMarkerAndBoxAndLineText (0.77, 0.230, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\% #/#it{pp}", 0.036);
    myMarkerAndBoxAndLineText (0.77, 0.180, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\% #/#it{pp}", 0.036);

    c2->SaveAs ("../Plots/FinalPlots/iaa_allptz_pTchOnly_onePlot.pdf");
  }




  {
    TCanvas* c3 = new TCanvas ("c3", "", 800, 1600);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    const double yPadMiddle = 0.5;

    TPad* uPad = nullptr;
    TPad* dPad = nullptr;

    uPad = new TPad ("uPad", "", 0, yPadMiddle, 1, 1);
    dPad = new TPad ("dPad", "", 0, 0, 1, yPadMiddle);

    uPad->SetLeftMargin (lMargin);
    uPad->SetRightMargin (rMargin);
    uPad->SetBottomMargin (bMargin);
    uPad->SetTopMargin (tMargin);
    dPad->SetLeftMargin (lMargin);
    dPad->SetRightMargin (rMargin);
    dPad->SetBottomMargin (bMargin);
    dPad->SetTopMargin (tMargin);

    uPad->Draw ();
    dPad->Draw ();

    uPad->SetLogx ();
    dPad->SetLogx ();
    uPad->SetLogy ();
    dPad->SetLogy ();

    {
      dPad->cd ();

      TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (30);
      //xax->SetTitleOffset (2.5);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      const double ymin = 0.10;
      const double ymax = 1000;
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (30);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (28);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (pTchBins[nPtZBins-1][0], 1, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1);
      l->DrawLine (pTchBins[nPtZBins-1][0], 10, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 10);
      l->DrawLine (pTchBins[nPtZBins-1][0], 100, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 100);
    } // end loop over iPtZ

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      dPad->cd ();

      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

        OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.040, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        g_syst->Draw ( "5P");

      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

        OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_stat);

        Style_t markerStyle = markerStyles[iCent-1];
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ


    {
      uPad->cd ();

      TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (30);
      //xax->SetTitleOffset (2.5);
      xax->SetLabelSize (0);

      yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
      const double ymin = 2e-3;
      const double ymax = 2e3;
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (30);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (28);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
    } // end loop over iPtZ

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      uPad->cd ();

      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

        OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.04, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

        Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

        OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_stat);

        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        if (iCent > 0) g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    uPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.19, 0.19, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->SetTextAngle (-35);
    tl->DrawLatex (6, 0.15, "#it{p}_{T}^{Z} = 15-30 GeV (#times 1)");
    tl->DrawLatex (9.5, 1.65, "#it{p}_{T}^{Z} = 30-60 GeV (#times 10)");
    tl->DrawLatex (15, 18, "#it{p}_{T}^{Z} = 60+ GeV (#times 10^{2})");
    tl->SetTextAngle (0);

    myMarkerAndBoxAndLineText (0.65, 0.780, 1.4, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.6, "#it{pp}", 0.032);
    myMarkerAndBoxAndLineText (0.65, 0.730, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.780, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.730, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\%", 0.032);

    dPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.19, 0.19, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->DrawLatex (13, 1.25, "#it{p}_{T}^{Z} = 15-30 GeV (#times 1)");
    tl->DrawLatex (13, 12.5, "#it{p}_{T}^{Z} = 30-60 GeV (#times 10)");
    tl->DrawLatex (13, 125, "#it{p}_{T}^{Z} = 60+ GeV (#times 10^{2})");

    myMarkerAndBoxAndLineText (0.77, 0.280, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\% #/#it{pp}", 0.036);
    myMarkerAndBoxAndLineText (0.77, 0.230, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\% #/#it{pp}", 0.036);
    myMarkerAndBoxAndLineText (0.77, 0.180, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\% #/#it{pp}", 0.036);

    c3->SaveAs ("../Plots/FinalPlots/yield_and_iaa_allptz_pTchOnly_onePlot.pdf");
  }




  {
    TCanvas* c4a = new TCanvas ("c4a", "", 1600, 800);
    double llMargin = 0.14;
    double lrMargin = 0.032;
    double rlMargin = 0.032;
    double rrMargin = 0.04;

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

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{hZ}");
      xax->SetRangeUser (xhZBins[iPtZ][0], 1.);
      xax->SetTitleFont (43);
      xax->SetTitleSize (36);
      xax->SetTitleOffset (1.5);
      xax->SetLabelSize (0);
      //xax->SetLabelFont (43);
      //xax->SetLabelSize (36);
      //xax->SetMoreLogLabels ();

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      yax->SetRangeUser (0.05, 10);
      yax->SetTitleFont (43);
      yax->SetTitleSize (36);
      yax->SetTitleOffset (1.30);
      yax->SetLabelFont (43);
      if (gPad == lPad) yax->SetLabelSize (36);
      else              yax->SetLabelSize (0);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      const double yoff = 0.034;
      if (iPtZ > 2) {
        if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
        else tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
      }
      tl->DrawLatex (1e-1,  yoff, "10^{-1}");
      tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
      tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
      tl->DrawLatex (1,     yoff, "1");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (xhZBins[iPtZ][0], 1, xhZBins[iPtZ][nXhZBins[iPtZ]], 1);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      //SetConstantXErrors (g_syst, 0.120, true);
      SetConstantXErrors (g_syst, 0.040, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 1 + 0.012*(iPtZ - 3.5), true);

      g_stat->SetMarkerStyle (kFullCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (3);
      //g_stat->SetMarkerColor (kBlack);
      //g_stat->SetLineColor (kBlack);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      g_stat->SetMarkerStyle (kOpenCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      gPad->SetLogx ();

      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_vitev_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (vitevColor, vitevAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    lPad->cd ();
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.32, 0.876, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.32, 0.820, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.32, 0.765, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    lPad->cd ();
    tl->SetTextSize (34);
    tl->DrawLatexNDC (0.22, 0.210, "#it{p}_{T}^{Z} = 30-60 GeV");
    rPad->cd ();
    tl->DrawLatexNDC (0.12, 0.210, "#it{p}_{T}^{Z} = 60+ GeV");

    rPad->cd ();
    myText             (-0.07+0.56,  0.05+0.855-0.0130, kBlack, "#it{p}_{T}^{Z} [GeV]", 0.044);
    myText             (-0.07+0.48,  0.05+0.852-0.0130, kBlack, "60+", 0.040);
    myText             (-0.07+0.375, 0.05+0.852-0.0130, kBlack, "30-60", 0.040);
    myText             (-0.07+0.56,  0.05+0.800-0.0130, kBlack, "ATLAS 0-10\% Pb+Pb", 0.048);
    myText             (-0.07+0.56,  0.05+0.745-0.0130, kBlack, "Hybrid Model", 0.048);
    //myText             (-0.07+0.56,  0.05+0.690-0.0130, kBlack, "Li & Vitev", 0.048);
    myText             (-0.07+0.56,  0.05+0.690-0.0130, kBlack, "SCET_{G}", 0.048);
    myText             (-0.07+0.56,  0.05+0.635-0.0130, kBlack, "JEWEL", 0.048);
    //myOnlyBoxText      (-0.07+0.550, 0.05+0.800, 1.2,   kGray, kBlack, 0, "", 0.050, 1001, 0.30);
    //myMarkerTextNoLine (-0.07+0.54,  0.05+0.800+0.0001, kBlack, kOpenCircle, "", 1.1 * 1.5, 0.040);
    myOnlyBoxText      (-0.070+0.477, 0.05+0.800, 1.4,   finalFillColors[2], finalColors[2], 0, "", 0.070, 1001, 0.30);
    myMarkerTextNoLine (-0.073+0.457, 0.05+0.800+0.0001, finalColors[2], kFullCircle, "", 1.2 * 1.5, 0.040);
    myMarkerTextNoLine (-0.073+0.457, 0.05+0.800+0.0001, kBlack,     kOpenCircle, "", 1.2 * 1.5, 0.040);
    myOnlyBoxText      (-0.070+0.565, 0.05+0.800, 1.4,   finalFillColors[3], finalColors[3], 0, "", 0.070, 1001, 0.30);
    myMarkerTextNoLine (-0.073+0.545, 0.05+0.800+0.0001, finalColors[3], kFullCircle, "", 1.2 * 1.5, 0.040);
    myMarkerTextNoLine (-0.073+0.545, 0.05+0.800+0.0001, kBlack,     kOpenCircle, "", 1.2 * 1.5, 0.040);

    MakeTheoryBox (0.495, 0.795, hybridColor, hybridAlpha);
    MakeTheoryBox (0.495, 0.740, vitevColor, vitevAlpha);
    MakeTheoryBox (0.495, 0.685, jewelColor, jewelAlpha);

    c4a->SaveAs ("../Plots/FinalPlots/iaa_xhz_theoryComp.pdf");
  }




  {
    TCanvas* c4b = new TCanvas ("c4b", "", 800, 1600);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c4b_uPad", "", 0, 0.5, 1, 1);
    TPad* dPad = new TPad ("c4b_dPad", "", 0, 0, 1, 0.5);

    uPad->SetLeftMargin (lMargin);
    dPad->SetLeftMargin (lMargin);
    uPad->SetRightMargin (rMargin);
    dPad->SetRightMargin (rMargin);
    uPad->SetBottomMargin (bMargin);
    dPad->SetBottomMargin (bMargin);
    uPad->SetTopMargin (tMargin);
    dPad->SetTopMargin (tMargin);

    uPad->Draw ();
    dPad->Draw ();

    TPad* pads[2] = {uPad, dPad};

    short iPtZ = 3;
    const short iCent = 3;

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();
      gPad->SetLogx ();
      //gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{hZ}");
      xax->SetRangeUser (xhZBins[iPtZ][0], 1.);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.0);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      //yax->SetRangeUser (0.05, 10);
      double ymin, ymax;
      if (iPtZ == 3) { ymin = -0.1; ymax = 4.0; }
      if (iPtZ == 4) { ymin = -0.1; ymax = 5.0; }
      else           { ymin = 0;    ymax = 4.0; }
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (36);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (36);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      //const double yoff = 0.034;

      const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      if (iPtZ > 2) {
        tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
        if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
      }
      tl->DrawLatex (1e-1,  yoff, "10^{-1}");
      tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
      tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
      tl->DrawLatex (1,     yoff, "1");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (xhZBins[iPtZ][0], 1, xhZBins[iPtZ][nXhZBins[iPtZ]], 1);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      //SetConstantXErrors (g_syst, 0.120, true);
      SetConstantXErrors (g_syst, 0.040, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 1 + 0.012*(iPtZ - 3.5), true);

      g_stat->SetMarkerStyle (kFullCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (3);
      //g_stat->SetMarkerColor (kBlack);
      //g_stat->SetLineColor (kBlack);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      g_stat->SetMarkerStyle (kOpenCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_vitev_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (vitevColor, vitevAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    uPad->cd ();
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.24, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.24, 0.840, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.24, 0.790, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    dPad->cd ();
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.24, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.24, 0.840, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.24, 0.790, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    uPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.42, 0.730, "#it{p}_{T}^{Z} = 30-60 GeV");
    dPad->cd ();
    tl->DrawLatexNDC (0.42, 0.730, "#it{p}_{T}^{Z} = 60+ GeV");

    uPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.49, 0.685-0.0130, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.49, 0.630-0.0130, "Hybrid Model");
    //tl->DrawLatexNDC (0.49, 0.575-0.0130, "Li & Vitev");
    tl->DrawLatexNDC (0.49, 0.575-0.0130, "SCET_{G}");
    tl->DrawLatexNDC (0.49, 0.520-0.0130, "JEWEL");

    dPad->cd ();
    tl->DrawLatexNDC (0.49, 0.685-0.0130, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.49, 0.630-0.0130, "Hybrid Model");
    //tl->DrawLatexNDC (0.49, 0.575-0.0130, "Li & Vitev");
    tl->DrawLatexNDC (0.49, 0.575-0.0130, "SCET_{G}");
    tl->DrawLatexNDC (0.49, 0.520-0.0130, "JEWEL");

    uPad->cd ();
    MakeDataBox   (0.50, 0.685, finalFillColors[2], 0.30, kFullCircle, 2.3);
    MakeTheoryBox (0.50, 0.630, hybridColor, hybridAlpha);
    MakeTheoryBox (0.50, 0.575, vitevColor, vitevAlpha);
    MakeTheoryBox (0.50, 0.520, jewelColor, jewelAlpha);
    dPad->cd ();
    MakeDataBox   (0.50, 0.685, finalFillColors[3], 0.30, kFullCircle, 2.3);
    MakeTheoryBox (0.50, 0.630, hybridColor, hybridAlpha);
    MakeTheoryBox (0.50, 0.575, vitevColor, vitevAlpha);
    MakeTheoryBox (0.50, 0.520, jewelColor, jewelAlpha);

    c4b->SaveAs ("../Plots/FinalPlots/iaa_xhz_theoryComp_vertical.pdf");
  }




  {
    TCanvas* c4c = new TCanvas ("c4c", "", 800, 1600);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c4c_uPad", "", 0, 0.5, 1, 1);
    TPad* dPad = new TPad ("c4c_dPad", "", 0, 0, 1, 0.5);

    uPad->SetLeftMargin (lMargin);
    dPad->SetLeftMargin (lMargin);
    uPad->SetRightMargin (rMargin);
    dPad->SetRightMargin (rMargin);
    uPad->SetBottomMargin (bMargin);
    dPad->SetBottomMargin (bMargin);
    uPad->SetTopMargin (tMargin);
    dPad->SetTopMargin (tMargin);

    uPad->Draw ();
    dPad->Draw ();

    TPad* pads[2] = {uPad, dPad};

    short iPtZ = 3;
    const short iCent = 3;

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();
      gPad->SetLogx ();
      //gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.2);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      //yax->SetRangeUser (0.05, 10);
      double ymin, ymax;
      if (iPtZ == 3) { ymin = 0; ymax = 4.0; }
      if (iPtZ == 4) { ymin = 0; ymax = 5.0; }
      else           { ymin = 0; ymax = 4.0; }
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (36);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (36);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      //const double yoff = 0.034;

      const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (kFullCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (3);
      //g_stat->SetMarkerColor (kBlack);
      //g_stat->SetLineColor (kBlack);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      g_stat->SetMarkerStyle (kOpenCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    //for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
    //  pads[iPtZ-3]->cd ();

    //  TGAE* g = (TGAE*) g_vitev_pth[iPtZ]->Clone ();
  
    //  g->SetFillColorAlpha (vitevColor, vitevAlpha);
    //  ((TGAE*) g->Clone ())->Draw ("3");
    //  SaferDelete (&g);
    //}

    for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
      pads[iPtZ-3]->cd ();

      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    uPad->cd ();
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.28, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.840, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.30, 0.790, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    dPad->cd ();
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.28, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.840, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.30, 0.790, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    uPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.42, 0.730, "#it{p}_{T}^{Z} = 30-60 GeV");
    dPad->cd ();
    tl->DrawLatexNDC (0.42, 0.730, "#it{p}_{T}^{Z} = 60+ GeV");

    uPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.49, 0.685-0.0130, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.49, 0.630-0.0130, "Hybrid Model");
    //tl->DrawLatexNDC (0.49, 0.575-0.0130, "Li & Vitev");
    tl->DrawLatexNDC (0.49, 0.575-0.0130, "SCET_{G}");
    tl->DrawLatexNDC (0.49, 0.520-0.0130, "JEWEL");

    dPad->cd ();
    tl->DrawLatexNDC (0.49, 0.685-0.0130, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.49, 0.630-0.0130, "Hybrid Model");
    //tl->DrawLatexNDC (0.49, 0.575-0.0130, "Li & Vitev");
    tl->DrawLatexNDC (0.49, 0.575-0.0130, "SCET_{G}");
    tl->DrawLatexNDC (0.49, 0.520-0.0130, "JEWEL");

    uPad->cd ();
    MakeDataBox   (0.50, 0.685, finalFillColors[2], 0.30, kFullCircle, 2.3);
    MakeTheoryBox (0.50, 0.630, hybridColor, hybridAlpha);
    MakeTheoryBox (0.50, 0.575, vitevColor, vitevAlpha);
    MakeTheoryBox (0.50, 0.520, jewelColor, jewelAlpha);
    dPad->cd ();
    MakeDataBox   (0.50, 0.685, finalFillColors[3], 0.30, kFullCircle, 2.3);
    MakeTheoryBox (0.50, 0.630, hybridColor, hybridAlpha);
    MakeTheoryBox (0.50, 0.575, vitevColor, vitevAlpha);
    MakeTheoryBox (0.50, 0.520, jewelColor, jewelAlpha);

    c4c->SaveAs ("../Plots/FinalPlots/iaa_pTch_theoryComp_vertical.pdf");
  }




  {
    TCanvas* c4d = new TCanvas ("c4d", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c4d->SetLeftMargin (lMargin);
    c4d->SetRightMargin (rMargin);
    c4d->SetBottomMargin (bMargin);
    c4d->SetTopMargin (tMargin);

    const short iPtZ = 4;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      //gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.2);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      //yax->SetRangeUser (0.05, 10);
      double ymin, ymax;
      if (iPtZ == 3) { ymin = 0; ymax = 4.0; }
      if (iPtZ == 4) { ymin = 0; ymax = 5.0; }
      else           { ymin = 0; ymax = 4.0; }
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (36);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (36);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      //const double yoff = 0.034;

      const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (kFullCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (3);
      //g_stat->SetMarkerColor (kBlack);
      //g_stat->SetLineColor (kBlack);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      g_stat->SetMarkerStyle (kOpenCircle);
      g_stat->SetMarkerSize (2.3);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }
    //{
    //  TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent]);

    //  RecenterGraph (g_stat);
    //  ResetXErrors (g_stat);
    //  deltaize (g_stat, 1.05, true);
    //  ResetXErrors (g_stat);

    //  g_stat->SetMarkerStyle (kFullCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (3);
    //  //g_stat->SetMarkerColor (kBlack);
    //  //g_stat->SetLineColor (kBlack);
    //  g_stat->SetMarkerColor (finalColors[2]);
    //  g_stat->SetLineColor (finalColors[2]);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  g_stat->SetMarkerStyle (kOpenCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (0);
    //  g_stat->SetMarkerColor (kBlack);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  SaferDelete (&g_stat);
    //}

    {
      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    //{
    //  TGAE* g = (TGAE*) g_vitev_pth[iPtZ]->Clone ();
  
    //  g->SetFillColorAlpha (vitevColor, vitevAlpha);
    //  ((TGAE*) g->Clone ())->Draw ("3");
    //  SaferDelete (&g);
    //}

    {
      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.28, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.840, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.30, 0.790, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.42, 0.730, "#it{p}_{T}^{Z} = 60+ GeV");

    tl->SetTextSize (32);
    //tl->DrawLatexNDC (0.49, 0.685-0.0130, "2018 Pb+Pb 0-10\%");
    //tl->DrawLatexNDC (0.49, 0.630-0.0130, "'18+'15 Pb+Pb (#it{proj.})");
    //tl->DrawLatexNDC (0.49, 0.575-0.0130, "Hybrid Model");
    //tl->DrawLatexNDC (0.49, 0.520-0.0130, "JEWEL");

    //MakeDataBox   (0.50, 0.685, finalFillColors[3], 0.30, kFullCircle, 2.3);
    //MakeDataBox   (0.50, 0.630, finalFillColors[2], 0.30, kFullCircle, 2.3);
    //MakeTheoryBox (0.50, 0.575, hybridColor, hybridAlpha);
    //MakeTheoryBox (0.50, 0.520, jewelColor, jewelAlpha);

    tl->DrawLatexNDC (0.49, 0.685-0.0130, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.49, 0.630-0.0130, "Hybrid Model");
    //tl->DrawLatexNDC (0.49, 0.575-0.0130, "Li & Vitev");
    tl->DrawLatexNDC (0.49, 0.575-0.0130, "SCET_{G}");
    tl->DrawLatexNDC (0.49, 0.520-0.0130, "JEWEL");

    MakeDataBox   (0.50, 0.685, finalFillColors[3], 0.30, kFullCircle, 2.3);
    MakeTheoryBox (0.50, 0.630, hybridColor, hybridAlpha);
    MakeTheoryBox (0.50, 0.575, vitevColor, vitevAlpha);
    MakeTheoryBox (0.50, 0.520, jewelColor, jewelAlpha);

    c4d->SaveAs ("../Plots/FinalPlots/iaa_pTch_theoryComp_iPtZ4.pdf");
  }




  {
    TCanvas* c5a = new TCanvas ("c5a", "", 1100,  1000);
    c5a->SetLogx ();
    c5a->SetLogy ();
    short iPtZ = nPtZBins-1;
    const short iCent = 3;

    TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{x}_{h,#gamma/Z} = #it{p}_{T}^{ch} / #it{p}_{T}^{#gamma/Z}");
    xax->SetRangeUser (xhZBins[iPtZ][0], 1.);
    xax->SetMoreLogLabels ();

    yax->SetTitle ("#it{I}_{AA}");
    yax->SetRangeUser (0.05, 10);
    //yax->SetRangeUser (0.0, 2.5);
    yax->SetMoreLogLabels ();

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (1./60., 1, 1, 1);

    TGE* g = nullptr;

    g = (TGE*) tg_CMS_IAA_syst->Clone ();
    g->SetMarkerSize (0);
    g->SetLineColor (cmsColor);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (cmsColor, 0.3);
    ((TGE*) g->Clone ())->Draw ("5P");

    SaferDelete (&g);

    g = (TGE*) tg_CMS_IAA_stat->Clone ();
    g->SetMarkerStyle (cmsMarker);
    g->SetMarkerSize (2.5);
    g->SetMarkerColor (cmsColor);
    g->SetLineColor (cmsColor);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    if (IsFullMarker (cmsMarker)) {
      g->SetMarkerStyle (FullToOpenMarker (cmsMarker));
      g->SetMarkerSize (2.5);
      g->SetLineWidth (0);
      g->SetMarkerColor (kBlack);
      ((TGE*) g->Clone ())->Draw ("P");
    }

    SaferDelete (&g);

    TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();
    RecenterGraph (g_syst);
    ResetXErrors (g_syst);
    //SetConstantXErrors (g_syst, 0.060, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (atlasColor);
    g_syst->SetLineColor (atlasColor);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (atlasFillColor, 0.3);

    g_syst->Draw ("5P");

    TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
    //deltaize (g_stat, 1 + 0.01*(iCent - (numCentBins-1)), true);

    g_stat->SetMarkerStyle (kFullCircle);
    g_stat->SetMarkerSize (2.0);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (atlasColor);
    g_stat->SetLineColor (atlasColor);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    g_stat->SetMarkerStyle (kOpenCircle);
    g_stat->SetMarkerSize (2.0);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    SaferDelete (&g_stat);

    myMarkerAndBoxAndLineText (0.31, 0.90-0.012, 1.8, 1001, atlasFillColor, 0.30, atlasColor, kFullCircle,  2.0, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 0-10% Pb+Pb", 0.04);
    myMarkerAndBoxAndLineText (0.31, 0.85-0.012, 1.8, 1001, cmsColor,       0.30, cmsColor,   cmsMarker, 2.5, "CMS #it{p}_{T}^{#gamma} > 60 GeV, #it{p}_{T}^{jet} > 30, 0-10\% Pb+Pb", 0.04);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (40);
    tl->DrawLatexNDC (0.24, 0.310, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.24, 0.260, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.24, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    c5a->SaveAs ("../Plots/FinalPlots/iaa_xhz_cmsComp.pdf");
  }




  {
    TCanvas* c5b = new TCanvas ("c5b", "", 1100,  1000);
    c5b->SetLogx ();
    c5b->SetLogy ();
    short iPtZ = nPtZBins-1;
    const short iCent = 3;

    TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{x}_{h,#gamma/Z} = #it{p}_{T}^{ch} / #it{p}_{T}^{#gamma/Z}");
    xax->SetRangeUser (xhZBins[iPtZ][0], 1.);
    xax->SetMoreLogLabels ();

    yax->SetTitle ("#it{I}_{AA}");
    yax->SetRangeUser (0.05, 10);
    //yax->SetRangeUser (0.0, 2.5);
    yax->SetMoreLogLabels ();

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (1./60., 1, 1, 1);

    TGE* g = nullptr;

    g = (TGE*) tg_PHENIX_IAA_syst->Clone ();
    g->SetMarkerSize (0);
    g->SetMarkerColor (phenixColor);
    g->SetLineColor (phenixColor);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (phenixColor, 0.3);
    ((TGE*) g->Clone ())->Draw ("5P");

    SaferDelete (&g);

    g = (TGE*) tg_PHENIX_IAA_stat->Clone ();
    g->SetMarkerStyle (phenixMarker);
    g->SetMarkerSize (2.0);
    g->SetMarkerColor (phenixColor);
    g->SetLineColor (phenixColor);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    if (IsFullMarker (phenixMarker)) {
      g->SetMarkerStyle (FullToOpenMarker (phenixMarker));
      g->SetMarkerSize (2.0);
      g->SetLineWidth (0);
      g->SetMarkerColor (kBlack);
      ((TGE*) g->Clone ())->Draw ("P");
    }

    SaferDelete (&g);

    g = (TGE*) tg_STAR_IAA_syst->Clone ();
    g->SetMarkerSize (0);
    g->SetLineColor (starColor);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (starColor, 0.2);
    ((TGE*) g->Clone ())->Draw ("5P");

    SaferDelete (&g);

    g = (TGE*) tg_STAR_IAA_stat->Clone ();
    g->SetMarkerStyle (starMarker);
    g->SetMarkerSize (2.4);
    g->SetMarkerColor (starColor);
    g->SetLineColor (starColor);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    if (IsFullMarker (starMarker)) {
      g->SetMarkerStyle (FullToOpenMarker (starMarker));
      g->SetMarkerSize (2.4);
      g->SetLineWidth (0);
      g->SetMarkerColor (kBlack);
      ((TGE*) g->Clone ())->Draw ("P");
    }

    SaferDelete (&g);

    TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();
    RecenterGraph (g_syst);
    ResetXErrors (g_syst);
    //SetConstantXErrors (g_syst, 0.040, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (atlasColor);
    g_syst->SetLineColor (atlasColor);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (atlasFillColor, 0.3);

    g_syst->Draw ("5P");

    TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
    //deltaize (g_stat, 1 + 0.01*(iCent - (numCentBins-1)), true);

    g_stat->SetMarkerStyle (kFullCircle);
    g_stat->SetMarkerSize (2.0);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (atlasColor);
    g_stat->SetLineColor (atlasColor);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    g_stat->SetMarkerStyle (kOpenCircle);
    g_stat->SetMarkerSize (2.0);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    SaferDelete (&g_stat);

    myMarkerAndBoxAndLineText (0.31, 0.90-0.012, 1.8, 1001, atlasFillColor, 0.30, atlasColor,   kFullCircle, 2.0, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 0-10% Pb+Pb", 0.04);
    myMarkerAndBoxAndLineText (0.31, 0.85-0.012, 1.8, 1001, phenixColor,    0.30, phenixColor,  phenixMarker, 2.0, "PHENIX 5 < #it{p}_{T}^{#gamma} < 9 GeV, 0-40\% Au+Au", 0.04);
    myMarkerAndBoxAndLineText (0.31, 0.80-0.012, 1.8, 1001, starColor,      0.20, starColor,    starMarker, 2.5, "STAR 12 < #it{p}_{T}^{#gamma} < 20 GeV, 0-12\% Au+Au", 0.04);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (40);
    tl->DrawLatexNDC (0.24, 0.310, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.24, 0.260, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.24, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    c5b->SaveAs ("../Plots/FinalPlots/iaa_xhz_rhicComp.pdf");
  }




  {
    TCanvas* c6 = new TCanvas ("c6", "", 1600, plotXhZ ? 1200 : 600);

    const double llMargin = 0.17;
    const double lrMargin = 0.032;
    const double clMargin = 0.032;
    const double crMargin = 0.032;
    const double rlMargin = 0.032;
    const double rrMargin = 0.040;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    const double deltaL = (1. - llMargin - lrMargin);
    const double deltaC = (1. - clMargin - crMargin);
    const double deltaR = (1. - rlMargin - rrMargin);

    const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
    const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

    const double xPadLCMiddle = a;
    const double xPadCRMiddle = a+b;

    const double yPadMiddle = plotXhZ ? 0.5 : 0.0;

    TPad* luPad = nullptr;
    TPad* cuPad = nullptr;
    TPad* ruPad = nullptr;
    TPad* ldPad = nullptr;
    TPad* cdPad = nullptr;
    TPad* rdPad = nullptr;

    luPad = new TPad ("luPad", "", 0, yPadMiddle, xPadLCMiddle, 1);
    cuPad = new TPad ("cuPad", "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
    ruPad = new TPad ("ruPad", "", xPadCRMiddle, yPadMiddle, 1, 1);
    if (plotXhZ) {
      ldPad = new TPad ("ldPad", "", 0, 0, xPadLCMiddle, yPadMiddle);
      cdPad = new TPad ("cdPad", "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
      rdPad = new TPad ("rdPad", "", xPadCRMiddle, 0, 1, yPadMiddle);
    }

    luPad->SetLeftMargin (llMargin);
    luPad->SetRightMargin (lrMargin);
    cuPad->SetLeftMargin (clMargin);
    cuPad->SetRightMargin (crMargin);
    ruPad->SetLeftMargin (rlMargin);
    ruPad->SetRightMargin (rrMargin);
    luPad->SetBottomMargin (bMargin);
    luPad->SetTopMargin (tMargin);
    cuPad->SetBottomMargin (bMargin);
    cuPad->SetTopMargin (tMargin);
    ruPad->SetBottomMargin (bMargin);
    ruPad->SetTopMargin (tMargin);
    if (plotXhZ) {
      ldPad->SetLeftMargin (llMargin);
      ldPad->SetRightMargin (lrMargin);
      cdPad->SetLeftMargin (clMargin);
      cdPad->SetRightMargin (crMargin);
      rdPad->SetLeftMargin (rlMargin);
      rdPad->SetRightMargin (rrMargin);
      ldPad->SetBottomMargin (bMargin);
      ldPad->SetTopMargin (tMargin);
      cdPad->SetBottomMargin (bMargin);
      cdPad->SetTopMargin (tMargin);
      rdPad->SetBottomMargin (bMargin);
      rdPad->SetTopMargin (tMargin);
    }

    TPad* uPads[3] = {luPad, cuPad, ruPad};
    TPad* dPads[3] = {ldPad, cdPad, rdPad};

    for (int i = 0; i < 3; i++) {
      uPads[i]->Draw ();
      if (plotXhZ)
        dPads[i]->Draw ();
    }
    for (int i = 0; i < 3; i++) {
      uPads[i]->SetLogx ();
      uPads[i]->SetLogy ();
      if (plotXhZ) {
        dPads[i]->SetLogx ();
        dPads[i]->SetLogy ();
      }
    }


    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      uPads[iPtZ-2]->cd ();

      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      xax->SetTitleFont (43);
      xax->SetTitleSize (30);
      xax->SetTitleOffset (plotXhZ ? 2.5 : 1.25);
      xax->SetLabelSize (0);

      yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
      yax->SetRangeUser (2e-3, 200);
      yax->SetTitleFont (43);
      yax->SetTitleSize (30);
      yax->SetTitleOffset (plotXhZ ? 2.6 : 1.30);
      yax->SetLabelFont (43);
      if (iPtZ == 2) yax->SetLabelSize (28);
      else yax->SetLabelSize (0);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (28);
      tl->SetTextAlign (21);
      const double yoff = 0.0009;
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
    } // end loop over iPtZ

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      uPads[iPtZ-2]->cd ();

      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.04, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

        Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        ResetXErrors (g_stat);

        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        if (iCent > 0) g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ


    if (plotXhZ) {
      for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
        dPads[iPtZ-2]->cd ();

        TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

        TAxis* xax = h->GetXaxis ();
        TAxis* yax = h->GetYaxis ();

        xax->SetTitle ("#it{x}_{hZ}");
        xax->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
        xax->SetTitleFont (43);
        xax->SetTitleSize (30);
        xax->SetTitleOffset (plotXhZ ? 2.5 : 1.25);
        xax->SetLabelSize (0);

        yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
        yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{x} d#Delta#phi)");
        yax->SetRangeUser (5e-3, 4e3);
        yax->SetTitleFont (43);
        yax->SetTitleSize (30);
        yax->SetTitleOffset (plotXhZ ? 2.6 : 1.30);
        yax->SetLabelFont (43);
        if (iPtZ == 2) yax->SetLabelSize (28);
        else yax->SetLabelSize (0);

        h->SetLineWidth (0);

        h->DrawCopy ("");
        SaferDelete (&h);

        tl->SetTextFont (43);
        tl->SetTextSize (28);
        tl->SetTextAlign (21);
        const double yoff = 1.95e-3;
        if (iPtZ > 2) {
          if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
          else tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
        }
        tl->DrawLatex (1e-1,  yoff, "10^{-1}");
        tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
        tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
        tl->DrawLatex (1,     yoff, "1");
      } // end loop over iPtZ

      for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
        dPads[iPtZ-2]->cd ();

        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_sub_syst[iPtZ][iCent]->Clone ();

          RecenterGraph (g_syst);
          ResetXErrors (g_syst);
          deltaize (g_syst, 0.05*(iCent - 0.5*(numCentBins-1)), true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
          ResetXErrors (g_syst);
          SetConstantXErrors (g_syst, 0.04, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

          g_syst->SetMarkerSize (0);
          g_syst->SetLineWidth (1);
          g_syst->SetMarkerColor (finalColors[iCent]);
          g_syst->SetLineColor (finalColors[iCent]);
          g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

          ((TGAE*) g_syst->Clone ())->Draw ("5P");

          SaferDelete (&g_syst);
        } // end loop over iCent

        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TGAE* g_stat = make_graph (h_trk_xhz_ptz_sub_stat[iPtZ][iCent]);

          Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
          float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

          RecenterGraph (g_stat);
          ResetXErrors (g_stat);
          deltaize (g_stat, 0.05*(iCent - 0.5*(numCentBins-1)), true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
          ResetXErrors (g_stat);

          g_stat->SetMarkerStyle (markerStyle);
          g_stat->SetMarkerSize (markerSize);
          g_stat->SetLineWidth (3);
          g_stat->SetMarkerColor (finalColors[iCent]);
          g_stat->SetLineColor (finalColors[iCent]);

          ((TGAE*) g_stat->Clone ())->Draw ("P");

          markerStyle = FullToOpenMarker (markerStyle);
          markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
          
          g_stat->SetMarkerStyle (markerStyle);
          g_stat->SetMarkerSize (markerSize);
          if (iCent > 0) g_stat->SetLineWidth (0);
          g_stat->SetMarkerColor (kBlack);

          ((TGAE*) g_stat->Clone ())->Draw ("P");

          SaferDelete (&g_stat);
        } // end loop over iCent
      } // end loop over iPtZ
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    luPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    if (plotXhZ) {
      ldPad->cd ();
      tl->SetTextSize (32);
      tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
      tl->DrawLatexNDC (0.22, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
    }

    tl->SetTextSize (30);
    luPad->cd ();
    tl->DrawLatexNDC (0.62, 0.86, "#it{p}_{T}^{Z} = 15-30 GeV");
    cuPad->cd ();
    tl->DrawLatexNDC (0.56, 0.86, "#it{p}_{T}^{Z} = 30-60 GeV");
    ruPad->cd ();
    tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} = 60+ GeV");

    if (plotXhZ) {
      ldPad->cd ();
      tl->DrawLatexNDC (0.62, 0.86, "#it{p}_{T}^{Z} = 15-30 GeV");
      cdPad->cd ();
      tl->DrawLatexNDC (0.56, 0.86, "#it{p}_{T}^{Z} = 30-60 GeV");
      rdPad->cd ();
      tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} = 60+ GeV");
    }

    cuPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.016 / (gPad->GetWNDC ()));
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Pb+Pb 30-80\%", 0.016 / (gPad->GetWNDC ()));
    ruPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "Pb+Pb 10-30\%", 0.016 / (gPad->GetWNDC ()));
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Pb+Pb 0-10\%", 0.016 / (gPad->GetWNDC ()));

    if (plotXhZ) {
      cdPad->cd ();
      myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.016 / (gPad->GetWNDC ()));
      myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Pb+Pb 30-80\%", 0.016 / (gPad->GetWNDC ()));
      rdPad->cd ();
      myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "Pb+Pb 10-30\%", 0.016 / (gPad->GetWNDC ()));
      myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Pb+Pb 0-10\%", 0.016 / (gPad->GetWNDC ()));
    }

    c6->SaveAs (plotXhZ ? "../Plots/FinalPlots/yield_allptz.pdf" : "../Plots/FinalPlots/yield_allptz_pTchOnly.pdf");
  }



  {
    TCanvas* c6b = new TCanvas ("c6b", "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c6b->SetLeftMargin (lMargin);
    c6b->SetRightMargin (rMargin);
    c6b->SetBottomMargin (bMargin);
    c6b->SetTopMargin (tMargin);

    c6b->SetLogx ();
    c6b->SetLogy ();

    TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
    //xax->SetTitleFont (43);
    //xax->SetTitleSize (30);
    //xax->SetTitleOffset (1.25);
    xax->SetLabelSize (0);

    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
    double ymin = 2e-3;
    double ymax = 2e3;
    yax->SetRangeUser (ymin, ymax);
    //yax->SetTitleFont (43);
    //yax->SetTitleSize (30);
    //yax->SetTitleOffset (1.3);
    //yax->SetLabelFont (43);
    //yax->SetLabelSize (28);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    tl->SetTextFont (43);
    tl->SetTextSize (36);
    tl->SetTextAlign (21);
    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

        OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.04, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

        Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

        OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_stat);

        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        if (iCent > 0) g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.19, 0.19, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->SetTextAngle (-35);
    tl->DrawLatex (6, 0.15, "#it{p}_{T}^{Z} = 15-30 GeV (#times 1)");
    tl->DrawLatex (9.5, 1.65, "#it{p}_{T}^{Z} = 30-60 GeV (#times 10)");
    tl->DrawLatex (15, 18, "#it{p}_{T}^{Z} = 60+ GeV (#times 10^{2})");
    tl->SetTextAngle (0);

    myMarkerAndBoxAndLineText (0.65, 0.780, 1.4, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.6, "#it{pp}", 0.032);
    myMarkerAndBoxAndLineText (0.65, 0.730, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.780, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.730, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\%", 0.032);

    c6b->SaveAs ("../Plots/FinalPlots/yield_allptz_pTchOnly_onePlot.pdf");
  }




  {
    TCanvas* c6c = new TCanvas ("c6c", "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c6c->SetLeftMargin (lMargin);
    c6c->SetRightMargin (rMargin);
    c6c->SetBottomMargin (0.12);
    c6c->SetTopMargin (0.04);

    c6c->SetLogx ();
    c6c->SetLogy ();

    TH1D* h = new TH1D ("", "", nXhZBins[nPtZBins-1], xhZBins[nPtZBins-1]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{x}_{hZ}");
    xax->SetRangeUser (xhZBins[nPtZBins-1][0], xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
    //xax->SetTitleFont (43);
    //xax->SetTitleSize (30);
    //xax->SetTitleOffset (1.25);
    xax->SetLabelSize (0);

    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{x} d#Delta#phi)");
    const double ymin = 5e-3;
    const double ymax = 4e7;
    yax->SetRangeUser (ymin, ymax);
    //yax->SetTitleFont (43);
    //yax->SetTitleSize (30);
    //yax->SetTitleOffset (1.30);
    //yax->SetLabelFont (43);
    //yax->SetLabelSize (28);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    tl->SetTextFont (43);
    tl->SetTextSize (36);
    tl->SetTextAlign (21);
    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
    tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
    tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
    tl->DrawLatex (1e-1,  yoff, "10^{-1}");
    tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
    tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
    tl->DrawLatex (1,     yoff, "1");

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_sub_syst[iPtZ][iCent]->Clone ();

        OffsetYAxis (g_syst, pow (100, iPtZ-2), true);
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.04, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (finalColors[iCent]);
        g_syst->SetLineColor (finalColors[iCent]);
        g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
      } // end loop over iCent

      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_stat = make_graph (h_trk_xhz_ptz_sub_stat[iPtZ][iCent]);

        Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

        OffsetYAxis (g_stat, pow (100, iPtZ-2), true);
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
        ResetXErrors (g_stat);

        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iCent]);
        g_stat->SetLineColor (finalColors[iCent]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        if (iCent > 0) g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPtZ

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.17, 0.17, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->SetTextAngle (-35);
    tl->DrawLatex (6, 0.15, "#it{p}_{T}^{Z} = 15-30 GeV (#times 1)");
    tl->DrawLatex (9.5, 1.65, "#it{p}_{T}^{Z} = 30-60 GeV (#times 10^{2})");
    tl->DrawLatex (15, 18, "#it{p}_{T}^{Z} = 60+ GeV (#times 10^{4})");
    tl->SetTextAngle (0);

    myMarkerAndBoxAndLineText (0.65, 0.780, 1.4, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.6, "#it{pp}", 0.032);
    myMarkerAndBoxAndLineText (0.65, 0.730, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.780, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.730, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\%", 0.032);

    c6c->SaveAs ("../Plots/FinalPlots/yield_allptz_xhZOnly_onePlot.pdf");
  }




  {
    TCanvas* c7 = new TCanvas ("c7", "", 1200, 1200);

    c7->SetLeftMargin (0.15);
    c7->SetRightMargin (0.04);
    c7->SetBottomMargin (0.15);
    c7->SetTopMargin (0.04);

    c7->SetLogx ();
    //c7->Divide (2, 1);
    short iPtZ = 4;
    const short iCent = 3;

    //const Color_t vitevColor = kCyan+3;

    TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    xax->SetMoreLogLabels ();

    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
    yax->SetRangeUser (0, 3.8);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      c7->cd (iPtZ-2+1);
      gPad->SetLogx ();

      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      ResetXErrors (g);
      deltaize (g, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);
      ResetXErrors (g);

      g->SetFillColorAlpha (finalModelFillColors[iPtZ-1], 0.6);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      //c7->cd (iPtZ-3+1);
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      deltaize (g_syst, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.04, true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      //c7->cd (iPtZ-3+1);
      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);
      ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 3.0 : 2.3);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 3.0 : 2.3);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (50);
    tl->DrawLatexNDC (0.32, 0.880, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (46);
    tl->DrawLatexNDC (0.32, 0.830, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.32, 0.780, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (40);
    tl->SetTextAlign (21);
    tl->DrawLatexNDC (0.378, 0.71, "15-30");
    tl->DrawLatexNDC (0.475, 0.71, "30-60");
    tl->DrawLatexNDC (0.570, 0.71, "60+");

    tl->SetTextAlign (11);
    tl->DrawLatexNDC (0.62, 0.71, "#it{p}_{T}^{Z} [GeV]");
    tl->DrawLatexNDC (0.62, 0.65, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.62, 0.59, "Hybrid Model");
    myMarkerAndBoxAndLineText (0.614, 0.65, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 1.3*2.0, "", 0.045);
    myMarkerAndBoxAndLineText (0.520, 0.65, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.3*2.0, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.65, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.3*2.0, "", 0.045);
    myMarkerAndBoxAndLineText (0.614, 0.59, 1.4, 1001, finalModelFillColors[3], 0.60, finalModelFillColors[3], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.520, 0.59, 1.4, 1001, finalModelFillColors[2], 0.60, finalModelFillColors[2], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.59, 1.4, 1001, finalModelFillColors[1], 0.60, finalModelFillColors[1], -1, 1, "", 0.045);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (pTchBins[nPtZBins-1][0], 1, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1);

    c7->SaveAs ("../Plots/FinalPlots/iaa_ptch_hybridComp.pdf");
  }





  {
    TCanvas* c8 = new TCanvas ("c8", "", 1200, 1200);

    c8->SetLeftMargin (0.15);
    c8->SetRightMargin (0.04);
    c8->SetBottomMargin (0.15);
    c8->SetTopMargin (0.04);

    c8->SetLogx ();
    //c8->Divide (2, 1);
    short iPtZ = 4;
    const short iCent = 3;

    //const Color_t vitevColor = kCyan+3;

    TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    xax->SetMoreLogLabels ();

    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
    yax->SetRangeUser (0, 3.8);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      c8->cd (iPtZ-3+1);
      gPad->SetLogx ();

      TGAE* g = g_jewel_pth[iPtZ];
      ResetXErrors (g);
      deltaize (g, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);
      ResetXErrors (g);

      g->SetFillColorAlpha (finalModelFillColors[iPtZ-1], 0.6);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      //c7->cd (iPtZ-3+1);
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      deltaize (g_syst, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.04, true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      //c7->cd (iPtZ-3+1);
      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      deltaize (g_stat, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);
      ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 3.0 : 2.3);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 3.0 : 2.3);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (50);
    tl->DrawLatexNDC (0.32, 0.880, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (46);
    tl->DrawLatexNDC (0.32, 0.830, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.32, 0.780, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (40);
    tl->SetTextAlign (21);
    tl->DrawLatexNDC (0.378, 0.71, "15-30");
    tl->DrawLatexNDC (0.475, 0.71, "30-60");
    tl->DrawLatexNDC (0.570, 0.71, "60+");

    tl->SetTextAlign (11);
    tl->DrawLatexNDC (0.62, 0.71, "#it{p}_{T}^{Z} [GeV]");
    tl->DrawLatexNDC (0.62, 0.65, "ATLAS 0-10\% Pb+Pb");
    tl->DrawLatexNDC (0.62, 0.59, "JEWEL");
    myMarkerAndBoxAndLineText (0.614, 0.65, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 1.3*2.0, "", 0.045);
    myMarkerAndBoxAndLineText (0.520, 0.65, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.3*2.0, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.65, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.3*2.0, "", 0.045);
    myMarkerAndBoxAndLineText (0.614, 0.59, 1.4, 1001, finalModelFillColors[3], 0.60, finalModelFillColors[3], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.520, 0.59, 1.4, 1001, finalModelFillColors[2], 0.60, finalModelFillColors[2], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.59, 1.4, 1001, finalModelFillColors[1], 0.60, finalModelFillColors[1], -1, 1, "", 0.045);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (pTchBins[nPtZBins-1][0], 1, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1);

    c8->SaveAs ("../Plots/FinalPlots/iaa_ptch_jewelComp.pdf");
  }





  {
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Plots mean track <pTch> vs. <pT^Z>
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* c9 = new TCanvas ("c_mean_pTch_comb", "", 800, 800);
    c9->cd ();

    TH1D* h = new TH1D ("htemp", "", 1, 5, 120);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#LT#it{p}_{T}^{Z}#GT [GeV]");
    xax->SetRangeUser (5, 120);
    xax->SetMoreLogLabels ();

    yax->SetTitle ("#LT#it{p}_{T}^{ch}#GT [GeV]");
    yax->SetRangeUser (0, 12);

    h->SetLineWidth (0);
    h->SetMarkerSize (0);

    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);
    tl->SetTextSize (40);
    tl->DrawLatexNDC (0.22, 0.87, "#bf{#it{ATLAS}} Internal");
    myMarkerAndBoxAndLineText (0.32, 0.810, 2.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.5, "#it{pp}", 0.045);
    myMarkerAndBoxAndLineText (0.32, 0.750, 2.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.5, "30-80\%", 0.045);
    myMarkerAndBoxAndLineText (0.57, 0.810, 2.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.5, "10-30\%", 0.045);
    myMarkerAndBoxAndLineText (0.57, 0.750, 2.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\%", 0.045);

    tl->SetTextSize (34);
    tl->DrawLatexNDC (0.22, 0.27, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
    tl->DrawLatexNDC (0.22, 0.21, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.68, "2 < #it{p}_{T}^{ch} < 240 GeV");

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TGAE* g = (TGAE*) g_trk_avg_pt_ptz_syst[iCent]->Clone ();

      g->SetMarkerSize (0);
      //g->SetLineWidth (0);
      g->SetLineWidth (1);
      g->SetLineColor (finalColors[iCent]);
      g->SetFillColorAlpha (finalFillColors[iCent], 0.3);
      ((TGAE*)g->Clone ())->Draw ("5P");
      g->Draw ("2P");
    } // end loop over iCent

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TGAE* g = (TGAE*) g_trk_avg_pt_ptz_stat[iCent]->Clone ();

      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (markerSize);
      g->SetMarkerColor (finalColors[iCent]);
      g->SetLineColor (finalColors[iCent]);
      g->SetLineWidth (3);
      ((TGAE*) g->Clone ())->Draw ("P");

      if (iCent != 0) {
        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
          
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetLineWidth (0);
        g->SetMarkerColor (kBlack);

        ((TGAE*) g->Clone ())->Draw ("P");
      }

      SaferDelete (&g);
    } // end loop over iCent

    c9->SaveAs ("../Plots/FinalPlots/mean_ptch.pdf");
  }





  {
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Plots mean track <xhZ> vs. <pT^Z>
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* c10 = new TCanvas ("c_mean_xhZ_comb", "", 800, 800);
    c10->cd ();

    TH1D* h = new TH1D ("htemp", "", 1, 5, 120);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#LT#it{p}_{T}^{Z}#GT [GeV]");
    xax->SetRangeUser (5, 120);
    xax->SetMoreLogLabels ();

    yax->SetTitle ("#LT#it{x}_{hZ}#GT");
    yax->SetRangeUser (0, 0.30);

    h->SetLineWidth (0);
    h->SetMarkerSize (0);

    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);
    tl->SetTextSize (40);
    tl->DrawLatexNDC (0.22, 0.87, "#bf{#it{ATLAS}} Internal");
    myMarkerAndBoxAndLineText (0.32, 0.810, 2.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.5, "#it{pp}", 0.045);
    myMarkerAndBoxAndLineText (0.32, 0.750, 2.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.5, "30-80\%", 0.045);
    myMarkerAndBoxAndLineText (0.57, 0.810, 2.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.5, "10-30\%", 0.045);
    myMarkerAndBoxAndLineText (0.57, 0.750, 2.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\%", 0.045);

    tl->SetTextSize (34);
    tl->DrawLatexNDC (0.22, 0.27, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
    tl->DrawLatexNDC (0.22, 0.21, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.68, "1/15 < #it{x}_{hZ} < 2");

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TGAE* g = (TGAE*) g_trk_avg_xhz_ptz_syst[iCent]->Clone ();

      g->SetMarkerSize (0);
      //g->SetLineWidth (0);
      g->SetLineWidth (1);
      g->SetLineColor (finalColors[iCent]);
      g->SetFillColorAlpha (finalFillColors[iCent], 0.3);
      ((TGAE*)g->Clone ())->Draw ("5P");
      g->Draw ("2P");
    } // end loop over iCent

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TGAE* g = (TGAE*) g_trk_avg_xhz_ptz_stat[iCent]->Clone ();

      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (markerSize);
      g->SetMarkerColor (finalColors[iCent]);
      g->SetLineColor (finalColors[iCent]);
      g->SetLineWidth (3);
      ((TGAE*) g->Clone ())->Draw ("P");

      if (iCent != 0) {
        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
          
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetLineWidth (0);
        g->SetMarkerColor (kBlack);

        ((TGAE*) g->Clone ())->Draw ("P");
      }

      SaferDelete (&g);
    } // end loop over iCent

    c10->SaveAs ("../Plots/FinalPlots/mean_xhz.pdf");
  }




  {
    TCanvas* c11 = new TCanvas ("c11", "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c11->SetLeftMargin (lMargin);
    c11->SetRightMargin (rMargin);
    c11->SetBottomMargin (bMargin);
    c11->SetTopMargin (tMargin);

    c11->SetLogx ();
    c11->SetLogy ();

    TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
    //xax->SetTitleFont (43);
    //xax->SetTitleSize (30);
    //xax->SetTitleOffset (1.25);
    xax->SetLabelSize (0);

    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
    double ymin = 2e-3;
    double ymax = 8e3;
    yax->SetRangeUser (ymin, ymax);
    //yax->SetTitleFont (43);
    //yax->SetTitleSize (30);
    //yax->SetTitleOffset (1.3);
    //yax->SetLabelFont (43);
    //yax->SetLabelSize (28);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    tl->SetTextFont (43);
    tl->SetTextSize (36);
    tl->SetTextAlign (21);
    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      const int iCent = 0;
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

      OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      //deltaize (g_syst, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
      //ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.05, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetMarkerColor (finalColors[iPtZ-1]);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      ((TGAE*) g_syst->Clone ())->Draw ("5P");

      SaferDelete (&g_syst);

      TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

      OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.05*(iCent - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
      //ResetXErrors (g_stat);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (kBlack);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = kDot;
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    } // end loop over iPtZ


    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      TGAE* g = g_pythia_pth_ptz[iPtZ];
      OffsetYAxis (g, pow (10, iPtZ-2), true);
      RecenterGraph (g);
      ResetXErrors (g);

      //g->SetMarkerSize (0);
      //g->SetLineWidth (1);
      //g->SetLineColor (finalColors[iPtZ-1]);
      g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g->Draw ("3");
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (40);
    tl->DrawLatexNDC (0.20, 0.265, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.20, 0.200, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.645, 0.89, "#it{p}_{T}^{Z} [GeV]");
    tl->DrawLatexNDC (0.385, 0.885, "15-30");
    tl->DrawLatexNDC (0.477, 0.885, "30-60");
    tl->DrawLatexNDC (0.578, 0.885, "60+");
    myMarkerAndBoxAndLineText (0.46, 0.790, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], -1, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.46, 0.840, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], kOpenCircle, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.55, 0.790, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], -1, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.55, 0.840, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], kOpenCircle, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.64, 0.790, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], -1, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.64, 0.840, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], kOpenCircle, 1.6, "", 0.036);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.645, 0.84, "ATLAS");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.645, 0.79, "Powheg + Pythia8");

    c11->SaveAs ("../Plots/FinalPlots/yield_allptz_pTch_pythiaComp_onePlot.pdf");
  }
}

#endif

