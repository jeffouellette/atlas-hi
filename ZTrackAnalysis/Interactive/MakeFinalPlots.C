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
#include <TPave.h>
#include <TMarker.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace atlashi;
using namespace std;

typedef TGraphErrors TGE;
typedef TGraphAsymmErrors TGAE;

const Color_t starColor     = (Color_t) tcolor->GetColor (50, 175, 50); //(Color_t) tcolor->GetColor (240, 197,  40);
const Style_t starMarker    = kFullCrossX;
//const Color_t cmsColor      = (Color_t) tcolor->GetColor (153,  55,  55);
const Color_t cmsColor      = (Color_t) tcolor->GetColor ( 71,  97, 130);
const Style_t cmsMarker     = kFullCircle;
const Color_t phenixColor   = kMagenta-6;
const Style_t phenixMarker  = kFullSquare;
const Color_t atlasColor    = finalColors[3];
const Color_t atlasFillColor = finalFillColors[3];

const Color_t hybridColor = (Color_t) tcolor->GetColor (82, 82, 173);
const double  hybridAlpha = 0.7;
const Color_t scetgColor  = (Color_t) tcolor->GetColor (39, 180, 66);
const double  scetgAlpha  = 0.6;
const Color_t jewelColor  = tcolor->GetColor (255, 170, 50);
const double  jewelAlpha  = 0.5;

const bool plotXhZ = true;
const double minModelUnc = 0.08;


void SetMinErrors (TGAE* g, const double minerr, const bool logy) {
  double x, y;
  for (int i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    if (y > 0 && g->GetErrorYhigh (i) / pow (y, logy ? 1 : 0) < minerr) g->SetPointEYhigh (i, minerr * pow (y, logy ? 1 : 0));
    if (y > 0 && g->GetErrorYlow (i) / pow (y, logy ? 1 : 0) < minerr) g->SetPointEYlow (i, minerr * pow (y, logy ? 1 : 0));
  }
}


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


  TGAE* tg_CMS_Zh_IAA_0_30_syst = new TGAE ();
  TGAE* tg_CMS_Zh_IAA_0_30_stat = new TGAE ();

  {
    const string fileName = "../DataThieving/CMS_Zh_IAA_HP2020/0_30centrality.dat";
    ifstream f;
    f.open (fileName.c_str ());

    double xlo, x, xhi, y, ystathi, ystatlo, ysysthi, ysystlo;
    f >> xlo >> x >> xhi >> y >> ystathi >> ystatlo >> ysysthi >> ysystlo;
    do {
      tg_CMS_Zh_IAA_0_30_syst->SetPoint (tg_CMS_Zh_IAA_0_30_syst->GetN (), x, y);
      tg_CMS_Zh_IAA_0_30_stat->SetPoint (tg_CMS_Zh_IAA_0_30_stat->GetN (), x, y);

      tg_CMS_Zh_IAA_0_30_syst->SetPointEXlow (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (x-xlo));
      tg_CMS_Zh_IAA_0_30_syst->SetPointEXhigh (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (xhi-x));
      tg_CMS_Zh_IAA_0_30_syst->SetPointEYlow (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (y-ysystlo));
      tg_CMS_Zh_IAA_0_30_syst->SetPointEYhigh (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (ysysthi-y));
      tg_CMS_Zh_IAA_0_30_stat->SetPointEYlow (tg_CMS_Zh_IAA_0_30_stat->GetN () - 1, fabs (y-ystatlo));
      tg_CMS_Zh_IAA_0_30_stat->SetPointEYhigh (tg_CMS_Zh_IAA_0_30_stat->GetN () - 1, fabs (ystathi-y));
      f >> xlo >> x >> xhi >> y >> ystathi >> ystatlo >> ysysthi >> ysystlo;
    } while (f);
  }


  TGAE** g_pythia_pth_ptz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_pythia_finepth_ptz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_pythia_xhz_ptz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_pythia_finexhz_ptz = Get1DArray <TGAE*> (nPtZBins);
  {
    TFile* pythiaFile = new TFile ("../rootFiles/TruthAnalysis/Nominal/pythiaCompare.root", "read");
    TH1D* h_z_counts = (TH1D*) pythiaFile->Get ("h_z_counts");
    TH1D* h = nullptr, *hSig = nullptr;
    TH2D* h2 = nullptr;
    TGAE* g = nullptr;
    double x, y;

    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      const float nTotalEvents = h_z_counts->GetBinContent (iPtZ+1);

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_pTch_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_pTch_cov_iPtZ%i", iPtZ));

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

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_pTch_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_pTch_cov_bkg_iPtZ%i", iPtZ));

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



      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FinePtch_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FinePtch_cov_iPtZ%i", iPtZ));

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

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FinePtch_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FinePtch_cov_bkg_iPtZ%i", iPtZ));

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
      g_pythia_finepth_ptz[iPtZ] = make_graph (hSig);

      g = g_pythia_finepth_ptz[iPtZ];

      
      h = (TH1D*) pythiaFile->Get (Form ("h_trk_xhZ_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_xhZ_cov_iPtZ%i", iPtZ));

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

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_xhZ_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_xhZ_cov_bkg_iPtZ%i", iPtZ));

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



      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FineXhZ_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FineXhZ_cov_iPtZ%i", iPtZ));

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

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FineXhZ_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FineXhZ_cov_bkg_iPtZ%i", iPtZ));

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
      g_pythia_finexhz_ptz[iPtZ] = make_graph (hSig);
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
    //string modelFileName = "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_1.dat"; // means full medium response, including also the negative contribution from the wake
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
    //modelFileName = "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_1.dat";
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



  TGAE** g_scetg_pth = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_scetg_xhz = Get1DArray <TGAE*> (nPtZBins);
  for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    g_scetg_pth[iPtZ] = new TGAE ();
    g_scetg_xhz[iPtZ] = new TGAE ();
  }

  {
    ifstream f;
    double dummy = 0, x = 0, y = 0;

    TGAE* g;
    vector<double> xarr (0);
    vector<double> yarr_g1_8 (0), yarr_g2_0 (0), yarr_g2_2 (0); 
    string modelFileName;

    g = g_scetg_pth[3];
    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASpt.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASpt.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASpt.Z30-60";
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



    g = g_scetg_xhz[3];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASnew.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASnew.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASnew.Z30-60";
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



    g = g_scetg_pth[4];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();
    
    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASpt.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASpt.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASpt.Z60-100";
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



    g = g_scetg_xhz[4];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();
    
    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASnew.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASnew.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASnew.Z60-100";
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

      h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_medium_iSpc2_iPtZ%i", iPtZ));
      h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_vacuum_iSpc2_iPtZ%i", iPtZ));
      h_num->Divide (h_den);
      g_jewel_xhz[iPtZ] = make_graph (h_num);
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

        xax->SetTitle ("#it{x}_{hZ} = #it{p}_{T}^{ch} #/#it{p}_{T}^{Z}");
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
          tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
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
    tl->DrawLatexNDC (0.22, 0.80, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.74, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    if (plotXhZ) {
      ldPad->cd ();
      tl->SetTextSize (32);
      tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.80, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
      tl->DrawLatexNDC (0.22, 0.74, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
    }

    tl->SetTextSize (30);
    luPad->cd ();
    tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{Z} < 30 GeV");
    cuPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{Z} < 60 GeV");
    ruPad->cd ();
    tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} > 60 GeV");

    if (plotXhZ) {
      ldPad->cd ();
      tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{Z} < 30 GeV");
      cdPad->cd ();
      tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{Z} < 60 GeV");
      rdPad->cd ();
      tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} > 60 GeV");
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
    c2->SetBottomMargin (bMargin);
    c2->SetTopMargin (tMargin);

    c2->SetLogx ();
    c2->SetLogy ();

    {
      TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (30);
      //xax->SetTitleOffset (1.25);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      const double ymin = 0.06;
      const double ymax = 2000;
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
    }

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

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.19, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->DrawLatex (13, 1.28, "15 < #it{p}_{T}^{Z} < 30 GeV (#times 1)");
    tl->DrawLatex (13, 12.8, "30 < #it{p}_{T}^{Z} < 60 GeV (#times 10)");
    tl->DrawLatex (13, 128, "#it{p}_{T}^{Z} > 60 GeV (#times 10^{2})");

    tl->SetTextSize (28);
    tl->SetTextAlign (11);
    tl->DrawLatexNDC (0.765, 0.320-0.011, "30-80\% #/#it{pp}");
    tl->DrawLatexNDC (0.765, 0.270-0.011, "10-30\% #/#it{pp}");
    tl->DrawLatexNDC (0.765, 0.220-0.011, "0-10\% #/#it{pp}");
    MakeDataBox   (0.78, 0.320, finalFillColors[1], 0.30, markerStyles[0], 1.6);
    MakeDataBox   (0.78, 0.270, finalFillColors[2], 0.30, markerStyles[1], 1.6);
    MakeDataBox   (0.78, 0.220, finalFillColors[3], 0.30, markerStyles[2], 2.2);
    //myMarkerAndBoxAndLineText (0.76, 0.280, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\% #/#it{pp}", 0.032);
    //myMarkerAndBoxAndLineText (0.76, 0.230, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\% #/#it{pp}", 0.032);
    //myMarkerAndBoxAndLineText (0.76, 0.180, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\% #/#it{pp}", 0.032);

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
        const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
        deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
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
        const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
        deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
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
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.26, 0.19, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->SetTextAngle (-35);
    tl->DrawLatex (1.8, 0.2, "15 < #it{p}_{T}^{Z} < 30 GeV (#times 1)");
    tl->DrawLatex (3, 2.2, "30 < #it{p}_{T}^{Z} < 60 GeV (#times 10)");
    tl->DrawLatex (5.4, 15, "#it{p}_{T}^{Z} > 60 GeV (#times 10^{2})");
    tl->SetTextAngle (0);

    TLine* dashedLines = new TLine ();
    dashedLines->SetLineStyle (2);
    dashedLines->SetLineColor (kBlack);

    dashedLines->DrawLine (1, 17., 50., 2e-3);
    dashedLines->DrawLine (1, 350., 60., exp (((log(2e-3)-log(17.))/log(50.))*log(60.) + log(350.)));

    myMarkerAndBoxAndLineText (0.65, 0.780, 1.4, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.6, "#it{pp}", 0.032);
    myMarkerAndBoxAndLineText (0.65, 0.730, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.780, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\%", 0.032);
    myMarkerAndBoxAndLineText (0.82, 0.730, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\%", 0.032);

    dPad->cd ();
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.26, 0.19, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->DrawLatex (13, 1.25, "15 < #it{p}_{T}^{Z} < 30 GeV (#times 1)");
    tl->DrawLatex (13, 12.5, "30 < #it{p}_{T}^{Z} < 60 GeV (#times 10)");
    tl->DrawLatex (13, 125, "#it{p}_{T}^{Z} > 60 GeV (#times 10^{2})");

    myMarkerAndBoxAndLineText (0.77, 0.280, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\% #/#it{pp}", 0.036);
    myMarkerAndBoxAndLineText (0.77, 0.230, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\% #/#it{pp}", 0.036);
    myMarkerAndBoxAndLineText (0.77, 0.180, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\% #/#it{pp}", 0.036);

    c3->SaveAs ("../Plots/FinalPlots/yield_and_iaa_allptz_pTchOnly_onePlot.pdf");
  }




  {
    TCanvas* c4 = new TCanvas ("c4", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c4->SetLeftMargin (lMargin);
    c4->SetRightMargin (rMargin);
    c4->SetBottomMargin (bMargin);
    c4->SetTopMargin (tMargin);

    short iPtZ = nPtZBins-1;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{h,#gamma/Z} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#gamma/Z}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.05;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.054 * (ymax-ymin) / (1.-tMargin-bMargin);
      if (iPtZ > 2) {
        if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
        tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
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
    g->SetMarkerSize (1.8);
    g->SetMarkerColor (cmsColor);
    g->SetLineColor (cmsColor);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    if (IsFullMarker (cmsMarker)) {
      g->SetMarkerStyle (FullToOpenMarker (cmsMarker));
      g->SetMarkerSize (1.8);
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

    g_stat->SetMarkerStyle (kFullDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (atlasColor);
    g_stat->SetLineColor (atlasColor);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    g_stat->SetMarkerStyle (kOpenDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    SaferDelete (&g_stat);

    myMarkerAndBoxAndLineText (0.33, 0.90-0.012, 1.4, 1001, atlasFillColor, 0.30, atlasColor, kFullDiamond,  2.4, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 0-10% Pb+Pb", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.85-0.012, 1.4, 1001, cmsColor,       0.30, cmsColor,   cmsMarker, 1.8, "CMS #it{p}_{T}^{#gamma} > 60 GeV, #it{p}_{T}^{jet} > 30 GeV, 0-10\% Pb+Pb", 0.032);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.24, 0.300, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.26, 0.255, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.26, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    c4->SaveAs ("../Plots/FinalPlots/iaa_xhz_cmsComp.pdf");
  }




  {
    TCanvas* c5 = new TCanvas ("c5", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c5->SetLeftMargin (lMargin);
    c5->SetRightMargin (rMargin);
    c5->SetBottomMargin (bMargin);
    c5->SetTopMargin (tMargin);

    short iPtZ = nPtZBins-1;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{h,#gamma/Z} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#gamma/Z}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.05;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.054 * (ymax-ymin) / (1.-tMargin-bMargin);
      if (iPtZ > 2) {
        if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
        tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
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
    g->SetMarkerSize (1.8);
    g->SetMarkerColor (phenixColor);
    g->SetLineColor (phenixColor);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    if (IsFullMarker (phenixMarker)) {
      g->SetMarkerStyle (FullToOpenMarker (phenixMarker));
      g->SetMarkerSize (1.8);
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
    g->SetMarkerSize (1.8);
    g->SetMarkerColor (starColor);
    g->SetLineColor (starColor);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    if (IsFullMarker (starMarker)) {
      g->SetMarkerStyle (FullToOpenMarker (starMarker));
      g->SetMarkerSize (1.8);
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

    g_stat->SetMarkerStyle (kFullDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (atlasColor);
    g_stat->SetLineColor (atlasColor);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    g_stat->SetMarkerStyle (kOpenDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);

    ((TGAE*) g_stat->Clone ())->Draw ("P");

    SaferDelete (&g_stat);

    myMarkerAndBoxAndLineText (0.33, 0.90-0.012, 1.4, 1001, atlasFillColor, 0.30, atlasColor, kFullDiamond,  2.4, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 0-10% Pb+Pb", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.85-0.012, 1.4, 1001, phenixColor,    0.30, phenixColor,  phenixMarker, 1.8, "PHENIX 5 < #it{p}_{T}^{#gamma} < 9 GeV, 0-40\% Au+Au", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.80-0.012, 1.4, 1001, starColor,      0.20, starColor,    starMarker, 1.8, "STAR 12 < #it{p}_{T}^{#gamma} < 20 GeV, 0-12\% Au+Au", 0.032);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.24, 0.310, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.26, 0.260, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.26, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    c5->SaveAs ("../Plots/FinalPlots/iaa_xhz_rhicComp.pdf");
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
        const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
        deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
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
        const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
        deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
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

        xax->SetTitle ("#it{x}_{hZ} = #it{p}_{T}^{ch} #/#it{p}_{T}^{Z}");
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
          const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
          deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
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
          const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
          deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
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
    tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{Z} < 30 GeV");
    cuPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{Z} < 60 GeV");
    ruPad->cd ();
    tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} > 60 GeV");

    if (plotXhZ) {
      ldPad->cd ();
      tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{Z} < 30 GeV");
      cdPad->cd ();
      tl->DrawLatexNDC (0.50, 0.86, "30 <#it{p}_{T}^{Z} < 60 GeV");
      rdPad->cd ();
      tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{Z} > 60 GeV");
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
    double ymin = 5e-4;
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
      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

        OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
        deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
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
        const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
        deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
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

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.19, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (24);
    tl->SetTextAngle (-35);
    tl->DrawLatex (4, 1, "15 < #it{p}_{T}^{Z} < 30 GeV (#times 1)");
    tl->DrawLatex (8, 6, "30 < #it{p}_{T}^{Z} < 60 GeV (#times 10)");
    tl->DrawLatex (15, 40, "#it{p}_{T}^{Z} > 60 GeV (#times 10^{2})");
    tl->SetTextAngle (0);

    //TLine* dashedLines = new TLine ();
    //dashedLines->SetLineStyle (2);
    //dashedLines->SetLineColor (kBlack);

    //dashedLines->DrawLine (1, 17., 50., ymin);
    //dashedLines->DrawLine (1, 350., 60., exp (((log(ymin)-log(17.))/log(50.))*log(60.) + log(350.)));

    myMarkerAndBoxAndLineText (0.30, 0.245, 1.4, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.6, "#it{pp}", 0.032);
    myMarkerAndBoxAndLineText (0.30, 0.200, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\%", 0.032);
    myMarkerAndBoxAndLineText (0.47, 0.245, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\%", 0.032);
    myMarkerAndBoxAndLineText (0.47, 0.200, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\%", 0.032);

    c6b->SaveAs ("../Plots/FinalPlots/yield_allptz_pTchOnly_onePlot.pdf");
  }




  //{
  //  TCanvas* c6c = new TCanvas ("c6c", "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c6c->SetLeftMargin (lMargin);
  //  c6c->SetRightMargin (rMargin);
  //  c6c->SetBottomMargin (0.12);
  //  c6c->SetTopMargin (0.04);

  //  c6c->SetLogx ();
  //  c6c->SetLogy ();

  //  TH1D* h = new TH1D ("", "", nXhZBins[nPtZBins-1], xhZBins[nPtZBins-1]);

  //  TAxis* xax = h->GetXaxis ();
  //  TAxis* yax = h->GetYaxis ();

  //  xax->SetTitle ("#it{x}_{hZ} = #it{p}_{T}^{ch} #/#it{p}_{T}^{Z}");
  //  xax->SetRangeUser (xhZBins[nPtZBins-1][0], xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
  //  //xax->SetTitleFont (43);
  //  //xax->SetTitleSize (30);
  //  //xax->SetTitleOffset (1.25);
  //  xax->SetLabelSize (0);

  //  yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{x} d#Delta#phi)");
  //  const double ymin = 5e-3;
  //  const double ymax = 4e7;
  //  yax->SetRangeUser (ymin, ymax);
  //  //yax->SetTitleFont (43);
  //  //yax->SetTitleSize (30);
  //  //yax->SetTitleOffset (1.30);
  //  //yax->SetLabelFont (43);
  //  //yax->SetLabelSize (28);

  //  h->SetLineWidth (0);

  //  h->DrawCopy ("");
  //  SaferDelete (&h);

  //  tl->SetTextFont (43);
  //  tl->SetTextSize (36);
  //  tl->SetTextAlign (21);
  //  const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //  tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
  //  tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
  //  tl->DrawLatex (1e-1,  yoff, "10^{-1}");
  //  tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
  //  tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
  //  tl->DrawLatex (1,     yoff, "1");

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_sub_syst[iPtZ][iCent]->Clone ();

  //      OffsetYAxis (g_syst, pow (100, iPtZ-2), true);
  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.04, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      ((TGAE*) g_syst->Clone ())->Draw ("5P");

  //      SaferDelete (&g_syst);
  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_xhz_ptz_sub_stat[iPtZ][iCent]);

  //      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

  //      OffsetYAxis (g_stat, pow (100, iPtZ-2), true);
  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //      ResetXErrors (g_stat);

  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
  //      
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      if (iCent > 0) g_stat->SetLineWidth (0);
  //      g_stat->SetMarkerColor (kBlack);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.17, 0.17, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (26);
  //  tl->DrawLatexNDC (0.45, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.45, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (24);
  //  tl->SetTextAngle (-35);
  //  tl->DrawLatex (6, 0.15, "15 < #it{p}_{T}^{Z} < 30 GeV (#times 1)");
  //  tl->DrawLatex (9.5, 1.65, "30 < #it{p}_{T}^{Z} < 60 GeV (#times 10^{2})");
  //  tl->DrawLatex (15, 18, "#it{p}_{T}^{Z} > 60 GeV (#times 10^{4})");
  //  tl->SetTextAngle (0);

  //  myMarkerAndBoxAndLineText (0.65, 0.780, 1.4, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.6, "#it{pp}", 0.032);
  //  myMarkerAndBoxAndLineText (0.65, 0.730, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\%", 0.032);
  //  myMarkerAndBoxAndLineText (0.82, 0.780, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\%", 0.032);
  //  myMarkerAndBoxAndLineText (0.82, 0.730, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\%", 0.032);

  //  c6c->SaveAs ("../Plots/FinalPlots/yield_allptz_xhZOnly_onePlot.pdf");
  //}




  {
    TCanvas* c7 = new TCanvas ("c7", "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c7->SetLeftMargin (lMargin);
    c7->SetRightMargin (rMargin);
    c7->SetBottomMargin (bMargin);
    c7->SetTopMargin (tMargin);

    c7->SetLogx ();
    //c7->Divide (2, 1);
    short iPtZ = 4;
    const short iCent = 3;

    {
      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.2);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      //yax->SetRangeUser (0.05, 10);
      const double ymin = 0.;
      const double ymax = 3.6;
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

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      c7->cd (iPtZ-2+1);
      gPad->SetLogx ();

      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, false);
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
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (30);
    tl->SetTextAlign (21);
    tl->DrawLatexNDC (0.378, 0.71, "15-30");
    tl->DrawLatexNDC (0.485, 0.71, "30-60");
    tl->DrawLatexNDC (0.590, 0.71, "60+");

    tl->SetTextAlign (11);
    tl->DrawLatexNDC (0.64, 0.71, "#it{p}_{T}^{Z} [GeV]");
    tl->DrawLatexNDC (0.64, 0.65, "ATLAS 0-10\%");
    tl->DrawLatexNDC (0.64, 0.59, "Hybrid Model");
    myMarkerAndBoxAndLineText (0.634, 0.65, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "", 0.045);
    myMarkerAndBoxAndLineText (0.530, 0.65, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.65, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "", 0.045);
    myMarkerAndBoxAndLineText (0.634, 0.59, 1.4, 1001, finalModelFillColors[3], 0.60, finalModelFillColors[3], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.530, 0.59, 1.4, 1001, finalModelFillColors[2], 0.60, finalModelFillColors[2], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.59, 1.4, 1001, finalModelFillColors[1], 0.60, finalModelFillColors[1], -1, 1, "", 0.045);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (pTchBins[nPtZBins-1][0], 1, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1);

    c7->SaveAs ("../Plots/FinalPlots/iaa_ptch_hybridComp.pdf");
  }





  {
    TCanvas* c8 = new TCanvas ("c8", "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c8->SetLeftMargin (lMargin);
    c8->SetRightMargin (rMargin);
    c8->SetBottomMargin (bMargin);
    c8->SetTopMargin (tMargin);

    c8->SetLogx ();
    //c8->Divide (2, 1);
    short iPtZ = 4;
    const short iCent = 3;

    {
      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.2);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      //yax->SetRangeUser (0.05, 10);
      const double ymin = 0.;
      const double ymax = 3.6;
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

    for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      c8->cd (iPtZ-3+1);
      gPad->SetLogx ();

      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, false);
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
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (30);
    tl->SetTextAlign (21);
    tl->DrawLatexNDC (0.378, 0.71, "15-30");
    tl->DrawLatexNDC (0.485, 0.71, "30-60");
    tl->DrawLatexNDC (0.590, 0.71, "60+");

    tl->SetTextAlign (11);
    tl->DrawLatexNDC (0.64, 0.71, "#it{p}_{T}^{Z} [GeV]");
    tl->DrawLatexNDC (0.64, 0.65, "ATLAS 0-10\%");
    tl->DrawLatexNDC (0.64, 0.59, "JEWEL");
    myMarkerAndBoxAndLineText (0.634, 0.65, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "", 0.045);
    myMarkerAndBoxAndLineText (0.530, 0.65, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "", 0.045);
    myMarkerAndBoxAndLineText (0.424, 0.65, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "", 0.045);
    myMarkerAndBoxAndLineText (0.634, 0.59, 1.4, 1001, finalModelFillColors[3], 0.60, finalModelFillColors[3], -1, 1, "", 0.045);
    myMarkerAndBoxAndLineText (0.530, 0.59, 1.4, 1001, finalModelFillColors[2], 0.60, finalModelFillColors[2], -1, 1, "", 0.045);
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
    TCanvas* c11 = new TCanvas ("c11", "", 800, 880);
    c11->Draw ();

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double utMargin = 0.04;
    const double ubMargin = 0;//0.01;
    const double dtMargin = 0;//0.02;
    const double dbMargin = 0.35;

    TPad* uPad = new TPad ("c11_uPad", "", 0, 0.3, 1, 1);
    TPad* dPad = new TPad ("c11_dPad", "", 0, 0, 1, 0.3);

    uPad->SetLeftMargin (lMargin);
    uPad->SetRightMargin (rMargin);
    dPad->SetLeftMargin (lMargin);
    dPad->SetRightMargin (rMargin);
    uPad->SetTopMargin (utMargin);
    uPad->SetBottomMargin (ubMargin);
    dPad->SetTopMargin (dtMargin);
    dPad->SetBottomMargin (dbMargin);

    uPad->Draw ();
    dPad->Draw ();

    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    {
      TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleSize (0);
      xax->SetTickLength (0.03 * (1.-utMargin-ubMargin) / uPad->GetHNDC ());
      xax->SetTickLength (0.03);
      xax->SetLabelSize (0);

      xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

      yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
      yax->SetTitleFont (43);
      yax->SetTitleSize (36);
      yax->SetTickLength (0.02 * (1.-utMargin-ubMargin) / uPad->GetHNDC ());
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      double ymin = 2e-3;
      double ymax = 8e3;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);
    }

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      TGAE* g_pyth = (TGAE*) g_pythia_finepth_ptz[iPtZ]->Clone ();
      SetMinErrors (g_pyth, 0.10, true);

      OffsetYAxis (g_pyth, pow (10, iPtZ-2), true);
      RecenterGraph (g_pyth);
      //ResetXErrors (g_pyth);

      //g_pyth->SetMarkerSize (0);
      //g_pyth->SetLineWidth (1);
      //g_pyth->SetLineColor (finalColors[iPtZ-1]);
      g_pyth->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      ((TGAE*) g_pyth->Clone ())->Draw ("3");

      SaferDelete (&g_pyth);
    }

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      const int iCent = 0;
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

      OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      //deltaize (g_syst, 0.05*(0.5*(numCentBins-1)-iCent), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
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
      //ResetXErrors (g_stat);
      //deltaize (g_stat, 0.05*(0.5*(numCentBins-1)-iCent), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
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

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.20, 0.14, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.20, 0.09, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.20, 0.04, "Powheg + Pythia 8.186");

    tl->SetTextSize (22);
    tl->DrawLatexNDC (0.730, 0.885, "#it{p}_{T}^{Z} [GeV]");
    tl->DrawLatexNDC (0.730, 0.835, "15-30 (#times 1)");
    tl->DrawLatexNDC (0.730, 0.785, "30-60 (#times 10)");
    tl->DrawLatexNDC (0.730, 0.735, "60+ (#times 10^{2})");
    myMarkerAndBoxAndLineText (0.72, 0.840, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], -1, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.63, 0.840, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], kOpenCircle, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.72, 0.790, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], -1, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.63, 0.790, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], kOpenCircle, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.72, 0.740, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], -1, 1.6, "", 0.036);
    myMarkerAndBoxAndLineText (0.63, 0.740, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], kOpenCircle, 1.6, "", 0.036);

    tl->SetTextSize (22);
    tl->DrawLatexNDC (0.565, 0.885, "Data");
    tl->SetTextSize (22);
    tl->DrawLatexNDC (0.663, 0.885, "MC");


    dPad->cd ();
    dPad->SetLogx ();

    {
      TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleFont (43);
      xax->SetTitleSize (36);
      xax->SetTitleOffset (3.6);
      xax->SetTickLength (0.03 * (1.-dtMargin-dbMargin) / dPad->GetHNDC ());
      xax->SetLabelSize (0);
      xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

      yax->SetTitle ("MC / data");
      yax->SetTitleFont (43);
      yax->SetTitleSize (36);
      yax->CenterTitle (true);
      yax->SetTickLength (0.02 * (1.-dtMargin-dbMargin) / dPad->GetHNDC ());
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      double ymin = 0.88;
      double ymax = 1.12;
      yax->SetRangeUser (ymin, ymax);
      yax->SetNdivisions (504);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      const double yoff = ymin - (0.12 * (ymax - ymin) / (1.-dtMargin-dbMargin));
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

      TLine* line = new TLine (pTchBins[nPtZBins-1][0], 1., pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1.);
      line->SetLineColor (kBlack);
      line->SetLineStyle (2);
      line->Draw ("same");
    }

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      const int iCent = 0;
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();
      TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);
      TGAE* g_pyth = (TGAE*) g_pythia_pth_ptz[iPtZ]->Clone ();

      RecenterGraph (g_syst);
      RecenterGraph (g_stat);
      RecenterGraph (g_pyth);

      TGAE* g_ratio_stat = new TGAE ();
      TGAE* g_ratio_syst = new TGAE ();

      double x, y_data, y_pyth;
      double y_stat_hi, y_stat_lo, y_syst_hi, y_syst_lo, y_pyth_hi, y_pyth_lo, x_syst_hi, x_syst_lo;
      for (int iX = 0; iX < g_syst->GetN (); iX++) {
        g_syst->GetPoint (iX, x, y_data);
        y_syst_hi = g_syst->GetErrorYhigh (iX);
        y_syst_lo = g_syst->GetErrorYlow (iX);
        x_syst_hi = g_syst->GetErrorXhigh (iX);
        x_syst_lo = g_syst->GetErrorXlow (iX);
        g_stat->GetPoint (iX, x, y_data);
        y_stat_hi = g_stat->GetErrorYhigh (iX);
        y_stat_lo = g_stat->GetErrorYlow (iX);
        g_pyth->GetPoint (iX, x, y_pyth);
        y_pyth_hi = g_pyth->GetErrorYhigh (iX);
        y_pyth_lo = g_pyth->GetErrorYlow (iX);

        const double ratio = y_pyth / y_data;

        g_ratio_syst->SetPoint (g_ratio_syst->GetN (), x, ratio);
        g_ratio_syst->SetPointEXhigh (g_ratio_syst->GetN () - 1, x_syst_hi);
        g_ratio_syst->SetPointEXlow (g_ratio_syst->GetN () - 1, x_syst_lo);
        g_ratio_syst->SetPointEYhigh (g_ratio_syst->GetN () - 1, fabs (ratio) * fabs (y_syst_hi / y_data));
        g_ratio_syst->SetPointEYlow (g_ratio_syst->GetN () - 1, fabs (ratio) * fabs (y_syst_lo / y_data));

        g_ratio_stat->SetPoint (g_ratio_stat->GetN (), x, ratio);
        g_ratio_stat->SetPointEXhigh (g_ratio_stat->GetN () - 1, 0);
        g_ratio_stat->SetPointEXlow (g_ratio_stat->GetN () - 1, 0);
        g_ratio_stat->SetPointEYhigh (g_ratio_stat->GetN () - 1, fabs (ratio) * sqrt (pow (y_pyth_hi / y_pyth, 2) + pow (y_stat_hi / y_data, 2)));
        g_ratio_stat->SetPointEYlow (g_ratio_stat->GetN () - 1, fabs (ratio) * sqrt (pow (y_pyth_lo / y_pyth, 2) + pow (y_stat_lo / y_data, 2)));
      }

      g_ratio_syst->SetMarkerSize (0);
      g_ratio_syst->SetLineWidth (1);
      g_ratio_syst->SetMarkerColor (finalColors[iPtZ-1]);
      g_ratio_syst->SetLineColor (finalColors[iPtZ-1]);
      g_ratio_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      ((TGAE*) g_ratio_syst->Clone ())->Draw ("5P");

      SaferDelete (&g_ratio_syst);

      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

      g_ratio_stat->SetMarkerStyle (markerStyle);
      g_ratio_stat->SetMarkerSize (markerSize);
      g_ratio_stat->SetLineWidth (3);
      g_ratio_stat->SetMarkerColor (kBlack);
      g_ratio_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_ratio_stat->Clone ())->Draw ("P");

      markerStyle = kDot;
      
      g_ratio_stat->SetMarkerStyle (markerStyle);
      g_ratio_stat->SetMarkerSize (markerSize);
      g_ratio_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_ratio_stat->Clone ())->Draw ("P");

      SaferDelete (&g_ratio_stat);

      SaferDelete (&g_syst);
      SaferDelete (&g_stat);
      SaferDelete (&g_pyth);
    }

    c11->SaveAs ("../Plots/FinalPlots/yield_allptz_pTch_pythiaComp_onePlot.pdf");
  }




  {
    TCanvas* c12 = new TCanvas ("c12", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c12->SetLeftMargin (lMargin);
    c12->SetRightMargin (rMargin);
    c12->SetBottomMargin (bMargin);
    c12->SetTopMargin (tMargin);

    const short iPtZ = 3;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.18;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
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
      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
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

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
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

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.50, 0.750, "30 < #it{p}_{T}^{Z} < 60 GeV");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.59, 0.705-0.011, "ATLAS 0-10\%");
    tl->DrawLatexNDC (0.59, 0.660-0.011, "Hybrid Model");
    tl->DrawLatexNDC (0.59, 0.615-0.011, "SCET_{G}");
    tl->DrawLatexNDC (0.59, 0.570-0.011, "JEWEL");

    MakeDataBox   (0.60, 0.705, finalFillColors[2], 0.30, markerStyles[1], 1.8);
    MakeTheoryBox (0.60, 0.660, hybridColor, hybridAlpha);
    MakeTheoryBox (0.60, 0.615, scetgColor, scetgAlpha);
    MakeTheoryBox (0.60, 0.570, jewelColor, jewelAlpha);

    c12->SaveAs ("../Plots/FinalPlots/iaa_pTch_theoryComp_iPtZ3.pdf");
  }




  {
    TCanvas* c13 = new TCanvas ("c13", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c13->SetLeftMargin (lMargin);
    c13->SetRightMargin (rMargin);
    c13->SetBottomMargin (bMargin);
    c13->SetTopMargin (tMargin);

    const short iPtZ = 4;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.18;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
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
      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
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

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
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

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.50, 0.750, "#it{p}_{T}^{Z} > 60 GeV");

    tl->SetTextSize (28);
    //tl->DrawLatexNDC (0.59, 0.705-0.011, "2018 Pb+Pb 0-10\%");
    //tl->DrawLatexNDC (0.59, 0.660-0.011, "'18+'15 Pb+Pb (#it{proj.})");
    //tl->DrawLatexNDC (0.59, 0.615-0.011, "Hybrid Model");
    //tl->DrawLatexNDC (0.59, 0.570-0.011, "JEWEL");

    //MakeDataBox   (0.60, 0.705, finalFillColors[3], 0.30, kFullCircle, 2.3);
    //MakeDataBox   (0.60, 0.660, finalFillColors[2], 0.30, kFullCircle, 2.3);
    //MakeTheoryBox (0.60, 0.615, hybridColor, hybridAlpha);
    //MakeTheoryBox (0.60, 0.570, jewelColor, jewelAlpha);

    tl->DrawLatexNDC (0.59, 0.705-0.011, "ATLAS 0-10\%");
    tl->DrawLatexNDC (0.59, 0.660-0.011, "Hybrid Model");
    tl->DrawLatexNDC (0.59, 0.615-0.011, "SCET_{G}");
    tl->DrawLatexNDC (0.59, 0.570-0.011, "JEWEL");

    MakeDataBox   (0.60, 0.705, finalFillColors[3], 0.30, markerStyles[2], 2.4);
    MakeTheoryBox (0.60, 0.660, hybridColor, hybridAlpha);
    MakeTheoryBox (0.60, 0.615, scetgColor, scetgAlpha);
    MakeTheoryBox (0.60, 0.570, jewelColor, jewelAlpha);

    c13->SaveAs ("../Plots/FinalPlots/iaa_pTch_theoryComp_iPtZ4.pdf");
  }




  {
    TCanvas* c14 = new TCanvas ("c14", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c14->SetLeftMargin (lMargin);
    c14->SetRightMargin (rMargin);
    c14->SetBottomMargin (bMargin);
    c14->SetTopMargin (tMargin);

    const short iPtZ = 3;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{hZ} = #it{p}_{T}^{ch} #/#it{p}_{T}^{Z}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.05;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.054 * (ymax-ymin) / (1.-tMargin-bMargin);
      if (iPtZ > 2) {
        if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
        tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
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

    {
      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.24, 0.290, "30 < #it{p}_{T}^{Z} < 60 GeV");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.33, 0.245-0.011, "ATLAS 0-10\%");
    tl->DrawLatexNDC (0.33, 0.200-0.011, "Hybrid Model");
    tl->DrawLatexNDC (0.65, 0.245-0.011, "SCET_{G}");
    tl->DrawLatexNDC (0.65, 0.200-0.011, "JEWEL");

    MakeDataBox   (0.34, 0.245, finalFillColors[2], 0.30, markerStyles[1], 1.8);
    MakeTheoryBox (0.34, 0.200, hybridColor, hybridAlpha);
    MakeTheoryBox (0.66, 0.245, scetgColor, scetgAlpha);
    MakeTheoryBox (0.66, 0.200, jewelColor, jewelAlpha);

    c14->SaveAs ("../Plots/FinalPlots/iaa_xhZ_theoryComp_iPtZ3.pdf");
  }




  {
    TCanvas* c15 = new TCanvas ("c15", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c15->SetLeftMargin (lMargin);
    c15->SetRightMargin (rMargin);
    c15->SetBottomMargin (bMargin);
    c15->SetTopMargin (tMargin);

    const short iPtZ = 4;
    const short iCent = 3;

    {
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{hZ} = #it{p}_{T}^{ch} #/#it{p}_{T}^{Z}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.05;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.054 * (ymax-ymin) / (1.-tMargin-bMargin);
      if (iPtZ > 2) {
        if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
        tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
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

    {
      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.24, 0.290, "#it{p}_{T}^{Z} > 60 GeV");

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.33, 0.245-0.011, "ATLAS 0-10\%");
    tl->DrawLatexNDC (0.33, 0.200-0.011, "Hybrid Model");
    tl->DrawLatexNDC (0.65, 0.245-0.011, "SCET_{G}");
    tl->DrawLatexNDC (0.65, 0.200-0.011, "JEWEL");

    MakeDataBox   (0.34, 0.245, finalFillColors[3], 0.30, markerStyles[2], 2.4);
    MakeTheoryBox (0.34, 0.200, hybridColor, hybridAlpha);
    MakeTheoryBox (0.66, 0.245, scetgColor, scetgAlpha);
    MakeTheoryBox (0.66, 0.200, jewelColor, jewelAlpha);

    c15->SaveAs ("../Plots/FinalPlots/iaa_xhZ_theoryComp_iPtZ4.pdf");
  }



  {
    TCanvas* c16 = new TCanvas ("c16", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;
  
    c16->SetLeftMargin (lMargin);
    c16->SetRightMargin (rMargin);
    c16->SetBottomMargin (bMargin);
    c16->SetTopMargin (tMargin);
  
    short iPtZ = nPtZBins-1;
  
    {
      gPad->SetLogx ();
      gPad->SetLogy ();
  
      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);
  
      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();
  
      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      xax->SetLabelSize (0);
  
      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.18;
      const double ymax = 7;
      //const double ymin = 0.;
      //const double ymax = 3.6;
      yax->SetRangeUser (ymin, ymax);
  
      h->SetLineWidth (0);
  
      h->DrawCopy ("");
      SaferDelete (&h);
  
      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
  
      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
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
  
    TGAE* g = nullptr;
  
    g = (TGAE*) tg_CMS_Zh_IAA_0_30_syst->Clone ();
    g->SetMarkerSize (0);
    g->SetMarkerColor (cmsColor);
    g->SetLineColor (kGray+2);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (kGray+2, 0.3);
    ((TGE*) g->Clone ())->Draw ("5P");
  
    SaferDelete (&g);
  
    g = (TGAE*) tg_CMS_Zh_IAA_0_30_stat->Clone ();
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    SaferDelete (&g);

  
    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][3]->Clone ();
    RecenterGraph (g_syst);
    ResetXErrors (g_syst);
    SetConstantXErrors (g_syst, 0.060, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
  
    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (finalColors[3]);
    g_syst->SetLineColor (finalColors[3]);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (finalFillColors[3], 0.3);
  
    ((TGAE*) g_syst->Clone ())->Draw ("5P");

    SaferDelete (&g_syst);
  
    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][3]);
  
    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
  
    g_stat->SetMarkerStyle (kFullDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (finalColors[3]);
    g_stat->SetLineColor (finalColors[3]);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    g_stat->SetMarkerStyle (kOpenDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    SaferDelete (&g_stat);


    g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][2]->Clone ();
    RecenterGraph (g_syst);
    ResetXErrors (g_syst);
    SetConstantXErrors (g_syst, 0.060, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
  
    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (finalColors[2]);
    g_syst->SetLineColor (finalColors[2]);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (finalFillColors[2], 0.3);
  
    ((TGAE*) g_syst->Clone ())->Draw ("5P");

    SaferDelete (&g_syst);
  
    g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][2]);
  
    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
  
    g_stat->SetMarkerStyle (kFullSquare);
    g_stat->SetMarkerSize (1.8);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (finalColors[2]);
    g_stat->SetLineColor (finalColors[2]);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    g_stat->SetMarkerStyle (kOpenSquare);
    g_stat->SetMarkerSize (1.8);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    SaferDelete (&g_stat);
  
    myMarkerAndBoxAndLineText (0.33, 0.90-0.012, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], kFullDiamond,  2.4, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 0-10%", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.85-0.012, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], kFullSquare,  1.8, "ATLAS #it{p}_{T}^{Z} > 60 GeV, 10-30%", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.80-0.012, 1.4, 1001, kGray+2, 0.30, kBlack,  kFullCircle, 1.8, "CMS #it{p}_{T}^{Z} > 30 GeV, 0-30% (#it{preliminary})", 0.032);
  
    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);
  
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.24, 0.310, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.26, 0.260, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.26, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    //TPad* p_hamburglar = new TPad ("p_hamburglar", "", 0.65, 0.50, 0.90, 0.75);
    //p_hamburglar->Draw ();
    //p_hamburglar->cd ();

    //TImage* image = TImage::Open ("../DataThieving/CMS_Zh_IAA_HP2020/Hamburglar.jpg");
    //image->Draw ();


    c16->SaveAs ("../Plots/FinalPlots/iaa_pTch_cmsComp.pdf");
  }

} // end of macro

#endif

