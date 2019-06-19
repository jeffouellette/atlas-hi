#ifndef __Analysis_h__
#define __Analysis_h__

#include "Params.h"
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <AtlasUtils.h>

#include <TEfficiency.h>
#include <TClass.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>

#include <iostream>
#include <string>

using namespace atlashi;
using namespace std;

class Analysis {

  protected:
  string name;
  string directory;
  bool backgroundSubtracted = false;

  vector<TH1*> drawnHists;
  vector<TGraphAsymmErrors*> drawnGraphs;

  bool iaaCalculated = false;
  bool icpCalculated = false;

  TFile* trkEffFile = nullptr;
  TFile* histFile   = nullptr;
  bool histsLoaded  = false;
  bool histsScaled  = false;

  public:
  bool plotFill     = false; // whether to plot as filled (bar) graph or points w/ errors
  bool plotSignal   = true; // whether to plot background subtracted plots
  bool useAltMarker = false; // whether to plot as open markers (instead of closed)

  // Event info distributions (for reweighting)
  TH1D* h_PbPbFCal_weights = nullptr;
  TH1D* h_PbPbQ2_weights = nullptr;
  TH1D* h_PbPbVZ_weights = nullptr;
  TH1D* h_ppVZ_weights = nullptr;
  //TH3D* h_PbPb_event_reweights  = nullptr;
  //TH1D* h_pp_event_reweights    = nullptr;

  // Efficiencies
  TEfficiency***  h_trk_effs  = nullptr;
  TEfficiency**  h2_trk_effs  = nullptr;
  TH2D** h2_num_trk_effs = nullptr;
  TH2D** h2_den_trk_effs = nullptr;

  // Analysis checks
  TH1D*   h_fcal_et               = nullptr;
  //TH2D*   h_fcal_et_q2            = nullptr;
  TH1D*   h_fcal_et_reweighted    = nullptr;
  //TH2D*   h_fcal_et_q2_reweighted = nullptr;

  TH1D*   h_q2               = nullptr;
  TH1D*   h_q2_reweighted    = nullptr;
  TH1D*   h_PbPb_vz               = nullptr;
  TH1D*   h_PbPb_vz_reweighted    = nullptr;
  TH1D*   h_pp_vz               = nullptr;
  TH1D*   h_pp_vz_reweighted    = nullptr;

  TH1D*** h_z_phi         = nullptr;
  TH1D*** h_z_pt          = nullptr;
  TH2D*** h_z_y_phi       = nullptr;
  TH1D*** h_z_m           = nullptr;
  TH1D*** h_z_m_ratio     = nullptr;
  TH1D*** h_z_lepton_dphi = nullptr;
  TH1D*** h_lepton_pt     = nullptr;
  TH1D*** h_lepton_trk_pt = nullptr;
  TH1D*** h_trk_pt        = nullptr;
  TH2D*** h_lepton_trk_dr = nullptr;

  // Correlations plots
  TH2D*****   h_z_trk_pt_phi  = nullptr;
  TH1D******  h_z_trk_phi     = nullptr;
  
  // Physics plots
  TH2D*****   h_z_missing_pt      = nullptr;
  TH1D*****   h_z_missing_pt_avg  = nullptr;
  TH1D*****   h_z_missing_pt_int  = nullptr;
  TH1D******  h_z_trk_pt          = nullptr;
  //TH1D*****   h_z_trk_pt          = nullptr;
  TH1D****    h_z_counts          = nullptr;

  TH1D****** h_z_trk_pt_sub         = nullptr;
  TH1D****** h_z_trk_pt_sig_to_bkg  = nullptr;
  TH1D****** h_z_trk_pt_iaa         = nullptr;
  TH1D****** h_z_trk_pt_icp         = nullptr;

  Analysis () {
    name = "";
    directory = "";
    backgroundSubtracted = false;

    plotFill     = false;
    plotSignal   = true;
    useAltMarker = false;

    // Reweighting histograms
    h_PbPbFCal_weights = nullptr;
    h_PbPbQ2_weights = nullptr;
    h_PbPbVZ_weights = nullptr;
    h_ppVZ_weights = nullptr;
    //h_PbPb_event_reweights = nullptr;
    //h_pp_event_reweights   = nullptr;

    // Efficiencies
    h_trk_effs = Get2DArray <TEfficiency*> (numFinerCentBins, numEtaTrkBins); // iCent, iEta
    h2_trk_effs = Get1DArray <TEfficiency*> (numFinerCentBins); // iCent, iEta
    h2_num_trk_effs = Get1DArray <TH2D*> (numFinerCentBins); // iCent
    h2_den_trk_effs = Get1DArray <TH2D*> (numFinerCentBins);

    // Analysis checks
    h_fcal_et               = nullptr;
    //h_fcal_et_q2            = nullptr;
    h_fcal_et_reweighted    = nullptr;
    //h_fcal_et_q2_reweighted = nullptr;

    h_q2                  = nullptr;
    h_q2_reweighted       = nullptr;
    h_PbPb_vz             = nullptr;
    h_PbPb_vz_reweighted  = nullptr;
    h_pp_vz               = nullptr;
    h_pp_vz_reweighted    = nullptr;

    h_z_phi         = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_z_pt          = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_z_y_phi       = Get2DArray <TH2D*> (numCentBins, 3);             // iCent, iSpc
    h_z_m           = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_z_m_ratio     = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_z_lepton_dphi = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_lepton_pt     = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_lepton_trk_pt = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_trk_pt        = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    h_lepton_trk_dr = Get2DArray <TH2D*> (numCentBins, 3);             // iCent, iSpc

    // Correlations plots
    h_z_trk_pt_phi  = Get4DArray <TH2D*> (nPtZBins, nXZTrkBins, numCentBins, 3);             // iPtZ, iXZTrk, iCent, iSpc (0=ee, 1=mumu, 2=combined)
    h_z_trk_phi     = Get5DArray <TH1D*> (nPtTrkBins, nPtZBins, nXZTrkBins, numCentBins, 3); // iPtTrk, iPtZ, iXZTrk, iCent, iSpc
    
    // Physics plots
    h_z_missing_pt      = Get4DArray <TH2D*> (3, nPtZBins, numPhiTrkBins, numCentBins);           // iSpc, iPtZ, iPhi, iCent
    h_z_missing_pt_avg  = Get4DArray <TH1D*> (3, nPtZBins, numPhiTrkBins, numCentBins);           // iSpc, iPtZ, iPhi, iCent
    h_z_missing_pt_int  = Get4DArray <TH1D*> (3, nPtZBins, numPhiTrkBins, numCentBins);           // iSpc, iPtZ, iPhi, iCent
    h_z_trk_pt          = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins);  // iSpc, iPtZ, iXZTrk, iPhi, iCent
    h_z_counts          = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);                          // iSpc, iPtZ, iCent

    h_z_trk_pt_sub        = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent
    h_z_trk_pt_sig_to_bkg = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent
    h_z_trk_pt_iaa        = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent
    h_z_trk_pt_icp        = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent

  }

  protected:
  void LabelTrackingEfficiencies (const short iCent, const short iEta);
  void LabelCorrelations (const short iPtZ, const short iPtTrk, const short iXZTrk, const short iCent, const bool diffXZTrk);
  void LabelZMassSpectra (const short iSpc, const short iCent);
  void LabelTrkYield (const short iCent, const short iPhi);
  void LabelIAARatios (const short iCent, const short iPhi);
  void LabelICPRatios (const short iCent, const short iPhi);

  void GetDrawnObjects ();
  void GetMinAndMax (double &min, double &max, const bool log = false);
  void SetMinAndMax (double min, double max);

  public:
  string Name () { return name; }
  void SetName (string _name) { name = _name; }

  virtual void CreateHists ();
  virtual void CopyAnalysis (Analysis* a, const bool copyBkgs = false);
  virtual void CombineHists ();
  virtual void LoadHists ();
  virtual void SaveHists ();
  virtual void ScaleHists ();
  virtual void Execute ();
  virtual void SubtractBackground (Analysis* a = nullptr);

  virtual void LoadTrackingEfficiencies ();
  virtual double GetTrackingEfficiency (const float fcal_et, const float trk_pt, const float trk_eta, const bool isPbPb = true);

  void PrintZYields ();

  void Plot3DDist ();
  void PlotFCalDists (const char* plotTag = "data");
  void PlotQ2Dists (const char* plotTag = "data");
  void PlotVZDists (const char* plotTag = "data");

  void PlotCorrelations (const bool diffXZTrk = false, const short pSpc = 2);
  void PlotLeptonPtSpectra ();
  void PlotLeptonTrackPtSpectra ();
  void PlotLeptonTrackDR ();
  void PlotZPtSpectra ();
  void PlotZYPhiMap ();
  void PlotZMassSpectra ();
  void PlotZPhiYield ();
  void PlotZLeptonDPhi ();
  void PlotZMissingPt ();
  void PlotTrackingEfficiencies ();
  void PlotTrackingEfficiencies2D ();

  void CalculateZMassSpectraRatio (Analysis* a);
  void CalculateIAA ();
  void CalculateICP ();

  virtual void PlotTrkYield  (const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAARatios (const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotICPRatios (const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);

};


void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}


void Analysis :: GetDrawnObjects () {
  TList* primitives = gPad->GetListOfPrimitives ();
  drawnHists.clear ();
  drawnGraphs.clear ();
  for (int i = 0; i < primitives->GetSize (); i++) {
    TObject *obj = primitives->At (i);
    if (obj->IsA()->InheritsFrom (TH1::Class ())) {
      drawnHists.push_back ((TH1*)obj);
    }
    else if (obj->IsA()->InheritsFrom (TGraphAsymmErrors::Class ())) {
      drawnGraphs.push_back ((TGraphAsymmErrors*)obj);
    }
  }
}


void Analysis :: GetMinAndMax (double &min, double &max, const bool log) {
  for (TH1* h : drawnHists) {
    const double _max = log ? h->GetMaximum (0) : h->GetMaximum ();
    const double _min = log ? h->GetMinimum (0) : h->GetMinimum ();

    if (_max > max) max = _max;
    if (_min < min) min = _min;
  }
  for (TGraphAsymmErrors* g : drawnGraphs) {
    const double* ys = g->GetY ();
    for (int n = 0; n < g->GetN (); n++) {
      if (log && ys[n] <= 0)
        continue;
      if (ys[n] > max) max = ys[n];
      if (ys[n] < min) min = ys[n];
    }
  }
  return;
}


void Analysis :: SetMinAndMax (double min, double max) {
  for (TH1* h : drawnHists) {
    h->GetYaxis ()->SetRangeUser (min, max);
  }
  for (TGraphAsymmErrors* g : drawnGraphs) {
    g->GetYaxis ()->SetRangeUser (min, max);
  }
  gPad->Update ();
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: CreateHists () {
  h_fcal_et = new TH1D (Form ("h_fcal_et_%s", name.c_str ()), "", 300, 0, 6000); 
  h_fcal_et->Sumw2 ();
  //h_fcal_et_q2 = new TH2D (Form ("h_fcal_et_q2_%s", name.c_str ()), "", 300, 0, 6000, 150, 0, 300);
  //h_fcal_et_q2->Sumw2 ();
  h_fcal_et_reweighted = new TH1D (Form ("h_fcal_et_reweighted_%s", name.c_str ()), "", 300, 0, 6000);
  h_fcal_et_reweighted->Sumw2 ();
  //h_fcal_et_q2_reweighted = new TH2D (Form ("h_fcal_et_q2_reweighted_%s", name.c_str ()), "", 300, 0, 6000, 150, 0, 300);
  //h_fcal_et_q2_reweighted->Sumw2 ();

  h_q2 = new TH1D (Form ("h_q2_%s", name.c_str ()), "", 50, 0, 1);
  h_q2_reweighted = new TH1D (Form ("h_q2_reweighted_%s", name.c_str ()), "", 50, 0, 1);
  h_PbPb_vz = new TH1D (Form ("h_PbPb_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_PbPb_vz_reweighted = new TH1D (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz = new TH1D (Form ("h_pp_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz_reweighted = new TH1D (Form ("h_pp_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = new TH1D (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 80, 0, pi);
      h_z_pt[iCent][iSpc]           = new TH1D (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 300, 0, 300);
      h_z_y_phi[iCent][iSpc]        = new TH2D (Form ("h_z_y_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 20, -2.5, 2.5, 20, 0, 2*pi);
      h_z_m[iCent][iSpc]            = new TH1D (Form ("h_z_m_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 40, 76, 106);
      h_z_lepton_dphi[iCent][iSpc]  = new TH1D (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 45, 0, pi);
      h_lepton_pt[iCent][iSpc]      = new TH1D (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 250, 0, 250);
      h_lepton_trk_pt[iCent][iSpc]  = new TH1D (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 250, 0, 250);
      h_trk_pt[iCent][iSpc]         = new TH1D (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 100, 0, 100);
      h_lepton_trk_dr[iCent][iSpc]  = new TH2D (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 40, 0., 1., 40, 0., 1.);
      
      h_z_phi[iCent][iSpc]->Sumw2 ();
      h_z_pt[iCent][iSpc]->Sumw2 ();
      h_z_y_phi[iCent][iSpc]->Sumw2 ();
      h_z_m[iCent][iSpc]->Sumw2 ();
      h_z_lepton_dphi[iCent][iSpc]->Sumw2 ();
      h_lepton_pt[iCent][iSpc]->Sumw2 ();
      h_lepton_trk_pt[iCent][iSpc]->Sumw2 ();
      h_trk_pt[iCent][iSpc]->Sumw2 ();
      h_lepton_trk_dr[iCent][iSpc]->Sumw2 ();

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc] = new TH2D (Form ("h_z_trk_pt_phi_%s_iPtZ%i_iXZTrk%i_iCent%i_%s", spc, iPtZ, iXZTrk, iCent, name.c_str ()), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Sumw2 ();
        }
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent] = new TH2D (Form ("h_z_missing_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        }
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent] = new TH1D (Form ("h_z_trk_pt_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()), "", nPtTrkBins, ptTrkBins);
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Sumw2 ();
          }
        }
        h_z_counts[iSpc][iPtZ][iCent] = new TH1D (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", 1, 0, 1);
        h_z_counts[iSpc][iPtZ][iCent]->Sumw2 ();
      }
    }
  }

  histsLoaded = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms as clones from another analysis, where appropriate
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: CopyAnalysis (Analysis* a, const bool copyBkgs) {
  if (name == "")
    cout << "Warning in Analysis :: CopyAnalysis: name of analysis not set!" << endl;

  // Don't need to clone these histograms

  // Reweighting histograms
  h_PbPbFCal_weights = nullptr;
  h_PbPbQ2_weights = nullptr;
  h_PbPbVZ_weights = nullptr;
  h_ppVZ_weights = nullptr;
  //h_PbPb_event_reweights = a->h_PbPb_event_reweights;
  //h_pp_event_reweights   = a->h_pp_event_reweights;

  // Efficiencies
  h_trk_effs = a->h_trk_effs;
  h2_trk_effs = a->h2_trk_effs;
  h2_num_trk_effs = a->h2_num_trk_effs;
  h2_den_trk_effs = a->h2_den_trk_effs;

  // Should clone these histograms
  h_fcal_et = (TH1D*) a->h_fcal_et->Clone (Form ("h_fcal_et_%s", name.c_str ()));
  //h_fcal_et_q2 = (TH2D*) a->h_fcal_et_q2->Clone (Form ("h_fcal_et_q2_%s", name.c_str ()));
  h_fcal_et_reweighted = (TH1D*) a->h_fcal_et_reweighted->Clone (Form ("h_fcal_et_reweighted_%s", name.c_str ()));
  //h_fcal_et_q2_reweighted = (TH2D*) a->h_fcal_et_q2_reweighted->Clone (Form ("h_fcal_et_q2_reweighted_%s", name.c_str ()));

  h_q2 = (TH1D*) a->h_q2->Clone (Form ("h_q2_%s", name.c_str ()));
  h_q2_reweighted = (TH1D*) a->h_q2_reweighted->Clone (Form ("h_q2_reweighted_%s", name.c_str ()));
  h_PbPb_vz = (TH1D*) a->h_PbPb_vz->Clone (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted = (TH1D*) a->h_PbPb_vz_reweighted->Clone (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz = (TH1D*) a->h_pp_vz->Clone (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted = (TH1D*) a->h_pp_vz_reweighted->Clone (Form ("h_pp_vz_reweighted_%s", name.c_str ()));

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = (TH1D*) a->h_z_phi[iCent][iSpc]->Clone (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt[iCent][iSpc]           = (TH1D*) a->h_z_pt[iCent][iSpc]->Clone (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_y_phi[iCent][iSpc]        = (TH2D*) a->h_z_y_phi[iCent][iSpc]->Clone (Form ("h_z_y_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_m[iCent][iSpc]            = (TH1D*) a->h_z_m[iCent][iSpc]->Clone (Form ("h_z_m_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_z_m_ratio[iCent][iSpc])
        h_z_m_ratio[iCent][iSpc]    = (TH1D*) a->h_z_m_ratio[iCent][iSpc]->Clone (Form ("h_z_m_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      h_z_lepton_dphi[iCent][iSpc]  = (TH1D*) a->h_z_lepton_dphi[iCent][iSpc]->Clone (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_pt[iCent][iSpc]      = (TH1D*) a->h_lepton_pt[iCent][iSpc]->Clone (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_pt[iCent][iSpc]  = (TH1D*) a->h_lepton_trk_pt[iCent][iSpc]->Clone (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_trk_pt[iCent][iSpc]         = (TH1D*) a->h_trk_pt[iCent][iSpc]->Clone (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_dr[iCent][iSpc]  = (TH2D*) a->h_lepton_trk_dr[iCent][iSpc]->Clone (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc] = (TH2D*) a->h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Clone (Form ("h_z_trk_pt_phi_%s_iPtZ%i_iXZTrk%i_iCent%i_%s", spc, iPtZ, iXZTrk, iCent, name.c_str ()));
        } // end loop over iXZTrk
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent] = (TH2D*) a->h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_missing_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        } // end loop over iPhi
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent] = (TH1D*) a->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Clone (Form ("h_z_trk_pt_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()));
          } // end loop over iPhi
        } // end loop over iXZTrk
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_counts[iSpc][iPtZ][iCent]->Clone (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  if (copyBkgs) {
    if (a->backgroundSubtracted) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
            for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
              for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
                if (a->h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent])
                  h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Clone (Form ("h_z_trk_pt_subtracted_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()));
                  h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Clone (Form ("h_z_trk_pt_sig_to_bkg_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()));
              } // end loop over phi
            } // end loop over xztrk bins
          } // end loop over pT^Z bins
        } // end loop over cents
      } // end loop over species
      backgroundSubtracted = true;
    }

    if (a->iaaCalculated) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            for (short iCent = 1; iCent < numCentBins; iCent++) {
              if (a->h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent])
                h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent]->Clone (Form ("h_z_trk_pt_iaa_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, 0, iPhi, iCent, name.c_str ()));
            } // end loop over cents
          } // end loop over phi
        } // end loop over pT^Z bins
      } // end loop over species
      iaaCalculated = true;
    }

    if (a->icpCalculated) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            for (short iCent = 2; iCent < numCentBins; iCent++) {
              if (a->h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent])
                h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent]->Clone (Form ("h_z_trk_pt_icp_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, 0, iPhi, iCent, name.c_str ()));
            } // end loop over cents
          } // end loop over phi
        } // end loop over pT^Z bins
      } // end loop over species
      icpCalculated = true;
    }
  }

  histsLoaded = true;
  histsScaled = true;
  return;    
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LoadHists () {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");
  if (histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = (TH1D*) histFile->Get (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt[iCent][iSpc]           = (TH1D*) histFile->Get (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_y_phi[iCent][iSpc]        = (TH2D*) histFile->Get (Form ("h_z_y_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_m[iCent][iSpc]            = (TH1D*) histFile->Get (Form ("h_z_m_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_m_ratio[iCent][iSpc]      = (TH1D*) histFile->Get (Form ("h_z_m_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_lepton_dphi[iCent][iSpc]  = (TH1D*) histFile->Get (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_pt[iCent][iSpc]      = (TH1D*) histFile->Get (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_pt[iCent][iSpc]  = (TH1D*) histFile->Get (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_trk_pt[iCent][iSpc]         = (TH1D*) histFile->Get (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_dr[iCent][iSpc]  = (TH2D*) histFile->Get (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc] = (TH2D*) histFile->Get (Form ("h_z_trk_pt_phi_%s_iPtZ%i_iXZTrk%i_iCent%i_%s", spc, iPtZ, iXZTrk, iCent, name.c_str ()));
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc] = (TH1D*) histFile->Get (Form ("h_z_trk_phi_iPtTrk%i_iPtZ%i_iXZTrk%i_iCent%i_%s", iPtTrk, iPtZ, iXZTrk, iCent, name.c_str ()));
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()));
          }
        }
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent] = (TH2D*) histFile->Get (Form ("h_z_missing_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        }
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      }
    }
  }
  h_fcal_et               = (TH1D*) histFile->Get (Form ("h_fcal_et_%s", name.c_str ()));
  //h_fcal_et_q2            = (TH2D*) histFile->Get (Form ("h_fcal_et_q2_%s", name.c_str ()));
  h_fcal_et_reweighted    = (TH1D*) histFile->Get (Form ("h_fcal_et_reweighted_%s", name.c_str ()));
  //h_fcal_et_q2_reweighted = (TH2D*) histFile->Get (Form ("h_fcal_et_q2_reweighted_%s", name.c_str ()));

  h_q2            = (TH1D*) histFile->Get (Form ("h_q2_%s", name.c_str ()));
  h_q2_reweighted = (TH1D*) histFile->Get (Form ("h_q2_reweighted_%s", name.c_str ()));
  h_PbPb_vz            = (TH1D*) histFile->Get (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted = (TH1D*) histFile->Get (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz            = (TH1D*) histFile->Get (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted = (TH1D*) histFile->Get (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  
  histsLoaded = true;
  histsScaled = true;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: SaveHists () {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");
  if (!histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "recreate");
  histFile->cd ();
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      SafeWrite (h_z_phi[iCent][iSpc]);
      SafeWrite (h_z_pt[iCent][iSpc]);
      SafeWrite (h_z_y_phi[iCent][iSpc]);
      SafeWrite (h_z_m[iCent][iSpc]);
      SafeWrite (h_z_m_ratio[iCent][iSpc]);
      SafeWrite (h_z_lepton_dphi[iCent][iSpc]);
      SafeWrite (h_lepton_pt[iCent][iSpc]);
      SafeWrite (h_lepton_trk_pt[iCent][iSpc]);
      SafeWrite (h_trk_pt[iCent][iSpc]);
      SafeWrite (h_lepton_trk_dr[iCent][iSpc]);

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          SafeWrite (h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]);
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            SafeWrite (h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc]);
          }
        }
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          SafeWrite (h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]);
        }
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            SafeWrite (h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);
          }
        }
        SafeWrite (h_z_counts[iSpc][iPtZ][iCent]);
      }
    }
    
  }
  SafeWrite (h_fcal_et);
  //SafeWrite (h_fcal_et_q2);
  SafeWrite (h_fcal_et_reweighted);
  //SafeWrite (h_fcal_et_q2_reweighted);

  SafeWrite (h_q2);
  SafeWrite (h_q2_reweighted);
  SafeWrite (h_PbPb_vz);
  SafeWrite (h_PbPb_vz_reweighted);
  SafeWrite (h_pp_vz);
  SafeWrite (h_pp_vz_reweighted);
  
  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      h_z_phi[iCent][2]->Add (h_z_phi[iCent][iSpc]);
      h_z_pt[iCent][2]->Add (h_z_pt[iCent][iSpc]);
      h_z_y_phi[iCent][2]->Add (h_z_y_phi[iCent][iSpc]);
      h_z_m[iCent][2]->Add (h_z_m[iCent][iSpc]);
      h_lepton_pt[iCent][2]->Add (h_lepton_pt[iCent][iSpc]);
      h_lepton_trk_pt[iCent][2]->Add (h_lepton_pt[iCent][iSpc]);
      h_trk_pt[iCent][2]->Add (h_trk_pt[iCent][iSpc]);
      h_lepton_trk_dr[iCent][2]->Add (h_lepton_trk_dr[iCent][iSpc]);

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][2]->Add (h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]);
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            h_z_trk_pt[2][iPtZ][iXZTrk][iPhi][iCent]->Add (h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);
          } // end loop over phi
        }
        for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          h_z_missing_pt[2][iPtZ][iPhi][iCent]->Add (h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over phi
        h_z_counts[2][iPtZ][iCent]->Add (h_z_counts[iSpc][iPtZ][iCent]);
      } // end loop over pT^Z
    } // end loop over species
  } // end loop over centralities
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* h = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
            const double yieldNormFactor = countsHist->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
            //const double yieldNormFactorError = countsHist->GetBinError (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);

            //RescaleWithError (h, yieldNormFactor, yieldNormFactorError);
            if (yieldNormFactor > 0)
              h->Scale (1. / yieldNormFactor);
          } // end loop over phi
        } // end loop over xZTrk
      } // end loop over pT^Z

      double normFactor = 0;
      for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++)
        normFactor += h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1);

      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
            h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc] = h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->ProjectionX (Form ("h_z_trk_phi_iPtTrk%i_iPtZ%i_iXZTrk%i_iCent%i_%s", iPtTrk, iPtZ, iXZTrk, iCent, name.c_str ()), iPtTrk+1, iPtTrk+1);
            TH1D* h = h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc];
            h->Rebin (2);
            if (iPtTrk > 3)
              h->Rebin (2);
            if (iCent != 0)
              h->Rebin (2);
            if (normFactor > 0)
              h->Scale (1. / normFactor, "width");
          }
        }
      }

      h_lepton_pt[iCent][iSpc]->Rebin (5);
      if (normFactor > 0)
        h_lepton_pt[iCent][iSpc]->Scale (1. / normFactor, "width");

      h_lepton_trk_pt[iCent][iSpc]->Rebin (5);
      if (normFactor > 0)
        h_lepton_trk_pt[iCent][iSpc]->Scale (1. / normFactor, "width");

      h_trk_pt[iCent][iSpc]->Rebin (5);
      if (normFactor > 0)
        h_trk_pt[iCent][iSpc]->Scale (1. / normFactor, "width");

      if (h_z_m[iCent][iSpc]->GetMaximum () > 0)
        h_z_m[iCent][iSpc]->Scale (1. / (h_z_m[iCent][iSpc]->GetMaximum ()));

      h_z_phi[iCent][iSpc]->Rebin (8);
      if (h_z_phi[iCent][iSpc]->Integral () > 0)
        h_z_phi[iCent][iSpc]->Scale (1. / h_z_phi[iCent][iSpc]->Integral (), "width");
    } // end loop over centralities
  } // end loop over species

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
// Designed to be overloaded. The default here is for analyzing data.
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: Execute () {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_pt = 0, z_y = 0, z_phi = 0, z_m = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr, *l_trk_pt = nullptr, *l_trk_eta = nullptr, *l_trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("isEE",      &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",   &fcal_et);
    PbPbTree->SetBranchAddress ("q2",        &q2);
    PbPbTree->SetBranchAddress ("psi2",      &psi2);
    PbPbTree->SetBranchAddress ("vz",        &vz);
    PbPbTree->SetBranchAddress ("z_pt",      &z_pt);
    PbPbTree->SetBranchAddress ("z_y",       &z_y);
    PbPbTree->SetBranchAddress ("z_phi",     &z_phi);
    PbPbTree->SetBranchAddress ("z_m",       &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",     &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",    &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",    &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge", &l1_charge);
    PbPbTree->SetBranchAddress ("l2_pt",     &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",    &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",    &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge", &l2_charge);
    PbPbTree->SetBranchAddress ("l_trk_pt",  &l_trk_pt);
    PbPbTree->SetBranchAddress ("l_trk_eta", &l_trk_eta);
    PbPbTree->SetBranchAddress ("l_trk_phi", &l_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",      &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",    &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",   &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",   &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      //if (fabs (vz) > 1.5)
      //  continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      short iCent = 0;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      h_fcal_et->Fill (fcal_et);
      //h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);
      //h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

      h_q2->Fill (q2);
      h_q2_reweighted->Fill (q2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc]->Fill (z_y, InTwoPi (z_phi), event_weight);
      if (z_pt > zPtBins[1]) {
        h_z_m[iCent][iSpc]->Fill (z_m, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);

        for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
          h_lepton_trk_pt[iCent][iSpc]->Fill (l_trk_pt->at (iLTrk), event_weight);
        }
      }

      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        {
          float mindr = pi;
          //float phidiff = 0;
          float ptdiff = 0;
          for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
            const float dr = DeltaR (trk_eta->at (iTrk), l_trk_eta->at (iLTrk), trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            if (dr < mindr) {
              mindr = dr;
              ptdiff = 2. * fabs (trkpt - l_trk_pt->at (iLTrk)) / (trkpt + l_trk_pt->at (iLTrk));
              //phidiff = DeltaPhi (trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            }
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
          //h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, phidiff);
        }

        const float xZTrk = trkpt / z_pt;
        const short iXZTrk = GetiXZTrk (xZTrk);
        if (iXZTrk < 0 || iXZTrk > nXZTrkBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), true);
        if (trkEff == 0)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight / trkEff);

        // Add to missing pT (requires dphi in +/-pi/2 to +/-pi)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        bool awaySide = false;
        if (dphi > pi/2) {
          dphi = pi-dphi;
          awaySide = true;
        }

        short iPtTrk = 0;
        while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
          iPtTrk++;
        // start at the 1st phi bin and integrate outwards until the track is no longer contained 
        // e.g. so 7pi/8->pi is a subset of pi/2->pi
        short iPhi = 0;
        while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
          if (awaySide)
            trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi) / trkEff;
          else
            trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi) / trkEff;
          iPhi++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++)
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);
      } // end loop over tracks

      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }
    } // end loop over Pb+Pb tree
    cout << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("isEE",      &isEE);
    ppTree->SetBranchAddress ("vz",        &vz);
    ppTree->SetBranchAddress ("z_pt",      &z_pt);
    ppTree->SetBranchAddress ("z_y",       &z_y);
    ppTree->SetBranchAddress ("z_phi",     &z_phi);
    ppTree->SetBranchAddress ("z_m",       &z_m);
    ppTree->SetBranchAddress ("l1_pt",     &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",    &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",    &l1_phi);
    ppTree->SetBranchAddress ("l1_charge", &l1_charge);
    ppTree->SetBranchAddress ("l2_pt",     &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",    &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",    &l2_phi);
    ppTree->SetBranchAddress ("l2_charge", &l2_charge);
    ppTree->SetBranchAddress ("l_trk_pt",  &l_trk_pt);
    ppTree->SetBranchAddress ("l_trk_eta", &l_trk_eta);
    ppTree->SetBranchAddress ("l_trk_phi", &l_trk_phi);
    ppTree->SetBranchAddress ("ntrk",      &ntrk);
    ppTree->SetBranchAddress ("trk_pt",    &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",   &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",   &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      //if (fabs (vz) > 1.5)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc]->Fill (z_y, InTwoPi (z_phi), event_weight);
      if (z_pt > zPtBins[1]) {
        h_z_m[iCent][iSpc]->Fill (z_m, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);

        for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
          h_lepton_trk_pt[iCent][iSpc]->Fill (l_trk_pt->at (iLTrk), event_weight);
        }
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        {
          float mindr = pi;
          //float phidiff = 0;
          float ptdiff = 0;
          for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
            const float dr = DeltaR (trk_eta->at (iTrk), l_trk_eta->at (iLTrk), trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            if (dr < mindr) {
              mindr = dr;
              ptdiff = 2. * fabs (trkpt - l_trk_pt->at (iLTrk)) / (trkpt + l_trk_pt->at (iLTrk));
              //phidiff = DeltaPhi (trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            }
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
          //h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, phidiff);
        }

        const float xZTrk = trkpt / z_pt;
        const short iXZTrk = GetiXZTrk (xZTrk);
        if (iXZTrk < 0 || iXZTrk > nXZTrkBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), false);
        if (trkEff == 0)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight / trkEff);

        // Add to missing pT (requires dphi in -pi/2 to pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        bool awaySide = false;
        if (dphi > pi/2) {
          dphi = pi-dphi;
          awaySide = true;
        }

        short iPtTrk = 0;
        while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
          iPtTrk++;
        // start at the 1st phi bin and integrate outwards until the track is no longer contained 
        // e.g. so 7pi/8->pi is a subset of pi/2->pi
        short iPhi = 0;
        while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
          if (awaySide)
            trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi) / trkEff;
          else
            trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi) / trkEff;
          iPhi++;
        }
        
        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++)
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);
      } // end loop over tracks

      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }
    } // end loop over pp tree
    cout << endl;
  }

  CombineHists ();
  ScaleHists ();
  
  SaveHists ();
  //LoadHists ();

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  Delete2DArray (trkPtProj, numPhiBins, nPtTrkBins);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking efficiencies into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LoadTrackingEfficiencies () {
  SetupDirectories ("TrackingEfficiencies/", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;

  trkEffFile = new TFile (Form ("%s/trackingEfficiencies.root", rootPath.Data ()), "read");

  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    h2_num_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_matched_reco_tracks_iCent%i", iCent));
    h2_den_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_tracks_iCent%i", iCent));

    h2_trk_effs[iCent] = new TEfficiency (*(h2_num_trk_effs[iCent]), *(h2_den_trk_effs[iCent]));

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h_trk_effs[iCent][iEta] = (TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta));
    }
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking efficiency for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double Analysis :: GetTrackingEfficiency (const float fcal_et, const float trk_pt, const float trk_eta, const bool isPbPb) {
  short iCent = 0;
  if (isPbPb) {
    while (iCent < numFinerCentBins) {
      if (fcal_et < finerCentBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numFinerCentBins-1)
      return 0;
  }

  //short iEta = 0;
  //while (iEta < numEtaTrkBins) {
  //  if (fabs (trk_eta) < etaTrkBins[iEta+1])
  //    break;
  //  else
  //    iEta++;
  //}
  //if (iEta < 0 || iEta >= numEtaTrkBins)
  //  return 0;

  //if (iCent != 0)
  //  return 1; // no efficiencies for PbPb yet...

  //const double eff = h_trk_effs[iCent][iEta]->GetEfficiency (h_trk_effs[iCent][iEta]->FindFixBin (trk_pt));


  const double eff = h2_trk_effs[iCent]->GetEfficiency (h2_trk_effs[iCent]->FindFixBin (trk_eta, trk_pt));
  if (eff == 0)
    return 1;

  return eff;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - xZTrk 3d distribution
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: Plot3DDist () {
  TCanvas* c = new TCanvas ("c_z_trk_pt_phi", "", 800, 600);
  c->cd ();

  c->SetLogy ();

  h_z_trk_pt_phi[0][0][numCentBins-1][2]->GetXaxis ()->SetTitle ("#phi_{Z} - #phi_{Trk}");
  h_z_trk_pt_phi[0][0][numCentBins-1][2]->GetYaxis ()->SetTitle ("#it{p}_{T}^{ ch} / #it{p}_{T}^{ Z}");

  //h_z_trk_pt_phi[0][numCentBins-1][2]->RebinY (2);
  h_z_trk_pt_phi[0][0][numCentBins-1][2]->Draw ("lego2");

  c->SaveAs (Form ("%s/ZTrackCorr.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotFCalDists (const char* plotTag) {
  const char* canvasName = "c_fcal_et";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    c->cd ();
  }

  c->cd ();

  c->SetLogy ();

  h_fcal_et->SetLineColor (kBlack);

  h_fcal_et->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
  h_fcal_et->GetYaxis ()->SetTitle ("Counts");

  h_fcal_et->Draw (canvasExists ? "same hist" : "hist");

  if (strcmp (plotTag, "data") != 0) {
    h_fcal_et_reweighted->SetLineColor (kBlue);

    h_fcal_et_reweighted->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
    h_fcal_et_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_fcal_et_reweighted->Draw ("same hist");

    myText (0.65, 0.88, kBlack, "Unweighted", 0.04);
    myText (0.65, 0.81, kBlue, "Reweighted", 0.04);
  }

  myText (0.22, 0.28, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  //myText (0.22, 0.21, kBlack, "Z-tagged data", 0.04);
  myText (0.22, 0.21, kBlack, "Minimum bias", 0.04);

  c->SaveAs (Form ("%s/FCalDist_%s.pdf", plotPath.Data (), plotTag));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Q2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotQ2Dists (const char* plotTag) {
  const char* canvasName = "c_q2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    c->cd ();
  }

  c->cd ();

  c->SetLogy ();

  h_q2->SetLineColor (kBlack);

  h_q2->GetXaxis ()->SetTitle ("#it{q}_{2}");
  h_q2->GetYaxis ()->SetTitle ("Counts");

  h_q2->Draw (canvasExists ? "same hist" : "hist");

  if (strcmp (plotTag, "data") != 0) {
    h_q2_reweighted->SetLineColor (kBlue);

    h_q2_reweighted->GetXaxis ()->SetTitle ("#it{q}_{2}");
    h_q2_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_q2_reweighted->Draw ("same hist");

    myText (0.65, 0.72, kBlack, "Unweighted", 0.04);
    myText (0.65, 0.65, kBlue, "Reweighted", 0.04);
  }

  myText (0.65, 0.88, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  //myText (0.65, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.65, 0.81, kBlack, "Minimum bias", 0.04);

  c->SaveAs (Form ("%s/Q2Dist_%s.pdf", plotPath.Data (), plotTag));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot VZ distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotVZDists (const char* plotTag) {
  const char* canvasName = "c_vz";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1600, 600);
    gDirectory->Add (c);
    c->Divide (2, 1);
    c->cd ();
  }
  c->cd (1);

  gPad->SetLogy ();

  h_PbPb_vz->SetLineColor (kBlack);

  h_PbPb_vz->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
  h_PbPb_vz->GetYaxis ()->SetTitle ("Counts");

  h_PbPb_vz->Draw (canvasExists ? "same hist" : "hist");

  c->cd (2);

  gPad->SetLogy ();

  h_pp_vz->SetLineColor (kBlack);

  h_pp_vz->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
  h_pp_vz->GetYaxis ()->SetTitle ("Counts");

  h_pp_vz->Draw (canvasExists ? "same hist" : "hist");

  c->cd (1);

  if (strcmp (plotTag, "data") != 0) {
    h_PbPb_vz_reweighted->SetLineColor (kBlue);

    h_PbPb_vz_reweighted->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
    h_PbPb_vz_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_PbPb_vz_reweighted->Draw ("same hist");
  }

  myText (0.18, 0.88, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  //myText (0.18, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.18, 0.81, kBlack, "Minimum bias", 0.04);

  c->cd (2);

  if (strcmp (plotTag, "data") != 0) {
    h_pp_vz_reweighted->SetLineColor (kBlue);

    h_pp_vz_reweighted->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
    h_pp_vz_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_pp_vz_reweighted->Draw ("same hist");

    myText (0.75, 0.88, kBlack, "Unweighted", 0.04);
    myText (0.75, 0.81, kBlue, "Reweighted", 0.04);
  }

  myText (0.18, 0.88, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
  //myText (0.18, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.18, 0.81, kBlack, "Minimum bias", 0.04);

  c->SaveAs (Form ("%s/VZDist_%s.pdf", plotPath.Data (), plotTag));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - pTTrk 2d projections
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotCorrelations (const bool diffXZTrk, const short pSpc) {
  //TF1****** Fits = Get5DArray <TF1*> (nPtTrkBins, nPtZBins, nXZTrkBins, numCentBins, 3); // iPtTrk, iCent, iSpc
  const short pPtZ[3] = {1, 2, 3};
  const short pXZTrk[4] = {0, 1, 2, 3};
  const short pPtTrk[4] = {0, 1, 2, 3};

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    short nx = nXZTrkBins, ny = nPtTrkBins; // x = binning dimension, y = x-axis dimension
    short iXZTrk, iPtTrk;
    if (diffXZTrk) {
      nx = nPtTrkBins;
      ny = nXZTrkBins;
    }

    for (short ix = 0; ix < nx; ix++) {
      bool skipBin = true;
      if (diffXZTrk) { // true = binned in XZTrk, overall plot is 1 pTTrk bin
        iPtTrk = ix;
        for (short i = 0; i < sizeof (pPtTrk) / sizeof (pPtTrk[0]); i++)
          skipBin = skipBin && !(iPtTrk == pPtTrk[i]);
      } else {
        iXZTrk = ix;
        for (short i = 0; i < sizeof (pXZTrk) / sizeof (pXZTrk[0]); i++)
          skipBin = skipBin && !(iXZTrk == pXZTrk[i]);
      }
      if (skipBin)
        continue;

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        skipBin = true;
        for (short i = 0; i < sizeof (pPtZ) / sizeof (pPtZ[0]); i++) 
          skipBin = skipBin && !(iPtZ == pPtZ[i]);
        if (skipBin)
          continue;

        const char* canvasName = (diffXZTrk ? Form ("c_z_trk_phi_xztrk_iPtZ%i_iPtTrk%i_%s", iPtZ, iPtTrk, spc) : Form ("c_z_trk_phi_pttrk_iPtZ%i_iXZTrk%i_%s", iPtZ, iXZTrk, spc));
        const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
        TCanvas* c = nullptr;
        if (canvasExists)
          c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
        else {
          c = new TCanvas (canvasName, "", 1600, 1200);
          gDirectory->Add (c);
          c->cd ();
        //  c->Divide (2, 2);
        }

        for (short iCent = 0; iCent < 1; iCent++) {
          c->cd ();
          gPad->SetLogy ();
        //for (short iCent = 0; iCent < numCentBins; iCent++) {
        //  c->cd (iCent+1);
          GetDrawnObjects ();
          gPad->SetLogy ();

          gPad->SetTopMargin (0.01);
          gPad->SetBottomMargin (0.12);
          gPad->SetRightMargin (0.01);
          gPad->SetLeftMargin (0.12);

          double min = 1e30, max = 0;
          GetMinAndMax (min, max, true);
          for (short iy = 0; iy < ny; iy++) {
            skipBin = true;
            if (diffXZTrk) {
              iXZTrk = iy;
              for (short i = 0; i < sizeof (pXZTrk) / sizeof (pXZTrk[0]); i++)
                skipBin = skipBin && !(iXZTrk == pXZTrk[i]);
            } else {
              iPtTrk = iy;
              for (short i = 0; i < sizeof (pPtTrk) / sizeof (pPtTrk[0]); i++)
                skipBin = skipBin && !(iPtTrk == pPtTrk[i]);
            }
            if (skipBin)
              continue;
            TH1D* h = h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc];
            if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
            if (h->GetMaximum () > max) max = h->GetMaximum ();
          } // end loop over pT^trk bins
          min *= 0.5;
          max = max <= 0 ? 1 : 14*max;
          SetMinAndMax (min, max);

          if (plotFill) {
            for (short iy = 0; iy < ny; iy++) {
              skipBin = true;
              if (diffXZTrk) {
                iXZTrk = iy;
                for (short i = 0; i < sizeof (pXZTrk) / sizeof (pXZTrk[0]); i++)
                  skipBin = skipBin && !(iXZTrk == pXZTrk[i]);
              } else {
                iPtTrk = iy;
                for (short i = 0; i < sizeof (pPtTrk) / sizeof (pPtTrk[0]); i++)
                  skipBin = skipBin && !(iPtTrk == pPtTrk[i]);
              }
              if (skipBin)
                continue;

              TH1D* h = (TH1D*)h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc]->Clone ();

              h->GetYaxis ()->SetRangeUser (min, max);

              h->SetLineColor (fillColors[iy]);
              h->SetLineWidth (4);
              //h->SetFillColorAlpha (fillColors[iy], fillAlpha);
              h->SetMarkerSize (0);

              h->GetXaxis ()->SetTitle ("#Delta#phi");
              //h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{ch}/d#Delta#phi (\"ZYAM\")");
              h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{ch}/d#Delta#phi");
              h->GetXaxis ()->SetTitleOffset (0.6);
              h->GetYaxis ()->SetTitleOffset (0.8);
              h->GetXaxis ()->SetTitleSize (0.08);
              h->GetYaxis ()->SetTitleSize (0.06);
              h->GetXaxis ()->SetLabelSize (0.06);
              h->GetYaxis ()->SetLabelSize (0.06);

              h->DrawCopy (!canvasExists && iy == 0 ? "hist" : "same hist");
              //h->SetLineWidth (1);
              //h->Draw ("hist same");

              //TLine* line = new TLine (h->GetBinLowEdge (1), 0, h->GetBinLowEdge (h->GetNbinsX ()), 0);
              //line->Draw ("same");

              LabelCorrelations (iPtZ, iPtTrk, iXZTrk, iCent, diffXZTrk);
            } // end loop over iy
            gPad->RedrawAxis ();
          }
          else {
            for (short iy = 0; iy < ny; iy++) {
              skipBin = true;
              if (diffXZTrk) {
                iXZTrk = iy;
                for (short i = 0; i < sizeof (pXZTrk) / sizeof (pXZTrk[0]); i++)
                  skipBin = skipBin && !(iXZTrk == pXZTrk[i]);
              } else {
                iPtTrk = iy;
                for (short i = 0; i < sizeof (pPtTrk) / sizeof (pPtTrk[0]); i++)
                  skipBin = skipBin && !(iPtTrk == pPtTrk[i]);
              }
              if (skipBin)
                continue;

              TGraphAsymmErrors* g = make_graph (h_z_trk_phi[iPtTrk][iPtZ][iXZTrk][iCent][iSpc]);
              ResetXErrors (g);

              const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
              g->GetYaxis ()->SetRangeUser (min, max);

              g->SetMarkerStyle (markerStyle);
              g->SetLineColor (colors[iy]);
              g->SetMarkerColor (colors[iy]);

              g->Draw (!canvasExists && iy == 0 ? "AP" : "P");

              LabelCorrelations (iPtZ, iPtTrk, iXZTrk, iCent, diffXZTrk);
            } // end loop over iy
          }
        } // end loop over centrality

        c->SaveAs (diffXZTrk ? Form ("%s/dPhi_Distributions/dPhi_xZTrk_iPtZ%i_iPtTrk%i_%s.pdf", plotPath.Data (), iPtZ, iPtTrk, spc) : Form ("%s/dPhi_Distributions/dPhi_pTtrk_iPtZ%i_iXZTrk%i_%s.pdf", plotPath.Data (), iPtZ, iXZTrk, spc));
      } // end loop over ix
    } // end loop over pT^Z
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track dPhi distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LabelCorrelations (const short iPtZ, const short iPtTrk, const short iXZTrk, const short iCent, const bool diffXZTrk) {
  if (iCent == 0) {
    myText (0.8, 0.20, kBlack, "#it{pp}", 0.06);
    //myText (0.16, 0.93, kBlack, "Truth", 0.05);
    //myText (0.30, 0.93, kBlack, "Reco.", 0.05);
    myText (0.16, 0.93, kBlack, "Data", 0.05);
    myText (0.30, 0.93, kBlack, "Minbias", 0.05);
    TVirtualPad* cPad = gPad; // store current pad
    TBox* b = nullptr;
    if (diffXZTrk) {
      const float pt_lo = xZTrkBins[iXZTrk];
      const float pt_hi = xZTrkBins[iXZTrk+1];
      if (iXZTrk == 0)
        myText (0.4, 0.93, kBlack, "#it{x}_{Z}^{trk}", 0.05);
      myText (0.4, 0.87-0.065*iXZTrk, colors[iXZTrk], Form ("(%.1f, %.1f)", pt_lo, pt_hi), 0.05);
      myText (0.4, 0.87-0.065*iXZTrk, colors[iXZTrk], Form ("(%.1f, %.1f)", pt_lo, pt_hi), 0.05);

      b = TBoxNDC (0.33-0.024, 0.87-0.065*iXZTrk-0.016, 0.33+0.024, 0.87-0.065*iXZTrk+0.016);
      b->SetFillColorAlpha (fillColors[iXZTrk], fillAlpha);
      myMarkerText (0.23, 0.87-0.065*iXZTrk, colors[iXZTrk], kFullCircle, "", 2, 0.05);
    } else {
      const float pt_lo = ptTrkBins[iPtTrk];
      const float pt_hi = ptTrkBins[iPtTrk+1];
      if (iPtTrk == 0)
        myText (0.4, 0.93, kBlack, "#it{p}_{T}^{ ch} [GeV]", 0.05);
      myText (0.4, 0.87-0.065*iPtTrk, colors[iPtTrk], Form ("(%.1f, %.1f)", pt_lo, pt_hi), 0.05);
      myText (0.4, 0.87-0.065*iPtTrk, colors[iPtTrk], Form ("(%.1f, %.1f)", pt_lo, pt_hi), 0.05);

      b = TBoxNDC (0.33-0.024, 0.87-0.065*iPtTrk-0.016, 0.33+0.024, 0.87-0.065*iPtTrk+0.016);
      b->SetFillColorAlpha (fillColors[iPtTrk], fillAlpha);
      myMarkerText (0.23, 0.87-0.065*iPtTrk, colors[iPtTrk], kFullCircle, "", 2, 0.05);
    }
    b->Draw ("l");
    cPad->cd ();
  }
  else {
    myText (0.80, 0.20, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
  }
  //if (iCent == 1) {
    myText (0.65, 0.93, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    if (!diffXZTrk)
      myText (0.65, 0.87, kBlack, Form ("%g < #it{p}_{T}^{ ch} / #it{p}_{T}^{Z} < %g", xZTrkBins[iXZTrk], xZTrkBins[iXZTrk+1]), 0.05);
    else
      myText (0.65, 0.87, kBlack, Form ("%.1f < #it{p}_{T}^{ ch} < %.1f", ptTrkBins[iPtTrk], ptTrkBins[iPtTrk+1]), 0.05);
  //}
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot lepton Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotLeptonPtSpectra () {
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* canvasName = Form ("c_%s_pt", iSpc == 0 ? "electron" : (iSpc == 1 ? "muon" : "lepton"));
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 600, 600);
      gDirectory->Add (c);
      c->cd ();
    }

    const char* spc = iSpc == 0 ? "e" : "#mu";

    c->cd (iSpc+1);
    gPad->SetLogy ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TH1D* h = (TH1D*)h_lepton_pt[iCent][iSpc]->Clone ();
      //h->Rebin (5);
      //h->Scale (1./h_z_counts[iCent][iSpc]->Integral (), "width");

      h->GetXaxis ()->SetTitle (Form ("#it{p}_{T}^{ %s} [GeV]", spc));
      h->GetYaxis ()->SetTitle (Form ("1/N_{Z#rightarrow%s%s} dN_{%s}/d#it{p}_{T} [GeV^{-1}]", spc, spc, spc));

      h->SetLineColor (colors[iCent]);
      h->SetMarkerColor (colors[iCent]);
      h->SetMarkerStyle (kFullCircle);
      h->SetMarkerSize (0.75);

      h->Draw (iCent == 0 ? "e1" : "e1 same");
    }
    myText (0.76, 0.88, colors[0], "#it{pp}", 0.04);
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      myText (0.76, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
    }
    c->SaveAs (Form ("%s/%sPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : "Muon"));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot track Pt spectra for each lepton species
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotLeptonTrackPtSpectra () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* canvasName = Form ("c_%s_trk_pt", iSpc == 0 ? "electron" : (iSpc == 1 ? "muon" : "lepton"));
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 800, 600);
      gDirectory->Add (c);
    }
    c->cd ();
    c->SetLogy ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      //TH1D* h = (TH1D*)h_lepton_trk_pt[iCent][iSpc]->Clone (Form ("h_lepton_trk_pt_iSpc%i_iCent%i", iSpc, iCent));
      TH1D* h = h_lepton_trk_pt[iCent][iSpc];

      //h->Rebin (5);
      //h->Scale (1./h_z_pt[iCent][iSpc]->Integral (h_z_pt[iCent][iSpc]->GetXaxis ()->FindBin (5), h_z_pt[iCent][iSpc]->GetNbinsX ()), "width");
      //cout << iSpc << ", " << iCent << ", " << h_z_pt[iCent][iSpc]->Integral (h_z_pt[iCent][iSpc]->GetXaxis ()->FindBin (5), h_z_pt[iCent][iSpc]->GetNbinsX ()) << endl;

      h->GetYaxis ()->SetRangeUser (6e-6, 450);

      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
      if (iSpc == 0) {
        h->GetYaxis ()->SetTitle ("1/N_{Z#rightarrowee} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }
      else if (iSpc == 1) {
        h->GetYaxis ()->SetTitle ("1/N_{Z#rightarrow#mu#mu} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }
      else {
        h->GetYaxis ()->SetTitle ("1/N_{Z#rightarrowll} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }

      h->SetLineColor (colors[iCent]);
      h->SetMarkerColor (colors[iCent]);
      h->SetMarkerStyle (kFullCircle);
      h->SetMarkerSize (0.75);

      h->Draw (iCent == 0 ? "e1" : "e1 same");
    }
    gPad->RedrawAxis ();

    myText (0.66, 0.82, colors[0], "#it{pp}", 0.04);
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      myText (0.66, 0.82-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
    }

    myText (0.66, 0.88, kBlack, iSpc == 0 ? "Z#rightarrowee Events" : (iSpc == 1 ?"Z#rightarrow#mu#mu Events" : "Z#rightarrowll Events"), 0.04);
    c->SaveAs (Form ("%s/%sTrackPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Comb")));
  }

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dR between leptons and tracks
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotLeptonTrackDR () {
  
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      const char* canvasName = Form ("c_%s_trk_dr", iSpc == 0 ? "electron" : (iSpc == 1 ? "muon" : "lepton"));
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 800, 600);
        gDirectory->Add (c);
        FormatTH2Canvas (c, true);
        c->SetLogz ();
      }
      c->cd ();

      TH2D* h2 = h_lepton_trk_dr[iCent][iSpc];
      h2->RebinX (2);
      h2->RebinY (2);
      h2->GetXaxis ()->SetTitle (Form ("min (#DeltaR (track, %s))", iSpc == 0 ? "electrons" : (iSpc == 1 ? "muons" : "leptons")));
      h2->GetYaxis ()->SetTitle ("|#Delta#it{p}_{T}| / <#it{p}_{T}>");
      //h2->GetYaxis ()->SetTitle ("#Delta#phi");
      h2->GetZaxis ()->SetTitle ("Counts");

      h2->GetYaxis ()->SetTitleOffset (1.1);

      h2->Draw ("colz");

      myText (0.56, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);

      if (iCent != 0) {
        myText (0.56, 0.82, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
        myText (0.56, 0.76, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
      }
      else
        myText (0.56, 0.82, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);

      c->SaveAs (Form ("%s/%sTrackDist_iCent%i.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Lepton"), iCent));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotZPtSpectra () {
  const char* canvasName = "c_z_pt";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
  }
  c->cd ();
  c->SetLogy ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TH1D* h = (TH1D*)h_z_pt[iCent][2]->Clone ();
    //TH1D* h = (TH1D*)h_z_pt[iCent][1]->Clone ();
    //h->Add (h_z_pt[0][iCent][1]);
    h->Rebin (5);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");
    h->GetYaxis ()->SetTitle ("N_{Z}");

    //TH1D* integral = new TH1D (Form ("h_z_pt_integral_iCent%i", iCent), "", 300, 0, 300);
    //integral->Rebin (5);

    //if (iCent == 0)
    //  cout << "Centrality bin: pp" << endl;
    //else
    //  cout << "Centrality bin: " << Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]) << endl;

    //for (int ix = 1; ix <= integral->GetNbinsX (); ix++) {
    //  double content = 0, varsq = 0;
    //  content = h->IntegralAndError (ix, h->GetNbinsX (), varsq);
    //  integral->SetBinContent (ix, content);
    //  integral->SetBinError (ix, sqrt (varsq));
    //}

    //cout << "\t#Z's >= 0 GeV: " << integral->GetBinContent (1) << endl;
    //cout << "\t#Z's >= 5 GeV: " << integral->GetBinContent (6) << endl;
    //cout << "\t#Z's >= 25 GeV: " << integral->GetBinContent (26) << endl;

    h->SetLineColor (colors[iCent]);
    h->SetMarkerColor (colors[iCent]);
    h->SetMarkerStyle (kFullCircle);
    h->SetMarkerSize (0.5);

    //integral->SetLineColor (colors[iCent]);
    //integral->SetMarkerColor (colors[iCent]);
    //integral->SetMarkerStyle (kOpenCircle);
    //integral->SetMarkerSize (0.5);

    //integral->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");

    //if (iCent == 0)
    //  integral->Draw ("e1");
    //else
    //  integral->Draw ("same e1");
    h->Draw (iCent == 0 ? "e1" : "same e1");
  }
  gPad->RedrawAxis ();

  //myMarkerText (0.75, 0.88, kBlack, kFullCircle, "n_{Z}(#it{p}_{T})", 1.25, 0.04);
  //myMarkerText (0.75, 0.80, kBlack, kOpenCircle, "N_{Z}(#it{p}_{T}) =  #int_{#it{p}_{T}}^{#infty}n(x)dx", 1.25, 0.04);
  myText (0.74, 0.88, colors[0], "#it{pp}", 0.04);
  for (short iCent = 1; iCent < numCentBins; iCent++) {
    myText (0.74, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
  }

  c->SaveAs (Form ("%s/z_pt_spectrum.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints yield of Z's that meet the event selection criteria in each centrality
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PrintZYields () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TH1D* h = (TH1D*)h_z_pt[iCent][iSpc];
      cout << "Lower integration bound: " << h->GetBinLowEdge (h->FindBin (zPtBins[2])) << endl;

      float yield = h->Integral (h->FindBin (zPtBins[2]), h->GetNbinsX ());

      const char* spc = (iSpc == 0 ? "Z->ee" : (iSpc == 1 ? "Z->mumu" : "Z->ll"));
      if (iCent == 0) 
        cout << "pp " << spc << " # Z's > 25 GeV  =  " << yield << endl;
      else
        cout << Form ("Pb+Pb %i-%i%% ", (int)centCuts[iCent], (int)centCuts[iCent-1]) << spc << " # Z's > 25 GeV  =  " << yield << endl;
    }
  }
}





////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the eta-phi distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotZYPhiMap () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      const char* canvasName = Form ("c_z_y_phi_%s_iCent%i", spc, iCent);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 800, 600);
        FormatTH2Canvas (c, true);
        gDirectory->Add (c);
      }
      c->cd ();

      TH2D* h = h_z_y_phi[iCent][iSpc];

      h->GetXaxis ()->SetTitle ("y");
      h->GetYaxis ()->SetTitle ("#phi");
      h->GetZaxis ()->SetTitle ("Counts");

      h->GetXaxis ()->SetTitleOffset (1.4);
      h->GetYaxis ()->SetTitleOffset (1.1);

      h->Draw ("colz");

      myText (0.18, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);

      const char* spcLabel = iSpc == 0 ? "Z #rightarrow e^{+}e^{-} Events" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-} Events" : "Z #rightarrow l^{+}l^{-} Events");
      myText (0.62, 0.90, kBlack, spcLabel, 0.04);
      myText (0.62, 0.84, kBlack, "#it{p}_{T}^{Z} > 25 GeV", 0.04);
      if (iCent == 0) {
        myText (0.18, 0.84, kBlack, Form ("#it{pp}, #sqrt{s} = 5.02 TeV"), 0.04);
      }
      else
        myText (0.18, 0.84, kBlack, Form ("Pb+Pb %i-%i%%, #sqrt{s_{NN}} = 5.02 TeV", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);

      c->SaveAs (Form ("%s/ZYPhiDists/z%s_y_phi_iCent%i.pdf", plotPath.Data (), spc, iCent));
    }
  }

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates ratio of Z mass spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: CalculateZMassSpectraRatio (Analysis* a) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      h_z_m_ratio[iCent][iSpc] = (TH1D*) h_z_m[iCent][iSpc]->Clone (Form ("h_z_m_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_m_ratio[iCent][iSpc]->Divide (a->h_z_m[iCent][iSpc]);
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z mass spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotZMassSpectra () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
      const char* canvasName = Form ("c_z_m_%s_iCent%i", spc, iCent);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      TPad* uPad = nullptr, *dPad = nullptr;
      if (canvasExists) {
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName)); 
        uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
        dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
      }
      else {
        c = new TCanvas (canvasName, "", 800, 800);
        c->cd ();
        uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
        dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
        uPad->SetBottomMargin (0);
        dPad->SetTopMargin (0);
        dPad->SetBottomMargin (0.25);
        uPad->Draw ();
        dPad->Draw ();
        gDirectory->Add (c);
        gDirectory->Add (uPad);
        gDirectory->Add (dPad);
      }

      uPad->cd ();
      TH1D* h = h_z_m[iCent][iSpc];
      if (plotFill) {
        h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);
        h->GetYaxis ()->SetRangeUser (0, 1.3);

        h->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
        h->GetYaxis ()->SetTitle ("Arb. Units");
        h->GetXaxis ()->SetTitleSize (0.04/0.6);
        h->GetYaxis ()->SetTitleSize (0.04/0.6);
        h->GetXaxis ()->SetLabelSize (0.04/0.6);
        h->GetYaxis ()->SetLabelSize (0.04/0.6);
        h->GetXaxis ()->SetTitleOffset (1.5*0.6);
        h->GetYaxis ()->SetTitleOffset (1.5*0.6);

        h->DrawCopy (!canvasExists ? "bar" : "bar same");
        h->SetLineWidth (1);
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      }
      else {
        TGraphAsymmErrors* g = make_graph (h);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);
        g->GetYaxis ()->SetRangeUser (0, 1.3);

        g->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
        g->GetYaxis ()->SetTitle ("Arb. Units");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (!canvasExists ? "AP" : "P");
      }
      LabelZMassSpectra (iSpc, iCent);
      

      dPad->cd ();
      h = h_z_m_ratio[iCent][iSpc];
      if (h) {
        TGraphAsymmErrors* g = make_graph (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);
        g->GetYaxis ()->SetRangeUser (0.5, 1.5);

        g->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
        g->GetYaxis ()->SetTitle ("Data / MC");
        g->GetXaxis ()->SetTitleSize (0.04/0.4);
        g->GetYaxis ()->SetTitleSize (0.04/0.4);
        g->GetXaxis ()->SetLabelSize (0.04/0.4);
        g->GetYaxis ()->SetLabelSize (0.04/0.4);
        g->GetXaxis ()->SetTitleOffset (2.5*0.4);
        g->GetYaxis ()->SetTitleOffset (1.5*0.4);
        g->GetYaxis ()->CenterTitle ();
        g->Draw (!canvasExists ? "AP" : "P");
      }
      else {
        cout << "Warning in Analysis :: PlotZMassSpectra: Z mass spectra ratio not stored, needs to be calculated!" << endl;
      }

      c->SaveAs (Form ("%s/ZMassSpectra/z%s_mass_spectrum_iCent%i.pdf", plotPath.Data (), spc, iCent));
    }
  }
}




void Analysis :: LabelZMassSpectra (const short iSpc, const short iCent) {
  myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/0.6);
  const char* spc = iSpc == 0 ? "Z #rightarrow e^{+}e^{-} Events" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-} Events" : "Z #rightarrow l^{+}l^{-} Events");
  myText (0.66, 0.85, kBlack, spc, 0.04/0.6);
  if (iCent == 0) {
    myText (0.22, 0.75, colors[0], Form ("#it{pp}, #sqrt{s} = 5.02 TeV"), 0.04/0.6);
  }
  else
    myText (0.22, 0.75, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotZPhiYield () {
  const char* canvasName = "c_z_phi";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
  }
  c->cd ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    TH1D* h = h_z_phi[iCent][2];
    //h->Rebin (8);
    //h->Scale (1. / h->Integral (), "width");

    TF1* fit = new TF1 ("fit", "[0]+[1]*cos(x)+[2]*cos(2*x)", -pi/2, 3*pi/2);
    h->Fit (fit, "RN0");
    delete fit;
    fit = nullptr;

    TGraphAsymmErrors* g = make_graph (h);
    deltaize (g, (1.5-iCent)*0.02, false);

    g->GetXaxis ()->SetTitle ("2|#phi_{Z} - #Psi_{2}|");
    g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/d#Delta#phi");

    g->GetYaxis ()->SetRangeUser (0.1, 0.7);

    g->SetLineColor (colors[iCent]);
    g->SetMarkerColor (colors[iCent]);
    g->SetMarkerSize (0.75);
    if (iCent == 1)
      g->Draw ("ap");
    else
      g->Draw ("p");

  } // end loop over cents

  //myText (0.66, 0.88, colors[0], "#it{pp}", 0.04);
  myText (0.25, 0.88, colors[0], "#bf{#it{ATLAS Internal}}", 0.04);
  myText (0.72, 0.88, colors[0], "Pb+Pb, #sqrt{s}_{NN} = 5.02 TeV", 0.04);
  for (short iCent = 1; iCent < numCentBins; iCent++) {
    myText (0.72, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
  }

  c->SaveAs (Form ("%s/ZPhiYields.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotZLeptonDPhi () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    const char* canvasName = Form ("c_z_lepton_dphi_iCent%i", iCent);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 600, 600);
      gDirectory->Add (c);
      c->cd ();
    }
    c->cd ();

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      TH1D* h = h_z_lepton_dphi[iCent][iSpc];

      h->SetLineColor (colors[iSpc]);
      h->SetMarkerColor (colors[iSpc]);

      h->GetXaxis ()->SetTitle ("#Delta#phi");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw (canvasExists || iSpc != 0 ? "e1 same" : "e1");
    }
    c->SaveAs (Form ("%s/ZLeptonDPhi_iCent%i.pdf", plotPath.Data (), iCent));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotZMissingPt () {
  const char* canvasName = "c_z_missing_pt";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600*numPhiTrkBins, 500);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (2, 1);
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        c->cd (iPhi+1);
        gPad->SetLogx ();
        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TH1D* h = new TH1D (Form ("h_z_missing_pt_avg_%s_iPtZ%i_iPhi%i_iCent%i", spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          TH1D* integral = new TH1D (Form ("h_z_missing_pt_int_%s_iPtZ%i_iPhi%i_iCent%i", spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          h->Sumw2 ();
          integral->Sumw2 ();
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            TH1D* px = h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->ProjectionX ("_px", iPtTrk+1, iPtTrk+1);
            TF1* fit = new TF1 ("fit", "gaus(0)", zMissingPtBins[0], zMissingPtBins[numZMissingPtBins]);
            px->Fit (fit, "RN0Q");
            h->SetBinContent (iPtTrk+1, fit->GetParameter (1));
            h->SetBinError (iPtTrk+1, fit->GetParError (1));
            if (px) delete px;

            px = h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->ProjectionX ("_px", 0, iPtTrk+1);
            px->Fit (fit, "RN0Q");
            integral->SetBinContent (iPtTrk+1, fit->GetParameter (1));
            integral->SetBinError (iPtTrk+1, fit->GetParError (1));
            if (px) delete px;

            if (fit) delete fit;
          }
          h_z_missing_pt_avg[iSpc][iPtZ][iPhi][iCent] = h;
          h_z_missing_pt_int[iSpc][iPtZ][iPhi][iCent] = integral;

          integral->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch}  [GeV]");
          integral->GetYaxis ()->SetTitle ("<#it{p}_{T}^{ ||}>  [GeV]");
          integral->GetXaxis ()->SetTitleSize (0.06);
          integral->GetYaxis ()->SetTitleSize (0.06);
          integral->GetXaxis ()->SetLabelSize (0.06);
          integral->GetYaxis ()->SetLabelSize (0.06);
          integral->GetXaxis ()->SetTitleOffset (1.2);
          integral->GetYaxis ()->SetTitleOffset (1.2);

          integral->GetYaxis ()->SetRangeUser (-32, 16);

          integral->SetFillColorAlpha (fillColors[iCent], fillAlpha);
          integral->SetLineColor (kBlack);
          integral->SetMarkerSize (0);
          integral->SetLineWidth (0);
          integral->DrawCopy (iCent == numCentBins-1 ? "b" : "b same");
          integral->SetLineWidth (1);
          integral->Draw ("hist same");
        }

        TLine* l = new TLine (ptTrkBins[0], 0, ptTrkBins[nPtTrkBins], 0);
        l->Draw ("same");

        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TGraphAsymmErrors* g = make_graph (h_z_missing_pt_avg[iSpc][iPtZ][iPhi][iCent]);
          RecenterGraph (g);
          ResetXErrors (g);
          deltaize (g, (1.5-iCent)*0.04, false);

          g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch}  [GeV]");
          g->GetYaxis ()->SetTitle ("<#it{p}_{T}^{ ||}>  [GeV]");

          g->SetLineColor (colors[iCent]);
          g->SetMarkerColor (colors[iCent]);
          g->SetMarkerSize (0.75);
          g->Draw ("P");
        } // end loop over cents
        if (iPhi == 0) {
          myText (0.25, 0.88, kBlack, "0 < #Delta#phi < #pi/8 or 7#pi/8 < #Delta#phi < #pi", 0.05);
          myText (0.25, 0.81, kBlack, Form ("%g < #it{p}_{T}^{ Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
        }
        else if (iPhi == 1)
          myText (0.25, 0.88, kBlack, "0 < #Delta#phi < #pi", 0.05);
      } // end loop over directions

      myText (0.66, 0.88, colors[0], "#it{pp}", 0.05);
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        myText (0.66, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
      }

      c->SaveAs (Form ("%s/ZMissingPt.pdf", plotPath.Data ()));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotTrackingEfficiencies () {
  const char* canvasName = "c_trk_effs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 812);
    gDirectory->Add (c);
    c->cd ();
    //c->Divide (4, 3);
  }
  c->cd ();

  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    //c->cd (iCent+1);
    if (iCent > 0) continue;

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TEfficiency* eff = h_trk_effs[iCent][iEta];

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);

      eff->SetTitle (";#it{p}_{T} [GeV];Reco. Eff.");

      eff->Draw (iEta == 0 ? "APL" : "LP same");

      gPad->Update ();

      eff->GetPaintedGraph ()->GetXaxis ()->SetRangeUser (0.5, 80);
      eff->GetPaintedGraph ()->GetYaxis ()->SetRangeUser (0.5, 1.08);

      LabelTrackingEfficiencies (iCent, iEta);
    }
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotTrackingEfficiencies2D () {
  const char* canvasName = "c_trk_effs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 812);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    //c->Divide (4, 3);
  }
  c->cd ();

  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    //c->cd (iCent+1);
    if (iCent > 0) continue;

    TEfficiency* eff = h2_trk_effs[iCent];

    eff->Draw ("colz");
    eff->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track efficiency plots
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LabelTrackingEfficiencies (const short iCent, const short iEta) {
  if (iCent == 0)
    myText (0.22, 0.88, kBlack, "#it{pp}", 0.1);
  else
    myText (0.22, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.1);

  if (iCent == 0) {
  //  myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
    myMarkerText (0.5, 0.54-0.08*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 3, 0.08);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Subtracts default (HM) background from track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: SubtractBackground (Analysis* a) {
  if (backgroundSubtracted)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* h = new TH1D (Form ("h_z_trk_pt_sub_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()), "", nPtTrkBins, ptTrkBins);
            h->Sumw2 ();

            h->Add (h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);
            TH1D* sub = nullptr;
            if (a != nullptr) {
              //cout << "Info in SubtractBackground: Using MBM background subtraction method." << endl;
              sub = h_z_trk_pt[iSpc][iPtZ][iXZTrk][0][iCent];
            }
            else {
              //cout << "Info in SubtractBackground: Using HM background subtraction method." << endl;
              sub = a->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            }
            if (sub == nullptr)
              cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;
            else
              h->Add (sub, -1);

            h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent] = h;

            h = new TH1D (Form ("h_z_trk_pt_sigToBkg_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name.c_str ()), "", nPtTrkBins, ptTrkBins);
            h->Sumw2 ();

            if (sub == nullptr)
              cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;
            else {
              h->Add (h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);
              h->Add (h_z_trk_pt[iSpc][iPtZ][iXZTrk][0][iCent], -1);
              h->Divide (h_z_trk_pt[iSpc][iPtZ][iXZTrk][0][iCent]);
            }

            h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent] = h;
          } // end loop over phi
        } // end loop over xZTrk
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  backgroundSubtracted = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTtrk binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotTrkYield (const bool plotAsSystematic, const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (2*padRatio+1.0);
  const double mPadY = padRatio / (2*padRatio+1.0);
  const double uPadY = 1.0 - mPadY - dPadY;
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_TrkYield_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 400*numCentBins, 1000);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* topPadName = Form ("p_top_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
        const char* middlePadName = Form ("p_middle_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
        const char* bottomPadName = Form ("p_bottom_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);

        TPad* topPad = nullptr, *middlePad = nullptr, *bottomPad = nullptr;
        if (!canvasExists) {
          topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY+mPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
          middlePad = new TPad (middlePadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, dPadY+mPadY);
          bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

          gDirectory->Add (topPad);
          gDirectory->Add (middlePad);
          gDirectory->Add (bottomPad);

          topPad->SetTopMargin (0.04);
          topPad->SetBottomMargin (0);
          topPad->SetLeftMargin (0.17);
          topPad->SetRightMargin (0.06);
          middlePad->SetTopMargin (0);
          middlePad->SetBottomMargin (0);
          middlePad->SetLeftMargin (0.17);
          middlePad->SetRightMargin (0.06);
          bottomPad->SetTopMargin (0);
          bottomPad->SetBottomMargin (0.20);
          bottomPad->SetLeftMargin (0.17);
          bottomPad->SetRightMargin (0.06);
          topPad->Draw ();
          middlePad->Draw ();
          bottomPad->Draw ();
        }
        else {
          topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
          middlePad = dynamic_cast<TPad*> (gDirectory->Get (middlePadName));
          bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
        }

        topPad->cd ();
        GetDrawnObjects ();
        bool plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        double min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent];
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max)  max = h->GetMaximum ();
        } // end loop over phi
        min = (min > 0 ? 0.5*min : 0.1);
        max = (max > 0 ? 2*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            h->GetYaxis ()->SetTitle ("Per-trigger yield Y (#it{p}_{T})");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          }
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGraphAsymmErrors* g = make_graph (h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerColor (colors[iPhi]);
              g->SetLineColor (colors[iPhi]);
              g->SetMarkerSize (1);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (0);
              g->SetFillColorAlpha (fillColors[iPhi], 0.8);
              g->SetFillStyle (3001);
            }

            g->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            g->GetYaxis ()->SetTitle ("Per-trigger yield Y (#it{p}_{T})");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

            string drawString = (plotNewAxes && iPhi == 0 ? "A" : "");
            if (plotAsSystematic) drawString = drawString + "2";
            drawString = drawString + " P";
            
            g->Draw (drawString.c_str ());
          } // end loop over phi
        }

        if (!canvasExists)
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++)
            LabelTrkYield (iCent, iPhi);

        if (!plotSignal)
          continue;

        middlePad->cd ();
        GetDrawnObjects ();
        plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent];
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max) max = h->GetMaximum ();
        } // end loop over phi
        float delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* h = h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            h->GetYaxis ()->SetTitle ("Signal Yield");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetXaxis ()->SetTitleOffset (2.5 * h->GetXaxis ()->GetTitleOffset ());
            h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));

            TGraphAsymmErrors* g = make_graph (h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerColor (colors[iPhi]);
              g->SetLineColor (colors[iPhi]);
              g->SetMarkerSize (1);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (0);
              g->SetFillColorAlpha (fillColors[iPhi], 0.8);
              g->SetFillStyle (3001);
            }

            g->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            g->GetYaxis ()->SetTitle ("Signal Yield");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetXaxis ()->SetTitleOffset (2.5 * g->GetXaxis ()->GetTitleOffset ());
            g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

            string drawString = (plotNewAxes && iPhi == numPhiBins-1 ? "A" : "");
            if (plotAsSystematic) drawString = drawString + "2";
            drawString = drawString + " P";
            g->Draw (drawString.c_str ());
          } // end loop over phi
        }

        bottomPad->cd ();
        GetDrawnObjects ();
        plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][0][iPhi][iCent];
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max) max = h->GetMaximum ();
        } // end loop over phi
        //delta = max - min;
        //min = min - 0.3*delta;
        //max = max + 0.3*delta;
        delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* h = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][0][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            h->GetYaxis ()->SetTitle ("Signal / Bkg.");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetXaxis ()->SetTitleOffset (2.5 * h->GetXaxis ()->GetTitleOffset ());
            h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGraphAsymmErrors* g = make_graph (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][0][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerColor (colors[iPhi]);
              g->SetLineColor (colors[iPhi]);
              g->SetMarkerSize (1);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (0);
              g->SetFillColorAlpha (fillColors[iPhi], 0.8);
              g->SetFillStyle (3001);
            }

            g->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            g->GetYaxis ()->SetTitle ("Signal / Bkg.");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetXaxis ()->SetTitleOffset (2.5 * g->GetXaxis ()->GetTitleOffset ());
            g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

            string drawString = (plotNewAxes && iPhi == numPhiBins-1 ? "AP" : "P");
            if (plotAsSystematic) drawString = drawString + "2";
            drawString = drawString + " P";
            g->Draw (drawString.c_str ());
          } // end loop over phi
        }
      } // end loop over cents
      
      c->SaveAs (Form ("%s/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track pT distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LabelTrkYield (const short iCent, const short iPhi) {
  const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);

  if (iCent == 0)
    myText (0.22, 0.06, kBlack, "#it{pp}", 0.06);
  else
    myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

  if (iCent == 0)
    myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
  else if (iCent == numCentBins-1) {
    if (iPhi == 0) {
      myText (0.35, 0.91, kBlack, "Z-tagged", 0.054);
      //myText (0.45, 0.91, kBlack, "Truth", 0.06);
      myText (0.544, 0.91, kBlack, "Minbias", 0.054);
      myText (0.703, 0.91, kBlack, "#Delta#phi", 0.054);
    }

    TVirtualPad* cPad = gPad; // store current pad
    const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
    const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
    TBox* b = TBoxNDC (0.612-0.032, 0.85-0.06*iPhi-0.016, 0.612+0.032, 0.85-0.06*iPhi+0.016);
    b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    b->Draw ("l");
    myMarkerText (0.462, 0.852-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.5, 0.054);
    myText (0.70, 0.85-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.054);
    cPad->cd ();
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: CalculateIAA () {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (iaaCalculated)
    return;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        TH1D* ppHist = h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          if (!h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent]) {
            TH1D* PbPbHist = (TH1D*)(h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent]->Clone ());
            PbPbHist->Divide (ppHist);
            h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent] = PbPbHist;
          } 
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species
  iaaCalculated = true;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotIAARatios (const bool plotAsSystematic, const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_z_trk_pt_iaa_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-1), 500);
        c->Divide (numCentBins-1, 1);
        gDirectory->Add (c);
      }

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            TH1D* h = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            h->GetYaxis ()->SetRangeUser (0, 1.6);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            h->GetYaxis ()->SetTitle ("I_{AA}");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetXaxis ()->SetTitleOffset (0.9 * h->GetXaxis ()->GetTitleOffset ());
            //h->GetYaxis ()->SetTitleOffset (0.9 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (!canvasExists && iPhi == 1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");

            LabelIAARatios (iCent, iPhi);
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

            TGraphAsymmErrors* g = make_graph (h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              deltaize (g, 1+((numPhiBins-1)*((int)useAltMarker)-iPhi)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
              g->SetLineColor (colors[iPhi]);
              g->SetMarkerColor (colors[iPhi]);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (1.2);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (0);
              g->SetFillColorAlpha (fillColors[iPhi], 0.8);
              g->SetFillStyle (3001);
            }

            g->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            g->GetYaxis ()->SetRangeUser (0, 1.6);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            g->GetYaxis ()->SetTitle ("I_{AA}");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetXaxis ()->SetTitleOffset (0.9 * g->GetXaxis ()->GetTitleOffset ());
            //g->GetYaxis ()->SetTitleOffset (0.9 * g->GetYaxis ()->GetTitleOffset ());

            string drawString = (!canvasExists && iPhi == 1 ? "A" : "");
            if (plotAsSystematic) drawString = drawString + "2";
            drawString = drawString + " P";
            g->Draw (drawString.c_str ());

            LabelIAARatios (iCent, iPhi);
          } // end loop over phi
        }
      } // end loop over cents
      c->cd (1);

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
      } // end loop over cents

      c->SaveAs (Form ("%s/iaa_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LabelIAARatios (const short iCent, const short iPhi) {

  if (iCent == 1)
    myText (0.24, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
  else if (iCent == numCentBins-1) {
    if (iPhi == 1) {
      //myText (0.44, 0.91, kBlack, "MBM", 0.05);
      //myText (0.577, 0.91, kBlack, "HM", 0.05);
      //myText (0.44, 0.91, kBlack, "Data", 0.05);
      //myText (0.577, 0.91, kBlack, "MC", 0.05);
      myText (0.685, 0.91, kBlack, "#Delta#phi", 0.05);
    }
    const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
    const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerText (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    //myMarkerText (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myMarkerText (0.63, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myText (0.68, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: CalculateICP () {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (icpCalculated)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        TH1D* periphHist = h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent]) {
            TH1D* centHist = (TH1D*)(h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent] = centHist;
          }
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species
  icpCalculated = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: PlotICPRatios (const bool plotAsSystematic, const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!icpCalculated)
    CalculateICP ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_z_trk_pt_icp_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-2), 500);
        c->Divide (numCentBins-2, 1);
        gDirectory->Add (c);
      }

      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            TH1D* h = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            h->GetYaxis ()->SetRangeUser (0, 1.6);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            h->GetYaxis ()->SetTitle ("I_{CP}");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetXaxis ()->SetTitleOffset (0.9 * h->GetXaxis ()->GetTitleOffset ());
            //h->GetYaxis ()->SetTitleOffset (0.9 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (!canvasExists && iPhi == 1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");

            LabelICPRatios (iCent, iPhi);
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
            
            TGraphAsymmErrors* g = make_graph (h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              deltaize (g, 1+((numPhiBins-1)*((int)useAltMarker)-iPhi)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
              g->SetLineColor (colors[iPhi]);
              g->SetMarkerColor (colors[iPhi]);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (1.2);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (0);
              g->SetFillColorAlpha (fillColors[iPhi], 0.8);
              g->SetFillStyle (3001);
            }

            g->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
            g->GetYaxis ()->SetRangeUser (0, 1.6);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            g->GetYaxis ()->SetTitle ("I_{CP}");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetXaxis ()->SetTitleOffset (0.9 * g->GetXaxis ()->GetTitleOffset ());
            //g->GetYaxis ()->SetTitleOffset (0.9 * g->GetYaxis ()->GetTitleOffset ());

            string drawString = (!canvasExists && iPhi == 1 ? "A" : "");
            if (plotAsSystematic) drawString = drawString + "2";
            drawString = drawString + " P";
            g->Draw (drawString.c_str ());

            LabelICPRatios (iCent, iPhi);
          } // end loop over phi
        }
      } // end loop over cents

      c->cd (1);
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);
      } // end loop over cents

      c->SaveAs (Form ("%s/icp_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis :: LabelICPRatios (const short iCent, const short iPhi) {
  if (iCent == 2)
    myText (0.24, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
  else if (iCent == numCentBins-1) {
    if (iPhi == 1) {
      //myText (0.55, 0.91, kBlack, "Data", 0.05);
      //myText (0.577, 0.91, kBlack, "MC", 0.05);
      myText (0.44, 0.91, kBlack, "MBM", 0.05);
      myText (0.577, 0.91, kBlack, "HM", 0.05);
      myText (0.685, 0.91, kBlack, "#Delta#phi", 0.05);
    }
    const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
    const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.71-0.024, 0.91-0.06*iPhi-0.016, 0.71+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    myMarkerText (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    myMarkerText (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myText (0.68, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
  }
}

#endif
