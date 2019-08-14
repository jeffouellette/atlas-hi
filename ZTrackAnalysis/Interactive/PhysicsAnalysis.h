#ifndef __PhysicsAnalysis_h__
#define __PhysicsAnalysis_h__

#include "Params.h"

#include <ArrayTemplates.h>

#include <AtlasUtils.h>

#include <TEfficiency.h>
#include <TClass.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TBox.h>
#include <TVirtualFitter.h>

#include <iostream>
#include <string>

using namespace atlashi;
using namespace std;

typedef TGraphAsymmErrors TGAE;

class PhysicsAnalysis {

  protected:
  string name = "";
  string directory = "";
  bool backgroundSubtracted = false;
  bool sameSignsSubtracted = false;

  vector<TH1*> drawnHists;
  vector<TGAE*> drawnGraphs;

  bool iaaCalculated = false;
  bool icpCalculated = false;

  TFile* trkEffFile = nullptr;
  TFile* histFile   = nullptr;
  bool histsLoaded  = false;
  bool histsScaled  = false;

  public:
  bool plotFill       = false; // whether to plot as filled (bar) graph or points w/ errors
  bool plotSignal     = true; // whether to plot background subtracted plots
  bool useAltMarker   = false; // whether to plot as open markers (instead of closed)
  bool useHITight     = false; // whether to use HITight tracking efficiencies
  bool use2015Effs    = false; // whether to use tracking efficiencies fromo 2015
  float trkEffNSigma  = 0; // how many sigma to vary the track efficiency by (-1,0,+1 suggested)
  float trkPurNSigma  = 0; // how many sigma to vary the track purity by (-1,0,+1 suggested)

  // Event info distributions (for reweighting)
  TH1D* h_PbPbFCal_weights  = nullptr;
  TH1D** h_PbPbQ2_weights   = Get1DArray <TH1D*> (numFinerCentBins);
  TH1D* h_ppNch_weights     = nullptr;

  // Efficiencies
  TH1D*** h_trk_effs    = Get2DArray <TH1D*> (numCentBins, numEtaTrkBins); // iCent, iEta
  //TF1**   f_trk_effs    = Get1DArray <TH1D*> (numCentBins); // iCent (eta dependence is extrapolated)
  TH2D**  h2_trk_effs   = Get1DArray <TH2D*> (numCentBins); // iCent
  //TEfficiency*** h_trk_effs   = Get2DArray <TEfficiency*> (numCentBins, numEtaTrkBins); // iCent, iEta
  TH2D** h2_num_trk_effs      = Get1DArray <TH2D*> (numCentBins); // iCent
  TH2D** h2_den_trk_effs      = Get1DArray <TH2D*> (numCentBins);

  // Tracking purities
  //TH1D*** h_trk_purs    = Get2DArray <TH1D*> (numCentBins, numEtaTrkBins); // iCent, iEta
  TH2D**  h2_trk_purs   = Get1DArray <TH2D*> (numCentBins); // iCent
  TEfficiency*** h_trk_purs   = Get2DArray <TEfficiency*> (numCentBins, numEtaTrkBins); // iCent, iEta
  TH2D** h2_num_trk_purs      = Get1DArray <TH2D*> (numCentBins); // iCent
  TH2D** h2_den_trk_purs      = Get1DArray <TH2D*> (numCentBins);

  //// Correlations plots
  //TH2D****   h_z_trk_pt_phi  = Get3DArray <TH2D*> (nPtZBins, numCentBins, 3);             // iPtZ, iCent, iSpc (0=ee, 1=mumu, 2=combined)
  //TH1D*****  h_z_trk_phi     = Get4DArray <TH1D*> (nPtTrkBins, nPtZBins, numCentBins, 3); // iPtTrk, iPtZ, iCent, iSpc
  
  // Physics plots
  TH1D*****   h_z_trk_phi         = Get4DArray <TH1D*> (3, nPtZBins, nPtTrkBins, numCentBins); // iSpc, iPtZ, iPtTrk, iCent
  TH1D*****   h_z_trk_phi_sub     = Get4DArray <TH1D*> (3, nPtZBins, nPtTrkBins, numCentBins); // iSpc, iPtZ, iPtTrk, iCent

  TH1D*****   h_z_trk_raw_pt      = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);   // iSpc, iPtZ, iPhi, iCent
  TH1D****    h_z_counts          = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);               // iSpc, iPtZ, iCent

  TH1D*****  h_z_trk_pt             = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_z_trk_pt_sub         = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_z_trk_pt_sig_to_bkg  = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_z_trk_pt_iaa         = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_z_trk_pt_icp         = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_z_trk_zpt            = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zpt_sub        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zpt_sig_to_bkg = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zpt_iaa        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zpt_icp        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  TH1D***** h_z_trk_xzh             = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_z_trk_xzh_sub         = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_z_trk_xzh_sig_to_bkg  = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_z_trk_xzh_iaa         = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_z_trk_xzh_icp         = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_z_trk_zxzh             = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zxzh_sub         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zxzh_sig_to_bkg  = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zxzh_iaa         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_z_trk_zxzh_icp         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  PhysicsAnalysis () { }

  PhysicsAnalysis (const char* _name, const char* subDir, const bool _useHITight = false) {
    name = _name;
    directory = Form ("DataAnalysis/%s/", subDir);
    plotFill = false;
    useHITight = _useHITight;
    LoadTrackingEfficiencies ();
    LoadTrackingPurities ();
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  virtual ~PhysicsAnalysis () {
    Delete1DArray (h_PbPbQ2_weights,  numFinerCentBins);

    Delete2DArray (h_trk_effs,      numCentBins, numEtaTrkBins);
    Delete1DArray (h2_trk_effs,     numCentBins);
    Delete1DArray (h2_num_trk_effs, numCentBins);
    Delete1DArray (h2_den_trk_effs, numCentBins);

    Delete2DArray (h_trk_purs,      numCentBins, numEtaTrkBins);
    Delete1DArray (h2_trk_purs,     numCentBins);
    Delete1DArray (h2_num_trk_purs, numCentBins);
    Delete1DArray (h2_den_trk_purs, numCentBins);

    Delete4DArray (h_z_trk_phi,     3, nPtZBins, nPtTrkBins, numCentBins);
    Delete4DArray (h_z_trk_phi_sub,         3, nPtZBins, nPtTrkBins, numCentBins);

    Delete4DArray (h_z_trk_raw_pt,  3, nPtZBins, numPhiBins, numCentBins);
    Delete3DArray (h_z_counts,      3, nPtZBins, numCentBins);

    Delete4DArray (h_z_trk_pt,              3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_pt_sub,          3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_pt_sig_to_bkg,   3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_pt_iaa,          3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_pt_icp,          3, nPtZBins, numPhiBins, numCentBins);

    Delete3DArray (h_z_trk_zpt,             3, nPtZBins, numCentBins);
    Delete4DArray (h_z_trk_zpt_sub,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zpt_sig_to_bkg,  3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zpt_iaa,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zpt_icp,         3, nPtZBins, numPhiBins, numCentBins);

    Delete4DArray (h_z_trk_xzh,             3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_xzh_sub,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_xzh_sig_to_bkg,  3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_xzh_iaa,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_xzh_icp,         3, nPtZBins, numPhiBins, numCentBins);

    Delete4DArray (h_z_trk_zxzh,            3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zxzh_sub,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zxzh_sig_to_bkg, 3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zxzh_iaa,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_z_trk_zxzh_icp,        3, nPtZBins, numPhiBins, numCentBins);
  }

  protected:
  void LabelTrackingEfficiencies (const short iCent, const short iEta);
  void LabelTrackingPurities (const short iCent, const short iEta);
  void LabelCorrelations (const short iPtZ, const short iPtTrk, const short iCent);
  void LabelTrkYield (const short iCent, const short iPhi, const short iPtZ);
  void LabelTrkYieldZPt (const short iCent, const short iPtZ);
  void LabelIAAdPhi (const short iCent, const short iPhi, const short iPtZ);
  void LabelIAAdCent (const short iCent, const short iPhi, const short iPtZ);
  void LabelIAAdPtZ (const short iCent, const short iPtZ);
  void LabelICPdPhi (const short iCent, const short iPhi, const short iPtZ);
  void LabelICPdCent (const short iCent, const short iPhi, const short iPtZ);

  void GetDrawnObjects ();
  void GetMinAndMax (double &min, double &max, const bool log = false);
  void SetMinAndMax (double min, double max);

  public:
  string Name () { return name; }
  void SetName (string _name) { name = _name; }

  virtual TGAE* GetTGAE (TH1D* h);

  virtual void CreateHists ();
  virtual void CopyAnalysis (PhysicsAnalysis* a, const bool copyBkgs = false);
  virtual void CombineHists ();
  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true);
  virtual void SaveHists (const char* histFileName = "savedHists.root");
  virtual void ScaleHists ();
  virtual void Execute (const char* inFileName = "outFile.root", const char* outFileName = "savedHists.root");
  virtual void SubtractBackground (PhysicsAnalysis* a = nullptr);
  virtual void SubtractSameSigns (PhysicsAnalysis* a);

  virtual void ConvertToStatVariation (const bool upVar = true, const float nSigma = 1.); // adds or subtracts nSigma of statistical errors to analysis

  virtual void LoadTrackingEfficiencies (); // defaults to HILoose
  virtual double GetTrackingEfficiency (const float fcal_et, const float trk_pt, const float trk_eta, const bool isPbPb = true);

  virtual void LoadTrackingPurities (); // defaults to HILoose
  virtual double GetTrackingPurity (const float fcal_et, const float trk_pt, const float trk_eta, const bool isPbPb = true);

  void PrintZYields (const int iPtZ = 2);

  //void Plot3DDist ();
  void PlotCorrelations (const short pSpc = 2, const short pPtZ = nPtZBins-1, const bool _subBkg = false);
  //void PlotZMissingPt ();
  void PlotTrackingEfficiencies ();
  void PlotTrackingEfficiencies2D ();
  void PlotTrackingPurities ();
  void PlotTrackingPurities2D ();

  void CalculateIAA ();
  void CalculateICP ();

  virtual void PlotRawTrkYield (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotTrkYield (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotTrkYieldZPt (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotIAAdPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAAdCent (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAAdPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotICPdPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotICPdCent (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);

};


void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}


void PhysicsAnalysis :: GetDrawnObjects () {
  TList* primitives = gPad->GetListOfPrimitives ();
  drawnHists.clear ();
  drawnGraphs.clear ();
  for (int i = 0; i < primitives->GetSize (); i++) {
    TObject *obj = primitives->At (i);
    if (obj->IsA()->InheritsFrom (TH1::Class ())) {
      drawnHists.push_back ((TH1*)obj);
    }
    else if (obj->IsA()->InheritsFrom (TGAE::Class ())) {
      drawnGraphs.push_back ((TGAE*)obj);
    }
  }
}


void PhysicsAnalysis :: GetMinAndMax (double &min, double &max, const bool log) {
  for (TH1* h : drawnHists) {
    const double _max = log ? h->GetMaximum (0) : h->GetMaximum ();
    const double _min = log ? h->GetMinimum (0) : h->GetMinimum ();

    if (_max > max) max = _max;
    if (_min < min) min = _min;
  }
  for (TGAE* g : drawnGraphs) {
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


void PhysicsAnalysis :: SetMinAndMax (double min, double max) {
  for (TH1* h : drawnHists) {
    h->GetYaxis ()->SetRangeUser (min, max);
  }
  for (TGAE* g : drawnGraphs) {
    g->GetYaxis ()->SetRangeUser (min, max);
  }
  gPad->Update ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create TGraphAsymmErrors from a histogram
////////////////////////////////////////////////////////////////////////////////////////////////
TGAE* PhysicsAnalysis :: GetTGAE (TH1D* h) {
  return make_graph (h);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CreateHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        //for (short iZH = 0; iZH < nZHBins; iZH++) {
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc] = new TH2D (Form ("h_z_trk_pt_phi_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc]->Sumw2 ();
        //}
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_raw_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins, ptTrkBins[iPtZ]);
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins, ptTrkBins[iPtZ]);
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_xzh_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nZHBins, zHBins);
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        }
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent] = new TH1D (Form ("h_z_trk_phi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()), "", 80, -pi/2, 3*pi/2);
          h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Sumw2 ();
        }
        h_z_counts[iSpc][iPtZ][iCent] = new TH1D (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", 2, 0, 2);
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
void PhysicsAnalysis :: CopyAnalysis (PhysicsAnalysis* a, const bool copyBkgs) {
  if (name == "")
    cout << "Warning in PhysicsAnalysis :: CopyAnalysis: name of analysis not set!" << endl;

  // Don't need to clone these histograms

  // Reweighting histograms
  h_PbPbFCal_weights  = nullptr;
  h_PbPbQ2_weights    = Get1DArray <TH1D*> (numFinerCentBins);

  // Efficiencies
  h_trk_effs      = a->h_trk_effs;
  h2_trk_effs     = a->h2_trk_effs;
  h2_num_trk_effs = a->h2_num_trk_effs;
  h2_den_trk_effs = a->h2_den_trk_effs;

  // Should clone these histograms
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        //for (short iZH = 0; iZH < nZHBins; iZH++) {
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc] = (TH2D*) a->h_z_trk_pt_phi[iPtZ][iCent][iSpc]->Clone (Form ("h_z_trk_pt_phi_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //}

        h_z_trk_zpt[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_z_trk_zxzh[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_raw_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_xzh_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        } // end loop over iPhi
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) a->h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Clone (Form ("h_z_trk_phi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
        }
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

            if (a->h_z_trk_pt_sub[iSpc][iPtZ][iCent]) {
              h_z_trk_zpt_sub[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zpt_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
              h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }
            if (a->h_z_trk_pt_sub[iSpc][iPtZ][iCent]) {
              h_z_trk_zxzh_sub[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
              h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }

            for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
              if (a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]) {
                h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
                h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
              }
              if (a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]) {
                h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_xzh_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
                h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_xzh_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
              }
            } // end loop over phi
            for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
              if (a->h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]) {
                h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) a->h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]->Clone (Form ("h_z_trk_phi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
              }
            } // end loop over pT^trk
          } // end loop over pT^Z bins
        } // end loop over cents
      } // end loop over species
      backgroundSubtracted = true;
    }

    if (a->iaaCalculated) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {

            if (a->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]) {
              h_z_trk_zpt_iaa[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }
            if (a->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]) {
              h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }

            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              if (a->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent])
                h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
              if (a->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent])
                h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_xzh_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
            } // end loop over phi
          } // end loop over cents
        } // end loop over pT^Z bins
      } // end loop over species
      iaaCalculated = true;
    }

    if (a->icpCalculated) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (short iCent = 2; iCent < numCentBins; iCent++) {

            if (a->h_z_trk_zpt_icp[iSpc][iPtZ][iCent]) {
              h_z_trk_zpt_icp[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zpt_icp[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_icp_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }
            if (a->h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]) {
              h_z_trk_zxzh_icp[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_icp_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }

            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              if (a->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent])
                h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_icp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
              if (a->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent])
                h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_xzh_icp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
            } // end loop over phi
          } // end loop over cents
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
void PhysicsAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");
  if (histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "read");

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        //for (short iZH = 0; iZH < nZHBins; iZH++) {
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc] = (TH2D*) histFile->Get (Form ("h_z_trk_pt_phi_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_phi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
        }
        //}
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_raw_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        }
        h_z_trk_zpt[iSpc][iPtZ][iCent] = new TH1D (Form ("h_z_trk_zpt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nPtTrkBins, ptTrkBins[iPtZ]);
        h_z_trk_zpt[iSpc][iPtZ][iCent]->Sumw2 ();
        h_z_trk_zxzh[iSpc][iPtZ][iCent] = new TH1D (Form ("h_z_trk_zxzh_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nZHBins, zHBins);
        h_z_trk_zxzh[iSpc][iPtZ][iCent]->Sumw2 ();
        
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      }
    }
  }

  histsLoaded = true;

  if (_finishHists) {
    PhysicsAnalysis :: CombineHists ();
    PhysicsAnalysis :: ScaleHists ();
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SaveHists (const char* histFileName) {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");
  if (!histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  if (!histFile) {
    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "recreate");
    histFile->cd ();
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        //for (short iZH = 0; iZH < nZHBins; iZH++) {
        //SafeWrite (h_z_trk_pt_phi[iPtZ][iCent][iSpc]);
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          SafeWrite (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
        }
        //}
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SafeWrite (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        }
        SafeWrite (h_z_counts[iSpc][iPtZ][iCent]);
      }
    }
  }
  
  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iSpc = 0; iSpc < 2; iSpc++) {
        //for (short iZH = 0; iZH < nZHBins; iZH++) {
        //h_z_trk_pt_phi[iPtZ][iCent][2]->Add (h_z_trk_pt_phi[iPtZ][iCent][iSpc]);
        //}
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          if (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]) h_z_trk_phi[2][iPtZ][iPtTrk][iCent]->Add (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]) h_z_trk_raw_pt[2][iPtZ][iPhi][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);

          if (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]) h_z_trk_pt[2][iPtZ][iPhi][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
          if (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]) h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);

          if (iPhi != 0) {
            if (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zpt[iSpc][iPtZ][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
            if (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zpt[2][iPtZ][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
            if (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zxzh[iSpc][iPtZ][iCent]->Add (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
            if (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zxzh[2][iPtZ][iCent]->Add (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
          }

          if (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]) h_z_trk_xzh[2][iPtZ][iPhi][iCent]->Add (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
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
void PhysicsAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          const double yieldNormFactor = countsHist->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);

          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1./ countsHist->GetBinContent (1));
          if (yieldNormFactor > 0) {
            h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1. / yieldNormFactor, "width");
            h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Scale (1. / yieldNormFactor, "width");
          }
        } // end loop over phi

        h_z_trk_zpt[iSpc][iPtZ][iCent]->Scale (1./ (countsHist->GetBinContent (1) * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
        h_z_trk_zxzh[iSpc][iPtZ][iCent]->Scale (1./ (countsHist->GetBinContent (1) * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
        
      } // end loop over pT^Z

      //double normFactor = 0;
      //for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++)
      //  normFactor += h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1);

      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          //for (short iZH = 0; iZH < nZHBins; iZH++) {
            const double normFactor = h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1);
            //h_z_trk_phi[iPtTrk][iPtZ][iCent][iSpc] = h_z_trk_pt_phi[iPtZ][iCent][iSpc]->ProjectionX (Form ("h_z_trk_phi_iPtTrk%i_iPtZ%i_iCent%i_%s", iPtTrk, iPtZ, iCent, name.c_str ()), iPtTrk+1, iPtTrk+1);

            TH1D* h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
            h->Rebin (2);
            if (iPtTrk > 3)
              h->Rebin (2);
            if (iCent != 0)
              h->Rebin (2);
            if (normFactor > 0)
              h->Scale (1. / normFactor);
          //}
        }
      }
    } // end loop over centralities
  } // end loop over species

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
// Designed to be overloaded. The default here is for analyzing data.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: Execute (const char* inFileName, const char* outFileName) {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("isEE",       &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchAddress ("psi2",       &psi2);
    PbPbTree->SetBranchAddress ("vz",         &vz);
    PbPbTree->SetBranchAddress ("z_pt",       &z_pt);
    PbPbTree->SetBranchAddress ("z_y",        &z_y);
    PbPbTree->SetBranchAddress ("z_phi",      &z_phi);
    PbPbTree->SetBranchAddress ("z_m",        &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",      &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",     &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",     &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",  &l1_charge);
    PbPbTree->SetBranchAddress ("l1_trk_pt",  &l1_trk_pt);
    PbPbTree->SetBranchAddress ("l1_trk_eta", &l1_trk_eta);
    PbPbTree->SetBranchAddress ("l1_trk_phi", &l1_trk_phi);
    PbPbTree->SetBranchAddress ("l2_pt",      &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",     &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",     &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",  &l2_charge);
    PbPbTree->SetBranchAddress ("l2_trk_pt",  &l2_trk_pt);
    PbPbTree->SetBranchAddress ("l2_trk_eta", &l2_trk_eta);
    PbPbTree->SetBranchAddress ("l2_trk_phi", &l2_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",       &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",     &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",    &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",    &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      //if (z_m < 86 || z_m > 96)
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

      short iFinerCent = 0;
      while (iFinerCent < numFinerCentBins) {
        if (fcal_et < finerCentBins[iFinerCent])
          break;
        else
          iFinerCent++;
      }
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const float zH = trkpt / z_pt;
        const short iZH = GetiZH (zH);
        if (iZH < 0 || iZH > nZHBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), true);
        if (trkEff == 0)
          continue;

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (zH, event_weight / trkEff);
          }
        }

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight / trkEff);
        }
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);dd
      } // end loop over tracks

    } // end loop over Pb+Pb tree
    cout << "Done primary Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("isEE",       &isEE);
    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchAddress ("z_pt",       &z_pt);
    ppTree->SetBranchAddress ("z_y",        &z_y);
    ppTree->SetBranchAddress ("z_phi",      &z_phi);
    ppTree->SetBranchAddress ("z_m",        &z_m);
    ppTree->SetBranchAddress ("l1_pt",      &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",     &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",     &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",  &l1_charge);
    ppTree->SetBranchAddress ("l1_trk_pt",  &l1_trk_pt);
    ppTree->SetBranchAddress ("l1_trk_eta", &l1_trk_eta);
    ppTree->SetBranchAddress ("l1_trk_phi", &l1_trk_phi);
    ppTree->SetBranchAddress ("l2_pt",      &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",     &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",     &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",  &l2_charge);
    ppTree->SetBranchAddress ("l2_trk_pt",  &l2_trk_pt);
    ppTree->SetBranchAddress ("l2_trk_eta", &l2_trk_eta);
    ppTree->SetBranchAddress ("l2_trk_phi", &l2_trk_phi);
    ppTree->SetBranchAddress ("ntrk",       &ntrk);
    ppTree->SetBranchAddress ("trk_pt",     &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",    &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",    &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      //if (z_m < 86 || z_m > 96)
      //  continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const float zH = trkpt / z_pt;
        const short iZH = GetiZH (zH);
        if (iZH < 0 || iZH > nZHBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), false);
        if (trkEff == 0)
          continue;

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (zH, event_weight / trkEff);
          }
        }

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight / trkEff);
        }
      } // end loop over tracks

    } // end loop over pp tree
    cout << "Done primary pp loop." << endl;
  }

  //CombineHists ();
  //ScaleHists ();
  
  SaveHists (outFileName);
  //LoadHists ();

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  //Delete2DArray (trkPtProj, numPhiBins, nPtTrkBins);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking efficiencies into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadTrackingEfficiencies () {
  SetupDirectories ("TrackingEfficiencies/", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;

  if (useHITight)
    trkEffFile = new TFile (Form ("%s/Variations/TrackHITightWPVariation/trackingEfficiencies.root", rootPath.Data ()), "read");
  else if (use2015Effs)
    trkEffFile = new TFile (Form ("%s/Nominal/PbPb_Hijing_15.root", rootPath.Data ()), "read");
  else if (!useHITight)
    trkEffFile = new TFile (Form ("%s/Nominal/trackingEfficiencies.root", rootPath.Data ()), "read");

  for (int iCent = 0; iCent < numCentBins; iCent++) {
  //for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    h2_num_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_matched_reco_tracks_iCent%i", iCent));
    h2_den_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_tracks_iCent%i", iCent));

    if (iCent > 0) {
      h2_num_trk_effs[iCent]->RebinX (2);
      h2_num_trk_effs[iCent]->RebinY (2);
      h2_den_trk_effs[iCent]->RebinX (2);
      h2_den_trk_effs[iCent]->RebinY (2);
    }

    //h2_trk_effs[iCent] = new TEfficiency (*(h2_num_trk_effs[iCent]), *(h2_den_trk_effs[iCent]));
    h2_trk_effs[iCent] = (TH2D*) h2_num_trk_effs[iCent]->Clone (Form ("h2_trk_eff_iCent%i", iCent));
    h2_trk_effs[iCent]->Divide (h2_den_trk_effs[iCent]);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TH1D* num = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta)))->GetCopyPassedHisto ();
      //TH1D* den = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta)))->GetCopyTotalHisto ();

      TH1D* num = (TH1D*) trkEffFile->Get (Form ("h_trk_eff_num_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkEffFile->Get (Form ("h_trk_eff_den_iCent%i_iEta%i", iCent, iEta));

      if (iCent > 0) {
        num->Rebin (2);
        den->Rebin (2);
      }

      //h_trk_effs[iCent][iEta] =  new TEfficiency (*num, *den);
      h_trk_effs[iCent][iEta] = (TH1D*) num->Clone (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta));
      h_trk_effs[iCent][iEta]->Divide (den);
      //h_trk_effs[iCent][iEta]->SetDirectory (_gDirectory);

      //delete num;
      //delete den;
      //h_trk_effs[iCent][iEta]->SetName ("h_trk_eff_iCent%i_iEta%i");
      
      //h_trk_effs[iCent][iEta] = (TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta));
    }
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking efficiency for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double PhysicsAnalysis :: GetTrackingEfficiency (const float fcal_et, const float trk_pt, const float trk_eta, const bool isPbPb) {
  short iCent = 0;
  if (isPbPb) {
    while (iCent < numCentBins) {
      if (fcal_et < centBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numCentBins-1)
      return 0;
  }

  //TEfficiency* t = h2_trk_effs[iCent];

  //double eff = t->GetEfficiency (t->FindFixBin (trk_eta, trk_pt));
  //if (trkEffNSigma > 0)
  //  eff += trkEffNSigma * t->GetEfficiencyErrorUp (t->FindFixBin (trk_eta, trk_pt));
  //else if (trkEffNSigma < 0)
  //  eff += trkEffNSigma * t->GetEfficiencyErrorLow (t->FindFixBin (trk_eta, trk_pt));

  TH2D* t = h2_trk_effs[iCent];
  double eff = t->GetBinContent (t->FindFixBin (trk_eta, trk_pt)) + trkEffNSigma * t->GetBinError (t->FindFixBin (trk_eta, trk_pt));

  //if (eff == 0)
  //  return 1;

  return eff;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking purities into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadTrackingPurities () {
  SetupDirectories ("TrackingPurities/", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;

  if (!useHITight)
    trkEffFile = new TFile (Form ("%s/Nominal/trackingPurities.root", rootPath.Data ()), "read");
  else
    trkEffFile = new TFile (Form ("%s/Variations/TrackHITightWPVariation/trackingPurities.root", rootPath.Data ()), "read");

  if (!trkEffFile || !trkEffFile->IsOpen ()) {
    cout << "Error in PhysicsAnalysis.h:: LoadTrackingPurities can not find file for " << name << endl;
    return;
  }

  for (int iCent = 0; iCent < numCentBins; iCent++) {
  //for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    h2_num_trk_purs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_primary_reco_tracks_iCent%i", iCent));
    h2_den_trk_purs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_reco_tracks_iCent%i", iCent));

    if (iCent > 0) {
      h2_num_trk_purs[iCent]->RebinX (2);
      h2_num_trk_purs[iCent]->RebinY (2);
      h2_den_trk_purs[iCent]->RebinX (2);
      h2_den_trk_purs[iCent]->RebinY (2);
    }

    //h2_trk_purs[iCent] = new TEfficiency (*(h2_num_trk_purs[iCent]), *(h2_den_trk_purs[iCent]));
    h2_trk_purs[iCent] = (TH2D*) h2_num_trk_purs[iCent]->Clone (Form ("h2_trk_pur_iCent%i", iCent));
    h2_trk_purs[iCent]->Divide (h2_den_trk_purs[iCent]);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TH1D* num = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta)))->GetCopyPassedHisto ();
      //TH1D* den = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta)))->GetCopyTotalHisto ();

      TH1D* num = (TH1D*) trkEffFile->Get (Form ("h_trk_pur_num_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkEffFile->Get (Form ("h_trk_pur_den_iCent%i_iEta%i", iCent, iEta));

      if (iCent > 0) {
        num->Rebin (2);
        den->Rebin (2);
      }

      h_trk_purs[iCent][iEta] =  new TEfficiency (*num, *den);
      h_trk_purs[iCent][iEta]->SetName (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta));
      //h_trk_purs[iCent][iEta] = (TH1D*) num->Clone (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta));
      //h_trk_purs[iCent][iEta]->Divide (den);
      //h_trk_purs[iCent][iEta]->SetDirectory (_gDirectory);

      //delete num;
      //delete den;
      //h_trk_purs[iCent][iEta]->SetName ("h_trk_pur_iCent%i_iEta%i");
      
      //h_trk_purs[iCent][iEta] = (TEfficiency*) trkEffFile->Get (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta));
    }
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking purity factor for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double PhysicsAnalysis :: GetTrackingPurity (const float fcal_et, const float trk_pt, const float trk_eta, const bool isPbPb) {
  short iCent = 0;
  if (isPbPb) {
    while (iCent < numCentBins) {
      if (fcal_et < centBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numCentBins-1)
      return 0;
  }

  //TEfficiency* t = h2_trk_effs[iCent];

  //double eff = t->GetEfficiency (t->FindFixBin (trk_eta, trk_pt));
  //if (trkEffNSigma > 0)
  //  eff += trkEffNSigma * t->GetEfficiencyErrorUp (t->FindFixBin (trk_eta, trk_pt));
  //else if (trkEffNSigma < 0)
  //  eff += trkEffNSigma * t->GetEfficiencyErrorLow (t->FindFixBin (trk_eta, trk_pt));

  TH2D* t = h2_trk_purs[iCent];
  double eff = t->GetBinContent (t->FindFixBin (trk_eta, trk_pt)) + trkPurNSigma * t->GetBinError (t->FindFixBin (trk_eta, trk_pt));

  //if (eff == 0)
  //  return 1;

  return eff;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints yield of Z's that meet the event selection criteria in each centrality
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PrintZYields (const int iPtZ) {
  cout << "\t\t\t\\multirow{4}{*}";
  if (iPtZ == nPtZBins)
    cout << Form ("{\\pt > \\SI{%g}{\\GeV}$}& ", zPtBins[iPtZ]);
  else
    cout << Form ("{$%g < \\pt < \\SI{%g}{\\GeV}$}& ", zPtBins[iPtZ], zPtBins[iPtZ+1]);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    if (iCent == 0)
      cout << "\\pp ";
    else
      cout << Form ("& Pb+Pb / %i-%i\\%% ", (int)centCuts[iCent], (int)centCuts[iCent-1]);

    for (short iSpc = 0; iSpc < 3; iSpc++)
      cout << Form ("& %g ", h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2));
    cout << "\\\\";

    if (iCent == numCentBins-1)
      cout << " \\hline";
    cout << endl << "\t\t\t";
  }
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Plot dPhi - zH 3d distribution
//////////////////////////////////////////////////////////////////////////////////////////////////
//void PhysicsAnalysis :: Plot3DDist () {
//  TCanvas* c = new TCanvas ("c_z_trk_pt_phi", "", 800, 600);
//  c->cd ();
//
//  c->SetLogy ();
//
//  h_z_trk_pt_phi[0][numCentBins-1][2]->GetXaxis ()->SetTitle ("#phi_{Z} - #phi_{Trk}");
//  h_z_trk_pt_phi[0][numCentBins-1][2]->GetYaxis ()->SetTitle ("#it{p}_{T}^{ ch} / #it{p}_{T}^{ Z}");
//
//  //h_z_trk_pt_phi[0][numCentBins-1][2]->RebinY (2);
//  h_z_trk_pt_phi[0][numCentBins-1][2]->Draw ("lego2");
//
//  c->SaveAs (Form ("%s/ZTrackCorr.pdf", plotPath.Data ()));
//}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - pTTrk 2d projections
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotCorrelations (const short pSpc, const short pPtZ, const bool _subBkg) {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue;

      const char* canvasName = Form ("c_z_trk_phi_pttrk_iPtZ%i_%s%s", iPtZ, spc, _subBkg ? "_subBkg":"");
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 1200, 900);
        gDirectory->Add (c);
        c->cd ();
        c->Divide (2, 2);
      }

      //for (short iCent = 0; iCent < 1; iCent++) {
      //  c->cd ();
      //  gPad->SetLogy ();
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd (iCent+1);
        GetDrawnObjects ();
        gPad->SetLogy ();

        gPad->SetTopMargin (0.01);
        gPad->SetBottomMargin (0.12);
        gPad->SetRightMargin (0.01);
        gPad->SetLeftMargin (0.12);

        double min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          TH1D* h = (!_subBkg ? h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent] : h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max) max = h->GetMaximum ();
        } // end loop over iPtTrk
        min *= 0.5;
        max = max <= 0 ? 1 : 14*max;
        SetMinAndMax (min, max);

        if (plotFill) {
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            TH1D* h = (TH1D*) (!_subBkg ? h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent] : h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent])->Clone ();

            h->GetYaxis ()->SetRangeUser (min, max);

            h->SetLineColor (fillColors[iPtTrk]);
            h->SetLineWidth (4);
            h->SetMarkerSize (0);

            h->GetXaxis ()->SetTitle ("#Delta#phi");
            //h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{ch}/d#Delta#phi (\"ZYAM\")");
            h->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

            h->GetXaxis ()->SetTitleOffset (0.6);
            h->GetYaxis ()->SetTitleOffset (0.8);
            h->GetXaxis ()->SetTitleSize (0.08);
            h->GetYaxis ()->SetTitleSize (0.06);
            h->GetXaxis ()->SetLabelSize (0.06);
            h->GetYaxis ()->SetLabelSize (0.06);

            h->DrawCopy (!canvasExists && iPtTrk == 0 ? "hist" : "same hist");
            //h->SetLineWidth (1);
            //h->Draw ("hist same");

            //TLine* line = new TLine (h->GetBinLowEdge (1), 0, h->GetBinLowEdge (h->GetNbinsX ()), 0);
            //line->Draw ("same");

            LabelCorrelations (iPtZ, iPtTrk, iCent);
          } // end loop over iPtTrk
          gPad->RedrawAxis ();
        } else {
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            TGAE* g = GetTGAE (!_subBkg ? h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent] : h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
            ResetXErrors (g);

            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetTitle ("#Delta#phi");
            //g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{ch}/d#Delta#phi (\"ZYAM\")");
            g->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

            g->SetMarkerStyle (markerStyle);
            g->SetLineColor (colors[iPtTrk]);
            g->SetMarkerColor (colors[iPtTrk]);

            g->GetXaxis ()->SetTitleOffset (0.6);
            g->GetYaxis ()->SetTitleOffset (0.8);
            g->GetXaxis ()->SetTitleSize (0.08);
            g->GetYaxis ()->SetTitleSize (0.06);
            g->GetXaxis ()->SetLabelSize (0.06);
            g->GetYaxis ()->SetLabelSize (0.06);

            g->Draw (!canvasExists && iPtTrk == 0 ? "AP" : "P");

            LabelCorrelations (iPtZ, iPtTrk, iCent);

            TLine* line1 = new TLine (phiLowBins[1], min, phiLowBins[1], max);
            TLine* line2 = new TLine (2*pi - phiLowBins[1], min, 2*pi - phiLowBins[1], max);
            TLine* line3 = new TLine (phiLowBins[2], min, phiLowBins[2], max);
            TLine* line4 = new TLine (2*pi - phiLowBins[2], min, 2*pi - phiLowBins[2], max);

            line1->SetLineStyle (2);
            line2->SetLineStyle (2);
            line3->SetLineStyle (2);
            line4->SetLineStyle (2);

            line1->Draw ("same");
            line2->Draw ("same");
            line3->Draw ("same");
            line4->Draw ("same");
          } // end loop over iPtTrk
        }
      } // end loop over centrality

      c->SaveAs (Form ("%s/dPhi_Distributions/dPhi_pTtrk_iPtZ%i_%s%s.pdf", plotPath.Data (), iPtZ, spc, _subBkg ? "_subBkg":""));
    } // end loop over pT^Z
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track dPhi distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelCorrelations (const short iPtZ, const short iPtTrk, const short iCent) {
  if (iCent == 0) {
    //myText (0.67, 0.26, kBlack, "#bf{#it{ATLAS}} Simulation", 0.05);
    myText (0.67, 0.26, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
    myText (0.67, 0.20, kBlack, "#it{pp}, 5.02 TeV", 0.05);
    //myText (0.16, 0.93, kBlack, "Data", 0.04);
    //myText (0.30, 0.93, kBlack, "Minbias", 0.05);
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = nullptr;

    const float pt_lo = ptTrkBins[iPtZ][iPtTrk];
    const float pt_hi = ptTrkBins[iPtZ][iPtTrk+1];
    //if (iPtTrk == 0)
    //  myText (0.3, 0.93, kBlack, "#it{p}_{T}^{ ch} [GeV]", 0.04);
    myText (0.2, 0.93-0.04*iPtTrk, colors[iPtTrk], Form ("%.1f < #it{p}_{T} < %.1f GeV", pt_lo, pt_hi), 0.04);
    //myText (0.3, 0.89-0.04*iPtTrk, colors[iPtTrk], Form ("(%.1f, %.1f)", pt_lo, pt_hi), 0.04);

    //b = TBoxNDC (0.33-0.024, 0.87-0.065*iPtTrk-0.016, 0.33+0.024, 0.87-0.065*iPtTrk+0.016);
    //b->SetFillColorAlpha (fillColors[iPtTrk], fillAlpha);
    //myMarkerTextNoLine (0.23, 0.89-0.04*iPtTrk, colors[iPtTrk], kFullCircle, "", 1.25, 0.04);

    //b->Draw ("l");
    //cPad->cd ();
  } else {
    myText (0.67, 0.20, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
  }

  if (iCent == 1) {
    if (iPtZ == nPtZBins-1) {
      myText (0.67, 0.88, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.67, 0.88, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies () {
  const char* canvasName = "c_trk_effs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 800);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (2, 2);
  }
  c->cd ();

  TGraphErrors* g0[numCentBins];
  TGraphErrors* g1[numCentBins];
  TGraphErrors* g2[numCentBins];
  TGraphErrors* g3[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
//    if (iCent > 0) continue;
    gPad->SetLogx ();

    g0[iCent] = new TGraphErrors (numEtaTrkBins);
    g1[iCent] = new TGraphErrors (numEtaTrkBins);
    g2[iCent] = new TGraphErrors (numEtaTrkBins);
    g3[iCent] = new TGraphErrors (numEtaTrkBins);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TEfficiency* eff = h_trk_effs[iCent][iEta];
      TGAE* eff = GetTGAE (h_trk_effs[iCent][iEta]);

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);
      eff->SetMarkerSize (0.5);

      eff->SetTitle (";#it{p}_{T} [GeV];Weighted Reco. Eff.");
      eff->GetXaxis ()->SetRangeUser (0.5, 65);
      eff->GetYaxis ()->SetRangeUser (0.3, 1.08);

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (iEta == 0 ? "AP" : "P");
      //eff->Draw (iEta == 0 ? "APL" : "LP same");

      //gPad->Update ();

      //eff->GetPaintedGraph ()->GetXaxis ()->SetRangeUser (0.5, 65);
      //eff->GetPaintedGraph ()->GetYaxis ()->SetRangeUser (0.3, 1.08);

      //TH1D* confInts = (TH1D*) h_trk_effs[iCent][iEta]->Clone ("confInts");
      TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2 + [3]*(log(x))^3", 1, 65);
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2", 2, 65);
      fit->SetParameter (0, 1);
      fit->SetParameter (1, 0);
      fit->SetParameter (2, 0);
      //fit->SetParameter (3, 0);
      h_trk_effs[iCent][iEta]->Fit (fit, "RN0Q");

      g0[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (0));
      g1[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (1));
      g2[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (2));
      g3[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (3));

      g0[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (0));
      g1[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (1));
      g2[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (2));
      g3[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (3));

      //fit->SetLineColor (colors[iEta]);
      //fit->SetFillColor (colors[iEta]);
      //fit->SetLineWidth (1);
      //fit->Draw ("same");
      //(TVirtualFitter::GetFitter ())->GetConfidenceIntervals (confInts, 0.68);
      //confInts->SetMarkerSize (0);
      //confInts->SetFillColorAlpha (colors[iEta], 0.3);
      //confInts->SetLineColor (colors[iEta]);
      //confInts->DrawCopy ("e3 same");
      //delete confInts;

      LabelTrackingEfficiencies (iCent, iEta);
    }
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies.pdf", plotPath.Data ()));

  //for (int iCent = 0; iCent < numCentBins; iCent++) {
  //  g0[iCent]->SetLineColor (colors[iCent]);
  //  g1[iCent]->SetLineColor (colors[iCent]);
  //  g2[iCent]->SetLineColor (colors[iCent]);
  //  g3[iCent]->SetLineColor (colors[iCent]);
  //  g0[iCent]->SetMarkerColor (colors[iCent]);
  //  g1[iCent]->SetMarkerColor (colors[iCent]);
  //  g2[iCent]->SetMarkerColor (colors[iCent]);
  //  g3[iCent]->SetMarkerColor (colors[iCent]);

  //  c->cd (1);
  //  gPad->SetLogx (false);
  //  g0[iCent]->Draw (iCent == 0 ? "AP":"P");
  //  c->cd (2);
  //  gPad->SetLogx (false);
  //  g1[iCent]->Draw (iCent == 0 ? "AP":"P");
  //  c->cd (3);
  //  gPad->SetLogx (false);
  //  g2[iCent]->Draw (iCent == 0 ? "AP":"P");
  //  c->cd (4);
  //  gPad->SetLogx (false);
  //  g3[iCent]->Draw (iCent == 0 ? "AP":"P");
  //}
  //c->SaveAs (Form ("%s/TrackingEfficienciesEtaDep.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies2D () {
  const char* canvasName = "c_trk_effs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 800);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (2, 2);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    //if (iCent > 0) continue;
    gPad->SetLogy ();

    //TEfficiency* eff = h2_trk_effs[iCent];
    TH2D* eff = h2_trk_effs[iCent];
    eff->GetZaxis ()->SetRangeUser (0.3, 1.00);

    eff->GetYaxis ()->SetMoreLogLabels ();

    eff->Draw ("colz");
    //eff->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();

    if (iCent == 0)
      myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
    else
      myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track efficiency plots
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrackingEfficiencies (const short iCent, const short iEta) {
  if (iCent == 0)
    myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
  else
    myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);

  if (iCent == 0) {
  //  myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
    myMarkerTextNoLine (0.5, 0.5-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.08);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities () {
  const char* canvasName = "c_trk_purs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 800);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (2, 2);
    c->Draw ();
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    gPad->SetLogx ();

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TGAE* pur = TEff2TGAE (h_trk_purs[iCent][iEta]);

      pur->SetLineColor (colors[iEta]);
      pur->SetMarkerColor (colors[iEta]);
      pur->SetMarkerSize (0.5);

      pur->SetTitle (";#it{p}_{T} [GeV];Primary Track Fraction");
      pur->GetXaxis ()->SetRangeUser (0.5, 65);
      pur->GetYaxis ()->SetRangeUser (0.97, 1.005);

      pur->GetXaxis ()->SetMoreLogLabels ();

      pur->Draw (iEta == 0 ? "AP" : "P");

      //TH1D* confInts = (TH1D*) h_trk_purs[iCent][iEta]->Clone ("confInts");
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2 + [3]*(log(x))^3", 1, 65);
      ////TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2", 2, 65);
      //fit->SetParameter (0, 1);
      //fit->SetParameter (1, 0);
      //fit->SetParameter (2, 0);
      ////fit->SetParameter (3, 0);
      //h_trk_purs[iCent][iEta]->Fit (fit, "RN0Q");
      //fit->SetLineColor (colors[iEta]);
      //fit->SetFillColor (colors[iEta]);
      //fit->SetLineWidth (1);
      //fit->Draw ("same");
      //(TVirtualFitter::GetFitter ())->GetConfidenceIntervals (confInts, 0.68);
      //confInts->SetMarkerSize (0);
      //confInts->SetFillColorAlpha (colors[iEta], 0.3);
      //confInts->SetLineColor (colors[iEta]);
      //confInts->DrawCopy ("e3 same");
      //delete confInts;

      LabelTrackingPurities (iCent, iEta);
    }
  }

  c->SaveAs (Form ("%s/TrackingPurities.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities2D () {
  const char* canvasName = "c_trk_purs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 800);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (2, 2);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    //if (iCent > 0) continue;
    gPad->SetLogy ();

    //TEfficiency* pur = h2_trk_purs[iCent];
    TH2D* pur = h2_trk_purs[iCent];
    pur->GetZaxis ()->SetRangeUser (0.96, 1.00);

    pur->GetYaxis ()->SetMoreLogLabels ();

    pur->Draw ("colz");
    //pur->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();

    if (iCent == 0)
      myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
    else
      myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
  }

  c->SaveAs (Form ("%s/TrackingPurities2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track purity plots
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrackingPurities (const short iCent, const short iEta) {
  if (iCent == 0)
    myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
  else
    myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);

  if (iCent == 0) {
  //  myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
    myMarkerTextNoLine (0.5, 0.5-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.08);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Subtracts default (HM) background from track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SubtractBackground (PhysicsAnalysis* a) {
  if (backgroundSubtracted)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        //******** Do subtraction of integrated dPhi plot ********//
        TH1D* h = (TH1D*) h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));

        TH1D* sub = a->h_z_trk_zpt[iSpc][iPtZ][iCent];
        if (sub != nullptr)
          h->Add (sub, -1);
        else
          cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;

        h_z_trk_zpt_sub[iSpc][iPtZ][iCent] = h;

        h = (TH1D*) h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //h->Sumw2 ();

        if (sub != nullptr) {
          h->Add (sub, -1);
          h->Divide (sub);
        }

        h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent] = h;


        h = (TH1D*) h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));

        sub = a->h_z_trk_zxzh[iSpc][iPtZ][iCent];
        if (sub != nullptr)
          h->Add (sub, -1);
        else
          cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;

        h_z_trk_zxzh_sub[iSpc][iPtZ][iCent] = h;

        h = (TH1D*) h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //h->Sumw2 ();

        if (sub != nullptr) {
          h->Add (sub, -1);
          h->Divide (sub);
        }

        h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent] = h;


        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          //******** Do subtraction of pT ********//
          h = (TH1D*) h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          //h->Sumw2 ();

          //h->Add (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
          sub = nullptr;
          if (a == nullptr) {
            //cout << "Info in SubtractBackground: Using HM background subtraction method." << endl;
            sub = h_z_trk_pt[iSpc][iPtZ][0][iCent];
          }
          else {
            //cout << "Info in SubtractBackground: Using MBM background subtraction method." << endl;
            sub = a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
          }
          if (sub != nullptr)
            h->Add (sub, -1);
          else
            cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;

          h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] = h;

          h = (TH1D*) h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          //h->Sumw2 ();

          if (sub != nullptr) {
            //h->Add (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
            h->Add (sub, -1);
            h->Divide (sub);
          }

          h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = h;


          //******** Do subtraction of z_h ********//
          h = new TH1D (Form ("h_z_trk_xzh_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nZHBins, zHBins);
          h->Sumw2 ();

          h->Add (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
          sub = nullptr;
          if (a == nullptr) {
            //cout << "Info in SubtractBackground: Using MBM background subtraction method." << endl;
            sub = h_z_trk_xzh[iSpc][iPtZ][0][iCent];
          }
          else {
            //cout << "Info in SubtractBackground: Using HM background subtraction method." << endl;
            sub = a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
          }
          if (sub != nullptr)
            h->Add (sub, -1);
          else
            cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;

          h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent] = h;

          h = new TH1D (Form ("h_z_trk_xzh_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nZHBins, zHBins);
          h->Sumw2 ();

          if (sub != nullptr) {
            h->Add (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
            h->Add (sub, -1);
            h->Divide (sub);
          }

          h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = h;
        } // end loop over phi


        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          //******** Do background subtraction of phi distributions ********//

          TH1D* h = (TH1D*) h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Clone (Form ("h_z_trk_phi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
          TH1D* sub = nullptr;
          if (a == nullptr) {
            //cout << "Info in SubtractBackground: Using HM background subtraction method." << endl;
            sub = h_z_trk_phi[iSpc][iPtZ][0][iCent];
          }
          else {
            //cout << "Info in SubtractBackground: Using MBM background subtraction method." << endl;
            sub = a->h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
            //for (int ix = 1; ix <= sub->GetNbinsX (); ix++) {
            //  sub->SetBinContent (ix, sub->GetBinContent (ix));
            //  sub->SetBinError (ix, sub->GetBinError (ix));
            //}
          }
          if (sub != nullptr) {
            while (sub->GetNbinsX () > h->GetNbinsX ()) {
              sub->Rebin (2);
            }
            h->Add (sub, -1);
          }
          else
            cout << "Error in SubtractBackground: Trying to subtract a null histogram!" << endl;

          h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent] = h;

        } // end loop over pT^trk
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  backgroundSubtracted = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Statistically subtracts same sign pairs
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SubtractSameSigns (PhysicsAnalysis* a) {
  if (sameSignsSubtracted)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {

          TH1D* h = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
          const float bkgCountsOverObsCounts = (a->h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1)) / (h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1));
          if (iCent == 0 && iPhi == 0 && iSpc == 2 && iPtZ == 2) {
            cout << bkgCountsOverObsCounts << endl;
            cout << "2nd to last bin before: " << h->GetBinContent (6) << endl;
          }

          h->Add (a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent], -bkgCountsOverObsCounts);
          h->Scale (1./ (1-bkgCountsOverObsCounts));
          //if (iCent == 0 && iPhi == 0 && iSpc == 2 && iPtZ == 2) {
          //  cout << "2nd to last bin after:  " << h->GetBinContent (6) << endl;
          //}

          h = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
          h->Add (a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent], -bkgCountsOverObsCounts);
          h->Scale (1./ (1-bkgCountsOverObsCounts));
        } // end loop over phi
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  sameSignsSubtracted = true;
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Converts analysis to a systematic variation by adding or subtracting statistical errors
// Only converts relevant histograms, e.g. for minbias will only do this to track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ConvertToStatVariation (const bool upVar, const float nSigma) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {

      // Hadron yield systematics, signal & signal+bkg levels
      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          AddStatVar (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);

          AddStatVar (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);

          AddStatVar (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        } // end loop over cents
      } // end loop over phi
      // IAA, ICP systematics
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          AddStatVar (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        } // end loop over cents
        for (short iCent = 2; iCent < numCentBins; iCent++) {
          AddStatVar (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTtrk binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotRawTrkYield (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_RawTrkYield_%s_%s_iPtZ%i", useTrkPt ? "pt" : "zh", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 400*numCentBins, 400);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* padName = Form ("p_RawTrkYield_%s_pad_%s_iPtZ%i_iCent%i", useTrkPt ? "pt" : "zh", spc, iPtZ, iCent);

        TPad* pad = nullptr;
        if (!canvasExists) {
          pad = new TPad (padName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, 1);

          gDirectory->Add (pad);

          pad->SetTopMargin (0.04);
          pad->SetBottomMargin (0.20);
          pad->SetLeftMargin (0.17);
          pad->SetRightMargin (0.06);
          pad->Draw ();
        }
        else {
          pad = dynamic_cast<TPad*> (gDirectory->Get (padName));
        }

        pad->cd ();
        GetDrawnObjects ();
        bool plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        double min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* h = useTrkPt ? h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over phi
        min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
        max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = 0; iPhi < 1; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = useTrkPt ? h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            if (useTrkPt) h->GetYaxis ()->SetTitle ("Y (#it{p}_{T}, #Delta#phi)");
            else h->GetYaxis ()->SetTitle ("Y (#it{x}_{zh}, #Delta#phi)");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleOffset (1.1);

            h->DrawCopy (plotNewAxes && iPhi == 0 ? "bar" : "bar same");
            //h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          }
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
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
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            if (useTrkPt) g->GetYaxis ()->SetTitle ("Y (#it{p}_{T}, #Delta#phi)");
            else g->GetYaxis ()->SetTitle ("Y (#it{x}_{zh}, #Delta#phi)");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleOffset (1.1);

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iPhi == 0 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == 0 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }
          } // end loop over phi
        }

        if (!canvasExists) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            if (iCent == 0)
              myText (0.22, 0.26, kBlack, "#it{pp}, 5.02 TeV", 0.06);
            else
              myText (0.22, 0.26, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

            if (iCent == 0)
              myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
            else if (iCent == numCentBins-1) {
              if (iPhi == 1) {
                TVirtualPad* cPad = gPad; // store current pad
                myText (0.653, 0.91, kBlack, "#Delta#phi", 0.054);
                TBox* b = TBoxNDC (0.598-0.018, 0.91-0.06*numPhiBins-0.016, 0.598+0.018, 0.91-0.06*numPhiBins+0.016);
                b->SetFillColorAlpha (fillColors[0], fillAlpha);
                b->Draw ("l");
                cPad->cd ();

                myText (0.65, 0.91-0.06*numPhiBins, kBlack, "MinBias", 0.054);
              }

              const char* lo = GetPiString (phiLowBins[iPhi]);
              const char* hi = GetPiString (phiHighBins[iPhi]);

              myMarkerTextNoLine (0.62, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.5, 0.054);
              myText (0.65, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.054);
            }
          }
        }
      } // end loop over cents
      
      c->SaveAs (Form ("%s/RawTrkYields/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrkYield (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (2*padRatio+1.0);
  const double mPadY = padRatio / (2*padRatio+1.0);
  //const double uPadY = 1.0 - mPadY - dPadY;
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_TrkYield_%s_%s_iPtZ%i", useTrkPt ? "pt" : "zh", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 395*numCentBins, 985);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* topPadName = Form ("p_TrkYield_%s_top_%s_iPtZ%i_iCent%i", useTrkPt ? "pt" : "zh", spc, iPtZ, iCent);
        const char* middlePadName = Form ("p_TrkYield_%s_middle_%s_iPtZ%i_iCent%i", useTrkPt ? "pt" : "zh", spc, iPtZ, iCent);
        const char* bottomPadName = Form ("p_TrkYield_%s_bottom_%s_iPtZ%i_iCent%i", useTrkPt ? "pt" : "zh", spc, iPtZ, iCent);

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
          TH1D* h = useTrkPt ? h_z_trk_pt[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over phi
        min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
        max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = 0; iPhi < 1; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = useTrkPt ? h_z_trk_pt[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (min, max);
            h->GetYaxis ()->SetRangeUser (10,100);

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            if (useTrkPt) h->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
            else h->GetYaxis ()->SetTitle ("dY/d#it{x}_{zh}d#Delta#phi");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (plotNewAxes && iPhi == 0 ? "bar" : "bar same");
            //h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          }
          gPad->RedrawAxis ();
        } else {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
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
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (min, max);
            g->GetYaxis ()->SetRangeUser (10,100);

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            if (useTrkPt) g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
            else g->GetYaxis ()->SetTitle ("dY/d#it{x}_{zh}d#Delta#phi");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iPhi == 1 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == 1 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }
          } // end loop over phi
        }

        if (!canvasExists || plotFill)
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++)
            LabelTrkYield (iCent, iPhi, iPtZ);

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
          TH1D* h = useTrkPt ? h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over phi
        float delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* h = useTrkPt ? h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            if (useTrkPt) h->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
            else h->GetYaxis ()->SetTitle ("dY/d#it{x}_{zh}d#Delta#phi");

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

            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
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
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            if (useTrkPt) g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
            else g->GetYaxis ()->SetTitle ("dY/d#it{x}_{zh}d#Delta#phi");

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

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iPhi == numPhiBins-1 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == numPhiBins-1 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }
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
          TH1D* h = useTrkPt ? h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
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
            TH1D* h = useTrkPt ? h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            h->GetYaxis ()->SetTitle ("(dY/d#Delta#phi) / (dY_{total}/d#Delta#phi)");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetXaxis ()->SetTitleOffset (2.6 * h->GetXaxis ()->GetTitleOffset ());
            h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

            h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
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
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            g->GetYaxis ()->SetTitle ("(dY/d#Delta#phi) / (dY_{bkg}/d#Delta#phi)");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
            g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iPhi == numPhiBins-1 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == numPhiBins-1 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }
          } // end loop over phi
        }
      } // end loop over cents
      
      c->SaveAs (Form ("%s/TrkYields/%s_dists_%s_iPtZ%i.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xZh", spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track pT distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrkYield (const short iCent, const short iPhi, const short iPtZ) {
  //const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);

  if (iCent == 0)
    myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
  else
    myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

  if (iCent == 0)
    myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
  else if (iCent == 1) {
    if (iPtZ == nPtZBins-1) {
      myText (0.485, 0.88, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.06);
    }
    else {
      myText (0.485, 0.88, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.06);
    }
  }
  else if (iCent == numCentBins-1) {
    if (iPhi == 1) {
      //myText (0.35, 0.91, kBlack, "Z-tagged", 0.054);
      //myText (0.44, 0.91, kBlack, "Reco.", 0.06);
      //myText (0.57, 0.91, kBlack, "Truth", 0.06);
      //myText (0.544, 0.91, kBlack, "Minbias", 0.054);
      //myText (0.703, 0.91, kBlack, "#Delta#phi", 0.054);

      TVirtualPad* cPad = gPad; // store current pad
      myText (0.653, 0.91, kBlack, "#Delta#phi", 0.054);
      TBox* b = TBoxNDC (0.598-0.018, 0.91-0.06*numPhiBins-0.016, 0.598+0.018, 0.91-0.06*numPhiBins+0.016);
      b->SetFillColorAlpha (fillColors[0], fillAlpha);
      b->Draw ("l");
      cPad->cd ();

      myText (0.65, 0.91-0.06*numPhiBins, kBlack, "Minimum Bias", 0.054);
    }

    const char* lo = GetPiString (phiLowBins[iPhi]);
    const char* hi = GetPiString (phiHighBins[iPhi]);

    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.612-0.032, 0.85-0.06*iPhi-0.016, 0.612+0.032, 0.85-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //myMarkerTextNoLine (0.462, 0.852-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.5, 0.054);
    //myText (0.70, 0.85-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.054);

    myMarkerTextNoLine (0.62, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.5, 0.054); // for plotting data vs bkg.

    //myMarkerTextNoLine (0.52, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.5, 0.054); // for plotting MC reco vs truth
    //myMarkerTextNoLine (0.62, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.5, 0.054);
    myText (0.65, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.054);
    //cPad->cd ();
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrkYieldZPt (const bool useTrkPt, const bool plotAsSystematic, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (2*padRatio+1.0);
  const double mPadY = padRatio / (2*padRatio+1.0);
  //const double uPadY = 1.0 - mPadY - dPadY;
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_TrkYield_pt_%s", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 395*numCentBins, 985);
      gDirectory->Add (c);
    }

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      c->cd ();

      const char* topPadName = Form ("p_TrkYield_pt_top_%s_iCent%i", spc, iCent);
      const char* middlePadName = Form ("p_TrkYield_pt_middle_%s_iCent%i", spc, iCent);
      const char* bottomPadName = Form ("p_TrkYield_pt_bottom_%s_iCent%i", spc, iCent);

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
      for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        TH1D* h = (useTrkPt ? h_z_trk_zpt[iSpc][iPtZ][iCent] : h_z_trk_zxzh[iSpc][iPtZ][iCent]);
        min = fmin (min, h->GetMinimum (0));
        max = fmax (max, h->GetMaximum ());
      } // end loop over phi
      min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
      max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
      SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_z_trk_zpt[iSpc][iPtZ][iCent] :  h_z_trk_zxzh[iSpc][iPtZ][iCent]);

          h->SetFillColorAlpha (fillColors[iPtZ], fillAlpha);
          h->SetMarkerSize (0);
          h->SetLineColor (kBlack);
          h->SetLineWidth (0);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, 65) : h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{Zh}");
          useTrkPt ? h->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("dY/d#it{x}d#Delta#phi");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

          h->DrawCopy (plotNewAxes && iPtZ == 0 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        }
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = (useAltMarker ? (iPtZ == 0 ? kOpenSquare : kOpenCircle) : (iPtZ == 0 ? kFullSquare : kFullCircle));
          TGAE* g = GetTGAE (useTrkPt ? h_z_trk_zpt[iSpc][iPtZ][iCent] : h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerColor (colors[iPtZ]);
            g->SetLineColor (colors[iPtZ]);
            g->SetMarkerSize (1);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ]);
            g->SetFillColorAlpha (fillColors[iPtZ], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, 65) : g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          g->GetYaxis ()->SetRangeUser (min, max);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
          useTrkPt ? g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("dY/d#it{x}d#Delta#phi");

          g->GetXaxis ()->SetTitleFont (43);
          g->GetXaxis ()->SetTitleSize (axisTextSize);
          g->GetXaxis ()->SetLabelFont (43);
          g->GetXaxis ()->SetLabelSize (axisTextSize);

          g->GetYaxis ()->SetTitleFont (43);
          g->GetYaxis ()->SetTitleSize (axisTextSize);
          g->GetYaxis ()->SetLabelFont (43);
          g->GetYaxis ()->SetLabelSize (axisTextSize);

          g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

          if (!plotAsSystematic) {
            string drawString = string (!canvasExists && iPtZ == 1 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ == 1 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }
        } // end loop over phi
      }

      if (!canvasExists || plotFill)
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          LabelTrkYieldZPt (iCent, iPtZ);
        }

      if (!plotSignal)
        continue;

      middlePad->cd ();
      GetDrawnObjects ();
      plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
      gPad->SetLogx ();
      gPad->SetLogy ();

      min = 1e30, max = 0;
      GetMinAndMax (min, max, true);
      for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        TH1D* h = h_z_trk_zpt_sub[iSpc][iPtZ][iCent];
        min = fmin (min, h->GetMinimum (0));
        max = fmax (max, h->GetMaximum ());
      } // end loop over phi
      float delta = log10 (max) - log10 (min);
      min = pow (10, log10 (min) - 0.1*delta);
      max = pow (10, log10 (max) + 0.1*delta);
      SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_z_trk_zpt_sub[iSpc][iPtZ][iCent] : h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);

          if (!h) continue;

          h->SetFillColorAlpha (fillColors[iPtZ], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          h->SetMarkerStyle (kFullCircle);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, 65) : h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{Zh}");
          useTrkPt ? h->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("dY/d#it{x}d#Delta#phi");

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

          h->DrawCopy (plotNewAxes && iPtZ == nPtZBins-1 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = (useAltMarker ? (iPtZ == 0 ? kOpenSquare : kOpenCircle) : (iPtZ == 0 ? kFullSquare : kFullCircle));

          TGAE* g = GetTGAE (useTrkPt ? h_z_trk_zpt_sub[iSpc][iPtZ][iCent] : h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerColor (colors[iPtZ]);
            g->SetLineColor (colors[iPtZ]);
            g->SetMarkerSize (1);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ]);
            g->SetFillColorAlpha (fillColors[iPtZ], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, 65) : g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          g->GetYaxis ()->SetRangeUser (min, max);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
          useTrkPt ? g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("dY/d#it{x}d#Delta#phi");

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

          if (!plotAsSystematic) {
            string drawString = string (!canvasExists && iPtZ == 0 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ == 0 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }
        } // end loop over phi
      }

      bottomPad->cd ();
      GetDrawnObjects ();
      plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
      gPad->SetLogx ();
      gPad->SetLogy ();

      min = 1e30, max = 0;
      GetMinAndMax (min, max, true);
      for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        TH1D* h = h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent];
        min = fmin (min, h->GetMinimum (0));
        max = fmax (max, h->GetMaximum ());
      } // end loop over phi
      //delta = max - min;
      //min = min - 0.3*delta;
      //max = max + 0.3*delta;
      delta = log10 (max) - log10 (min);
      min = pow (10, log10 (min) - 0.1*delta);
      max = pow (10, log10 (max) + 0.1*delta);
      SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent] : h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);

          if (!h) continue;

          h->SetFillColorAlpha (fillColors[iPtZ], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          h->SetMarkerStyle (kFullCircle);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, 65) : h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{Zh}");
          useTrkPt ? h->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("dY/d#it{x}d#Delta#phi");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (2.6 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

          h->DrawCopy (plotNewAxes && iPtZ == 0 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = (useAltMarker ? (iPtZ == 0 ? kOpenSquare : kOpenCircle) : (iPtZ == 0 ? kFullSquare : kFullCircle));
          TGAE* g = GetTGAE (useTrkPt ? h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent] : h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerColor (colors[iPtZ]);
            g->SetLineColor (colors[iPtZ]);
            g->SetMarkerSize (1);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ]);
            g->SetFillColorAlpha (fillColors[iPtZ], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, 65) : g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          g->GetYaxis ()->SetRangeUser (min, max);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
          useTrkPt ? g->GetYaxis ()->SetTitle ("dY/d#it{p}_{T}d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("dY/d#it{x}d#Delta#phi");

          g->GetXaxis ()->SetTitleFont (43);
          g->GetXaxis ()->SetTitleSize (axisTextSize);
          g->GetXaxis ()->SetLabelFont (43);
          g->GetXaxis ()->SetLabelSize (axisTextSize);

          g->GetYaxis ()->SetTitleFont (43);
          g->GetYaxis ()->SetTitleSize (axisTextSize);
          g->GetYaxis ()->SetLabelFont (43);
          g->GetYaxis ()->SetLabelSize (axisTextSize);

          g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
          g->GetYaxis ()->SetTitleOffset (2.2 * g->GetYaxis ()->GetTitleOffset ());

          if (!plotAsSystematic) {
            string drawString = string (!canvasExists && iPtZ == 0 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ == 0 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }
        } // end loop over phi
      }
    } // end loop over cents
    
    c->SaveAs (Form ("%s/TrkYields/pTTrk_dists_%s.pdf", plotPath.Data (), spc));
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track pT distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrkYieldZPt (const short iCent, const short iPtZ) {
  //const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);

  if (iCent == 0)
    myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
  else
    myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

  if (iCent == 0)
    myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
  else if (iCent == 1) {
    const char* lo = GetPiString (phiLowBins[1]);
    const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
    myText (0.485, 0.88, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.06);
  }
  else if (iCent == numCentBins-1) {
    if (iPtZ == 1) {
      //myText (0.35, 0.91, kBlack, "Z-tagged", 0.054);
      //myText (0.44, 0.91, kBlack, "Reco.", 0.06);
      //myText (0.57, 0.91, kBlack, "Truth", 0.06);
      //myText (0.544, 0.91, kBlack, "Minbias", 0.054);
      //myText (0.703, 0.91, kBlack, "#Delta#phi", 0.054);

      TVirtualPad* cPad = gPad; // store current pad
      //myText (0.653, 0.91, kBlack, "#Delta#phi", 0.054);
      TBox* b = TBoxNDC (0.598-0.018, 0.91-0.06*nPtZBins-0.016, 0.598+0.018, 0.91-0.06*nPtZBins+0.016);
      b->SetFillColorAlpha (fillColors[0], fillAlpha);
      b->Draw ("l");
      cPad->cd ();

      myText (0.65, 0.91-0.06*nPtZBins, kBlack, "Minimum Bias", 0.054);
    }

    myMarkerTextNoLine (0.56, 0.912-0.06*iPtZ, colors[iPtZ], kFullCircle, "", 1.5, 0.054); // for plotting data vs bkg.

    if (iPtZ == nPtZBins-1) {
      myText (0.56, 0.91-0.06*iPtZ, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.054);
    }
    else {
      myText (0.56, 0.91-0.06*iPtZ, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
    }
  }
}





////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CalculateIAA () {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (iaaCalculated)
    return;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      TH1D* ppHist = nullptr, *PbPbHist = nullptr;

      ppHist = h_z_trk_zpt_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        PbPbHist = (TH1D*)(h_z_trk_zpt_sub[iSpc][iPtZ][iCent]->Clone ());
        PbPbHist->Divide (ppHist);
        h_z_trk_zpt_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      }

      ppHist = h_z_trk_zxzh_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        PbPbHist = (TH1D*)(h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]->Clone ());
        PbPbHist->Divide (ppHist);
        h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      }

      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        TH1D* ppHist = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          if (!h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]) {
            TH1D* PbPbHist = (TH1D*)(h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            PbPbHist->Divide (ppHist);
            h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
          } 
        } // end loop over cents
        ppHist = h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          if (!h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]) {
            TH1D* PbPbHist = (TH1D*)(h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            PbPbHist->Divide (ppHist);
            h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
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
void PhysicsAnalysis :: PlotIAAdPhi (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_z_trk_%s_iaa_dPhi_%s_iPtZ%i", useTrkPt ? "pt" : "xzh", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-1), 502);
        c->Divide (numCentBins-1, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            TH1D* h = useTrkPt ? h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (0, 1.4);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            LabelIAAdPhi (iCent, iPhi, iPtZ);
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
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
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (0, 1.4);

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
            g->GetYaxis ()->SetTitle ("I_{AA}");
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

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

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iPhi == 1 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == 1 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }

            LabelIAAdPhi (iCent, iPhi, iPtZ);
          } // end loop over phi
        }
      } // end loop over cents
      c->cd (1);

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over cents

      c->SaveAs (Form ("%s/IAA/iaa_dPhi_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAAdPhi (const short iCent, const short iPhi, const short iPtZ) {

  if (iCent == 1)
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
  else if (iCent == 2) {
    if (iPtZ == nPtZBins-1) {
      myText (0.50, 0.80, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.50, 0.80, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
  else if (iCent == numCentBins-1) {
    if (iPhi == 1) {
      //myText (0.44, 0.91, kBlack, "MBM", 0.05);
      //myText (0.577, 0.91, kBlack, "HM", 0.05);
      //myText (0.44, 0.91, kBlack, "Data", 0.05);
      //myText (0.577, 0.91, kBlack, "MC", 0.05);
      myText (0.655, 0.91, kBlack, "#Delta#phi", 0.05);
    }
    const char* lo = GetPiString (phiLowBins[iPhi]);
    const char* hi = GetPiString (phiHighBins[iPhi]);
    //const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/16", phiLowBins[iPhi]*16/pi)) : "0";
    //const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/16", phiHighBins[iPhi]*16/pi)) : "#pi";
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerTextNoLine (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    //myMarkerTextNoLine (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myMarkerTextNoLine (0.63, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myText (0.65, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAAdCent (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_z_trk_%s_iaa_dCent_%s_iPtZ%i", useTrkPt ? "pt" : "xzh", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numPhiBins-1), 500);
        c->Divide (numPhiBins-1, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iCent = 1; iCent < numCentBins; iCent++) {
            TH1D* h = useTrkPt ? h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (0, 1.4);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            string drawString = string ("bar") + string (!canvasExists && iCent == 1 ? "" : " same");

            h->DrawCopy (drawString.c_str ());
            h->SetLineWidth (1);
            h->Draw ("hist same");

            LabelIAAdCent (iCent, iPhi, iPtZ);
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iCent = 1; iCent < numCentBins; iCent++) {
            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              deltaize (g, 1+((numCentBins-1)*((int)useAltMarker)-iCent)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
              g->SetLineColor (colors[iCent]);
              g->SetMarkerColor (colors[iCent]);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (1.2);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (1);
              g->SetLineColor (colors[iCent]);
              g->SetFillColorAlpha (fillColors[iCent], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (0, 1.4);
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iCent == 1 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iCent == 1 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }

            LabelIAAdCent (iCent, iPhi, iPtZ);
          } // end loop over phi
        }
      } // end loop over cents
      c->cd (1);

      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        const char* lo = GetPiString (phiLowBins[iPhi]);
        const char* hi = GetPiString (phiHighBins[iPhi]);
        myText (0.22, 0.24, kBlack, Form ("%s < #Delta#phi < %s", lo, hi), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over cents

      c->SaveAs (Form ("%s/IAA/iaa_dCent_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAAdCent (const short iCent, const short iPhi, const short iPtZ) {

  if (iPhi == 1) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
    if (iPtZ == nPtZBins-1) {
      myText (0.50, 0.80, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.50, 0.80, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
  else if (iPhi == numPhiBins-1) {
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerTextNoLine (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    //myMarkerTextNoLine (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myMarkerTextNoLine (0.63, 0.912-0.06*iCent, colors[iCent], kFullCircle, "", 1.4, 0.05);
    myText (0.65, 0.91-0.06*iCent, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAAdPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_z_trk_zpt_iaa_%s", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 500*(numCentBins-1), 502);
      c->Divide (numCentBins-1, 1);
      gDirectory->Add (c);
    }

    double xmin = 0, xmax = 0;
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      c->cd (iCent);
      gPad->SetLogx ();

      if (plotFill) {
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_z_trk_zpt_iaa[iSpc][iPtZ][iCent] : h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);

          h->SetFillColorAlpha (fillColors[iPtZ], fillAlpha);
          h->SetMarkerSize (0);
          h->SetLineColor (kBlack);
          h->SetLineWidth (0);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins]) : h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          h->GetYaxis ()->SetRangeUser (0, 1.4);
          xmin = h->GetXaxis ()->GetXmin ();
          xmax = h->GetXaxis ()->GetXmax ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{Zh}");
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

          h->DrawCopy (!canvasExists && iPtZ == 0 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");

          LabelIAAdPtZ (iCent, iPtZ);
        } // end loop over pT^Z bins
        gPad->RedrawAxis ();
      } else {
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

          TGAE* g = GetTGAE (useTrkPt ? h_z_trk_zpt_iaa[iSpc][iPtZ][iCent] : h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            deltaize (g, 1+((nPtZBins-1)*((int)useAltMarker)-iPtZ)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
            g->SetLineColor (colors[iPtZ]);
            g->SetMarkerColor (colors[iPtZ]);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (1.2);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ]);
            g->SetFillColorAlpha (fillColors[iPtZ], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins]) : g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
          g->GetYaxis ()->SetRangeUser (0, 1.4);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{Zh}");
          g->GetYaxis ()->SetTitle ("I_{AA}");
          xmin = g->GetXaxis ()->GetXmin ();
          xmax = g->GetXaxis ()->GetXmax ();

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

          if (!plotAsSystematic) {
            string drawString = string (!canvasExists && iPtZ == 0 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ == 0 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }

          LabelIAAdPtZ (iCent, iPtZ);
        } // end loop over pT^Z bins
      }
    } // end loop over cents
    c->cd (1);

    for (short iCent = 1; iCent < numCentBins; iCent++) {
      c->cd (iCent);
      myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    } // end loop over cents

    c->SaveAs (Form ("%s/IAA/iaa_dPtZ_%s.pdf", plotPath.Data (), spc));
    // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAAdPtZ (const short iCent, const short iPtZ) {

  if (iCent == 1)
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
  else if (iCent == 2) {
    const char* lo = GetPiString (phiLowBins[1]);
    const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
    myText (0.50, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.05);
  }
  else if (iCent == numCentBins-1) {
    myMarkerTextNoLine (0.63, 0.912-0.06*iPtZ, colors[iPtZ], kFullCircle, "", 1.4, 0.05);
    if (iPtZ == nPtZBins-1) {
      myText (0.65, 0.91-0.06*iPtZ, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.65, 0.91-0.06*iPtZ, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CalculateICP () {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (icpCalculated)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      TH1D* periphHist = nullptr, *centHist = nullptr;

      periphHist = h_z_trk_zpt_sub[iSpc][iPtZ][0];
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        centHist = (TH1D*)(h_z_trk_zpt_sub[iSpc][iPtZ][iCent]->Clone ());
        centHist->Divide (periphHist);
        h_z_trk_zpt_icp[iSpc][iPtZ][iCent] = centHist;
      }
      periphHist = h_z_trk_zxzh_sub[iSpc][iPtZ][0];
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        centHist = (TH1D*)(h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]->Clone ());
        centHist->Divide (periphHist);
        h_z_trk_zxzh_icp[iSpc][iPtZ][iCent] = centHist;
      }

      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        periphHist = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]) {
            centHist = (TH1D*)(h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent] = centHist;
          }
        } // end loop over cents

        periphHist = h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]) {
            centHist = (TH1D*)(h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent] = centHist;
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
void PhysicsAnalysis :: PlotICPdPhi (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_z_trk_%s_icp_dPhi_%s_iPtZ%i", useTrkPt ? "pt" : "xzh", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-2), 500);
        c->Divide (numCentBins-2, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            TH1D* h = useTrkPt ? h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (0, 2.4);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            LabelICPdPhi (iCent, iPhi, iPtZ);
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
            
            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
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
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (0, 2.4);
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iPhi == 1 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == 1 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }

            LabelICPdPhi (iCent, iPhi, iPtZ);
          } // end loop over phi
        }
      } // end loop over cents

      c->cd (1);
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over cents

      c->SaveAs (Form ("%s/ICP/icp_dPhi_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelICPdPhi (const short iCent, const short iPhi, const short iPtZ) {
  if (iCent == 2) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
    if (iPtZ == nPtZBins-1) {
      myText (0.50, 0.80, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.50, 0.80, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
  else if (iCent == numCentBins-1) {
    if (iPhi == 1) {
      //myText (0.55, 0.91, kBlack, "Data", 0.05);
      //myText (0.577, 0.91, kBlack, "MC", 0.05);
      //myText (0.44, 0.91, kBlack, "MBM", 0.05);
      //myText (0.577, 0.91, kBlack, "HM", 0.05);
      myText (0.655, 0.91, kBlack, "#Delta#phi", 0.05);
    }
    const char* lo = GetPiString (phiLowBins[iPhi]);
    const char* hi = GetPiString (phiHighBins[iPhi]);
    //const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/16", phiLowBins[iPhi]*16/pi)) : "0";
    //const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/16", phiHighBins[iPhi]*16/pi)) : "#pi";
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.71-0.024, 0.91-0.06*iPhi-0.016, 0.71+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerTextNoLine (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    myMarkerTextNoLine (0.612, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myText (0.65, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotICPdCent (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_z_trk_%s_icp_dCent_%s_iPtZ%i", useTrkPt ? "pt" : "xzh", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numPhiBins-1), 500);
        c->Divide (numPhiBins-1, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iCent = 2; iCent < numCentBins; iCent++) {
            TH1D* h = useTrkPt ? h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            if (useTrkPt) h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else h->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            h->GetYaxis ()->SetRangeUser (0, 2.4);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else h->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            h->DrawCopy (!canvasExists && iCent == 2 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");

            LabelICPdCent (iCent, iPhi, iPtZ);
          } // end loop over phi
          gPad->RedrawAxis ();
        } else {
          for (int iCent = 2; iCent < numCentBins; iCent++) {
            const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
            
            TGAE* g = GetTGAE (useTrkPt ? h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent] : h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              deltaize (g, 1+((numCentBins-2)*((int)useAltMarker)-iCent)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
              g->SetLineColor (colors[iCent-1]);
              g->SetMarkerColor (colors[iCent-1]);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (1.2);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (1);
              g->SetLineColor (colors[iCent-1]);
              g->SetFillColorAlpha (fillColors[iCent-1], 0.3);
            }

            if (useTrkPt) g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins]);
            else g->GetXaxis ()->SetLimits (zHBins[0], zHBins[nZHBins]);
            g->GetYaxis ()->SetRangeUser (0, 2.4);
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

            g->GetXaxis ()->SetMoreLogLabels ();

            if (useTrkPt) g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            else g->GetXaxis ()->SetTitle ("#it{x}_{zh}");
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

            if (!plotAsSystematic) {
              string drawString = string (!canvasExists && iCent == 2 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iCent == 2 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }

            LabelICPdCent (iCent, iPhi, iPtZ);
          } // end loop over phi
        }
      } // end loop over cents

      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        const char* lo = GetPiString (phiLowBins[iPhi]);
        const char* hi = GetPiString (phiHighBins[iPhi]);
        myText (0.22, 0.24, kBlack, Form ("%s < #Delta#phi < %s", lo, hi), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over cents

      c->SaveAs (Form ("%s/ICP/icp_dCent_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelICPdCent (const short iCent, const short iPhi, const short iPtZ) {
  if (iPhi == 1) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
    if (iPtZ == nPtZBins-1) {
      myText (0.50, 0.80, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.50, 0.80, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
  else if (iPhi == numPhiBins-1) {
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.71-0.024, 0.91-0.06*iPhi-0.016, 0.71+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerTextNoLine (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    myMarkerTextNoLine (0.582, 0.912-0.06*(iCent-1), colors[iCent-1], kFullCircle, "", 1.4, 0.05);
    myText (0.6, 0.91-0.06*(iCent-1), kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.05);
  }
}

#endif
