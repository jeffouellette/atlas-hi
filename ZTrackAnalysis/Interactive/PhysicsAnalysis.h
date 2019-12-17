#ifndef __PhysicsAnalysis_h__
#define __PhysicsAnalysis_h__

#include "Params.h"

#include <ArrayTemplates.h>

#include <AtlasUtils.h>

#include <TEfficiency.h>
#include <TClass.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TChain.h>
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
  bool backgroundSubtracted = false;
  bool sameSignsSubtracted = false;
  bool hasBkg = true;

  vector<TH1*> drawnHists;
  vector<TGAE*> drawnGraphs;

  bool iaaCalculated = false;
  //bool icpCalculated = false;

  TFile* eventWeightsFile = nullptr;
  string eventWeightsFileName = "DataAnalysis/Nominal/eventWeightsFile.root";
  bool eventWeightsLoaded = false;
  TFile* trkEffFile = nullptr;
  bool effsLoaded   = false;
  TFile* trkPurFile = nullptr;
  bool pursLoaded   = false;
  TFile* trigEffFile = nullptr;
  bool trigEffsLoaded = false;
  TFile* histFile   = nullptr;
  bool histsLoaded  = false;
  bool histsScaled  = false;

  public:
  bool histsUnfolded = false;
  bool plotFill       = false; // whether to plot as filled (bar) graph or points w/ errors
  bool plotSignal     = true; // whether to plot background subtracted plots
  bool useAltMarker   = false; // whether to plot as open markers (instead of closed)

  bool isMC           = false;
  bool is2015Conds    = false; // whether this analysis uses 2015 data (different conditions)
  bool useHITight     = false; // whether to use HITight tracking efficiencies
  bool useHijingEffs  = false; // whether to use tracking efficiencies derived from Hijing
  bool doLeptonRejVar = false; // whether to impose an additional dR cut on tracks away from the leptons
  bool doTrackPurVar  = false; // whether to impose an additional correction based on tracking purity
  bool doTrackEffVar  = false; // whether to use pions-only tracking efficiency variation
  float trkEffNSigma  = 0; // how many sigma to vary the track efficiency by (-1,0,+1 suggested)
  float trkPurNSigma  = 0; // how many sigma to vary the track purity by (-1,0,+1 suggested)

  //// Event info distributions (for reweighting)
  TH1D*** h_PbPbFCal_weights   = Get2DArray <TH1D*> (3, nPtZBins+1);
  TH1D**** h_PbPbQ2_weights    = Get3DArray <TH1D*> (3, numFineCentBins, nPtZBins+1);
  TH1D**** h_PbPbPsi2_weights  = Get3DArray <TH1D*> (3, numFineCentBins, nPtZBins+1);
  //TH1D* h_ppNch_weights        = nullptr;

  // Efficiencies
  TH1D*** h_trk_effs          = Get2DArray <TH1D*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  //TF1**   f_trk_effs          = Get1DArray <TH1D*> (numTrkCorrCentBins); // iCent (eta dependence is extrapolated)
  TH2D**  h2_trk_effs         = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  //TEfficiency*** h_trk_effs   = Get2DArray <TEfficiency*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  TH2D** h2_num_trk_effs      = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  TH2D** h2_den_trk_effs      = Get1DArray <TH2D*> (numTrkCorrCentBins);

  // Tracking purities
  TH1D*** h_trk_purs          = Get2DArray <TH1D*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  TH2D**  h2_trk_purs         = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  //TEfficiency*** h_trk_purs   = Get2DArray <TEfficiency*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  TH2D** h2_num_trk_purs      = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  TH2D** h2_den_trk_purs      = Get1DArray <TH2D*> (numTrkCorrCentBins);

  // Bin migration factors
  TF1**** f_z_trk_zpt_binMigration  = Get3DArray <TF1*> (2, nPtZBins, numBBBCorrCentBins); // iSpc, iPtZ, iCent
  TF1**** f_z_trk_zxzh_binMigration = Get3DArray <TF1*> (2, nPtZBins, numBBBCorrCentBins); // iSpc, iPtZ, iCent

  // Physics plots
  TH1D*****   h_trk_dphi        = Get4DArray <TH1D*> (3, nPtZBins, maxNPtTrkBins, numCentBins); // iSpc, iPtZ, iPtTrk, iCent
  TH1D*****   h_trk_dphi_sub    = Get4DArray <TH1D*> (3, nPtZBins, maxNPtTrkBins, numCentBins); // iSpc, iPtZ, iPtTrk, iCent

  TH1D*****   h_trk_pt_dphi_raw = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);   // iSpc, iPtZ, iPhi, iCent
  TH1D****    h_z_counts        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);               // iSpc, iPtZ, iCent

  TH1D*****  h_trk_pt_dphi            = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_sub        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_sig_to_bkg = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_iaa        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  //TH1D*****  h_trk_pt_dphi_icp        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_trk_pt_ptz             = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_sub         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_sig_to_bkg  = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_iaa         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  //TH1D****   h_trk_pt_ptz_icp         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  TH1D***** h_trk_xhz_dphi            = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_sub        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_sig_to_bkg = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_iaa        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  //TH1D***** h_trk_xhz_dphi_icp        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_trk_xhz_ptz            = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_xhz_ptz_sub        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_xhz_ptz_sig_to_bkg = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_xhz_ptz_iaa        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  //TH1D****   h_trk_xhz_ptz_icp        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  PhysicsAnalysis () { }

  PhysicsAnalysis (const char* _name) {
    name = _name;
    plotFill = false;
  }

  virtual ~PhysicsAnalysis () {
    Delete2DArray (h_PbPbFCal_weights,  3, nPtZBins+1);
    Delete3DArray (h_PbPbQ2_weights,    3, numFineCentBins, nPtZBins+1);
    Delete3DArray (h_PbPbPsi2_weights,  3, numFineCentBins, nPtZBins+1);

    ClearHists ();

    Delete2DArray (h_trk_effs,      numTrkCorrCentBins, numEtaTrkBins);
    Delete1DArray (h2_trk_effs,     numTrkCorrCentBins);
    Delete1DArray (h2_num_trk_effs, numTrkCorrCentBins);
    Delete1DArray (h2_den_trk_effs, numTrkCorrCentBins);

    Delete2DArray (h_trk_purs,      numTrkCorrCentBins, numEtaTrkBins);
    Delete1DArray (h2_trk_purs,     numTrkCorrCentBins);
    Delete1DArray (h2_num_trk_purs, numTrkCorrCentBins);
    Delete1DArray (h2_den_trk_purs, numTrkCorrCentBins);

    Delete3DArray (f_z_trk_zpt_binMigration,  3, nPtZBins, numBBBCorrCentBins);
    Delete3DArray (f_z_trk_zxzh_binMigration, 3, nPtZBins, numBBBCorrCentBins);

    Delete4DArray (h_trk_dphi,                3, nPtZBins, maxNPtTrkBins, numCentBins);
    Delete4DArray (h_trk_dphi_sub,            3, nPtZBins, maxNPtTrkBins, numCentBins);

    Delete4DArray (h_trk_pt_dphi_raw,         3, nPtZBins, numPhiBins, numCentBins);
    Delete3DArray (h_z_counts,                3, nPtZBins, numCentBins);

    Delete4DArray (h_trk_pt_dphi,             3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_dphi_sub,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_dphi_sig_to_bkg,  3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_dphi_iaa,         3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (h_trk_pt_dphi_icp,         3, nPtZBins, numPhiBins, numCentBins);

    Delete3DArray (h_trk_pt_ptz,              3, nPtZBins, numCentBins);
    Delete4DArray (h_trk_pt_ptz_sub,          3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_ptz_sig_to_bkg,   3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_ptz_iaa,          3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (h_trk_pt_ptz_icp,          3, nPtZBins, numPhiBins, numCentBins);

    Delete4DArray (h_trk_xhz_dphi,            3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_dphi_sub,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_dphi_sig_to_bkg, 3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_dphi_iaa,        3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (h_trk_xhz_dphi_icp,        3, nPtZBins, numPhiBins, numCentBins);

    Delete4DArray (h_trk_xhz_ptz,             3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_ptz_sub,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_ptz_sig_to_bkg,  3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_ptz_iaa,         3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (h_trk_xhz_ptz_icp,         3, nPtZBins, numPhiBins, numCentBins);
  }

  protected:
  void LabelTrackingEfficiencies (const short iCent, const short iEta);
  void LabelTrackingPurities (const short iCent, const short iEta);
  void LabelCorrelations (const short iPtZ, const short iPtTrk, const short iCent, const bool subBkg);
  void LabelTrkYield (const short iCent, const short iPhi, const short iPtZ, const short iSpc);
  void LabelTrkYieldZPt (const short iCent, const short iPtZ, const short iSpc);
  void LabelIAAdPhi (const short iCent, const short iPhi, const short iPtZ);
  void LabelIAAdCent (const short iCent, const short iPhi, const short iPtZ);
  void LabelIAAdPtZ (const short iCent, const short iPtZ);
  //void LabelICPdPhi (const short iCent, const short iPhi, const short iPtZ);
  //void LabelICPdCent (const short iCent, const short iPhi, const short iPtZ);
  //void LabelICPdPtZ (const short iCent, const short iPtZ);

  void GetDrawnObjects ();
  void GetMinAndMax (double &min, double &max, const bool log = false);
  void SetMinAndMax (double min, double max);

  public:
  string Name () { return name; }
  void SetName (string _name) { name = _name; }

  virtual TGAE* GetTGAE (TH1D* h);

  virtual void CreateHists ();
  virtual void CopyAnalysis (PhysicsAnalysis* a, const bool copyBkgs = false);
  virtual void ClearHists ();
  virtual void CombineHists ();
  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true);
  virtual void SaveHists (const char* histFileName = "savedHists.root");
  virtual void SaveResults (const char* saveFileName = "resultsHists.root");
  virtual void ScaleHists ();
  virtual void Execute (const char* inFileName, const char* outFileName);
  virtual void GenerateWeights (const char* weightedSampleInFileName, const char* matchedSampleInFileName, const char* outFileName);
  virtual void LoadEventWeights ();
  virtual void SubtractBackground (PhysicsAnalysis* a = nullptr);
  virtual void UnfoldSubtractedYield ();
  virtual void InflateStatUnc (const float amount);
  virtual void SubtractSameSigns (PhysicsAnalysis* a);

  virtual void ApplyRelativeVariation (float**** relVar, const bool upVar = true); // multiplies yield results by relErr in each bin (or divides if not upVar)
  virtual void ConvertToStatVariation (const bool upVar = true, const float nSigma = 1.); // adds or subtracts nSigma of statistical errors to analysis

  virtual void LoadTrackingEfficiencies (const bool doRebin = false); // defaults to HILoose
  virtual double GetTrackingEfficiency (const float fcal_et, float trk_pt, const float trk_eta, const bool isPbPb = true);

  virtual void LoadTrackingPurities (const bool doRebin = false); // defaults to HILoose
  virtual double GetTrackingPurity (const float fcal_et, float trk_pt, const float trk_eta, const bool isPbPb = true);

  void PrintZYields (const int iPtZ = 2);

  void PlotCorrelations (const short pSpc = 2, const short pPtZ = nPtZBins-1, const bool _subBkg = false);
  void PlotTrackingEfficiencies (PhysicsAnalysis* a = nullptr);
  void PlotTrackingEfficiencies2D ();
  void PlotTrackingPurities (PhysicsAnalysis* a = nullptr);
  void PlotTrackingPurities2D ();

  void CalculateIAA ();
  //void CalculateICP ();

  virtual void PlotRawTrkYield (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotTrkYield (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotTrkYieldZPt (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotTrkYieldZPtSpcComp (const bool useTrkPt = true, const bool plotAsSystematic = false);
  virtual void PlotIAAdPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAAdCent (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAAdPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotSingleIAAdPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pPtZ = -1, const short iCent = numCentBins-1, const short pSpc = 2);
  virtual void PlotIAASpcComp (const bool useTrkPt = true, const bool plotAsSystematic = false);
  //virtual void PlotICPdPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //virtual void PlotICPdCent (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //virtual void PlotICPdPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotSignalToBkg (const bool useTrkPt = true, const short iSpc = 2);
  //virtual void PlotPullDist (const bool useTrkPt = true);

  virtual void WriteIAAs ();
  virtual void PrintIAA (const bool printErrs, const bool useTrkPt = true, const short iCent = numCentBins-1, const short iPtZ = nPtZBins-1, const short iSpc = 2);
};


void SafeWrite (TObject* tobj, TString nname = "") {
  if (tobj)
    tobj->Write (nname);
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
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]); // old name: h_z_trk_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]); // old name: h_z_trk_xzh
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        }
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent] = new TH1D (Form ("h_trk_dphi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()), "", 80, -pi/2, 3*pi/2); // old name: h_z_trk_phi
          h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Sumw2 ();
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

  ClearHists ();

  // Should clone these histograms
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        //for (short iZH = 0; iZH < nXHZBins[iPtZ]; iZH++) {
        //h_trk_pt_dphi_phi[iPtZ][iCent][iSpc] = (TH2D*) a->h_trk_pt_dphi_phi[iPtZ][iCent][iSpc]->Clone (Form ("h_trk_pt_dphi_phi_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //}
        h_trk_pt_ptz[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h_trk_xhz_ptz[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
        } // end loop over iPhi
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Clone (Form ("h_trk_dphi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ())); // old name: h_z_trk_phi
        }
      } // end loop over iPtZ
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_counts[iSpc][iPtZ][iCent]->Clone (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  if (copyBkgs) {
    if (a->backgroundSubtracted) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

            if (a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]) {
              h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_sub
              h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_sig_to_bkg
            }
            if (a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]) {
              h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_sub
              h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_sig_to_bkg
            }

            for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
              if (a->h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]) {
                h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_sub
                h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_sig_to_bkg
              }
              if (a->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]) {
                h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_sub
                h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_sig_to_bkg
              }
            } // end loop over phi

            for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
              if (a->h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]) {
                h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) a->h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ())); // old name: h_z_trk_dphi
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
        for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {

            if (a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]) {
              if (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]) SaferDelete (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]);
              h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_iaa
            }
            if (a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]) {
              if (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]) SaferDelete (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]);
              h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_iaa
            }

            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              if (a->h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]) {
                if (h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
                h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_iaa
              }
              if (a->h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]) {
                if (h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
                h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_iaa
              }
            } // end loop over phi
          } // end loop over cents
        } // end loop over pT^Z bins
      } // end loop over species
      iaaCalculated = true;
    }

    /*
    if (a->icpCalculated) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          for (short iCent = 2; iCent < numCentBins; iCent++) {

            if (a->h_trk_pt_ptz_icp[iSpc][iPtZ][iCent]) {
              if (h_trk_pt_ptz_icp[iSpc][iPtZ][iCent]) SaferDelete (h_trk_pt_ptz_icp[iSpc][iPtZ][iCent]);
              h_trk_pt_ptz_icp[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz_icp[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_icp_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }
            if (a->h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]) {
              if (h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]) SaferDelete (h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]);
              h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_icp_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
            }

            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              if (a->h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]) {
                if (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]);
                h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_icp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
              }
              if (a->h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent]) {
                if (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]);
                h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_icp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
              }
            } // end loop over phi
          } // end loop over cents
        } // end loop over pT^Z bins
      } // end loop over species
      icpCalculated = true;
    }
    */
  }

  histsLoaded = true;
  histsScaled = true;
  return;    
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Clears histograms from memory (allows them to be overwritten).
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ClearHists () {
  // Should clone these histograms
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        //for (short iZH = 0; iZH < nXHZBins[iPtZ]; iZH++) {
        //h_trk_pt_dphi_phi[iPtZ][iCent][iSpc] = (TH2D*) a->h_trk_pt_dphi_phi[iPtZ][iCent][iSpc]->Clone (Form ("h_trk_pt_dphi_phi_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //}
        if (h_trk_pt_ptz[iSpc][iPtZ][iCent])   SaferDelete (h_trk_pt_ptz[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz[iSpc][iPtZ][iCent])  SaferDelete (h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent])      SaferDelete (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent])     SaferDelete (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]) SaferDelete (h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]);
        }
      } // end loop over iPtZ
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        if (h_z_counts[iSpc][iPtZ][iCent])  SaferDelete (h_z_counts[iSpc][iPtZ][iCent]);
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

        if (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent])         SaferDelete (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]);
        if (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent])  SaferDelete (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent])        SaferDelete (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]) SaferDelete (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]);

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent])          SaferDelete (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent])   SaferDelete (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent])         SaferDelete (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over phi
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]) SaferDelete (h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]);
        } // end loop over pT^trk
      } // end loop over pT^Z bins
    } // end loop over cents
  } // end loop over species
  backgroundSubtracted = false;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        if (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent])   SaferDelete (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent])  SaferDelete (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]);

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over phi
      } // end loop over cents
    } // end loop over pT^Z bins
  } // end loop over species
  iaaCalculated = false;
  

  /*
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        if (h_trk_pt_ptz_icp[iSpc][iPtZ][iCent])   SaferDelete (h_trk_pt_ptz_icp[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent])  SaferDelete (h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]);

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over phi
      } // end loop over cents
    } // end loop over pT^Z bins
  } // end loop over species
  icpCalculated = false;
  */

  histsLoaded = false;
  histsScaled = false;

  if (histFile)
    SaferDelete (histFile);
  return;    

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  SetupDirectories ("", "ZTrackAnalysis/");
  //if (histsLoaded)
  ClearHists ();

  TDirectory* _gDirectory = gDirectory;
  histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "read");

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) histFile->Get (Form ("h_trk_dphi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ())); // old name: h_z_trk_phi
        }

        h_trk_pt_ptz[iSpc][iPtZ][iCent]  = new TH1D (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]); // old name: h_z_trk_zpt
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Sumw2 ();
        h_trk_xhz_ptz[iSpc][iPtZ][iCent] = new TH1D (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]); // old name: h_z_trk_zxzh
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Sumw2 ();

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]     = (TH1D*) histFile->Get (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]    = (TH1D*) histFile->Get (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
          if (iSpc != 2 && iPhi != 0) {
            if (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent])  h_trk_pt_ptz[iSpc][iPtZ][iCent]->Add   (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
            if (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]) h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Add  (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          }
        }
      }  
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      }
    }
  }

  histsLoaded = true;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h_z_counts[2][iPtZ][iCent]->Add (h_z_counts[iSpc][iPtZ][iCent]);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]) h_trk_pt_dphi_raw[2][iPtZ][iPhi][iCent]->Add (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
        }
      } // end loop over centralities
    } // end loop over pT^Z
  } // end loop over species

  if (_finishHists) {
    //PhysicsAnalysis :: CombineHists (); // deprecated function call
    PhysicsAnalysis :: ScaleHists ();
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SaveHists (const char* histFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");
  if (!histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  if (!histFile) {
    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "recreate");
    histFile->cd ();
  }

  cout << "Saving histograms to " << Form ("%s/%s", rootPath.Data (), histFileName) << endl;

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          SafeWrite (h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]);
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SafeWrite (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
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
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SaveResults (const char* saveFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");
  if (!histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  TFile* saveFile = new TFile (Form ("%s/%s", rootPath.Data (), saveFileName), "recreate");
  saveFile->cd ();

  cout << "Saving results to " << Form ("%s/%s", rootPath.Data (), saveFileName) << endl;

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          SafeWrite (h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent],      Form ("h_ztrk_phi_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
          SafeWrite (h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent],  Form ("h_ztrk_phi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
        }
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  SafeWrite (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
        //  SafeWrite (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
        //  SafeWrite (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
        //}

        SafeWrite (h_trk_pt_ptz[iSpc][iPtZ][iCent],                Form ("h_ztrk_pt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent],            Form ("h_ztrk_pt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent],     Form ("h_ztrk_pt_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz[iSpc][iPtZ][iCent],               Form ("h_ztrk_xhz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent],           Form ("h_ztrk_xhz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent],    Form ("h_ztrk_xhz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));

        SafeWrite (h_z_counts[iSpc][iPtZ][iCent],                 Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      }
    }
  }

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent],            Form ("h_ztrk_pt_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent],           Form ("h_ztrk_xhz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      }
    }
  }
  
  saveFile->Close ();
  saveFile = nullptr;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iSpc = 0; iSpc < 2; iSpc++) {

        // Gets the weighting factor needed for this species.
        // E.g. if there are 2 muon events and 1 electron event,
        // the per-Z yield should be weighted by 2/3 in the muon
        // channel and 1/3 in the electron channel.
        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
        float spcWeight = countsHist->GetBinContent (1);
        countsHist = h_z_counts[2][iPtZ][iCent];
        if (countsHist->GetBinContent (1) > 0)
          spcWeight = spcWeight / countsHist->GetBinContent (1);
        else {
          cout << "Warning: In PhysicsAnalysis :: CombineHists: Found 0 total Z bosons in this bin, weight is set to 0!" << endl;
          spcWeight = 0;
        }

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          while (h_trk_dphi[2][iPtZ][iPtTrk][iCent]->GetNbinsX () > h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->GetNbinsX ())
            h_trk_dphi[2][iPtZ][iPtTrk][iCent]->Rebin (2);
          if (h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent])     h_trk_dphi[2][iPtZ][iPtTrk][iCent]->Add      (h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent], spcWeight);

          if (hasBkg) {
            while (h_trk_dphi_sub[2][iPtZ][iPtTrk][iCent]->GetNbinsX () > h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]->GetNbinsX ())
              h_trk_dphi_sub[2][iPtZ][iPtTrk][iCent]->Rebin (2);
            if (h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]) h_trk_dphi_sub[2][iPtZ][iPtTrk][iCent]->Add  (h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent], spcWeight);
          }
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent])              h_trk_pt_dphi[2][iPtZ][iPhi][iCent]->Add             (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent], spcWeight);
          if (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent])             h_trk_xhz_dphi[2][iPtZ][iPhi][iCent]->Add            (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], spcWeight);
          if (hasBkg) {
            if (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent])          h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]->Add         (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent], spcWeight);
            if (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent])   h_trk_pt_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Add  (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], spcWeight);
            if (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent])         h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]->Add        (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent], spcWeight);
            if (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent])  h_trk_xhz_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Add (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], spcWeight);
          }
        } // end loop over phi

        if (h_trk_pt_ptz[iSpc][iPtZ][iCent])             h_trk_pt_ptz[2][iPtZ][iCent]->Add              (h_trk_pt_ptz[iSpc][iPtZ][iCent], spcWeight);
        if (h_trk_xhz_ptz[iSpc][iPtZ][iCent])            h_trk_xhz_ptz[2][iPtZ][iCent]->Add             (h_trk_xhz_ptz[iSpc][iPtZ][iCent], spcWeight);
        if (hasBkg) {
          if (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent])         h_trk_pt_ptz_sub[2][iPtZ][iCent]->Add          (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], spcWeight);
          if (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent])  h_trk_pt_ptz_sig_to_bkg[2][iPtZ][iCent]->Add   (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent], spcWeight);
          if (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent])        h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Add         (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], spcWeight);
          if (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]) h_trk_xhz_ptz_sig_to_bkg[2][iPtZ][iCent]->Add  (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent], spcWeight);
        }
      } // end loop over species
    } // end loop over pT^Z
  } // end loop over centralities

  InflateStatUnc (0.54);
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
        const float counts = countsHist->GetBinContent (1);
        //const float countsErr = countsHist->GetBinError (1);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          const double countsdPhi = counts * (phiHighBins[iPhi]-phiLowBins[iPhi]);

          if (countsdPhi > 0) {
            h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Scale  (1. / countsdPhi, "width");
            h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
          }
        } // end loop over phi

        if (counts > 0) {
          h_trk_pt_ptz[iSpc][iPtZ][iCent]->Scale   (1./ (counts * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
          h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Scale  (1./ (counts * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
        }
        
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          TH1D* h = h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent];
          h->Rebin (2);
          if (iPtTrk > 3)
            h->Rebin (2);
          if (iCent != 0)
            h->Rebin (2);
          if (counts > 0)
            h->Scale (1. / counts);
        }
      } // end loop over pT^Z
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

  SetupDirectories ("", "ZTrackAnalysis/");

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
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);
    PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
    PbPbTree->SetBranchAddress ("z_y",          &z_y);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("z_m",          &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",        &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",       &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",       &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",    &l1_charge);
    PbPbTree->SetBranchAddress ("l1_trk_pt",    &l1_trk_pt);
    PbPbTree->SetBranchAddress ("l1_trk_eta",   &l1_trk_eta);
    PbPbTree->SetBranchAddress ("l1_trk_phi",   &l1_trk_phi);
    PbPbTree->SetBranchAddress ("l2_pt",        &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",       &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",       &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",    &l2_charge);
    PbPbTree->SetBranchAddress ("l2_trk_pt",    &l2_trk_pt);
    PbPbTree->SetBranchAddress ("l2_trk_eta",   &l2_trk_eta);
    PbPbTree->SetBranchAddress ("l2_trk_phi",   &l2_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",         &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      if (fabs (vz) > 150)
        continue;

      if (event_weight == 0)
        continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03))
          continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const double trkWeight = event_weight * trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
          }
        }
      } // end loop over tracks

    } // end loop over Pb+Pb tree
    cout << "Done primary Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("isEE",         &isEE);
    ppTree->SetBranchAddress ("vz",           &vz);
    ppTree->SetBranchAddress ("z_pt",         &z_pt);
    ppTree->SetBranchAddress ("z_y",          &z_y);
    ppTree->SetBranchAddress ("z_phi",        &z_phi);
    ppTree->SetBranchAddress ("z_m",          &z_m);
    ppTree->SetBranchAddress ("l1_pt",        &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",       &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",       &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",    &l1_charge);
    ppTree->SetBranchAddress ("l1_trk_pt",    &l1_trk_pt);
    ppTree->SetBranchAddress ("l1_trk_eta",   &l1_trk_eta);
    ppTree->SetBranchAddress ("l1_trk_phi",   &l1_trk_phi);
    ppTree->SetBranchAddress ("l2_pt",        &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",       &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",       &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",    &l2_charge);
    ppTree->SetBranchAddress ("l2_trk_pt",    &l2_trk_pt);
    ppTree->SetBranchAddress ("l2_trk_eta",   &l2_trk_eta);
    ppTree->SetBranchAddress ("l2_trk_phi",   &l2_trk_phi);
    ppTree->SetBranchAddress ("ntrk",         &ntrk);
    ppTree->SetBranchAddress ("trk_pt",       &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",      &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      if (fabs (vz) > 150)
        continue;

      if (event_weight == 0)
        continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03))
          continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const double trkWeight = event_weight * trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
          }
        }
      } // end loop over tracks

    } // end loop over pp tree
    cout << "Done primary pp loop." << endl;
  }

  SaveHists (outFileName);

  inFile->Close ();
  SaferDelete (inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Generates weights between a weighted sample and a matched sample.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: GenerateWeights (const char* weightedSampleInFilePattern, const char* matchedSampleInFilePattern, const char* outFileName) {

  const int nQ2Bins = 20;
  const double* q2Bins = linspace (0, 0.3, nQ2Bins);
  const int nPsi2Bins = 8;
  const double* psi2Bins = linspace (0, pi/2, nPsi2Bins);

  TH1D* h_fcal_et_dist[2][3][nPtZBins];
  TH1D* h_q2_dist[2][3][numFineCentBins][nPtZBins];
  TH1D* h_psi2_dist[2][3][numFineCentBins][nPtZBins];

  for (int iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      h_fcal_et_dist[0][iSpc][iPtZ] = new TH1D (Form ("h_weighted_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
      h_fcal_et_dist[0][iSpc][iPtZ]->Sumw2 ();
      h_fcal_et_dist[1][iSpc][iPtZ] = new TH1D (Form ("h_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
      h_fcal_et_dist[1][iSpc][iPtZ]->Sumw2 ();

      for (int iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        h_q2_dist[0][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_weighted_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_q2_dist[0][iSpc][iFineCent][iPtZ]->Sumw2 ();
        h_q2_dist[1][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_q2_dist[1][iSpc][iFineCent][iPtZ]->Sumw2 ();

        h_psi2_dist[0][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_weighted_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_psi2_dist[0][iSpc][iFineCent][iPtZ]->Sumw2 ();
        h_psi2_dist[1][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_psi2_dist[1][iSpc][iFineCent][iPtZ]->Sumw2 ();
      }
    }
  }


  bool isEE = false;
  float event_weight = 0, fcal_et = 0, z_pt = 0, z_phi = 0, q2 = 0, psi2 = 0;

  TChain* PbPbTree = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  PbPbTree->Add (weightedSampleInFilePattern);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  int nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;

    h_fcal_et_dist[0][iSpc][iPtZ]->Fill (fcal_et);
  }
  cout << "Done 1st Pb+Pb loop over weighted sample." << endl;

  PbPbTree->Reset ();
  PbPbTree->Add (matchedSampleInFilePattern);

  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;

    h_fcal_et_dist[1][iSpc][iPtZ]->Fill (fcal_et);
  }
  cout << "Done 1st Pb+Pb loop over matched events." << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Finalize FCal weighting histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////
  for (short iSample : {0, 1})
    for (short iSpc = 0; iSpc < 2; iSpc++)
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
        h_fcal_et_dist[iSample][2][iPtZ]->Add (h_fcal_et_dist[iSample][iSpc][iPtZ]);

  for (short iSample : {0, 1})
    for (short iSpc = 0; iSpc < 3; iSpc++)
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
        h_fcal_et_dist[iSample][iSpc][iPtZ]->Scale (1./h_fcal_et_dist[iSample][iSpc][iPtZ]->Integral ());

  for (short iSpc = 0; iSpc < 3; iSpc++)
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
      h_fcal_et_dist[1][iSpc][iPtZ]->Divide (h_fcal_et_dist[0][iSpc][iPtZ]);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree for additional q2, psi2 weights
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->Reset ();
  PbPbTree->Add (weightedSampleInFilePattern);

  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;
    const short iFineCent = GetFineCentBin (fcal_et);
    if (iFineCent < 1 || iFineCent > numFineCentBins-1)
      continue;

    event_weight = h_fcal_et_dist[1][iSpc][iPtZ]->GetBinContent (h_fcal_et_dist[1][iSpc][iPtZ]->FindFixBin (fcal_et));

    h_q2_dist[0][iSpc][iFineCent][iPtZ]->Fill (q2, event_weight);

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    h_psi2_dist[0][iSpc][iFineCent][iPtZ]->Fill (dphi);
  }
  cout << "Done 2nd Pb+Pb loop over weighted sample." << endl;

  PbPbTree->Reset ();
  PbPbTree->Add (matchedSampleInFilePattern);

  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;
    const short iFineCent = GetFineCentBin (fcal_et);
    if (iFineCent < 1 || iFineCent > numFineCentBins-1)
      continue;

    h_q2_dist[1][iSpc][iFineCent][iPtZ]->Fill (q2);

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    h_psi2_dist[1][iSpc][iFineCent][iPtZ]->Fill (dphi);
  }
  cout << "Done 2nd Pb+Pb loop over matched events." << endl;


  delete PbPbTree;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Finalize weighting histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////
  for (short iSample : {0, 1}) {
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
          h_q2_dist[iSample][2][iFineCent][iPtZ]->Add (h_q2_dist[iSample][iSpc][iFineCent][iPtZ]);
          h_psi2_dist[iSample][2][iFineCent][iPtZ]->Add (h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]);
        }
      }
    }
  }

  for (short iSample : {0, 1}) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
          h_q2_dist[iSample][iSpc][iFineCent][iPtZ]->Scale (1./h_q2_dist[iSample][iSpc][iFineCent][iPtZ]->Integral ());
          h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]->Scale (1./h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]->Integral ());
        }
      }
    }
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        h_q2_dist[1][iSpc][iFineCent][iPtZ]->Divide (h_q2_dist[0][iSpc][iFineCent][iPtZ]);
        h_psi2_dist[1][iSpc][iFineCent][iPtZ]->Divide (h_psi2_dist[0][iSpc][iFineCent][iPtZ]);
      }
    }
  }


  if (eventWeightsFile && eventWeightsFile->IsOpen ()) {
    eventWeightsFile->Close ();
    eventWeightsLoaded = false;
  }

  eventWeightsFile = new TFile (outFileName, "recreate");
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      SafeWrite (h_fcal_et_dist[1][iSpc][iPtZ]);
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        SafeWrite (h_q2_dist[1][iSpc][iFineCent][iPtZ]);
        SafeWrite (h_psi2_dist[1][iSpc][iFineCent][iPtZ]);
      }
    }
  }
  eventWeightsFile->Close ();
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load the event weights into memory
//////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadEventWeights () {
  if (eventWeightsLoaded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory;

  eventWeightsFile = new TFile (Form ("%s/%s", rootPath.Data (), eventWeightsFileName.c_str ()), "read");

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      h_PbPbFCal_weights[iSpc][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, name.c_str ()));
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        h_PbPbQ2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()));
        h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()));
      }
    }
  }

  eventWeightsLoaded = true;

  _gDirectory-> cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking efficiencies into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadTrackingEfficiencies (const bool doRebin) {
  if (effsLoaded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory;

  TString _effDir = "Nominal";
  if (useHITight)
    _effDir = "Variations/TrackHITightWPVariation";
  else if (doTrackEffVar)
    _effDir = "Variations/TrackEffPionsVariation";

  trkEffFile = new TFile (Form ("%s/TrackingEfficiencies/%s/trackingEfficiencies_%s.root", rootPath.Data (), _effDir.Data (), is2015Conds ? (useHijingEffs ? "Hijing_15":"15") : (useHijingEffs ? "Hijing_18":"18")), "read");

  for (int iCent = 0; iCent < numTrkCorrCentBins; iCent++) {
  //for (int iCent = 0; iCent < numFineCentBins; iCent++) {
    h2_num_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_matched_reco_tracks_iCent%i", iCent));
    h2_den_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_tracks_iCent%i", iCent));

    //if (iCent > 0) {
    //  h2_num_trk_effs[iCent]->RebinX (2);
    //  h2_num_trk_effs[iCent]->RebinY (2);
    //  h2_den_trk_effs[iCent]->RebinX (2);
    //  h2_den_trk_effs[iCent]->RebinY (2);
    //}

    //h2_trk_effs[iCent] = new TEfficiency (*(h2_num_trk_effs[iCent]), *(h2_den_trk_effs[iCent]));
    h2_trk_effs[iCent] = (TH2D*) h2_num_trk_effs[iCent]->Clone (Form ("h2_trk_eff_iCent%i", iCent));
    h2_trk_effs[iCent]->Divide (h2_den_trk_effs[iCent]);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TH1D* num = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta)))->GetCopyPassedHisto ();
      //TH1D* den = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta)))->GetCopyTotalHisto ();

      TH1D* num = (TH1D*) trkEffFile->Get (Form ("h_trk_eff_num_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkEffFile->Get (Form ("h_trk_eff_den_iCent%i_iEta%i", iCent, iEta));

      if (doRebin) {
        RebinSomeBins (num, maxNPtTrkBins, allPtTrkBins);
        RebinSomeBins (den, maxNPtTrkBins, allPtTrkBins);
      }

      //if (iCent > 0) {
      //  num->Rebin (2);
      //  den->Rebin (2);
      //}

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

  effsLoaded = true;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking efficiency for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double PhysicsAnalysis :: GetTrackingEfficiency (const float fcal_et, float trk_pt, const float trk_eta, const bool isPbPb) {
  if (!effsLoaded)
    LoadTrackingEfficiencies ();

  short iCent = 0;
  if (isPbPb) {
    while (iCent < numTrkCorrCentBins) {
      if (fcal_et < trkCorrCentBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }

  TH2D* t = h2_trk_effs[iCent];

  //const int bin = t->FindFixBin (trk_eta, trk_pt);
  //if (bin < 1 || t->GetNbinsX () * t->GetNbinsY () < bin)
  //  return 0;

  const int xbin = t->GetXaxis ()->FindFixBin (trk_eta);
  const int ybin = t->GetYaxis ()->FindFixBin (trk_pt);
  if (xbin < 1 || t->GetXaxis ()->GetNbins () < xbin)
    return 0;
  if (ybin < 1)
    return 0;
  else if (t->GetYaxis ()->GetNbins () < ybin)
    trk_pt = t->GetYaxis ()->GetBinCenter (t->GetYaxis ()->GetNbins ());

  return t->GetBinContent (t->FindFixBin (trk_eta, trk_pt)) + trkEffNSigma * t->GetBinError (t->FindFixBin (trk_eta, trk_pt));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking purities into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadTrackingPurities (const bool doRebin) {
  if (pursLoaded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory;

  TString _purDir = "Nominal";
  if (useHITight)
    _purDir = "Variations/TrackHITightWPVariation";

  trkPurFile = new TFile (Form ("%s/TrackingPurities/%s/trackingPurities_%s.root", rootPath.Data (), _purDir.Data (), is2015Conds ? "15" : "18"), "read");

  if (!trkPurFile || !trkPurFile->IsOpen ()) {
    cout << "Error in PhysicsAnalysis.h:: LoadTrackingPurities can not find file for " << name << endl;
    return;
  }

  for (int iCent = 0; iCent < numTrkCorrCentBins; iCent++) {
    h2_num_trk_purs[iCent] = (TH2D*) trkPurFile->Get (Form ("h2_primary_reco_tracks_iCent%i", iCent));
    h2_den_trk_purs[iCent] = (TH2D*) trkPurFile->Get (Form ("h2_reco_tracks_iCent%i", iCent));

    h2_num_trk_purs[iCent]->RebinX (4);
    h2_den_trk_purs[iCent]->RebinX (4);

    h2_trk_purs[iCent] = (TH2D*) h2_num_trk_purs[iCent]->Clone (Form ("h2_trk_pur_iCent%i", iCent));
    h2_trk_purs[iCent]->Divide (h2_den_trk_purs[iCent]);

    if (doTrackPurVar) {
      for (int ix = 1; ix <= h2_trk_purs[iCent]->GetNbinsX (); ix++) {
        for (int iy = 1; iy <= h2_trk_purs[iCent]->GetNbinsY (); iy++) {
          float fakeRate = 1. - h2_trk_purs[iCent]->GetBinContent (ix, iy);
          fakeRate = fakeRate + trkPurNSigma * 0.25 * fakeRate;
          h2_trk_purs[iCent]->SetBinContent (ix, iy, 1.-fakeRate);
        }
      }
    }

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* num = (TH1D*) trkPurFile->Get (Form ("h_primary_reco_tracks_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkPurFile->Get (Form ("h_reco_tracks_iCent%i_iEta%i", iCent, iEta));

      if (doRebin) {
        RebinSomeBins (num, maxNPtTrkBins, allPtTrkBins);
        RebinSomeBins (den, maxNPtTrkBins, allPtTrkBins);
      }

      h_trk_purs[iCent][iEta] = (TH1D*) num->Clone (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta));
      //h_trk_purs[iCent][iEta]->Divide (den);

      for (int ix = 1; ix <= h_trk_purs[iCent][iEta]->GetNbinsX (); ix++) {
        const float passes = num->GetBinContent (ix);
        const float trials = den->GetBinContent (ix);
        h_trk_purs[iCent][iEta]->SetBinContent (ix, passes/trials);
        h_trk_purs[iCent][iEta]->SetBinError (ix, sqrt ((passes/trials)*(1.-(passes/trials)) * (*den->GetSumw2 ())[ix] / pow (trials, 2)));
        //h_trk_purs[iCent][iEta]->SetBinError (ix, sqrt ((passes+1)*(passes+2) / ((trials+2)*(trials+3)) - pow (passes+1, 2) / pow (trials+2, 2)));
      }

      if (doTrackPurVar) {
        for (int ix = 1; ix <= h_trk_purs[iCent][iEta]->GetNbinsX (); ix++) {
          float fakeRate = 1. - h_trk_purs[iCent][iEta]->GetBinContent (ix);
          fakeRate = fakeRate + trkPurNSigma * 0.25 * fakeRate;
          h_trk_purs[iCent][iEta]->SetBinContent (ix, 1.-fakeRate);
        }
      }
    }
  }

  pursLoaded = true;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking purity factor for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double PhysicsAnalysis :: GetTrackingPurity (const float fcal_et, float trk_pt, const float trk_eta, const bool isPbPb) {
  if (!pursLoaded)
    LoadTrackingPurities ();

  short iCent = 0;
  if (isPbPb) {
    while (iCent < numTrkCorrCentBins) {
      if (fcal_et < trkCorrCentBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }

  TH2D* t = h2_trk_purs[iCent];

  //const int bin = t->FindFixBin (trk_eta, trk_pt);
  //if (bin < 1 || t->GetNbinsX () * t->GetNbinsY () < bin)
  //  return 0;

  const int xbin = t->GetXaxis ()->FindFixBin (trk_eta);
  const int ybin = t->GetYaxis ()->FindFixBin (trk_pt);
  if (xbin < 1 || t->GetXaxis ()->GetNbins () < xbin)
    return 0;
  if (ybin < 1)
    return 0;
  else if (t->GetYaxis ()->GetNbins () < ybin)
    trk_pt = t->GetYaxis ()->GetBinCenter (t->GetYaxis ()->GetNbins ());

  const double pur = t->GetBinContent (t->FindFixBin (trk_eta, trk_pt));

  return pur;
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

      //const char* canvasName = Form ("c_z_trk_phi_pttrk_iPtZ%i_%s%s", iPtZ, spc, _subBkg ? "_subBkg":"");
      const char* canvasName = Form ("c_z_trk_phi_pttrk_iPtZ%i_%s", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 900, 450*numCentBins);
        gDirectory->Add (c);
        c->cd ();
        c->Divide (2, numCentBins);
      }

      //for (short iCent = 0; iCent < 1; iCent++) {
      //  c->cd ();
      //  gPad->SetLogy ();
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (bool subBkg : {false, true}) {
          if (subBkg && !_subBkg)
            continue;

          c->cd ((2*iCent)+1 + (int)subBkg);
          GetDrawnObjects ();
          //gPad->SetLogy ();

          gPad->SetTopMargin (0.01);
          gPad->SetBottomMargin (0.12);
          gPad->SetRightMargin (0.01);
          gPad->SetLeftMargin (0.12);

          double min = -1.5, max = 3.5;
          if (!subBkg) {
            GetMinAndMax (min, max, true);
            for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
              TH1D* h = (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]);
              if (h->GetMinimum () < min) min = h->GetMinimum ();
              if (h->GetMaximum () > max) max = h->GetMaximum ();
            } // end loop over iPtTrk
            //min *= 0.5;
            if (max != 3.5)
              max = max <= 0 ? 1 : 1.2*max;
            SetMinAndMax (min, max);
          }

          if (plotFill) {
            for (short iPtTrk = 0; iPtTrk < std::min (3, nPtTrkBins[iPtZ]); iPtTrk++) {
              TH1D* h = (TH1D*) (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent])->Clone ();

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
              h->GetYaxis ()->SetTitleSize (0.07);
              h->GetXaxis ()->SetLabelSize (0.06);
              h->GetYaxis ()->SetLabelSize (0.06);

              h->DrawCopy ((subBkg || !canvasExists) && iPtTrk == 0 ? "hist" : "same hist");
              //h->SetLineWidth (1);
              //h->Draw ("hist same");

              //TLine* line = new TLine (h->GetBinLowEdge (1), 0, h->GetBinLowEdge (h->GetNbinsX ()), 0);
              //line->Draw ("same");

              if (_subBkg) LabelCorrelations (iPtZ, iPtTrk, iCent, subBkg);
            } // end loop over iPtTrk
            gPad->RedrawAxis ();
          } else {
            for (short iPtTrk = 0; iPtTrk < std::min (3, nPtTrkBins[iPtZ]); iPtTrk++) {
              TGAE* g = GetTGAE (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent]);
              ResetXErrors (g);

              //const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
              const Style_t markerStyle = markerStyles[iPtTrk];
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
              g->GetYaxis ()->SetTitleSize (0.07);
              g->GetXaxis ()->SetLabelSize (0.06);
              g->GetYaxis ()->SetLabelSize (0.06);

              g->Draw ((subBkg || !canvasExists) && iPtTrk == 0 ? "AP" : "P");

              if (_subBkg) LabelCorrelations (iPtZ, iPtTrk, iCent, subBkg);

              TLine* line1 = new TLine (phiLowBins[1], min, phiLowBins[1], max);
              TLine* line2 = new TLine (2*pi - phiLowBins[1], min, 2*pi - phiLowBins[1], max);
              //TLine* line3 = new TLine (phiLowBins[2], min, phiLowBins[2], max);
              //TLine* line4 = new TLine (2*pi - phiLowBins[2], min, 2*pi - phiLowBins[2], max);

              line1->SetLineStyle (2);
              line2->SetLineStyle (2);
              //line3->SetLineStyle (2);
              //line4->SetLineStyle (2);

              line1->SetLineWidth (2);
              line2->SetLineWidth (2);
              //line3->SetLineWidth (2);
              //line4->SetLineWidth (2);

              line1->SetLineColor (kBlack);
              line2->SetLineColor (kBlack);
              //line3->SetLineColor (kGray);
              //line4->SetLineColor (kGray);

              line1->Draw ("same");
              line2->Draw ("same");
              //line3->Draw ("same");
              //line4->Draw ("same");
            } // end loop over iPtTrk
          }
        }
      } // end loop over centrality

      //c->SaveAs (Form ("%s/dPhi_Distributions/dPhi_pTtrk_iPtZ%i_%s%s.pdf", plotPath.Data (), iPtZ, spc, _subBkg ? "_subBkg":""));
      c->SaveAs (Form ("%s/dPhi_Distributions/dPhi_pTtrk_iPtZ%i_%s.pdf", plotPath.Data (), iPtZ, spc));
    } // end loop over pT^Z
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track dPhi distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelCorrelations (const short iPtZ, const short iPtTrk, const short iCent, const bool subBkg) {
  if (iCent == 0 && iPtTrk == 0) {
    if (!subBkg) {
      //myText (0.2, 0.91, kBlack, "#bf{#it{ATLAS}} Simulation", 0.07);
      myText (0.2, 0.91, kBlack, "#bf{#it{ATLAS}} Internal", 0.077);
      myText (0.2, 0.82, kBlack, "#it{pp}, 5.02 TeV", 0.07);
      //myText (0.16, 0.93, kBlack, "Data", 0.04);
      //myText (0.30, 0.93, kBlack, "Minbias", 0.05);
      //TVirtualPad* cPad = gPad; // store current pad
      //TBox* b = nullptr;
      if (iPtZ == nPtZBins-1)
        myText (0.2, 0.73, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.07);
      else
        myText (0.2, 0.73, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.07);
      myMarkerTextNoLine (0.24, 0.64, kBlack, kFullCircle, "Z-tagged Data", 1.25, 0.07);
      myOnlyBoxText (0.24, 0.55, 1.2, fillColors[0], kBlack, 1, "Minimum Bias", 0.07, 1001, 1);
    }
    else {
      myText (0.2, 0.91, kBlack, "After UE subtraction", 0.07);
    }
  }
  if (iCent == 0 && subBkg) {
    const float pt_lo = ptTrkBins[iPtZ][iPtTrk];
    const float pt_hi = ptTrkBins[iPtZ][iPtTrk+1];
    //if (iPtTrk == 0)
    //  myText (0.3, 0.93, kBlack, "#it{p}_{T}^{ ch} [GeV]", 0.04);
    myMarkerTextNoLine (0.23, 0.85-0.065*iPtTrk, colors[iPtTrk], markerStyles[iPtTrk], "", 1.4, 0.06);
    myText             (0.24, 0.83-0.065*iPtTrk,  kBlack,                               Form ("%.0f < #it{p}_{T} < %.0f GeV", pt_lo, pt_hi), 0.06);
    //myText (0.3, 0.89-0.04*iPtTrk, colors[iPtTrk], Form ("(%.1f, %.1f)", pt_lo, pt_hi), 0.04);

    //b = TBoxNDC (0.33-0.024, 0.87-0.065*iPtTrk-0.016, 0.33+0.024, 0.87-0.065*iPtTrk+0.016);
    //b->SetFillColorAlpha (fillColors[iPtTrk], fillAlpha);
    //myMarkerTextNoLine (0.23, 0.89-0.04*iPtTrk, colors[iPtTrk], kFullCircle, "", 1.25, 0.04);

    //b->Draw ("l");
    //cPad->cd ();
  }
  else if (!subBkg && iPtTrk == 0 && iCent != 0) {
    myText (0.2, 0.90, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.07);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies (PhysicsAnalysis* a) {
  if (!effsLoaded)
    LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 800);
    //c = new TCanvas (canvasName, "", 1800, 400);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 2);
    //c->Divide (4, 1);
  }
  c->cd ();

  //TGraphErrors* g0[numCentBins];
  //TGraphErrors* g1[numCentBins];
  //TGraphErrors* g2[numCentBins];
  //TGraphErrors* g3[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
//    if (iCent > 0) continue;
    gPad->SetLogx ();

    //g0[iCent] = new TGraphErrors (numEtaTrkBins);
    //g1[iCent] = new TGraphErrors (numEtaTrkBins);
    //g2[iCent] = new TGraphErrors (numEtaTrkBins);
    //g3[iCent] = new TGraphErrors (numEtaTrkBins);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TEfficiency* eff = h_trk_effs[iCent][iEta];
      TGAE* eff = GetTGAE (h_trk_effs[iCent][iEta]);

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);
      eff->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      eff->SetMarkerSize (0.5);

      eff->SetTitle (";#it{p}_{T} [GeV];Weighted Reco. Eff.");
      eff->GetXaxis ()->SetRangeUser (0.5, 60);
      eff->GetYaxis ()->SetRangeUser (0.3, 1.12);

      eff->GetXaxis ()->SetTitleSize (0.07);
      eff->GetYaxis ()->SetTitleSize (0.07);
      eff->GetXaxis ()->SetTitleOffset (0.7 * eff->GetXaxis ()->GetTitleOffset ());
      eff->GetYaxis ()->SetTitleOffset (0.7 * eff->GetYaxis ()->GetTitleOffset ());

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (!canvasExists && iEta == 0 ? "AP" : "P");
      //eff->Draw (iEta == 0 ? "APL" : "LP same");

      //gPad->Update ();

      //eff->GetPaintedGraph ()->GetXaxis ()->SetRangeUser (0.5, 65);
      //eff->GetPaintedGraph ()->GetYaxis ()->SetRangeUser (0.3, 1.08);

      ////TH1D* confInts = (TH1D*) h_trk_effs[iCent][iEta]->Clone ("confInts");
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2 + [3]*(log(x))^3", 1, 65);
      ////TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2", 2, 65);
      //fit->SetParameter (0, 1);
      //fit->SetParameter (1, 0);
      //fit->SetParameter (2, 0);
      ////fit->SetParameter (3, 0);
      //h_trk_effs[iCent][iEta]->Fit (fit, "RN0Q");

      //g0[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (0));
      //g1[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (1));
      //g2[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (2));
      //g3[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (3));

      //g0[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (0));
      //g1[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (1));
      //g2[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (2));
      //g3[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (3));

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
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");

    if (!a)
      continue;
    else
      a->LoadTrackingEfficiencies ();

    gPad->SetBottomMargin (0);
    c->cd (iCent+5);
    gPad->SetLogx ();
    gPad->SetTopMargin (0);
    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* h = (TH1D*) a->h_trk_effs[iCent][iEta]->Clone ("temp");
      h->Divide (h_trk_effs[iCent][iEta]);

      TGAE* eff = GetTGAE (h);
      delete h;

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);
      eff->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      eff->SetMarkerSize (0.5);

      //eff->SetTitle (";#it{p}_{T} [GeV];Pions / Inclusive hadrons");
      eff->SetTitle (";#it{p}_{T} [GeV];HITight / HILoose");
      eff->GetXaxis ()->SetRangeUser (0.5, 60);
      //eff->GetYaxis ()->SetRangeUser (0.89, 1.11);
      eff->GetYaxis ()->SetRangeUser (0.7, 1.11);

      eff->GetXaxis ()->SetTitleSize (0.07);
      eff->GetYaxis ()->SetTitleSize (0.07);
      eff->GetXaxis ()->SetTitleOffset (0.7 * eff->GetXaxis ()->GetTitleOffset ());
      eff->GetYaxis ()->SetTitleOffset (0.7 * eff->GetYaxis ()->GetTitleOffset ());

      eff->GetYaxis ()->CenterTitle ();

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (!canvasExists && iEta == 0 ? "AP" : "P");

      //TF1* fit = new TF1 ("fit", "[0]+[1]*log(x)+[2]*(log(x))^2", 0.7, 15);//+[3]*(log(x))^3+[4]*(log(x))^4", 0.7, 15);
      //fit->SetParameter (0, 1);
      //fit->SetParameter (1, 0);
      //fit->SetParameter (2, 0);
      ////fit->SetParameter (3, 0);
      ////fit->SetParameter (4, 0);
      //eff->Fit (fit, "RN0Q");
      //fit->SetLineColor (colors[iEta]);
      //fit->SetLineStyle (2);
      //fit->Draw ("same");
    }

    l->Draw ("same");
  }

  if (a) {
    bool temp = a->useAltMarker;
    a->useAltMarker = true;
    a->PlotTrackingEfficiencies ();
    a->useAltMarker = temp;
  }
  else
    c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficiencies.pdf", plotPath.Data ()));

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
  if (!effsLoaded)
    LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    //if (iCent > 0) continue;
    gPad->SetLogy ();

    //TEfficiency* eff = h2_trk_effs[iCent];
    TH2D* eff = h2_trk_effs[iCent];
    eff->GetZaxis ()->SetRangeUser (0.3, 1.00);
    eff->GetZaxis ()->SetTitle ("Weighted Reco. Eff.");
    eff->GetZaxis ()->SetTitleOffset (1.4 * eff->GetZaxis ()->GetTitleOffset ());

    eff->GetYaxis ()->SetMoreLogLabels ();

    eff->Draw ("colz");
    //eff->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();

    if (iCent == 0)
      myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
    else
      myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficiencies2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track efficiency plots
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrackingEfficiencies (const short iCent, const short iEta) {
  if (iEta == 0) {
    if (iCent == 0)
      myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
    else
      myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
  }

  if (iCent == 0) {
  //  myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
    myMarkerTextNoLine (0.5, 0.34-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
  //  myMarkerTextNoLine (0.5, 0.50-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
  }
  else if (iCent == 1 && iEta == 0) {
    myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "HILoose tracks", 1.2, 0.06);
    myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "HITight tracks", 1.2, 0.06);
    //myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "Inclusive hadrons", 1.2, 0.06);
    //myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "Pions only", 1.2, 0.06);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities (PhysicsAnalysis* a) {
  if (!pursLoaded)
    LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 800);
    //c = new TCanvas (canvasName, "", 1800, 400);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 2);
    c->Draw ();
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    gPad->SetLogx ();

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TGAE* pur = GetTGAE (h_trk_purs[iCent][iEta]);

      pur->SetLineColor (colors[iEta]);
      pur->SetMarkerColor (colors[iEta]);
      pur->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      pur->SetMarkerSize (0.5);

      pur->SetTitle (";#it{p}_{T} [GeV];Primary Track Fraction");
      pur->GetXaxis ()->SetRangeUser (0.5, 60);
      pur->GetYaxis ()->SetRangeUser (0.94, 1.010);

      pur->GetXaxis ()->SetTitleSize (0.07);
      pur->GetYaxis ()->SetTitleSize (0.07);
      pur->GetXaxis ()->SetTitleOffset (0.7 * pur->GetXaxis ()->GetTitleOffset ());
      pur->GetYaxis ()->SetTitleOffset (0.7 * pur->GetYaxis ()->GetTitleOffset ());

      pur->GetXaxis ()->SetMoreLogLabels ();

      pur->Draw (!canvasExists && iEta == 0 ? "AP" : "P");

      //TH1D* confInts = (TH1D*) h_trk_purs[iCent][iEta]->Clone ("confInts");
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2 + [3]*(log(x))^3 + [4]*(log(x))^4", 1, 65);
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2", 2, 65);
      //fit->SetParameter (0, 1);
      //fit->SetParameter (1, 0);
      //fit->SetParameter (2, 0);
      //fit->SetParameter (3, 0);
      //fit->FixParameter (4, 0);
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
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");

    if (!a)
      continue;
    else
      a->LoadTrackingPurities ();

    gPad->SetBottomMargin (0);
    c->cd (iCent+5);
    gPad->SetLogx ();
    gPad->SetTopMargin (0);
    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* h = (TH1D*) a->h_trk_purs[iCent][iEta]->Clone ("temp");
      h->Divide (h_trk_purs[iCent][iEta]);

      TGAE* pur = GetTGAE (h);
      delete h;

      pur->SetLineColor (colors[iEta]);
      pur->SetMarkerColor (colors[iEta]);
      pur->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      pur->SetMarkerSize (0.5);

      pur->SetTitle (";#it{p}_{T} [GeV];HITight / HILoose");
      pur->GetXaxis ()->SetRangeUser (0.5, 60);
      pur->GetYaxis ()->SetRangeUser (0.94, 1.06);

      pur->GetXaxis ()->SetTitleSize (0.07);
      pur->GetYaxis ()->SetTitleSize (0.07);
      pur->GetXaxis ()->SetTitleOffset (0.7 * pur->GetXaxis ()->GetTitleOffset ());
      pur->GetYaxis ()->SetTitleOffset (0.7 * pur->GetYaxis ()->GetTitleOffset ());

      pur->GetYaxis ()->CenterTitle ();

      pur->GetXaxis ()->SetMoreLogLabels ();

      pur->Draw (!canvasExists && iEta == 0 ? "AP" : "P");
    }

    l->Draw ("same");
  }

  if (a) {
    bool temp = a->useAltMarker;
    a->useAltMarker = true;
    a->PlotTrackingPurities ();
    a->useAltMarker = temp;
  }
  else
    c->SaveAs (Form ("%s/TrackingPurities/TrackingPurities.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities2D () {
  if (!pursLoaded)
    LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    //if (iCent > 0) continue;
    gPad->SetLogy ();

    //TEfficiency* pur = h2_trk_purs[iCent];
    TH2D* pur = h2_trk_purs[iCent];
    pur->GetZaxis ()->SetRangeUser (0.96, 1.00);
    pur->GetZaxis ()->SetTitle ("Primary Track Fraction");
    pur->GetZaxis ()->SetTitleOffset (1.4 * pur->GetZaxis ()->GetTitleOffset ());

    pur->GetYaxis ()->SetMoreLogLabels ();

    pur->Draw ("colz");
    //pur->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();

    if (iCent == 0)
      myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
    else
      myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
  }

  c->SaveAs (Form ("%s/TrackingPurities/TrackingPurities2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track purity plots
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrackingPurities (const short iCent, const short iEta) {
  if (iEta == 0) {
    if (iCent == 0)
      myText (0.22, 0.88, kBlack, "#it{pp}", 0.08);
    else
      myText (0.22, 0.88, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
  }

  if (iCent == 0) {
  //  myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
    myMarkerTextNoLine (0.5, 0.34-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
  //  myMarkerTextNoLine (0.5, 0.50-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
  }
  else if (iCent == 1 && iEta == 0) {
    myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "HILoose tracks", 1.2, 0.06);
    myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "HITight tracks", 1.2, 0.06);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Subtracts mixed event background from track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SubtractBackground (PhysicsAnalysis* a) {
  if (backgroundSubtracted) {
    cout << "Background already subtracted for " << name << endl;
    return;
  }

  //cout << "Subtracting bkg. " << a->Name () << " from " << Name () << endl;

  //**** Create empty subtracted histograms for combined channel yield measurement ****//
  if (!backgroundSubtracted) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
        h_trk_pt_ptz[2][iPtZ][iCent]->Reset ();
        h_trk_pt_ptz_sub[2][iPtZ][iCent]        = new TH1D (Form ("h_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
        h_trk_pt_ptz_sub[2][iPtZ][iCent]->Sumw2 ();
        h_trk_pt_ptz_sig_to_bkg[2][iPtZ][iCent] = new TH1D (Form ("h_trk_pt_ptz_sigToBkg_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
        h_trk_pt_ptz_sig_to_bkg[2][iPtZ][iCent]->Sumw2 ();
        h_trk_xhz_ptz[2][iPtZ][iCent]->Reset ();
        h_trk_xhz_ptz_sub[2][iPtZ][iCent]         = new TH1D (Form ("h_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
        h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Sumw2 ();
        h_trk_xhz_ptz_sig_to_bkg[2][iPtZ][iCent]  = new TH1D (Form ("h_trk_xhz_ptz_sigToBkg_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
        h_trk_xhz_ptz_sig_to_bkg[2][iPtZ][iCent]->Sumw2 ();
        for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi[2][iPtZ][iPhi][iCent]->Reset ();
          h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]         = new TH1D (Form ("h_trk_pt_dphi_sub_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
          h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_pt_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]  = new TH1D (Form ("h_trk_pt_dphi_sigToBkg_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
          h_trk_pt_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_xhz_dphi[2][iPtZ][iPhi][iCent]->Reset ();
          h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]        = new TH1D (Form ("h_trk_xhz_dphi_sub_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
          h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_xhz_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_xhz_dphi_sigToBkg_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
          h_trk_xhz_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Sumw2 ();
        }
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          h_trk_dphi[2][iPtZ][iPtTrk][iCent]->Reset ();
          h_trk_dphi_sub[2][iPtZ][iPtTrk][iCent] = new TH1D (Form ("h_trk_dphi_sub_comb_iPtZ%i_iPtTrk%i_iCent%i_%s", iPtZ, iPtTrk, iCent, name.c_str ()), "", 80, -pi/2, 3*pi/2);
          h_trk_dphi_sub[2][iPtZ][iPtTrk][iCent]->Sumw2 ();
        }
      } // end loop over iPtZ
    } // end loop over iCent
  }
  
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        //******** Do subtraction of integrated dPhi plot ********//
        TH1D* h = (TH1D*) h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        TH1D* sub = nullptr;
        if (a != nullptr) {
          sub = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
          h->Add (sub, -1);
        }
        h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] = h;

        h = (TH1D*) h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        if (a != nullptr && sub != nullptr) {
          h->Add (sub, -1);
          h->Divide (sub);
        }
        else
          h->Divide (h);
        h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = h;

        h = (TH1D*) h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        if (a != nullptr) {
          sub = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
          h->Add (sub, -1);
        }
        h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent] = h;

        h = (TH1D*) h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        if (a != nullptr && sub != nullptr) {
          h->Add (sub, -1);
          h->Divide (sub);
        }
        else
          h->Divide (h);
        h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = h;


        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {

          //******** Do subtraction of pT ********//
          h = (TH1D*) h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          if (a != nullptr) {
            sub = a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
            h->Add (sub, -1);
          }
          h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] = h;

          h = (TH1D*) h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          if (a != nullptr && sub != nullptr) {
            h->Add (sub, -1);
            h->Divide (sub);
          }
          h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = h;


          //******** Do subtraction of z_h ********//
          h = new TH1D (Form ("h_trk_xhz_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
          h->Sumw2 ();
          h->Add (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (a != nullptr) {
            sub = a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
            h->Add (sub, -1);
          }
          h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent] = h;

          h = new TH1D (Form ("h_trk_xhz_dphi_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
          h->Sumw2 ();
          h->Add (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (a != nullptr && sub != nullptr) {
            h->Add (sub, -1);
            h->Divide (sub);
          }
          h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = h;
        } // end loop over phi


        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          //******** Do background subtraction of phi distributions ********//
          TH1D* h = (TH1D*) h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
          if (a != nullptr) {
            TH1D* sub = a->h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent];
            while (sub->GetNbinsX () > h->GetNbinsX ())
              sub->Rebin (2);
            h->Add (sub, -1);
          }
          h_trk_dphi_sub[iSpc][iPtZ][iPtTrk][iCent] = h;

        } // end loop over pT^trk
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  UnfoldSubtractedYield ();

  PhysicsAnalysis :: CombineHists ();

  backgroundSubtracted = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Applies bin migration factors to subtracted yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: UnfoldSubtractedYield () {

  if (histsUnfolded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TFile* f_binMigrationFile = new TFile (Form ("%s/BinMigrationFactors/binmigration_corrfactors_master.root", rootPath.Data ()), "read");

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        const short iBBBCent = (iCent == 0 ? 0 : GetBBBCorrCentBin (0.5 * (centBins[iCent-1] + centBins[iCent])));

        const char* cent = (iBBBCent == 0 ? "pp" : Form ("CENT%i", numBBBCorrCentBins-iBBBCent-1));

        TF1* f = (TF1*) f_binMigrationFile->Get (Form ("tf1_%s_pt_ZPT%i_%s", spc, iPtZ-2, cent));
        TH1D* h = h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];

        for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
          const float x = h->GetBinCenter (ix);
          float y = h->GetBinContent (ix);
          float yerr = h->GetBinError (ix);
          y = y / f->Eval (x);
          yerr = yerr / f->Eval (x);
          h->SetBinContent (ix, y);
          h->SetBinError (ix, yerr);
        }

        f = f_z_trk_zxzh_binMigration[iSpc][iPtZ][iCent] = (TF1*) f_binMigrationFile->Get (Form ("tf1_%s_xh_ZPT%i_%s", spc, iPtZ-2, cent));
        h = h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
        for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
          const float x = h->GetBinCenter (ix);
          float y = h->GetBinContent (ix);
          float yerr = h->GetBinError (ix);
          y = y / f->Eval (x);
          yerr = yerr / f->Eval (x);
          h->SetBinContent (ix, y);
          h->SetBinError (ix, yerr);
        }
      }
    }
  }

  f_binMigrationFile->Close ();
  histsUnfolded = true;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Inflates statistical uncertainty in final results histograms by some amount
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: InflateStatUnc (const float amount) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (TH1D* h : {h_trk_pt_ptz[iSpc][iPtZ][iCent], h_trk_xhz_ptz[iSpc][iPtZ][iCent], h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent], h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]}) {
          if (!h)
            continue;
          for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
            h->SetBinError (ix, h->GetBinError (ix) * (1+amount));
          }
        }
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species
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

          TH1D* h = h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
          const float bkgCountsOverObsCounts = (a->h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1)) / (h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1));
          if (iCent == 0 && iPhi == 0 && iSpc == 2 && iPtZ == 2) {
            cout << bkgCountsOverObsCounts << endl;
            cout << "2nd to last bin before: " << h->GetBinContent (6) << endl;
          }

          h->Add (a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent], -bkgCountsOverObsCounts);
          h->Scale (1./ (1-bkgCountsOverObsCounts));
          //if (iCent == 0 && iPhi == 0 && iSpc == 2 && iPtZ == 2) {
          //  cout << "2nd to last bin after:  " << h->GetBinContent (6) << endl;
          //}

          h = h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          h->Add (a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], -bkgCountsOverObsCounts);
          h->Scale (1./ (1-bkgCountsOverObsCounts));
        } // end loop over phi
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  sameSignsSubtracted = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Applies a constant (x-axis independent) systematic to the analysis, e.g. a luminosity unc.
// Only applies to pT dependent histograms, i.e. will need to subtract background still.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ApplyRelativeVariation (float**** relVar, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      // Hadron yield systematics, signal & signal+bkg levels
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1./relVar[iSpc][iPtZ][iPhi][iCent]);
        //  h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1./relVar[iSpc][iPtZ][iPhi][iCent]);
        //} // end loop over phi

        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1./relVar[iSpc][iPtZ][numPhiBins][iCent]);
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1./relVar[iSpc][iPtZ][numPhiBins][iCent]);
      } // end loop over cents
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Converts analysis to a systematic variation by adding or subtracting statistical errors
// Only converts relevant histograms, e.g. for minbias will only do this to track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ConvertToStatVariation (const bool upVar, const float nSigma) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      // Hadron yield systematics, signal & signal+bkg levels
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          AddStatVar (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);

          AddStatVar (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);

          AddStatVar (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
          AddStatVar (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        } // end loop over phi

        AddStatVar (h_trk_pt_ptz[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent], upVar, nSigma);
      } // end loop over cents
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

      const char* canvasName = Form ("c_RawTrkYield_%s_%s_iPtZ%i", useTrkPt ? "pttrk" : "zh", spc, iPtZ);
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
          TH1D* h = useTrkPt ? h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over phi
        min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
        max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = 0; iPhi < 1; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = useTrkPt ? h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            useTrkPt ? h->GetYaxis ()->SetTitle ("Y (#it{p}_{T}, #Delta#phi)") : h->GetYaxis ()->SetTitle ("Y (#it{x}_{hZ}, #Delta#phi)");

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
            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
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

            useTrkPt ? g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            useTrkPt ? g->GetYaxis ()->SetTitle ("Y (#it{p}_{T}, #Delta#phi)") : g->GetYaxis ()->SetTitle ("Y (#it{x}_{hZ}, #Delta#phi)");

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

      const char* canvasName = Form ("c_TrkYield_%s_%s_iPtZ%i", useTrkPt ? "pttrk" : "zh", spc, iPtZ);
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
          TH1D* h = useTrkPt ? h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over phi
        min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
        max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = useTrkPt ? h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];

            //h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            //h->SetMarkerSize (0);
            h->SetLineColor (colors[iPhi]);
            h->SetLineStyle (iPhi+1);
            h->SetLineWidth (1);

            useTrkPt ? h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            useTrkPt ? h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x}_{hZ} d#Delta#phi");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

            h->Draw (plotNewAxes && iPhi == 1 ? "hist ][" : "hist ][ same");
            //h->DrawCopy (plotNewAxes && iPhi == 1 ? "bar" : "bar same");
            //h->SetLineWidth (1);
            //h->Draw ("hist same");
          }
          gPad->RedrawAxis ();
        } else {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
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

            useTrkPt ? g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            useTrkPt ? g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x}_{hZ}d#Delta#phi");

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
            LabelTrkYield (iCent, iPhi, iPtZ, iSpc);

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
          TH1D* h = useTrkPt ? h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over phi
        float delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* h = useTrkPt ? h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            useTrkPt ? h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            useTrkPt ? h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x}_{hZ}d#Delta#phi");

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

            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]);
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

            useTrkPt ? g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            useTrkPt ? g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x}_{hZ}d#Delta#phi");

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

        if (plotAsSystematic)
          continue;

        bottomPad->cd ();
        GetDrawnObjects ();
        plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = (useTrkPt ? h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
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
            TH1D* h = useTrkPt ? h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            useTrkPt ? h->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            h->GetYaxis ()->SetTitle ("Y / Y_{bkg}");
            h->GetYaxis ()->CenterTitle ();

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
            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
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

            useTrkPt ? g->GetXaxis ()->SetLimits (ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
            g->GetYaxis ()->SetTitle ("Y / Y_{bkg}");
            g->GetYaxis ()->CenterTitle ();

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
void PhysicsAnalysis :: LabelTrkYield (const short iCent, const short iPhi, const short iPtZ, const short iSpc) {
  //const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);

  if (iCent == 0)
    myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
  else
    myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

  if (iCent == 2 && iSpc == 0)
    myText (0.72, 0.88, kBlack, "Z#rightarrow ee", 0.06);
  else if (iCent == 2 && iSpc == 1)
    myText (0.72, 0.88, kBlack, "Z#rightarrow#mu#mu", 0.06);

  if (iCent == 0)
    myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
  else if (iCent == 1) {
    if (iPtZ == nPtZBins-1) {
      myText (0.485, 0.88, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.06);
    }
    else {
      myText (0.485, 0.88, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.06);
    }
  }
  else if (iCent == numCentBins-1) {
    myText (0.43, 0.90, kBlack, "MB", 0.06);
    myText (0.52, 0.90, kBlack, "Z-tag", 0.06);
    myLineText (0.52, 0.85-0.06*(iPhi-1), colors[iPhi], iPhi+1, "", 2.0, 0.054) ;
    myMarkerTextNoLine (0.59, 0.85-0.06*(iPhi-1), colors[iPhi], kFullCircle, "", 1.5, 0.054); // for plotting data vs bkg.

    const char* lo = GetPiString (phiLowBins[iPhi]);
    const char* hi = GetPiString (phiHighBins[iPhi]);

    //if (iPtZ == 2) {
    //  myText (0.34, 0.91, kBlack, "Reco.", 0.06);
    //  myText (0.47, 0.91, kBlack, "Truth", 0.06);
    //}
    //myMarkerTextNoLine (0.42, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kFullCircle, "", 1.5, 0.054); // for plotting MC reco vs truth
    //myMarkerTextNoLine (0.52, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kOpenCircle, "", 1.5, 0.054);
    myText (0.61, 0.84-0.06*(iPhi-1), kBlack, Form ("(%s, %s)", lo, hi), 0.054);

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

    const char* canvasName = Form ("c_TrkYield_zpt_%s", useTrkPt ? "pttrk" : "xhz");
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

      const char* topPadName = Form ("p_TrkYield_pt_top_iCent%i", iCent);
      const char* middlePadName = Form ("p_TrkYield_pt_middle_iCent%i", iCent);
      const char* bottomPadName = Form ("p_TrkYield_pt_bottom_iCent%i", iCent);

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
      if (useTrkPt) {
        if (iCent == 0) { min = 1e-7; max = 2e1; }
        else if (iCent == 1) { min = 1e-6; max = 2e2; }
        else if (iCent == 2) { min = 5e-5; max = 1e3; }
        else if (iCent == 3) { min = 1e-4; max = 2e3; }
      }
      else {
        if (iCent == 0) { min = 2e-5; max = 2e3; }
        else if (iCent == 1) { min = 1e-3; max = 2e4; }
        else if (iCent == 2) { min = 1e-3; max = 5e4; }
        else if (iCent == 3) { min = 1e-3; max = 4e5; }
      }
      //GetMinAndMax (min, max, true);
      //for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
      //  TH1D* h = (useTrkPt ? h_trk_pt_ptz[iSpc][iPtZ][iCent] : h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
      //  min = fmin (min, h->GetMinimum (0));
      //  max = fmax (max, h->GetMaximum ());
      //} // end loop over phi
      //min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
      //max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
      //SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
          TH1D* h = (useTrkPt ? h_trk_pt_ptz[iSpc][iPtZ][iCent] :  h_trk_xhz_ptz[iSpc][iPtZ][iCent]);

          //h->SetFillColorAlpha (fillColors[iPtZ-2], fillAlpha);
          //h->SetMarkerSize (0);
          h->SetLineColor (colors[iPtZ-2]);
          h->SetLineStyle (iPtZ);
          h->SetLineWidth (1);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
          useTrkPt ? h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x} d#Delta#phi");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

          h->Draw (plotNewAxes && iPtZ == nPtZBins-1 ? "hist ][" : "hist ][ same");
          //h->DrawCopy (plotNewAxes && iPtZ == nPtZBins-1 ? "bar" : "bar same");
          //h->SetLineWidth (1);
          //h->Draw ("hist same");
        }
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = markerStyles[iPtZ-3];
          TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz[iSpc][iPtZ][iCent] : h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerColor (colors[iPtZ-2]);
            g->SetLineColor (colors[iPtZ-2]);
            g->SetMarkerSize (1);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ-2]);
            g->SetFillColorAlpha (fillColors[iPtZ-2], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
          g->GetYaxis ()->SetRangeUser (min, max);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
          useTrkPt ? g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x} d#Delta#phi");

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
            string drawString = string (!canvasExists && iPtZ-3 == 0 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ-3 == 0 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }
        } // end loop over phi
      }

      if (!canvasExists)
        for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
          LabelTrkYieldZPt (iCent, iPtZ, iSpc);
        }

      if (!plotSignal)
        continue;

      middlePad->cd ();
      GetDrawnObjects ();
      plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
      gPad->SetLogx ();
      gPad->SetLogy ();

      if (useTrkPt) {
        min = 3e-3;
        max = 20;
      }
      else {
        min = 1e-2;
        max = 8e2;
      }
      //min = 1e30, max = 0;
      //GetMinAndMax (min, max, true);
      //for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
      //  TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]);
      //  min = fmin (min, h->GetMinimum (0));
      //  max = fmax (max, h->GetMaximum ());
      //} // end loop over phi
      //float delta = log10 (max) - log10 (min);
      //min = pow (10, log10 (min) - 0.1*delta);
      //max = pow (10, log10 (max) + 0.1*delta);
      //SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]);

          if (!h) continue;

          h->SetFillColorAlpha (fillColors[iPtZ-2], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          h->SetMarkerStyle (kFullCircle);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
          useTrkPt ? h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x} d#Delta#phi");

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

          h->DrawCopy (plotNewAxes && iPtZ-3 == 0 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = markerStyles[iPtZ-3];

          TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]);
          //if (iCent == 3)
          //  g->Print ("ALL");
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerColor (colors[iPtZ-2]);
            g->SetLineColor (colors[iPtZ-2]);
            g->SetMarkerSize (1);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ-2]);
            g->SetFillColorAlpha (fillColors[iPtZ-2], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
          g->GetYaxis ()->SetRangeUser (min, max);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
          useTrkPt ? g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x} d#Delta#phi");

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
            string drawString = string (!canvasExists && iPtZ-3 == 0 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ-3 == 0 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }
        } // end loop over phi
      }

      if (!hasBkg)
        continue;

      if (plotAsSystematic)
        continue;

      bottomPad->cd ();
      GetDrawnObjects ();
      plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
      gPad->SetLogx ();
      gPad->SetLogy ();

      min = 1e30, max = 0;
      GetMinAndMax (min, max, true);
      for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
        TH1D* h = (useTrkPt ? h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]);
        min = fmin (min, h->GetMinimum (0));
        max = fmax (max, h->GetMaximum ());
      } // end loop over phi
      //delta = max - min;
      //min = min - 0.3*delta;
      //max = max + 0.3*delta;
      float delta = log10 (max) - log10 (min);
      min = pow (10, log10 (min) - 0.1*delta);
      max = pow (10, log10 (max) + 0.1*delta);
      SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]);

          if (!h) continue;

          h->SetFillColorAlpha (fillColors[iPtZ-2], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          h->SetMarkerStyle (kFullCircle);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
          h->GetYaxis ()->SetTitle ("Y / Y_{bkg}");
          h->GetYaxis ()->CenterTitle ();

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

          h->DrawCopy (plotNewAxes && iPtZ-3 == 0 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = markerStyles[iPtZ-3];
          TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerColor (colors[iPtZ-2]);
            g->SetLineColor (colors[iPtZ-2]);
            g->SetMarkerSize (1);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ-2]);
            g->SetFillColorAlpha (fillColors[iPtZ-2], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
          g->GetYaxis ()->SetRangeUser (min, max);

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
          g->GetYaxis ()->SetTitle ("Y / Y_{bkg}");
          g->GetYaxis ()->CenterTitle ();

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
            string drawString = string (!canvasExists && iPtZ-3 == 0 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ-3 == 0 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }
        } // end loop over phi
      }
    } // end loop over cents
    
    c->SaveAs (Form ("%s/TrkYields/%s_dists_zPt_%s.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh", spc));
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track pT distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelTrkYieldZPt (const short iCent, const short iPtZ, const short iSpc) {
  //const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);

  if (iCent == 0)
    myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
  else
    myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

  if (iCent == 2 && iSpc == 0)
    myText (0.72, 0.88, kBlack, "Z#rightarrow ee", 0.06);
  else if (iCent == 2 && iSpc == 1)
    myText (0.72, 0.88, kBlack, "Z#rightarrow#mu#mu", 0.06);

  if (iCent == 0)
    myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
  else if (iCent == 1) {
    const char* lo = GetPiString (phiLowBins[1]);
    const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
    myText (0.485, 0.88, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.06);
  }
  else if (iCent == numCentBins-1) {
    myText (0.38, 0.90, kBlack, "MB", 0.06);
    myText (0.47, 0.90, kBlack, "Z-tag", 0.06);
    myLineText (0.47, 0.85-0.06*(iPtZ-3), colors[iPtZ-2], iPtZ, "", 2.0, 0.054) ;
    myMarkerTextNoLine (0.54, 0.85-0.06*(iPtZ-3), colors[iPtZ-2], markerStyles[iPtZ-3], "", 1.5, 0.054); // for plotting data vs bkg.

    if (iPtZ == nPtZBins-1)
      myText (0.56, 0.84-0.06*(iPtZ-3), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.054);
    else
      myText (0.56, 0.84-0.06*(iPtZ-3), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrkYieldZPtSpcComp (const bool useTrkPt, const bool plotAsSystematic) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  const char* canvasName = Form ("c_TrkYield_zpt_%s", useTrkPt ? "pttrk" : "xhz");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 395*numCentBins, 685);
    gDirectory->Add (c);
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* topPadName = Form ("p_TrkYield_pt_top_iCent%i", iCent);
    const char* bottomPadName = Form ("p_TrkYield_pt_bottom_iCent%i", iCent);

    TPad* topPad = nullptr, *bottomPad = nullptr;
    if (!canvasExists) {
      topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
      bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

      gDirectory->Add (topPad);
      gDirectory->Add (bottomPad);

      topPad->SetTopMargin (0.04);
      topPad->SetBottomMargin (0);
      topPad->SetLeftMargin (0.17);
      topPad->SetRightMargin (0.06);
      bottomPad->SetTopMargin (0);
      bottomPad->SetBottomMargin (0.20);
      bottomPad->SetLeftMargin (0.17);
      bottomPad->SetRightMargin (0.06);
      topPad->Draw ();
      bottomPad->Draw ();
    }
    else {
      topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
      bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
    }

    topPad->cd ();
    GetDrawnObjects ();
    //bool plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
    gPad->SetLogx ();
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetMinAndMax (min, max, true);
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      for (int iPtZ = 2; iPtZ < 3; iPtZ++) {
        TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]);
        min = fmin (min, h->GetMinimum (0));
        max = fmax (max, h->GetMaximum ());
      } // end loop over phi
    } // end loop over species
    float delta = log10 (max) - log10 (min);
    min = pow (10, log10 (min) - 0.1*delta);
    max = pow (10, log10 (max) + 0.1*delta);
    SetMinAndMax (min, max);

    for (short iSpc = 0; iSpc < 2; iSpc++) {
      for (int iPtZ = 2; iPtZ < 3; iPtZ++) {
        const Style_t markerStyle = (iSpc == 1 ? kOpenCircle : kFullCircle);

        TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]);
        //if (iCent == 3)
        //  g->Print ("ALL");
        RecenterGraph (g);

        if (!plotAsSystematic) {
          ResetXErrors (g);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerColor (colors[iPtZ-1]);
          g->SetLineColor (colors[iPtZ-1]);
          g->SetMarkerSize (1);
          g->SetLineWidth (2);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (1);
          g->SetLineColor (colors[iPtZ-1]);
          g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
        g->GetYaxis ()->SetRangeUser (min, max);

        g->GetXaxis ()->SetMoreLogLabels ();

        useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
        useTrkPt ? g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]") : g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{x} d#Delta#phi");

        g->GetXaxis ()->SetTitleFont (43);
        g->GetXaxis ()->SetTitleSize (axisTextSize);
        g->GetXaxis ()->SetLabelFont (43);
        g->GetXaxis ()->SetLabelSize (axisTextSize);

        g->GetYaxis ()->SetTitleFont (43);
        g->GetYaxis ()->SetTitleSize (axisTextSize);
        g->GetYaxis ()->SetLabelFont (43);
        g->GetYaxis ()->SetLabelSize (axisTextSize);

        g->GetXaxis ()->SetTitleOffset (1.5 * g->GetXaxis ()->GetTitleOffset ());
        g->GetYaxis ()->SetTitleOffset (1.2 * g->GetYaxis ()->GetTitleOffset ());

        if (!plotAsSystematic) {
          string drawString = string (!canvasExists && iSpc == 0 && iPtZ-2 == 0 ? "AP" : "P");
          g->Draw (drawString.c_str ());
        } else {
          string drawString = string (!canvasExists && iSpc == 0 && iPtZ-2 == 0 ? "A5P" : "5P");
          ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
          g->Draw ("2P");
        }
      } // end loop over phi
    } // end loop over species

    bottomPad->cd ();
    GetDrawnObjects ();
    //plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
    gPad->SetLogx ();

    for (int iPtZ = 2; iPtZ < 3; iPtZ++) {

      TH1D* h = useTrkPt ? h_trk_pt_ptz_sub[0][iPtZ][iCent] : h_trk_xhz_ptz_sub[0][iPtZ][iCent];
      h->Divide (useTrkPt ? h_trk_pt_ptz_sub[1][iPtZ][iCent] : h_trk_xhz_ptz_sub[1][iPtZ][iCent]);

      TGAE* g = GetTGAE (h);
      delete h;
      //if (iCent == 3)
      //  g->Print("ALL");
      RecenterGraph (g);

      if (!plotAsSystematic) {
        ResetXErrors (g);
        g->SetMarkerColor (colors[iPtZ-1]);
        g->SetLineColor (colors[iPtZ-1]);
        g->SetMarkerSize (1);
        g->SetLineWidth (2);
      } else {
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (colors[iPtZ-1]);
        g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
      }

      useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
      g->GetYaxis ()->SetRangeUser (0, 2);

      g->GetXaxis ()->SetMoreLogLabels ();

      useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
      g->GetYaxis ()->SetTitle ("Electrons / Muons");

      g->GetXaxis ()->SetTitleFont (43);
      g->GetXaxis ()->SetTitleSize (axisTextSize);
      g->GetXaxis ()->SetLabelFont (43);
      g->GetXaxis ()->SetLabelSize (axisTextSize);

      g->GetYaxis ()->SetTitleFont (43);
      g->GetYaxis ()->SetTitleSize (axisTextSize);
      g->GetYaxis ()->SetLabelFont (43);
      g->GetYaxis ()->SetLabelSize (axisTextSize);

      g->GetXaxis ()->SetTitleOffset (1.5 * g->GetXaxis ()->GetTitleOffset ());
      g->GetYaxis ()->SetTitleOffset (1.2 * g->GetYaxis ()->GetTitleOffset ());

      if (!plotAsSystematic) {
        string drawString = string (!canvasExists && iPtZ-2 == 0 ? "AP" : "P");
        g->Draw (drawString.c_str ());
      } else {
        string drawString = string (!canvasExists && iPtZ-2 == 0 ? "A5P" : "5P");
        ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
        g->Draw ("2P");
      }
    } // end loop over phi

    topPad->cd ();
    if (iCent == 0) {
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
      myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    }
    else
      myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

    if (iCent == 1) {
      const char* lo = GetPiString (phiLowBins[1]);
      const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
      myText (0.485, 0.88, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.06);
    }

    for (short iPtZ = 2; iPtZ < 3; iPtZ++) {
    //for (short iPtZ = 2; iPtZ < 3; iPtZ++) {
      if (iPtZ == 3) {
        myText (0.34, 0.91, kBlack, "ee", 0.06);
        myText (0.47, 0.91, kBlack, "#mu#mu", 0.06);
      }
      myMarkerTextNoLine (0.42, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kFullCircle, "", 1.5, 0.054); // for plotting MC reco vs truth
      myMarkerTextNoLine (0.52, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kOpenCircle, "", 1.5, 0.054);

      if (iPtZ == nPtZBins-1)
        myText (0.56, 0.85-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.054);
      else
        myText (0.56, 0.85-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
    }
  } // end loop over cents
  
  c->SaveAs (Form ("%s/TrkYields/%s_dists_zPt_SpcComp.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh"));
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
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* ppHist = nullptr, *PbPbHist = nullptr;

      ppHist = h_trk_pt_ptz_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        PbPbHist = (TH1D*)(h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())));
        PbPbHist->Divide (ppHist);
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      }

      ppHist = h_trk_xhz_ptz_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        PbPbHist = (TH1D*)(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())));
        PbPbHist->Divide (ppHist);
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      }

      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        TH1D* ppHist = h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* PbPbHist = (TH1D*)(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())));
          PbPbHist->Divide (ppHist);
          h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
        } // end loop over cents
        ppHist = h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* PbPbHist = (TH1D*)(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())));
          PbPbHist->Divide (ppHist);
          h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
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
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_z_trk_%s_iaa_dPhi_%s_iPtZ%i", useTrkPt ? "pttrk" : "xzh", spc, iPtZ);
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
            TH1D* h = useTrkPt ? h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (0, 1.4);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
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

            useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (0, max_iaa);

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

      c->SaveAs (Form ("%s/IAA/iaa_%s_dPhi_%s_iPtZ%i.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh", spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAAdPhi (const short iCent, const short iPhi, const short iPtZ) {

  if (iCent == 1)
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
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
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_z_trk_%s_iaa_dCent_%s_iPtZ%i", useTrkPt ? "pttrk" : "xzh", spc, iPtZ);
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
            TH1D* h = useTrkPt ? h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (0, max_iaa);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              deltaize (g, 1+(iCent-2)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
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

            useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (0, max_iaa);
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

      c->SaveAs (Form ("%s/IAA/iaa_%s_dCent_%s_iPtZ%i.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh", spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAAdCent (const short iCent, const short iPhi, const short iPtZ) {

  if (iPhi == 1) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
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

    const char* canvasName = Form ("c_z_trk_z%s_iaa_%s", useTrkPt ? "zpttrk" : "zxhz", spc);
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

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        const Style_t markerStyle = markerStyles[iPtZ-2];

        TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]);
        RecenterGraph (g);

        if (!plotAsSystematic) {
          ResetXErrors (g);
          //deltaize (g, 1+(iPtZ-3)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
          g->SetLineColor (colors[iPtZ-1]);
          g->SetMarkerColor (colors[iPtZ-1]);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1.2);
          g->SetLineWidth (2);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (1);
          g->SetLineColor (colors[iPtZ-1]);
          g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
        g->GetYaxis ()->SetRangeUser (0, max_iaa);

        g->GetXaxis ()->SetMoreLogLabels ();

        useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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
          string drawString = string (!canvasExists && iPtZ-2 == 0 ? "AP" : "P");
          g->Draw (drawString.c_str ());
        } else {
          string drawString = string (!canvasExists && iPtZ-2 == 0 ? "A5P" : "5P");
          ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
          g->Draw ("2P");
        }

        LabelIAAdPtZ (iCent, iPtZ);
      } // end loop over pT^Z bins
    } // end loop over cents

    for (short iCent = 1; iCent < numCentBins; iCent++) {
      c->cd (iCent);
      myText (0.22, 0.22, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    } // end loop over cents

    c->SaveAs (Form ("%s/IAA/iaa_%s_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xhz", spc));
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAAdPtZ (const short iCent, const short iPtZ) {

  if (iCent == 1) {
    myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.065);
    myText (0.22, 0.78, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  }
  else if (iCent == 2) {
    const char* lo = GetPiString (phiLowBins[1]);
    const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
    myText (0.50, 0.88, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.05);
  }
  else if (iCent == numCentBins-1) {
    myMarkerTextNoLine (0.55, 0.894-0.06*(iPtZ-2), colors[iPtZ-1], markerStyles[iPtZ-2], "", 1.4, 0.05);
    if (iPtZ == nPtZBins-1) {
      myText (0.565, 0.89-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    }
    else {
      myText (0.565, 0.89-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between some Pb+Pb centrality and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSingleIAAdPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pPtZ, const short iCent, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  short iPtZLo = 2;
  short iPtZHi = nPtZBins;
  if (pPtZ != -1) {
    iPtZLo = pPtZ;
    iPtZHi = pPtZ+1;
  }

  const int axisTextSize = 35;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    //const char* canvasName = Form ("c_z_trk_zpt_%s_iaa_%s_iCent%i", useTrkPt ? "pttrk" : "xhz", spc, iCent);
    const char* canvasName = Form ("c_z_trk_zpt_%s_iaa_iCent%i", useTrkPt ? "pttrk" : "xhz", iCent);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 800, 800);
      gDirectory->Add (c);
    }

    double xmin = 0, xmax = 0;
    gPad->SetLogx ();
    gPad->SetLogy ();

    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      const Style_t markerStyle = markerStyles[iPtZ-3];

      TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]);
      RecenterGraph (g);

      if (!plotAsSystematic) {
        //deltaize (g, 1+((nPtZBins-2)*((int)useAltMarker)-iPtZ)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
        ResetXErrors (g);
        //deltaize (g, 1+(iPtZ-3)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
        g->SetLineColor (colors[iPtZ-iPtZLo+1]);
        g->SetMarkerColor (colors[iPtZ-iPtZLo+1]);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1.2);
        g->SetLineWidth (2);
      } else {
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (colors[iPtZ-iPtZLo+1]);
        g->SetFillColorAlpha (fillColors[iPtZ-iPtZLo+1], 0.3);
      }

      useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
      //g->GetYaxis ()->SetRangeUser (0.1, max_iaa);
      g->GetYaxis ()->SetRangeUser (0.15, max_iaa);

      g->GetXaxis ()->SetMoreLogLabels ();
      g->GetYaxis ()->SetMoreLogLabels ();

      useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

      g->GetXaxis ()->SetTitleOffset (1.2 * g->GetXaxis ()->GetTitleOffset ());
      //g->GetYaxis ()->SetTitleOffset (0.9 * g->GetYaxis ()->GetTitleOffset ());

      if (!plotAsSystematic) {
        string drawString = string (!canvasExists && iPtZ-iPtZLo == 0 ? "AP" : "P");
        g->Draw (drawString.c_str ());
      } else {
        string drawString = string (!canvasExists && iPtZ-iPtZLo == 0 ? "A5P" : "5P");
        ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
        g->Draw ("2P");
      }
    } // end loop over pT^Z bins

    if (!canvasExists) {
      myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
      myText (0.20, 0.83, kBlack, "Pb+Pb, 5.02 TeV", 0.045);
      const char* lo = GetPiString (phiLowBins[1]);
      const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
      myText (0.20, 0.22, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
      myText (0.20, 0.77, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.04);
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        myMarkerTextNoLine (0.674, 0.868-0.05*(iPtZ-iPtZLo), colors[iPtZ-iPtZLo+1], markerStyles[iPtZ-3], "", 2.0, 0.04);
        if (iPtZ == nPtZBins-1)
          myText (0.690, 0.86-0.05*(iPtZ-iPtZLo), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.032);
        else
          myText (0.690, 0.86-0.05*(iPtZ-iPtZLo), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.032);
      }

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    }

    //c->SaveAs (Form ("%s/IAA/iaa_%s_iCent%i_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xhz", iCent, spc));
    c->SaveAs (Form ("%s/IAA/iaa_%s_iCent%i_dPtZ.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xhz", iCent));
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between some Pb+Pb centrality and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAASpcComp (const bool useTrkPt, const bool plotAsSystematic) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  short iPtZLo = 3;
  short iPtZHi = nPtZBins;
  short iCentLo = 1;
  short iCentHi = numCentBins;

  const int axisTextSize = 22;

  const char* canvasName = Form ("c_z_trk_zpt_%s_iaa_spccomp", useTrkPt ? "pttrk" : "xhz");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 400*(numCentBins-1), 800);
    c->Divide (numCentBins-1, 2);
    gDirectory->Add (c);
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    double xmin = 0, xmax = 0;

    const Style_t markerStyle = (iSpc == 0 ? kOpenCircle : kFullCircle);

    for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);

        gPad->SetLogx ();

        TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]);
        RecenterGraph (g);

        if (!plotAsSystematic) {
          //deltaize (g, 1+((nPtZBins-2)*((int)useAltMarker)-iPtZ)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
          ResetXErrors (g);
          //deltaize (g, 1+(iPtZ-3)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
          g->SetLineColor (colors[iPtZ-iPtZLo+1]);
          g->SetMarkerColor (colors[iPtZ-iPtZLo+1]);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1.2);
          g->SetLineWidth (2);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (1);
          g->SetLineColor (colors[iPtZ-iPtZLo+1]);
          g->SetFillColorAlpha (fillColors[iPtZ-iPtZLo+1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
        //g->GetYaxis ()->SetRangeUser (0.1, max_iaa);
        g->GetYaxis ()->SetRangeUser (0, max_iaa);

        g->GetXaxis ()->SetMoreLogLabels ();
        g->GetYaxis ()->SetMoreLogLabels ();

        useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

        g->GetXaxis ()->SetTitleOffset (2.0 * g->GetXaxis ()->GetTitleOffset ());
        g->GetYaxis ()->SetTitleOffset (2.0 * g->GetYaxis ()->GetTitleOffset ());

        if (!plotAsSystematic) {
          string drawString = string (!canvasExists && iSpc == 0 ? "AP" : "P");
          g->Draw (drawString.c_str ());
        } else {
          string drawString = string (!canvasExists && iSpc == 0 ? "A5P" : "5P");
          ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
          g->Draw ("2P");
        }
      } // end loop over iPtZ
    } // end loop over iCent

    if (!canvasExists) {
      for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
        for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
          c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);
          if (iPtZ == iPtZHi-1)
            myText (0.65, 0.85, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
          if (iCent == iCentLo && iPtZ == iPtZLo) {
            myText (0.25, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.08);
            myText (0.25, 0.75, kBlack, "Pb+Pb, 5.02 TeV", 0.07);
          }
          if (iCent == iCentLo+1 && iPtZ == iPtZLo) {
            const char* lo = GetPiString (phiLowBins[1]);
            const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
            myText (0.25, 0.85, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.06);
          }
        }
        for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
          c->cd (iCent-iCentLo + 1);
          if (iCent == iCentLo+2) {
            if (iPtZ == iPtZLo) {
              myText (0.356, 0.875, kBlack, "ee", 0.06);
              myText (0.436, 0.875, kBlack, "#mu#mu", 0.06);
            }
            myMarkerTextNoLine (0.490, 0.810-0.10*(iPtZ-iPtZLo), colors[iPtZ-iPtZLo+1], kFullCircle, "", 1.8, 0.06);
            myMarkerTextNoLine (0.410, 0.810-0.10*(iPtZ-iPtZLo), colors[iPtZ-iPtZLo+1], kOpenCircle, "", 1.8, 0.06);
            if (iPtZ == nPtZBins-1)
              myText (0.500, 0.80-0.10*(iPtZ-iPtZLo), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.06);
            else
              myText (0.500, 0.80-0.10*(iPtZ-iPtZLo), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.06);
          }

        } // end loop over iPtZ
      } // end loop over iCent
    }
  } // end loop over species
  for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);
      const float xmin = useTrkPt ? trk_min_pt : allXHZBins[0];
      const float xmax = useTrkPt ? ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]] : allXHZBins[maxNXHZBins];
      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    }
  } // end loop over iCent

  //if (!plotAsSystematic) {
  if (false) {
    for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);

        for (short iSpc = 0; iSpc < 2; iSpc++) {
          TGAE* g1 = GetTGAE (useTrkPt ? h_trk_pt_ptz_iaa[0][iPtZ][iCent] : h_trk_xhz_ptz_iaa[0][iPtZ][iCent]);
          TGAE* g2 = GetTGAE (useTrkPt ? h_trk_pt_ptz_iaa[1][iPtZ][iCent] : h_trk_xhz_ptz_iaa[1][iPtZ][iCent]);

          double chisq = 0;
          double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
          for (int ix = 0; ix < g1->GetN (); ix++) {
            g1->GetPoint (ix, x1, y1);
            g2->GetPoint (ix, x2, y2);

            if (y1 > y2) {
              chisq += pow (y2-y1, 2)/ (pow (g2->GetErrorYhigh (ix), 2) + pow (g1->GetErrorYlow (ix), 2));
            
            }
            else if (y1 < y2) {
              chisq += pow (y2-y1, 2)/ (pow (g1->GetErrorYhigh (ix), 2) + pow (g2->GetErrorYlow (ix), 2));
            }
          }
          if (iPtZ == 4)
            myText (0.50, 0.61, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chisq, 6), 0.06);
          else
            myText (0.50, 0.61, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chisq, 5), 0.06);
        } // end loop over species
      } // end loop over iPtZ
    } // end loop over iCent
  }


  c->SaveAs (Form ("%s/IAASpcComp/iaa_%s_dPtZ.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xhz"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: CalculateICP () {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (icpCalculated)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* periphHist = nullptr, *centHist = nullptr;

      periphHist = h_trk_pt_ptz_sub[iSpc][iPtZ][1];
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        centHist = (TH1D*)(h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone ());
        centHist->Divide (periphHist);
        h_trk_pt_ptz_icp[iSpc][iPtZ][iCent] = centHist;
      }
      periphHist = h_trk_xhz_ptz_sub[iSpc][iPtZ][1];
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        centHist = (TH1D*)(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone ());
        centHist->Divide (periphHist);
        h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent] = centHist;
      }

      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        periphHist = h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]) {
            centHist = (TH1D*)(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent] = centHist;
          }
        } // end loop over cents

        periphHist = h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent]) {
            centHist = (TH1D*)(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent] = centHist;
          }
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species
  icpCalculated = true;
  return;
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
/*
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
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
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
            TH1D* h = useTrkPt ? h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (0, max_icp);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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
            
            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSystematic) {
              ResetXErrors (g);
              deltaize (g, 1+(iPhi-1.5)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
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

            useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (0, max_icp);
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

      c->SaveAs (Form ("%s/ICP/icp_%s_dPhi_%s_iPtZ%i.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh", spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: LabelICPdPhi (const short iCent, const short iPhi, const short iPtZ) {
  if (iCent == 2) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
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
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
/*
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
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
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
            TH1D* h = useTrkPt ? h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (0, max_icp);
            xmin = h->GetXaxis ()->GetXmin ();
            xmax = h->GetXaxis ()->GetXmax ();

            h->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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
            
            TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent] : h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent]);
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

            useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (0, max_icp);
            xmin = g->GetXaxis ()->GetXmin ();
            xmax = g->GetXaxis ()->GetXmax ();

            g->GetXaxis ()->SetMoreLogLabels ();

            useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

      c->SaveAs (Form ("%s/ICP/icp_%s_dCent_%s_iPtZ%i.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh", spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: LabelICPdCent (const short iCent, const short iPhi, const short iPtZ) {
  if (iPhi == 1) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
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
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: PlotICPdPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!icpCalculated)
    CalculateICP ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_z_trk_%s_icp_dPtZ_%s", useTrkPt ? "pttrk" : "xzh", spc);
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
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = useTrkPt ? h_trk_pt_ptz_icp[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent];

          h->SetFillColorAlpha (fillColors[iPtZ-1], fillAlpha);
          h->SetMarkerSize (0);
          h->SetLineColor (kBlack);
          h->SetLineWidth (0);

          useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : h->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
          h->GetYaxis ()->SetRangeUser (0, max_icp);
          xmin = h->GetXaxis ()->GetXmin ();
          xmax = h->GetXaxis ()->GetXmax ();

          h->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : h->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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

          h->DrawCopy (!canvasExists && iPtZ == 2 ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");

          LabelICPdPtZ (iCent, iPtZ);
        } // end loop over phi
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
          
          TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_icp[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]);
          RecenterGraph (g);

          if (!plotAsSystematic) {
            ResetXErrors (g);
            deltaize (g, 1+((numCentBins-2)*((int)useAltMarker)-iCent)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
            g->SetLineColor (colors[iPtZ-1]);
            g->SetMarkerColor (colors[iPtZ-1]);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (1.2);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPtZ-1]);
            g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, ptTrkBins[nPtZBins-1][nPtTrkBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[nPtZBins-1]]);
          g->GetYaxis ()->SetRangeUser (0, max_icp);
          xmin = g->GetXaxis ()->GetXmin ();
          xmax = g->GetXaxis ()->GetXmax ();

          g->GetXaxis ()->SetMoreLogLabels ();

          useTrkPt ? g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]") : g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
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
            string drawString = string (!canvasExists && iPtZ == 2 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iPtZ == 2 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }

          LabelICPdPtZ (iCent, iPtZ);
        } // end loop over phi
      }
    } // end loop over cents

    for (short iCent = 2; iCent < numCentBins; iCent++) {
      c->cd (iCent-1);
      myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    } // end loop over cents

    c->SaveAs (Form ("%s/ICP/icp_%s_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTTrk":"xzh", spc));
  } // end loop over species
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for I_AA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: LabelICPdPtZ (const short iCent, const short iPtZ) {
  if (iCent == 2) {
    myText (0.50, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
    myText (0.50, 0.82, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
    const char* lo = GetPiString (phiLowBins[1]);
    const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
    myText (0.50, 0.76, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.05);
  }
  else if (iCent == numCentBins-1) {
    myMarkerTextNoLine (0.582, 0.912-0.06*(iPtZ-2), colors[iPtZ-1], kFullCircle, "", 1.4, 0.05);
    if (iPtZ == nPtZBins-1)
      myText (0.60, 0.89-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    else
      myText (0.60, 0.89-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
  }
}
*/



void PhysicsAnalysis :: WriteIAAs () {
  SetupDirectories ("", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;

  const char* outFileName = "DataAnalysis/Nominal/data18hi_iaa_fits.root"; 
  TFile* outFile = new TFile (Form ("%s/%s", rootPath.Data (), outFileName), "recreate");

  TF1**** f_z_trk_zpt_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);
  TF1**** f_z_trk_zxzh_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_z_trk_zpt_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Fit (f_z_trk_zpt_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_z_trk_zxzh_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[iPtZ]]);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Fit (f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->Write ();
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->Write ();
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
      }
    }
  }

  outFile->Close ();

  Delete3DArray (f_z_trk_zpt_iaa, 3, nPtZBins, numCentBins);
  Delete3DArray (f_z_trk_zxzh_iaa, 3, nPtZBins, numCentBins);

  _gDirectory->cd ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots signal-to-background panels for data
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSignalToBkg (const bool useTrkPt, const short iSpc) {
  TCanvas* c = new TCanvas ("c_stb", "", 1500, 800);
  c->Divide (2, 1);

  for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    c->cd (iPtZ-3+1);
    gPad->SetLogx ();
    gPad->SetLogy ();

    TH1D* htemp = new TH1D (Form ("htemp_iPtZ%i", iPtZ), "", useTrkPt ? nPtTrkBins[iPtZ] : nXHZBins[iPtZ], useTrkPt ? ptTrkBins[iPtZ] : xHZBins[iPtZ]);
    htemp->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    htemp->GetYaxis ()->SetTitle ("Y / Y_{bkg}");

    htemp->GetXaxis ()->SetMoreLogLabels ();

    htemp->GetYaxis ()->SetRangeUser (8e-4, 2e4);
    htemp->Draw ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TGAE* g = GetTGAE (useTrkPt ? h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] : h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]);
      RecenterGraph (g);
      ResetXErrors (g);
      deltaize (g, 1 + 0.01*(iCent - (numCentBins-1)), true);
      
      g->SetMarkerStyle (markerStyles[iCent]);
      g->SetMarkerSize (markerStyles[iCent] == kFullDiamond ? 1.9 : 1.4);
      g->SetMarkerColor (colors[iCent]);
      g->SetLineColor (colors[iCent]);
      g->SetLineWidth (2);
      g->Draw ("P");
    }

    if (iPtZ == 3) {
      myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
      myText (0.22, 0.83, kBlack, "Pb+Pb, 5.02 TeV, 1.7 nb^{-1}", 0.04);
      myText (0.65, 0.24, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.04);
    }
    else if (iPtZ == 4) {
      myMarkerTextNoLine (0.25, 0.88, colors[0], markerStyles[0], "#it{pp}, 258 pb^{-1}", 1.4, 0.04);
      myMarkerTextNoLine (0.25, 0.83, colors[1], markerStyles[1], "30-80%", 1.4, 0.04);
      myMarkerTextNoLine (0.25, 0.78, colors[2], markerStyles[2], "10-30%", 1.9, 0.04);
      myMarkerTextNoLine (0.25, 0.73, colors[3], markerStyles[3], "0-10%", 1.4, 0.04);
      myText (0.65, 0.24, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.04);
    }
  }

  c->SaveAs (Form ("%s/TrkYields/sigToBkg_%s_dPtZ.pdf", plotPath.Data (), useTrkPt ? "pTTrk" : "xhz"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots pull between electron & muon channels
////////////////////////////////////////////////////////////////////////////////////////////////
//void PhysicsAnalysis :: PlotPullDist (const bool useTrkPt = true) {
//  TCanvas* c = new TCanvas ("c_pull", "", 800, 600);
//  TH1D* h_pull = new TH1D ("h_pull", ";#Delta#it{I}_{AA}
//}



void PhysicsAnalysis :: PrintIAA (const bool printErrs, const bool useTrkPt, const short iCent, const short iPtZ, const short iSpc) {
  TH1D* h = nullptr;
  if (useTrkPt)
    h = h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent];
  else
    h = h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent];
  cout << name << "\t";
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    cout << h->GetBinContent (ix) << "\t";
    if (printErrs)
      cout << h->GetBinError (ix) << "\t";
    cout << endl;
  }
}

#endif

