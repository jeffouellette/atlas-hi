#ifndef __PhysicsAnalysis_h__
#define __PhysicsAnalysis_h__

#include "Params.h"
#include "EventPlaneCalibrator.h"

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
  TFile* histFile   = nullptr;
  bool histsLoaded  = false;
  bool histsScaled  = false;

  public:
  bool histsUnfolded = false;
  bool plotFill       = false; // whether to plot as filled (bar) graph or points w/ errors
  bool plotSignal     = true; // whether to plot background subtracted plots
  bool useAltMarker   = false; // whether to plot as open markers (instead of closed)

  bool isMC           = false;
  bool useCentWgts    = false; // whether to reweight this analysis to the data centrality distribution (important for getting correct tracking efficiencies)
  bool useQ2Wgts      = false; // whether to reweight this analysis to the data |q2| distribution
  bool usePsi2Wgts    = false; // whether to reweight this analysis to the data psi2 distribution
  bool is2015Conds    = false; // whether this analysis uses 2015 data (different conditions)
  bool useHITight     = false; // whether to use HITight tracking efficiencies
  bool useHijingEffs  = false; // whether to use tracking efficiencies derived from Hijing
  bool doLeptonRejVar = false; // whether to impose an additional dR cut on tracks away from the leptons
  bool doTrackPurVar  = false; // whether to impose an additional correction based on tracking purity
  bool doTrackEffVar  = false; // whether to use pions-only tracking efficiency variation
  float trkEffNSigma  = 0; // how many sigma to vary the track efficiency by (-1,0,+1 suggested)
  float trkPurNSigma  = 0; // how many sigma to vary the track purity by (-1,0,+1 suggested)
  bool doPPTransMinMixing = false; // by default analyses are not performing trans-min mixing. Only really applies to pp bkg.

  // Analysis checks
  TH1D*   h_fcal_et               = nullptr;
  TH1D*   h_fcal_et_reweighted    = nullptr;
  TH1D**  h_centrality            = Get1DArray <TH1D*> (3);
  TH1D**  h_centrality_reweighted = Get1DArray <TH1D*> (3);

  TH1D**  h_q2                  = Get1DArray <TH1D*> (numFineCentBins);
  TH1D**  h_q2_reweighted       = Get1DArray <TH1D*> (numFineCentBins);
  TH1D**  h_psi2                = Get1DArray <TH1D*> (numFineCentBins);
  TH1D**  h_psi2_reweighted     = Get1DArray <TH1D*> (numFineCentBins);
  TH1D*   h_PbPb_vz             = nullptr;
  TH1D*   h_PbPb_vz_reweighted  = nullptr;
  TH1D*   h_pp_vz               = nullptr;
  TH1D*   h_pp_vz_reweighted    = nullptr;
  TH1D*   h_pp_nch              = nullptr;
  TH1D*   h_pp_nch_reweighted   = nullptr;

  // Event info distributions (for reweighting)
  TH1D*** h_PbPbFCal_weights   = Get2DArray <TH1D*> (3, nPtZBins+1);                  // iSpc, iPtZ
  TH1D**** h_PbPbQ2_weights    = Get3DArray <TH1D*> (3, numFineCentBins, nPtZBins+1); // iSpc, iFineCent, iPtZ
  TH1D**** h_PbPbPsi2_weights  = Get3DArray <TH1D*> (3, numFineCentBins, nPtZBins+1); // iSpc, iFineCent, iPtZ

  // Event plane calibration information
  EventPlaneCalibrator eventPlaneCalibrator;

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
  TF1**** f_trk_pt_ptz_binMigration  = Get3DArray <TF1*> (2, nPtZBins, numBBBCorrCentBins); // iSpc, iPtZ, iCent
  TF1**** f_trk_xhz_ptz_binMigration = Get3DArray <TF1*> (2, nPtZBins, numBBBCorrCentBins); // iSpc, iPtZ, iCent

  // Physics plots
  TH1D*****   h_trk_dphi        = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins, numCentBins); // iSpc, iPtZ, iPtch, iCent
  TH2D*****   h2_trk_dphi_cov   = Get4DArray <TH2D*> (3, nPtZBins, maxNPtchBins, numCentBins); // iSpc, iPtZ, iPtch, iCent
  TH1D*****   h_trk_dphi_sub    = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins, numCentBins); // iSpc, iPtZ, iPtch, iCent

  TH1D*****   h_trk_pt_dphi_raw = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);   // iSpc, iPtZ, iPhi, iCent
  TH1D****    h_z_counts        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);               // iSpc, iPtZ, iCent

  TH1D*****  h_trk_pt_dphi            = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH2D*****  h2_trk_pt_dphi_cov       = Get4DArray <TH2D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_sub        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_sig_to_bkg = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_iaa        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  //TH1D*****  h_trk_pt_dphi_icp        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_trk_pt_ptz             = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH2D****   h2_trk_pt_ptz_cov        = Get3DArray <TH2D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_sub         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_sig_to_bkg  = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_iaa         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  //TH1D****   h_trk_pt_ptz_icp         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  TH1D***** h_trk_xhz_dphi            = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH2D***** h2_trk_xhz_dphi_cov       = Get4DArray <TH2D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_sub        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_sig_to_bkg = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_iaa        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  //TH1D***** h_trk_xhz_dphi_icp        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_trk_xhz_ptz            = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH2D****   h2_trk_xhz_ptz_cov       = Get3DArray <TH2D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
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
    Delete1DArray (h_centrality,            3);
    Delete1DArray (h_centrality_reweighted, 3);
    Delete1DArray (h_q2,                    numFineCentBins);
    Delete1DArray (h_q2_reweighted,         numFineCentBins);
    Delete1DArray (h_psi2,                  numFineCentBins);
    Delete1DArray (h_psi2_reweighted,       numFineCentBins);

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

    Delete3DArray (f_trk_pt_ptz_binMigration,  2, nPtZBins, numBBBCorrCentBins);
    Delete3DArray (f_trk_xhz_ptz_binMigration, 2, nPtZBins, numBBBCorrCentBins);

    Delete4DArray (h_trk_dphi,                3, nPtZBins, maxNPtchBins, numCentBins);
    Delete4DArray (h2_trk_dphi_cov,           3, nPtZBins, maxNPtchBins, numCentBins);
    Delete4DArray (h_trk_dphi_sub,            3, nPtZBins, maxNPtchBins, numCentBins);

    Delete4DArray (h_trk_pt_dphi_raw,         3, nPtZBins, numPhiBins, numCentBins);
    Delete3DArray (h_z_counts,                3, nPtZBins, numCentBins);

    Delete4DArray (h_trk_pt_dphi,             3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h2_trk_pt_dphi_cov,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_dphi_sub,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_dphi_sig_to_bkg,  3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_pt_dphi_iaa,         3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (h_trk_pt_dphi_icp,         3, nPtZBins, numPhiBins, numCentBins);

    Delete3DArray (h_trk_pt_ptz,              3, nPtZBins, numCentBins);
    Delete3DArray (h2_trk_pt_ptz_cov,         3, nPtZBins, numCentBins);
    Delete3DArray (h_trk_pt_ptz_sub,          3, nPtZBins, numCentBins);
    Delete3DArray (h_trk_pt_ptz_sig_to_bkg,   3, nPtZBins, numCentBins);
    Delete3DArray (h_trk_pt_ptz_iaa,          3, nPtZBins, numCentBins);
    //Delete3DArray (h_trk_pt_ptz_icp,          3, nPtZBins, numCentBins);

    Delete4DArray (h_trk_xhz_dphi,            3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h2_trk_xhz_dphi_cov,       3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_dphi_sub,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_dphi_sig_to_bkg, 3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (h_trk_xhz_dphi_iaa,        3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (h_trk_xhz_dphi_icp,        3, nPtZBins, numPhiBins, numCentBins);

    Delete3DArray (h_trk_xhz_ptz,             3, nPtZBins, numCentBins);
    Delete3DArray (h2_trk_xhz_ptz_cov,        3, nPtZBins, numCentBins);
    Delete3DArray (h_trk_xhz_ptz_sub,         3, nPtZBins, numCentBins);
    Delete3DArray (h_trk_xhz_ptz_sig_to_bkg,  3, nPtZBins, numCentBins);
    Delete3DArray (h_trk_xhz_ptz_iaa,         3, nPtZBins, numCentBins);
    //Delete3DArray (h_trk_xhz_ptz_icp,         3, nPtZBins, numCentBins);

    if (eventWeightsFile && eventWeightsFile->IsOpen ()) {
      eventWeightsFile->Close ();
      SaferDelete (eventWeightsFile);
    }
    if (trkEffFile && trkEffFile->IsOpen ()) {
      trkEffFile->Close ();
      SaferDelete (trkEffFile);
    }
    if (trkPurFile && trkPurFile->IsOpen ()) {
      trkPurFile->Close ();
      SaferDelete (trkPurFile);
    }
    if (histFile && histFile->IsOpen ()) {
      histFile->Close ();
      SaferDelete (histFile);
    }
  }


  protected:
  void LabelCorrelations (const short iPtZ, const short iPtch, const short iCent, const bool subBkg);
  void LabelIAA_dPhi (const short iCent, const short iPhi, const short iPtZ);
  void LabelIAA_dCent (const short iCent, const short iPhi, const short iPtZ);
  void LabelIAA_dPtZ (const short iCent, const short iPtZ);
  //void LabelICP_dPhi (const short iCent, const short iPhi, const short iPtZ);
  //void LabelICP_dCent (const short iCent, const short iPhi, const short iPtZ);
  //void LabelICP_dPtZ (const short iCent, const short iPtZ);

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
  virtual void GenerateEventWeights (const char* weightedSampleInFileName, const char* matchedSampleInFileName, const char* outFileName);
  virtual void LoadEventWeights ();
  virtual void SubtractBackground (PhysicsAnalysis* a = nullptr);
  virtual void UnfoldSubtractedYield ();
  virtual void InflateStatUnc (const float amount);
  virtual void SubtractSameSigns (PhysicsAnalysis* a);

  virtual void ApplyRelativeVariation (float**** relVar, const bool upVar = true); // multiplies yield results by relErr in each bin (or divides if not upVar)
  virtual void ConvertToStatVariation (const bool upVar = true, const float nSigma = 1); // adds or subtracts nSigma of statistical errors to analysis

  virtual void LoadTrackingEfficiencies (const bool doRebin = false); // defaults to HILoose
  virtual double GetTrackingEfficiency (const float fcal_et, float trk_pt, const float trk_eta, const bool isPbPb = true);

  virtual void LoadTrackingPurities (const bool doRebin = false); // defaults to HILoose
  virtual double GetTrackingPurity (const float fcal_et, float trk_pt, const float trk_eta, const bool isPbPb = true);

  void CorrectQ2Vector (float& q2x_a, float& q2y_a, float& q2x_c, float& q2y_c);

  void PrintZYields (const int iPtZ = 2);

  void PlotFCalDists (const bool _treatAsData = true);
  void PlotQ2Dists (const bool _treatAsData = true);
  void PlotQ2Weights (PhysicsAnalysis* a);
  void PlotPsi2Dists (const bool _treatAsData = true);
  void PlotPsi2Weights (PhysicsAnalysis* a);
  void PlotVZDists (const bool _treatAsData = true);
  void PlotNchDists (const bool _treatAsData = true);

  void PlotCorrelations (const short pSpc = 2, const short pPtZ = nPtZBins-1, const bool _subBkg = false);
  void PlotTrackingEfficiencies ();
  void PlotTrackingEfficienciesComparison (PhysicsAnalysis* a = nullptr);
  void PlotTrackingEfficiencies2D ();
  void PlotTrackingPurities ();
  void PlotTrackingPuritiesComparison (PhysicsAnalysis* a = nullptr);
  void PlotTrackingPurities2D ();

  void CalculateIAA ();
  //void CalculateICP ();

  virtual void PlotUnweightedTrkYields (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotAllYields_dPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotAllYields_dPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotAllYields_dPtZ_SpcComp (const bool useTrkPt = true, const bool plotAsSystematic = false);
  virtual void PlotAllCovMatrices ();
  virtual void PlotCovMatrix (const bool useTrkPt = true, const short pSpc = 0, const short pPtZ = nPtZBins-1, const short pCent = numCentBins-1);
  virtual void PlotIAA_dPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAA_dCent (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAA_dPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotSingleIAA_dPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pPtZ = -1, const short iCent = numCentBins-1, const short pSpc = 2);
  virtual void PlotIAA_dPtZ_SpcComp (const bool useTrkPt = true, const bool plotAsSystematic = false);
  //virtual void PlotICP_dPhi (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //virtual void PlotICP_dCent (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //virtual void PlotICP_dPtZ (const bool useTrkPt = true, const bool plotAsSystematic = false, const short pSpc = 2);
  virtual void PlotSignalToBkg (const bool useTrkPt = true, const short iSpc = 2);

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

  h_fcal_et = new TH1D (Form ("h_fcal_et_%s", name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins); 
  h_fcal_et->Sumw2 ();
  h_fcal_et_reweighted = new TH1D (Form ("h_fcal_et_reweighted_%s", name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
  h_fcal_et_reweighted->Sumw2 ();
  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig] = new TH1D (Form ("h_centrality_trig%i_%s", iMBTrig, name.c_str ()), "", 80, 0, 80);
    h_centrality[iMBTrig]->Sumw2 ();
    h_centrality_reweighted[iMBTrig] = new TH1D (Form ("h_centrality_reweighted_trig%i_%s", iMBTrig, name.c_str ()), "", 80, 0, 80);
    h_centrality_reweighted[iMBTrig]->Sumw2 ();
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    h_q2[iFineCent]               = new TH1D (Form ("h_q2_iCent%i_%s", iFineCent, name.c_str ()), "", 20, 0, 0.3);
    h_q2[iFineCent]->Sumw2 ();
    h_q2_reweighted[iFineCent]    = new TH1D (Form ("h_q2_reweighted_iCent%i_%s", iFineCent, name.c_str ()), "", 20, 0, 0.3);
    h_q2_reweighted[iFineCent]->Sumw2 ();
    h_psi2[iFineCent]             = new TH1D (Form ("h_psi2_iCent%i_%s", iFineCent, name.c_str ()), "", 8, -pi/2, pi/2);
    h_psi2[iFineCent]->Sumw2 ();
    h_psi2_reweighted[iFineCent]  = new TH1D (Form ("h_psi2_reweighted_iCent%i_%s", iFineCent, name.c_str ()), "", 8, -pi/2, pi/2);
    h_psi2_reweighted[iFineCent]->Sumw2 ();
  } // end loop over iFineCent
  h_PbPb_vz = new TH1D (Form ("h_PbPb_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_PbPb_vz_reweighted = new TH1D (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz = new TH1D (Form ("h_pp_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz_reweighted = new TH1D (Form ("h_pp_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_nch = new TH1D (Form ("h_pp_nch_%s", name.c_str ()), "", 80, -0.5, 160.5);
  h_pp_nch_reweighted = new TH1D (Form ("h_pp_nch_reweighted_%s", name.c_str ()), "", 80, -0.5, 160.5);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]); // old name: h_z_trk_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent] = new TH2D (Form ("h2_trk_pt_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ], nPtchBins[iPtZ], pTchBins[iPtZ]);
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]); // old name: h_z_trk_xzh
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent] = new TH2D (Form ("h2_trk_xhz_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ], nXhZBins[iPtZ], xhZBins[iPtZ]);
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        } // end loop over iPhi
        for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (iPtch, iCent), -pi/2, 3*pi/2); // old name: h_z_trk_phi
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Sumw2 ();
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent] = new TH2D (Form ("h2_trk_dphi_cov_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (iPtch, iCent), -pi/2, 3*pi/2, GetNdPhiBins (iPtch, iCent), -pi/2, 3*pi/2);
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->Sumw2 ();
        } // end loop over iPtch
        h_trk_pt_ptz[iSpc][iPtZ][iCent] = new TH1D (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]); // old name: h_z_trk_zpt
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Sumw2 ();
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent] = new TH2D (Form ("h2_trk_pt_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ], nPtchBins[iPtZ], pTchBins[iPtZ]);
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->Sumw2 ();
        h_trk_xhz_ptz[iSpc][iPtZ][iCent] = new TH1D (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]); // old name: h_z_trk_zxzh
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Sumw2 ();
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent] = new TH2D (Form ("h2_trk_xhz_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ], nXhZBins[iPtZ], xhZBins[iPtZ]);
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->Sumw2 ();

        h_z_counts[iSpc][iPtZ][iCent] = new TH1D (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", 3, 0, 3);
        h_z_counts[iSpc][iPtZ][iCent]->Sumw2 ();
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

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
  h_fcal_et               = (TH1D*) a->h_fcal_et->Clone (Form ("h_fcal_et_%s", name.c_str ()));
  h_fcal_et_reweighted    = (TH1D*) a->h_fcal_et_reweighted->Clone (Form ("h_fcal_et_reweighted_%s", name.c_str ()));
  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig]             = (TH1D*) a->h_centrality[iMBTrig]->Clone (Form ("h_centrality_trig%i_%s", iMBTrig, name.c_str ()));
    h_centrality_reweighted[iMBTrig]  = (TH1D*) a->h_centrality[iMBTrig]->Clone (Form ("h_centrality_reweighted_trig%i_%s", iMBTrig, name.c_str ()));
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    h_q2[iFineCent]               = (TH1D*) a->h_q2[iFineCent]->Clone (Form ("h_q2_iCent%i_%s", iFineCent, name.c_str ()));
    h_q2_reweighted[iFineCent]    = (TH1D*) a->h_q2_reweighted[iFineCent]->Clone (Form ("h_q2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
    h_psi2[iFineCent]             = (TH1D*) a->h_psi2[iFineCent]->Clone (Form ("h_psi2_iCent%i_%s", iFineCent, name.c_str ()));
    h_psi2_reweighted[iFineCent]  = (TH1D*) a->h_psi2_reweighted[iFineCent]->Clone (Form ("h_psi2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
  } // end loop over iFineCent
  h_PbPb_vz             = (TH1D*) a->h_PbPb_vz->Clone (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted  = (TH1D*) a->h_PbPb_vz_reweighted->Clone (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz               = (TH1D*) a->h_pp_vz->Clone (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted    = (TH1D*) a->h_pp_vz_reweighted->Clone (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  h_pp_nch              = (TH1D*) a->h_pp_nch->Clone (Form ("h_pp_nch_%s", name.c_str ()));
  h_pp_nch_reweighted   = (TH1D*) a->h_pp_nch_reweighted->Clone (Form ("h_pp_nch_reweighted_%s", name.c_str ()));

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]       = (TH1D*) a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]  = (TH2D*) a->h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->Clone (Form ("h2_trk_pt_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]      = (TH1D*) a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent] = (TH2D*) a->h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->Clone (Form ("h2_trk_xhz_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]    = (TH1D*) a->h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]        = (TH1D*) a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]   = (TH2D*) a->h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h2_trk_pt_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]       = (TH1D*) a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]  = (TH2D*) a->h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h2_trk_xhz_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
        } // end loop over iPhi
        for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]       = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]  = (TH2D*) a->h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h2_trk_dphi_cov_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
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
              h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]         = (TH1D*) a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_sub
              h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]  = (TH1D*) a->h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_sig_to_bkg
            }
            if (a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]) {
              h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]        = (TH1D*) a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_sub
              h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_sig_to_bkg
            }

            for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
              if (a->h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]) {
                h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]        = (TH1D*) a->h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_sub
                h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_sig_to_bkg
              }
              if (a->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]) {
                h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]         = (TH1D*) a->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_sub
                h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]  = (TH1D*) a->h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_sig_to_bkg
              }
            } // end loop over iPhi

            for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
              if (a->h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]) {
                h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent] = (TH1D*) a->h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_dphi
              }
            } // end loop over iPtch
          } // end loop over iPtZ
        } // end loop over iCent
      } // end loop over iSpc
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
            } // end loop over iPhi
          } // end loop over iCent
        } // end loop over iPtZ
      } // end loop over iSpc
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
            } // end loop over iPhi
          } // end loop over iCent
        } // end loop over iPtZ
      } // end loop over iSpc
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
        if (h_trk_pt_ptz[iSpc][iPtZ][iCent])        SaferDelete (h_trk_pt_ptz[iSpc][iPtZ][iCent]);
        if (h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent])   SaferDelete (h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz[iSpc][iPtZ][iCent])       SaferDelete (h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
        if (h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent])  SaferDelete (h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent])   SaferDelete (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent])       SaferDelete (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent])      SaferDelete (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi
        for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          if (h_trk_dphi[iSpc][iPtZ][iPtch][iCent])      SaferDelete (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]);
          if (h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]) SaferDelete (h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]);
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
        } // end loop over iPhi
        for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          if (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]) SaferDelete (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]);
        } // end loop over iPtch
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc
  backgroundSubtracted = false;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        if (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent])   SaferDelete (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]);
        if (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent])  SaferDelete (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]);

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent])  SaferDelete (h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
          if (h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]) SaferDelete (h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
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
        } // end loop over iPhi
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
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
        for (short iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) histFile->Get (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent] = (TH2D*) histFile->Get (Form ("h2_trk_dphi_cov_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
        } // end loop over iPtch

        h_trk_pt_ptz[iSpc][iPtZ][iCent]       = (TH1D*) histFile->Get (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]  = (TH2D*) histFile->Get (Form ("h2_trk_pt_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]      = (TH1D*) histFile->Get (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent] = (TH2D*) histFile->Get (Form ("h2_trk_xhz_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]    = (TH1D*) histFile->Get (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]        = (TH1D*) histFile->Get (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]   = (TH2D*) histFile->Get (Form ("h2_trk_pt_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]       = (TH1D*) histFile->Get (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]  = (TH2D*) histFile->Get (Form ("h2_trk_xhz_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
        } // end loop over iPhi
      } // end loop over iPtZ
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  histsLoaded = true;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h_z_counts[2][iPtZ][iCent]->Add (h_z_counts[iSpc][iPtZ][iCent]);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]) h_trk_pt_dphi_raw[2][iPtZ][iPhi][iCent]->Add (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  double newbins[numSuperFineCentBins];
  for (int i = 0; i < numSuperFineCentBins; i++)
    newbins[i] = superFineCentBins[i]/1000.;

  h_fcal_et             = new TH1D (Form ("h_fcal_et_fixedBins_%s", name.c_str ()), "", sizeof (newbins)/sizeof (newbins[0]) - 1, newbins);
  h_fcal_et_reweighted  = new TH1D (Form ("h_fcal_et_reweighted_fixedBins_%s", name.c_str ()), "", sizeof (newbins)/sizeof (newbins[0]) - 1, newbins);

  TH1D* temp = (TH1D*) histFile->Get (Form ("h_fcal_et_%s", name.c_str ()));
  assert (h_fcal_et->GetNbinsX () == temp->GetNbinsX ());
  for (int iX = 1; iX <= temp->GetNbinsX (); iX++)
    h_fcal_et->SetBinContent (iX, temp->GetBinContent (iX));

  temp = (TH1D*) histFile->Get (Form ("h_fcal_et_reweighted_%s", name.c_str ()));
  assert (h_fcal_et_reweighted->GetNbinsX () == temp->GetNbinsX ());
  for (int iX = 1; iX <= temp->GetNbinsX (); iX++)
    h_fcal_et_reweighted->SetBinContent (iX, temp->GetBinContent (iX));

  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig]             = (TH1D*) histFile->Get (Form ("h_centrality_trig%i_%s", iMBTrig, name.c_str ()));
    h_centrality_reweighted[iMBTrig]  = (TH1D*) histFile->Get (Form ("h_centrality_reweighted_trig%i_%s", iMBTrig, name.c_str ()));
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    h_q2[iFineCent]               = (TH1D*) histFile->Get (Form ("h_q2_iCent%i_%s", iFineCent, name.c_str ()));
    h_q2_reweighted[iFineCent]    = (TH1D*) histFile->Get (Form ("h_q2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
    h_psi2[iFineCent]             = (TH1D*) histFile->Get (Form ("h_psi2_iCent%i_%s", iFineCent, name.c_str ()));
    h_psi2_reweighted[iFineCent]  = (TH1D*) histFile->Get (Form ("h_psi2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
  } // end loop over iFineCent
  h_PbPb_vz             = (TH1D*) histFile->Get (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted  = (TH1D*) histFile->Get (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz               = (TH1D*) histFile->Get (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted    = (TH1D*) histFile->Get (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  h_pp_nch               = (TH1D*) histFile->Get (Form ("h_pp_nch_%s", name.c_str ()));
  h_pp_nch_reweighted    = (TH1D*) histFile->Get (Form ("h_pp_nch_reweighted_%s", name.c_str ()));
  

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

        for (short iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          SafeWrite (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]);
          SafeWrite (h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]);
        } // end loop over iPtch

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SafeWrite (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi

        SafeWrite (h_trk_pt_ptz[iSpc][iPtZ][iCent]);
        SafeWrite (h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]);
        SafeWrite (h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
        SafeWrite (h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]);

        SafeWrite (h_z_counts[iSpc][iPtZ][iCent]);
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  SafeWrite (h_fcal_et);
  SafeWrite (h_fcal_et_reweighted);
  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    SafeWrite (h_centrality[iMBTrig]);
    SafeWrite (h_centrality_reweighted[iMBTrig]);
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    SafeWrite (h_q2[iFineCent]);
    SafeWrite (h_q2_reweighted[iFineCent]);
    SafeWrite (h_psi2[iFineCent]);
    SafeWrite (h_psi2_reweighted[iFineCent]);
  } // end loop over iFineCent

  SafeWrite (h_PbPb_vz);
  SafeWrite (h_PbPb_vz_reweighted);
  SafeWrite (h_pp_vz);
  SafeWrite (h_pp_vz_reweighted);
  SafeWrite (h_pp_nch);
  SafeWrite (h_pp_nch_reweighted);
  
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
        for (short iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          SafeWrite (h_trk_dphi[iSpc][iPtZ][iPtch][iCent],      Form ("h_ztrk_phi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
          SafeWrite (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent],  Form ("h_ztrk_phi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
        } // end loop over iPtZ

        SafeWrite (h_trk_pt_ptz[iSpc][iPtZ][iCent],                Form ("h_ztrk_pt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent],            Form ("h_ztrk_pt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent],     Form ("h_ztrk_pt_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz[iSpc][iPtZ][iCent],               Form ("h_ztrk_xhz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent],           Form ("h_ztrk_xhz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent],    Form ("h_ztrk_xhz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));

        SafeWrite (h_z_counts[iSpc][iPtZ][iCent],                 Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
   } // end loop over iCent

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent],            Form ("h_ztrk_pt_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent],           Form ("h_ztrk_xhz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent
  
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
        float spcWeight = countsHist->GetBinContent (2);
        countsHist = h_z_counts[2][iPtZ][iCent];
        if (countsHist->GetBinContent (2) > 0)
          spcWeight = spcWeight / countsHist->GetBinContent (2);
        else {
          cout << "Warning: In PhysicsAnalysis :: CombineHists: Found 0 total Z bosons in this bin, iCent = " << iCent << ", iPtZ = " << iPtZ << ", iSpc = " << iSpc << "; weight is set to 0!" << endl;
          spcWeight = 0;
        }

        for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          while (h_trk_dphi[2][iPtZ][iPtch][iCent]->GetNbinsX () > h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetNbinsX ())
            h_trk_dphi[2][iPtZ][iPtch][iCent]->Rebin (2);
          if (h_trk_dphi[iSpc][iPtZ][iPtch][iCent])     h_trk_dphi[2][iPtZ][iPtch][iCent]->Add      (h_trk_dphi[iSpc][iPtZ][iPtch][iCent], spcWeight);

          if (hasBkg) {
            while (h_trk_dphi_sub[2][iPtZ][iPtch][iCent]->GetNbinsX () > h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]->GetNbinsX ())
              h_trk_dphi_sub[2][iPtZ][iPtch][iCent]->Rebin (2);
            if (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]) h_trk_dphi_sub[2][iPtZ][iPtch][iCent]->Add  (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent], spcWeight);
          }
        } // end loop over iPtch

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent])              h_trk_pt_dphi[2][iPtZ][iPhi][iCent]->Add             (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent], spcWeight);
          if (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent])             h_trk_xhz_dphi[2][iPtZ][iPhi][iCent]->Add            (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], spcWeight);
          if (hasBkg) {
            if (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent])          h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]->Add         (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent], spcWeight);
            if (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent])   h_trk_pt_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Add  (h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], spcWeight);
            if (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent])         h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]->Add        (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent], spcWeight);
            if (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent])  h_trk_xhz_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Add (h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent], spcWeight);
          }
        } // end loop over iPhi

        if (h_trk_pt_ptz[iSpc][iPtZ][iCent])  h_trk_pt_ptz[2][iPtZ][iCent]->Add   (h_trk_pt_ptz[iSpc][iPtZ][iCent], spcWeight);
        if (h_trk_xhz_ptz[iSpc][iPtZ][iCent]) h_trk_xhz_ptz[2][iPtZ][iCent]->Add  (h_trk_xhz_ptz[iSpc][iPtZ][iCent], spcWeight);

        if (hasBkg) {
          if (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent])          h_trk_pt_ptz_sub[2][iPtZ][iCent]->Add         (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], spcWeight);
          if (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent])         h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Add        (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], spcWeight);
          if (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent])   h_trk_pt_ptz_sig_to_bkg[2][iPtZ][iCent]->Add  (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent], spcWeight);
          if (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent])  h_trk_xhz_ptz_sig_to_bkg[2][iPtZ][iCent]->Add (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent], spcWeight);
        }
      } // end loop over iSpc
    } // end loop over iPtZ
  } // end loop over iCent

  //InflateStatUnc (0.54);
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
        const float counts = countsHist->GetBinContent (2);
        if (counts <= 0)  continue;

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          // finalize covariance calculation by normalizing to bin widths and subtracting off product of means
          // then use diagonals of the covariance matrix to determine statistical uncertainties
          TH1D* h = h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
          h->Scale (1/counts);
          TH2D* h2 = h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent];
          assert (h2->GetNbinsX () == h->GetNbinsX ());
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (counts)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
          h2->Scale (1 / (countsHist->GetBinContent (1)*(counts-1)));

          for (int iX = 1; iX <= h->GetNbinsX (); iX++)
            h->SetBinError (iX, sqrt (h2->GetBinContent (iX, iX)));

          h->Scale (1. / (doPPTransMinMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi])), "width");
          h2->Scale (1. / pow (doPPTransMinMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi]), 2), "width");

          h = h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          h->Scale (1/counts);
          h2 = h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent];
          assert (h2->GetNbinsX () == h->GetNbinsX ());
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (counts)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
          h2->Scale (1 / (countsHist->GetBinContent (1)*(counts-1)));

          for (int iX = 1; iX <= h->GetNbinsX (); iX++)
            h->SetBinError (iX, sqrt (h2->GetBinContent (iX, iX)));

          h->Scale (1/ (doPPTransMinMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi])), "width");
          h2->Scale (1/ pow (doPPTransMinMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi]), 2), "width");
        } // end loop over iPhi


        // finalize covariance calculation by normalizing to bin widths and subtracting off product of means
        // then use diagonals of the covariance matrix to determine statistical uncertainties
        TH1D* h = h_trk_pt_ptz[iSpc][iPtZ][iCent];
        h->Scale (1/counts);
        TH2D* h2 = h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent];
        assert (h2->GetNbinsX () == h->GetNbinsX ());
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (counts)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
        h2->Scale (1 / (countsHist->GetBinContent (1)*(counts-1)));

        for (int iX = 1; iX <= h->GetNbinsX (); iX++)
          h->SetBinError (iX, sqrt (h2->GetBinContent (iX, iX)));

        h->Scale (1/ (doPPTransMinMixing && iCent == 0 ? pi/8. : pi/4.), "width");
        h2->Scale (1/ pow (doPPTransMinMixing && iCent == 0 ? pi/8. : pi/4., 2), "width");

        h = h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        h->Scale (1/counts);
        h2 = h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent];
        assert (h2->GetNbinsX () == h->GetNbinsX ());
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (counts)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
        h2->Scale (1 / (countsHist->GetBinContent (1)*(counts-1)));

        for (int iX = 1; iX <= h->GetNbinsX (); iX++)
          h->SetBinError (iX, sqrt (h2->GetBinContent (iX, iX)));

        h->Scale (1/ (doPPTransMinMixing && iCent == 0 ? pi/8. : pi/4.), "width");
        h2->Scale (1/ pow (doPPTransMinMixing && iCent == 0 ? pi/8. : pi/4., 2), "width");
        
        
        for (short iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          TH1D* h = h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          h->Scale (1/counts);
          TH2D* h2 = h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent];
          assert (h2->GetNbinsX () == h->GetNbinsX ());
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (counts)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
          h2->Scale (1 / (countsHist->GetBinContent (1)*(counts-1)));

          for (int iX = 1; iX <= h->GetNbinsX (); iX++)
            h->SetBinError (iX, sqrt (h2->GetBinContent (iX, iX)));

          h->Scale (1, "width");
          h2->Scale (1, "width");
        } // end loop over iPtch
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
// Designed to be overloaded. The default here is for analyzing data.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  cout << "Arguments provided: " << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "outFileName = " << outFileName << endl;

  LoadEventWeights ();
  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_weight = 1, q2_weight = 1, psi2_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  //float q2x_a = 0, q2y_a = 0, q2x_c = 0, q2y_c = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  int***    trks_counts   = Get3DArray <int> (2, 6, numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, 6, numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, 6, numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (6, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (6, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (6, 40);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    //PbPbTree->SetBranchAddress ("q2x_a",         &q2x_a);
    //PbPbTree->SetBranchAddress ("q2y_a",         &q2y_a);
    //PbPbTree->SetBranchAddress ("q2x_c",         &q2x_c);
    //PbPbTree->SetBranchAddress ("q2y_c",         &q2y_c);
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

      if (fabs (vz) > 150) continue;

      //{
      //  CorrectQ2Vector (q2x_a, q2y_a, q2x_c, q2y_c);
      //  const float q2x = q2x_a + q2x_c;
      //  const float q2y = q2y_a + q2y_c;
      //  q2 = sqrt (q2x*q2x + q2y*q2y) / fcal_et;
      //  psi2 = 0.5 * atan2 (q2y, q2x);
      //}

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1) continue;

      const short iFineCent = GetFineCentBin (fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1) continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      // do a reweighting procedure
      {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)  dphi = pi - dphi;
        if (useCentWgts)  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
        if (useQ2Wgts)    q2_weight   = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
        if (usePsi2Wgts)  psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

        event_weight *= fcal_weight * q2_weight * psi2_weight;
      }

      if (event_weight == 0) continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFineCent]->Fill (q2);
      h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
      h_psi2[iFineCent]->Fill (psi2);
      h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03)) continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (pTchBins[iPtZ][0] <= trkpt) {
          iPtch = 0;
          while (iPtch < nPtchBins[iPtZ] && pTchBins[iPtZ][iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < 6) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (xhZBins[iPtZ][0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < nXhZBins[iPtZ] && xhZBins[iPtZ][iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < 6) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < 6) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < 6) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < 6) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < 6; iPtch++) {
        for (int idPhi = 0; idPhi < GetNdPhiBins (iPtch, iCent); idPhi++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi+1) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi]));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < nPtchBins[iPtZ]; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < nPtchBins[iPtZ]; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < nXhZBins[iPtZ]; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < nXhZBins[iPtZ]; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < nPtchBins[iPtZ]; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < nPtchBins[iPtZ]; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < nXhZBins[iPtZ]; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + event_weight*(trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < nXhZBins[iPtZ]; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

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

      if (fabs (vz) > 150) continue;

      if (event_weight == 0) continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03)) continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (pTchBins[iPtZ][0] <= trkpt) {
          iPtch = 0;
          while (iPtch < nPtchBins[iPtZ] && pTchBins[iPtZ][iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < 6) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (xhZBins[iPtZ][0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < nXhZBins[iPtZ] && xhZBins[iPtZ][iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < 6) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < 6) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < 6) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < 6) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < 6; iPtch++) {
        for (int idPhi = 0; idPhi < GetNdPhiBins (iPtch, iCent); idPhi++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi+1) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi]));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < nPtchBins[iPtZ]; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < nPtchBins[iPtZ]; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < nXhZBins[iPtZ]; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < nXhZBins[iPtZ]; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < nPtchBins[iPtZ]; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < nPtchBins[iPtZ]; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < nXhZBins[iPtZ]; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + event_weight*(trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < nXhZBins[iPtZ]; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over pp tree
    cout << "Done primary pp loop." << endl;
  }

  Delete3DArray (trks_counts, 2, 6, numPhiBins+1);
  Delete3DArray (trks_weights1, 2, 6, numPhiBins+1);
  Delete3DArray (trks_weights2, 2, 6, numPhiBins+1);
  Delete2DArray (trks_counts_inPhi, 6, 40);
  Delete2DArray (trks_weights1_inPhi, 6, 40);
  Delete2DArray (trks_weights2_inPhi, 6, 40);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Generates weights between a weighted sample and a matched sample.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: GenerateEventWeights (const char* weightedSampleInFilePattern, const char* matchedSampleInFilePattern, const char* outFileName) {

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
        h_fcal_et_dist[iSample][iSpc][iPtZ]->Scale (1/h_fcal_et_dist[iSample][iSpc][iPtZ]->Integral ());

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
          h_q2_dist[iSample][iSpc][iFineCent][iPtZ]->Scale (1/h_q2_dist[iSample][iSpc][iFineCent][iPtZ]->Integral ());
          h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]->Scale (1/h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]->Integral ());
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

  if (!useCentWgts && !useQ2Wgts && !usePsi2Wgts)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory;

  cout << "Loading event weights from " << rootPath.Data () << "/" << eventWeightsFileName.c_str () << endl;
  eventWeightsFile = new TFile (Form ("%s/%s", rootPath.Data (), eventWeightsFileName.c_str ()), "read");

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (useCentWgts)  h_PbPbFCal_weights[iSpc][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, name.c_str ()));
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        if (useQ2Wgts)    h_PbPbQ2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()));
        if (usePsi2Wgts)  h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()));
      }
    }
  }

  cout << "Event weights loaded." << endl;
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

  trkEffFile = new TFile (Form ("%s/TrackingEfficiencies/%s/trackingEfficiencies_%s.root", rootPath.Data (), _effDir.Data (), is2015Conds ? (useHijingEffs ? "Hijing_15" : "15") : (useHijingEffs ? "Hijing_18" : "18")), "read");

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
        RebinSomeBins (num, maxNPtchBins, allPtchBins);
        RebinSomeBins (den, maxNPtchBins, allPtchBins);
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
          float fakeRate = 1 - h2_trk_purs[iCent]->GetBinContent (ix, iy);
          fakeRate = fakeRate + trkPurNSigma * 0.25 * fakeRate;
          h2_trk_purs[iCent]->SetBinContent (ix, iy, 1-fakeRate);
        }
      }
    }

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* num = (TH1D*) trkPurFile->Get (Form ("h_primary_reco_tracks_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkPurFile->Get (Form ("h_reco_tracks_iCent%i_iEta%i", iCent, iEta));

      if (doRebin) {
        RebinSomeBins (num, maxNPtchBins, allPtchBins);
        RebinSomeBins (den, maxNPtchBins, allPtchBins);
      }

      h_trk_purs[iCent][iEta] = (TH1D*) num->Clone (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta));
      //h_trk_purs[iCent][iEta]->Divide (den);

      for (int ix = 1; ix <= h_trk_purs[iCent][iEta]->GetNbinsX (); ix++) {
        const float passes = num->GetBinContent (ix);
        const float trials = den->GetBinContent (ix);
        h_trk_purs[iCent][iEta]->SetBinContent (ix, passes/trials);
        h_trk_purs[iCent][iEta]->SetBinError (ix, sqrt ((passes/trials)*(1-(passes/trials)) * (*den->GetSumw2 ())[ix] / pow (trials, 2)));
        //h_trk_purs[iCent][iEta]->SetBinError (ix, sqrt ((passes+1)*(passes+2) / ((trials+2)*(trials+3)) - pow (passes+1, 2) / pow (trials+2, 2)));
      }

      if (doTrackPurVar) {
        for (int ix = 1; ix <= h_trk_purs[iCent][iEta]->GetNbinsX (); ix++) {
          float fakeRate = 1 - h_trk_purs[iCent][iEta]->GetBinContent (ix);
          fakeRate = fakeRate + trkPurNSigma * 0.25 * fakeRate;
          h_trk_purs[iCent][iEta]->SetBinContent (ix, 1-fakeRate);
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
// Takes references to the q2 vector entries in each detector side and corrects them.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CorrectQ2Vector (float& q2x_a, float& q2y_a, float& q2x_c, float& q2y_c) {
  float temp = eventPlaneCalibrator.CalibrateQ2XA (q2x_a, q2y_a);
  q2y_a = eventPlaneCalibrator.CalibrateQ2YA (q2x_a, q2y_a);
  q2x_a = temp;

  temp = eventPlaneCalibrator.CalibrateQ2XC (q2x_c, q2y_c);
  q2y_c = eventPlaneCalibrator.CalibrateQ2YC (q2x_c, q2y_c);
  q2x_c = temp;
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
      cout << Form ("& %g ", h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1));
    cout << "\\\\";

    if (iCent == numCentBins-1)
      cout << " \\hline";
    cout << endl << "\t\t\t";
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotFCalDists (const bool _treatAsData) {
  const char* canvasName = "c_fcal_et";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600, 600);
    gDirectory->Add (c);
    c->cd ();
  }

  c->cd ();

  c->SetLogy ();

  h_fcal_et_reweighted->Scale (1, "width");
  if (!_treatAsData) {
    h_fcal_et_reweighted->SetLineColor (kAzure-1);
  }
  else {
    h_fcal_et_reweighted->SetLineColor (kBlack);
  }

  h_fcal_et_reweighted->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
  h_fcal_et_reweighted->GetYaxis ()->SetTitle ("dN_{evt} / d#Sigma#it{E}_{T} [GeV^{-1}]");

  h_fcal_et_reweighted->GetYaxis ()->SetRangeUser (5e1, 8e6);

  h_fcal_et_reweighted->Draw (canvasExists ? "same hist" : "hist");

  if (!canvasExists) {
    myText (0.22, 0.27, kBlack, "Z-tagged events", 0.04);
    myText (0.22, 0.21, kAzure-1, "Mixed minimum bias", 0.04);

    myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
    myText (0.22, 0.81, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  }

  c->SaveAs (Form ("%s/FCalDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Q2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotQ2Dists (const bool _treatAsData) {
  const char* canvasName = "c_q2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (3, 3);
  }

  c->cd ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    c->cd (numCentBins-iCent);
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetDrawnObjects ();
    GetMinAndMax (min, max, true);
    SetMinAndMax (min, max);

    TH1D* h = h_q2[iCent];
    if (h->GetMinimum (0) < min)
      min = h->GetMinimum (0);
    if (h->GetMaximum () > max)
      max = h->GetMaximum ();

    h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    if (!_treatAsData) {
      h->SetLineColor (colors[iCent]);
    }
    else {
      h->SetLineColor (kBlack);
    }

    h->GetXaxis ()->SetTitle ("#it{q}_{2}");
    h->GetYaxis ()->SetTitle ("Counts");

    h->Draw (!canvasExists ? "hist" : "same hist");

    myText (0.51, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

    myText (0.51, 0.74, colors[iCent], "Mixed minimum bias", 0.05);
    myText (0.51, 0.80, kBlack, "Z-tagged Data", 0.05);

    SetMinAndMax (0.5*min, 2*max);

  }

  myText (0.52, 0.31, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.52, 0.23, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Q2Dists.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Q2 weights
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotQ2Weights (PhysicsAnalysis* a) {
  const char* canvasName = "c_q2_weights";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (2,2);
  }

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    c->cd (numCentBins-iCent);

    TH1D* h = (TH1D*) h_q2_reweighted[iCent]->Clone ();

    const float hint = h->Integral ();
    const float aint = a->h_q2_reweighted[iCent]->Integral ();
    h->Divide (a->h_q2_reweighted[iCent]);

    if (hint != 0)
      h->Scale (aint / hint);

    h->GetYaxis ()->SetRangeUser (0.5, 1.5);

    h->SetMarkerColor (colors[iCent]);
    h->SetLineColor (colors[iCent]);

    h->GetXaxis ()->SetTitle ("#left|#it{q}_{2}#right|");
    h->GetYaxis ()->SetTitle ("Event Weight");

    h->Draw (!canvasExists ? "e1" : "same e1");

    myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
  }

  c->cd (1);

  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.22, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Q2Weights.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Psi2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotPsi2Dists (const bool _treatAsData) {
  const char* canvasName = "c_psi2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (3, 3);
  }

  c->cd ();

  for (short iCent = 1; iCent < numFineCentBins; iCent++) {
    c->cd (numFineCentBins-iCent);
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetDrawnObjects ();
    GetMinAndMax (min, max, true);
    SetMinAndMax (min, max);

    TH1D* h = h_psi2[iCent];
    if (h->GetMinimum (0) < min)
      min = h->GetMinimum (0);
    if (h->GetMaximum () > max)
      max = h->GetMaximum ();

    h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    if (!_treatAsData) {
      h->SetLineColor (colors[iCent]);
    }
    else {
      h->SetLineColor (kBlack);
    }

    h->GetXaxis ()->SetTitle ("#psi_{2}");
    h->GetYaxis ()->SetTitle ("Counts");

    h->Draw (!canvasExists ? "hist" : "same hist");

    myText (0.51, 0.88, kBlack, Form ("%i-%i%%", (int)fineCentCuts[iCent], (int)fineCentCuts[iCent-1]), 0.06);

    myText (0.51, 0.74, colors[iCent], "Mixed minimum bias", 0.05);
    myText (0.51, 0.80, kBlack, "Z-tagged Data", 0.05);

    SetMinAndMax (0.5*min, 2*max);

  }

  myText (0.52, 0.31, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.52, 0.23, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Psi2Dists.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Psi2 weights
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotPsi2Weights (PhysicsAnalysis* a) {
  const char* canvasName = "c_psi2_weights";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (2, 2);
  }

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    c->cd (numCentBins-iCent);

    TH1D* h = (TH1D*) h_psi2_reweighted[iCent]->Clone ();

    const float hint = h->Integral ();
    const float aint = a->h_psi2_reweighted[iCent]->Integral ();
    h->Divide (a->h_psi2_reweighted[iCent]);

    if (hint != 0)
      h->Scale (aint / hint);

    h->GetYaxis ()->SetRangeUser (0.5, 1.5);

    h->SetMarkerColor (colors[iCent]);
    h->SetLineColor (colors[iCent]);

    h->GetXaxis ()->SetTitle ("#phi_{Z} - #Psi_{2}");
    h->GetYaxis ()->SetTitle ("Event Weight");

    h->Draw (!canvasExists ? "e1" : "same e1");

    myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
  }

  c->cd (1);

  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.22, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Psi2Weights.pdf", plotPath.Data ()));

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot VZ distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotVZDists (const bool _treatAsData) {
  const char* canvasName = "c_vz";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
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

  if (!_treatAsData) {
    h_PbPb_vz_reweighted->SetLineColor (kBlue);

    h_PbPb_vz_reweighted->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
    h_PbPb_vz_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_PbPb_vz_reweighted->Draw ("same hist");
  }

  myText (0.18, 0.88, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  //myText (0.18, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.18, 0.81, kBlack, "Minimum bias", 0.04);

  c->cd (2);

  if (!_treatAsData) {
    h_pp_vz_reweighted->SetLineColor (kBlue);

    h_pp_vz_reweighted->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
    h_pp_vz_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_pp_vz_reweighted->Draw ("same hist");

    myText (0.75, 0.88, kBlack, "Unweighted", 0.04);
    myText (0.75, 0.81, kBlue, "Reweighted", 0.04);
  }

  myText (0.18, 0.88, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  //myText (0.18, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.18, 0.81, kBlack, "Minimum bias", 0.04);

  c->SaveAs (Form ("%s/VZDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot VZ distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotNchDists (const bool _treatAsData) {
  const char* canvasName = "c_nch";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    c->cd ();
  }

  gPad->SetLogy ();

  if (_treatAsData)
    h_pp_nch->SetLineColor (kBlack);
  else
    h_pp_nch->SetLineColor (kGray+1);

  h_pp_nch->Scale (1/h_pp_nch->Integral ());
  h_pp_nch->GetXaxis ()->SetTitle ("N_{ch}");
  h_pp_nch->GetYaxis ()->SetTitle ("dN_{evt} / N_{evt}");

  h_pp_nch->Draw (canvasExists ? "same hist" : "hist");

  if (!_treatAsData) {
    h_pp_nch_reweighted->Scale (1/h_pp_nch_reweighted->Integral ());
    h_pp_nch_reweighted->SetLineColor (kRed+1);

    h_pp_nch_reweighted->GetXaxis ()->SetTitle ("N_{ch}");
    h_pp_nch_reweighted->GetYaxis ()->SetTitle ("dN_{evt} / N_{evt");

    h_pp_nch_reweighted->Draw ("same hist");
  }

  myText (0.68, 0.88, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  myText (0.68, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.68, 0.74, kGray+1, "Minimum bias", 0.04);
  myText (0.68, 0.67, kRed+1, "MB Reweighted", 0.04);

  c->SaveAs (Form ("%s/NchDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - pTch 2d projections
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotCorrelations (const short pSpc, const short pPtZ, const bool _subBkg) {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue;

      const char* canvasName = Form ("c_deltaPhi_pTch_iPtZ%i_%s", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 900, 300*numCentBins);
        gDirectory->Add (c);
        c->cd ();
        c->Divide (2, numCentBins);
      }

      //for (short iCent = 0; iCent < 1; iCent++) {
      //  c->cd ();
      //  gPad->SetLogy ();
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (bool subBkg : {false, true}) {
          if (subBkg && !_subBkg) continue;

          c->cd ((2*iCent)+1 + (int)subBkg);

          gPad->SetTopMargin (0.01);
          gPad->SetBottomMargin (0.12);
          gPad->SetRightMargin (0.01);
          gPad->SetLeftMargin (0.12);

          if (!canvasExists) {
            GetDrawnObjects ();

            double min = -1.5, max = 2.7;
            if (!subBkg) {
              GetMinAndMax (min, max, true);
              for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++) {
                TH1D* h = (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtch][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]);
                min = fmin (min, h->GetMinimum ());
                max = fmax (max, h->GetMaximum ());
              } // end loop over iPtch
              if (max != 2.7) max = max <= 0 ? 1 : 1.2*max;
              SetMinAndMax (min, max);
            }

            TH1D* h = new TH1D (Form ("_hplot_iCent%i_subBkg%i", iCent, subBkg ? 1 : 0), "", 1, -pi/2, 3*pi/2);

            h->GetYaxis ()->SetRangeUser (min, max);

            h->SetLineWidth (0);
            h->SetMarkerSize (0);

            h->GetXaxis ()->SetTitle ("#Delta#phi");
            h->GetYaxis ()->SetTitle ("dY / d#Delta#phi");

            h->GetXaxis ()->SetTitleOffset (0.6);
            h->GetYaxis ()->SetTitleOffset (0.8);
            h->GetXaxis ()->SetTitleSize (0.08);
            h->GetYaxis ()->SetTitleSize (0.07);
            h->GetXaxis ()->SetLabelSize (0.06);
            h->GetYaxis ()->SetLabelSize (0.06);

            h->DrawCopy ();
            delete h;

            TLine* line1 = new TLine (3*pi/4, min, 3*pi/4, (iCent == 0 ? 0.8:1.)*max);
            TLine* line2 = new TLine (5*pi/4, min, 5*pi/4, (iCent == 0 ? 0.8:1.)*max);

            line1->SetLineStyle (2);
            line2->SetLineStyle (2);
            line1->SetLineWidth (2);
            line2->SetLineWidth (2);
            line1->SetLineColor (kBlack);
            line2->SetLineColor (kBlack);

            line1->Draw ("same");
            line2->Draw ("same");
 
            for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++)
              LabelCorrelations (iPtZ, iPtch, iCent, subBkg);
          }

          if (plotFill) {
            for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++) {
              TH1D* h = (TH1D*) (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtch][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent])->Clone ();

              //h->GetYaxis ()->SetRangeUser (min, max);

              h->SetLineColorAlpha (fillColors[iPtch], 0.8);
              h->SetLineWidth (4);
              h->SetMarkerSize (0);

              //h->GetXaxis ()->SetTitle ("#Delta#phi");
              //h->GetYaxis ()->SetTitle ("dY / d#Delta#phi");

              //h->GetXaxis ()->SetTitleOffset (0.6);
              //h->GetYaxis ()->SetTitleOffset (0.8);
              //h->GetXaxis ()->SetTitleSize (0.08);
              //h->GetYaxis ()->SetTitleSize (0.07);
              //h->GetXaxis ()->SetLabelSize (0.06);
              //h->GetYaxis ()->SetLabelSize (0.06);

              //h->DrawCopy ((subBkg || !canvasExists) && iPtch == 0 ? "hist" : "same hist");
              h->DrawCopy ("same hist");
              //h->SetLineWidth (1);
              //h->Draw ("hist same");

              //TLine* line = new TLine (h->GetBinLowEdge (1), 0, h->GetBinLowEdge (h->GetNbinsX ()), 0);
              //line->Draw ("same");
            } // end loop over iPtch
            gPad->RedrawAxis ();
          } else {
            for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++) {
              TGAE* g = GetTGAE (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtch][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]);
              ResetXErrors (g);

              //const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
              const Style_t markerStyle = markerStyles[iPtch];
              //g->GetYaxis ()->SetRangeUser (min, max);

              g->SetMarkerStyle (markerStyle);
              g->SetLineColor (colors[iPtch]);
              g->SetMarkerColor (colors[iPtch]);

              //g->GetXaxis ()->SetTitle ("#Delta#phi");
              //g->GetYaxis ()->SetTitle ("dY / d#Delta#phi");

              //g->GetXaxis ()->SetTitleOffset (0.6);
              //g->GetYaxis ()->SetTitleOffset (0.8);
              //g->GetXaxis ()->SetTitleSize (0.08);
              //g->GetYaxis ()->SetTitleSize (0.07);
              //g->GetXaxis ()->SetLabelSize (0.06);
              //g->GetYaxis ()->SetLabelSize (0.06);

              //g->Draw ((subBkg || !canvasExists) && iPtch == 0 ? "AP" : "P");
              g->Draw ("P");

              //if (_subBkg) LabelCorrelations (iPtZ, iPtch, iCent, subBkg);

              //TLine* line1 = new TLine (phiLowBins[1], min, phiLowBins[1], max);
              //TLine* line2 = new TLine (2*pi - phiLowBins[1], min, 2*pi - phiLowBins[1], max);

              //line1->SetLineStyle (2);
              //line2->SetLineStyle (2);
              //line1->SetLineWidth (2);
              //line2->SetLineWidth (2);
              //line1->SetLineColor (kBlack);
              //line2->SetLineColor (kBlack);

              //line1->Draw ("same");
              //line2->Draw ("same");
            } // end loop over iPtch
          }
        }
      } // end loop over iCent

      c->SaveAs (Form ("%s/DeltaPhiCorrelations/deltaPhi_pTch_iPtZ%i_%s.pdf", plotPath.Data (), iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track dPhi distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelCorrelations (const short iPtZ, const short iPtch, const short iCent, const bool subBkg) {
  if (iCent == 0) {
    if (iPtch == 0) {
      if (!subBkg) {
        //myText (0.2, 0.91, kBlack, "#bf{#it{ATLAS}} Simulation", 0.07);
        myText (0.2, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.077);
        myText (0.2, 0.79, kBlack, "#it{pp}, 5.02 TeV", 0.07);
        if (iPtZ == nPtZBins-1) myText (0.2, 0.70, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.07);
        else                    myText (0.2, 0.70, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.07);
        //myMarkerTextNoLine (0.24, 0.64, kBlack, kFullCircle, "Z-tagged Data", 1.25, 0.07);
        //myOnlyBoxText (0.24, 0.55, 1.2, fillColors[0], kBlack, 1, "Minimum Bias", 0.07, 1001, 1);
        myText             (0.61, 0.90, kBlack, "Before subtraction", 0.07);
      }
      else {
        myText             (0.17, 0.88, kBlack, "MB", 0.07);
        myText             (0.26, 0.88, kBlack, "#it{Z}-tagged", 0.07);
        myText             (0.63, 0.90, kBlack, "After subtraction", 0.07);
      }
    }
    if (subBkg) {
      const float pt_lo = pTchBins[iPtZ][iPtch];
      const float pt_hi = pTchBins[iPtZ][iPtch+1];
      //if (iPtch == 0)
      //  myText (0.3, 0.93, kBlack, "#it{p}_{T}^{ ch} [GeV]", 0.04);
      myOnlyBoxText      (0.24, 0.81-0.075*(iPtch), 1.2, fillColors[iPtch], kBlack, 1, "", 0.07, 1001, 0.8);
      myMarkerTextNoLine (0.36, 0.81-0.075*(iPtch), colors[iPtch], markerStyles[iPtch], "", 1.4, 0.06);
      myText             (0.40, 0.79-0.075*(iPtch),  kBlack, Form ("%.0f < #it{p}_{T}^{ch} < %.0f GeV", pt_lo, pt_hi), 0.06);
    }
  }
  else if (!subBkg && iPtch == 0 && iCent != 0) {
    myText (0.2, 0.90, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.07);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies () {
  if (!effsLoaded) LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  //TGraphErrors* g0[numCentBins];
  //TGraphErrors* g1[numCentBins];
  //TGraphErrors* g2[numCentBins];
  //TGraphErrors* g3[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
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

      if (iCent == 0) {
        myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        myMarkerTextNoLine (0.5, 0.50-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
      }
      else myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficiencies.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies with comparison to some alternative
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficienciesComparison (PhysicsAnalysis* a) {
  if (!effsLoaded) LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs_comp";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 800);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 2);
  }
  c->cd ();

  //TGraphErrors* g0[numCentBins];
  //TGraphErrors* g1[numCentBins];
  //TGraphErrors* g2[numCentBins];
  //TGraphErrors* g3[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
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

      if (iEta == 0) {
        if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
        if (iCent == 1) {
          //myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "HILoose tracks", 1.2, 0.06);
          //myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "HITight tracks", 1.2, 0.06);
          myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "Inclusive hadrons", 1.2, 0.06);
          myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "Pions only", 1.2, 0.06);
        }
      }
      if (iCent == 0) myMarkerTextNoLine (0.5, 0.34-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");

    if (!a) continue;
    else    a->LoadTrackingEfficiencies ();

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

      eff->SetTitle (";#it{p}_{T} [GeV];Pions / Inclusive hadrons");
      //eff->SetTitle (";#it{p}_{T} [GeV];HITight / HILoose");
      eff->GetXaxis ()->SetRangeUser (0.5, 60);
      eff->GetYaxis ()->SetRangeUser (0.89, 1.11);
      //eff->GetYaxis ()->SetRangeUser (0.7, 1.11);

      eff->GetXaxis ()->SetTitleSize (0.07);
      eff->GetYaxis ()->SetTitleSize (0.07);
      eff->GetXaxis ()->SetTitleOffset (0.7 * eff->GetXaxis ()->GetTitleOffset ());
      eff->GetYaxis ()->SetTitleOffset (0.7 * eff->GetYaxis ()->GetTitleOffset ());

      eff->GetYaxis ()->CenterTitle ();

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (!canvasExists && iEta == 0 ? "AP" : "P");
    }

    l->Draw ("same");
  }

  if (a) {
    bool temp = a->useAltMarker;
    a->useAltMarker = true;
    a->PlotTrackingEfficienciesComparison ();
    a->useAltMarker = temp;
  }
  else c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficienciesComparison.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies2D () {
  if (!effsLoaded) LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
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

    if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
    else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficiencies2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities () {
  if (!pursLoaded) LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
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

      if (iCent == 0) {
        myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        myMarkerTextNoLine (0.5, 0.50-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
      }
      else myText (0.22, 0.86, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");
  }

  c->SaveAs (Form ("%s/TrackingPurities/TrackingPurities.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities with comparison to some alternative
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPuritiesComparison (PhysicsAnalysis* a) {
  if (!pursLoaded) LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs_comp";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 800);
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

      if (iEta == 0) {
        if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
        if (iCent == 1) {
          myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "HILoose tracks", 1.2, 0.06);
          myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "HITight tracks", 1.2, 0.06);
        }
      }
      if (iCent == 0) myMarkerTextNoLine (0.5, 0.34-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");

    if (!a) continue;
    else    a->LoadTrackingPurities ();

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
    a->PlotTrackingPuritiesComparison ();
    a->useAltMarker = temp;
  }
  else c->SaveAs (Form ("%s/TrackingPurities/TrackingPuritiesComparison.pdf", plotPath.Data ()));
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
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
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

    if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
    else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
  }

  c->SaveAs (Form ("%s/TrackingPurities/TrackingPurities2D.pdf", plotPath.Data ()));
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
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      h_trk_pt_ptz[2][iPtZ][iCent]->Reset ();
      h_trk_pt_ptz_sub[2][iPtZ][iCent]        = new TH1D (Form ("h_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]);
      h_trk_pt_ptz_sub[2][iPtZ][iCent]->Sumw2 ();
      h_trk_pt_ptz_sig_to_bkg[2][iPtZ][iCent] = new TH1D (Form ("h_trk_pt_ptz_sigToBkg_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]);
      h_trk_pt_ptz_sig_to_bkg[2][iPtZ][iCent]->Sumw2 ();
      h_trk_xhz_ptz[2][iPtZ][iCent]->Reset ();
      h_trk_xhz_ptz_sub[2][iPtZ][iCent]         = new TH1D (Form ("h_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]);
      h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Sumw2 ();
      h_trk_xhz_ptz_sig_to_bkg[2][iPtZ][iCent]  = new TH1D (Form ("h_trk_xhz_ptz_sigToBkg_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]);
      h_trk_xhz_ptz_sig_to_bkg[2][iPtZ][iCent]->Sumw2 ();
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        h_trk_pt_dphi[2][iPtZ][iPhi][iCent]->Reset ();
        h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]         = new TH1D (Form ("h_trk_pt_dphi_sub_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]);
        h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]->Sumw2 ();
        h_trk_pt_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]  = new TH1D (Form ("h_trk_pt_dphi_sigToBkg_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nPtchBins[iPtZ], pTchBins[iPtZ]);
        h_trk_pt_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Sumw2 ();
        h_trk_xhz_dphi[2][iPtZ][iPhi][iCent]->Reset ();
        h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]        = new TH1D (Form ("h_trk_xhz_dphi_sub_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]);
        h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]->Sumw2 ();
        h_trk_xhz_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_xhz_dphi_sigToBkg_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]);
        h_trk_xhz_dphi_sig_to_bkg[2][iPtZ][iPhi][iCent]->Sumw2 ();
      } // end loop over iPhi
      for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
        h_trk_dphi[2][iPtZ][iPtch][iCent]->Reset ();
        h_trk_dphi_sub[2][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_sub_comb_iPtZ%i_iPtch%i_iCent%i_%s", iPtZ, iPtch, iCent, name.c_str ()), "", 80, -pi/2, 3*pi/2);
        h_trk_dphi_sub[2][iPtZ][iPtch][iCent]->Sumw2 ();
      } // end loop over iPtch
    } // end loop over iPtZ
  } // end loop over iCent
  
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        //******** Do subtraction of integrated dPhi plot ********//
        TH1D* h = (TH1D*) h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] = h;
        TH1D* sub = nullptr;
        if (a != nullptr) {
          sub = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
          AddNoErrors (h, sub, -1);
        }

        h = (TH1D*) h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = h;
        if (a != nullptr && sub != nullptr) {
          AddNoErrors (h, sub, -1);
          h->Divide (sub);
          //MultiplyNoErrors (h, sub, -1);
        }

        h = (TH1D*) h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent] = h;
        if (a != nullptr) {
          sub = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
          AddNoErrors (h, sub, -1);
        }

        h = (TH1D*) h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = h;
        if (a != nullptr && sub != nullptr) {
          AddNoErrors (h, sub, -1);
          h->Divide (sub);
          //MultiplyNoErrors (h, sub, -1);
        }


        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {

          //******** Do subtraction of pT ********//
          h = (TH1D*) h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] = h;
          if (a != nullptr) {
            sub = a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
            AddNoErrors (h, sub, -1);
          }

          h = (TH1D*) h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = h;
          if (a != nullptr && sub != nullptr) {
            AddNoErrors (h, sub, -1);
            h->Divide (sub);
            //MultiplyNoErrors (h, sub, -1);
          }


          //******** Do subtraction of z_h ********//
          h = new TH1D (Form ("h_trk_xhz_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]);
          h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent] = h;
          h->Sumw2 ();
          h->Add (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (a != nullptr) {
            sub = a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
            AddNoErrors (h, sub, -1);
          }

          h = new TH1D (Form ("h_trk_xhz_dphi_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXhZBins[iPtZ], xhZBins[iPtZ]);
          h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = h;
          h->Sumw2 ();
          h->Add (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          if (a != nullptr && sub != nullptr) {
            AddNoErrors (h, sub, -1);
            h->Divide (sub);
            //MultiplyNoErrors (h, sub, -1);
          }
        } // end loop over iPhi


        for (int iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          //******** Do background subtraction of phi distributions ********//
          TH1D* h = (TH1D*) h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
          h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent] = h;
          if (a != nullptr) {
            TH1D* sub = a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
            //while (sub->GetNbinsX () > h->GetNbinsX ())
            //  sub->Rebin (2);
            AddNoErrors (h, sub, -1);
          }

        } // end loop over iPtch

      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

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

        f = f_trk_xhz_ptz_binMigration[iSpc][iPtZ][iCent] = (TF1*) f_binMigrationFile->Get (Form ("tf1_%s_xh_ZPT%i_%s", spc, iPtZ-2, cent));
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
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc
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
          const float bkgCountsOverObsCounts = (a->h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2)) / (h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2));
          if (iCent == 0 && iPhi == 0 && iSpc == 2 && iPtZ == 2) {
            cout << bkgCountsOverObsCounts << endl;
            cout << "2nd to last bin before: " << h->GetBinContent (6) << endl;
          }

          h->Add (a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent], -bkgCountsOverObsCounts);
          h->Scale (1/ (1-bkgCountsOverObsCounts));
          //if (iCent == 0 && iPhi == 0 && iSpc == 2 && iPtZ == 2) {
          //  cout << "2nd to last bin after:  " << h->GetBinContent (6) << endl;
          //}

          h = h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          h->Add (a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], -bkgCountsOverObsCounts);
          h->Scale (1/ (1-bkgCountsOverObsCounts));
        } // end loop over iPhi
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

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
        //  h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1/relVar[iSpc][iPtZ][iPhi][iCent]);
        //  h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1/relVar[iSpc][iPtZ][iPhi][iCent]);
        //} // end loop over iPhi

        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1/relVar[iSpc][iPtZ][numPhiBins][iCent]);
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1/relVar[iSpc][iPtZ][numPhiBins][iCent]);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
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
        } // end loop over iPhi

        AddStatVar (h_trk_pt_ptz[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent], upVar, nSigma);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTch distributions binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotUnweightedTrkYields (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_UnweightedTrkYields_%s_dPhi_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 400*numCentBins, 400);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* padName = Form ("p_UnweightedTrkYields_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, iCent, spc);

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
          TH1D* h = (useTrkPt ? h_trk_pt_dphi_raw : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over iPhi
        min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
        max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = 0; iPhi < 1; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = (useTrkPt ? h_trk_pt_dphi_raw : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
            h->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}, #Delta#phi)" : "Y (#it{x}_{hZ}, #Delta#phi)");

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
            TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_raw : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent]);
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
              //g->SetLineWidth (1);
              //g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
            g->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}, #Delta#phi)" : "Y (#it{x}_{hZ}, #Delta#phi)");

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
          } // end loop over iPhi
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
      } // end loop over iCent
      
      c->SaveAs (Form ("%s/UnweightedTrkYields/allYields_pTch_iPtZ%i_%s.pdf", plotPath.Data (), iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllYields_dPhi (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_AllYields_%s_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 395*numCentBins, 685);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* topPadName = Form ("p_AllYields_top_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, iCent, spc);
        const char* bottomPadName = Form ("p_AllYields_bottom_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, iCent, spc);

        TPad* topPad = nullptr, *bottomPad = nullptr;
        if (!canvasExists) {
          topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
          bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

          gDirectory->Add (topPad);
          gDirectory->Add (bottomPad);

          topPad->SetTopMargin (0.04);
          topPad->SetBottomMargin (0);
          topPad->SetLeftMargin (0.20);
          topPad->SetRightMargin (0.06);
          bottomPad->SetTopMargin (0);
          bottomPad->SetBottomMargin (0.20);
          bottomPad->SetLeftMargin (0.20);
          bottomPad->SetRightMargin (0.06);
          topPad->Draw ();
          bottomPad->Draw ();
        }
        else {
          topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
          bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
        }

        topPad->cd ();
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

        if (!canvasExists) {
          TH1D* h = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);
          useTrkPt ? h->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->SetLineWidth (0);
          h->SetMarkerSize (0);

          h->GetXaxis ()->SetMoreLogLabels ();

          h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
          h->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ} d#Delta#phi");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (1.5 * h->GetYaxis ()->GetTitleOffset ());

          h->DrawCopy ("hist ][");
          delete h;
        }

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = (useTrkPt ? h_trk_pt_dphi : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent];

            //h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            //h->SetMarkerSize (0);
            h->SetLineColor (colors[iPhi]);
            h->SetLineStyle (iPhi+1);
            h->SetLineWidth (1);

            h->Draw ("hist ][ same");
          }
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent]);
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
              //g->SetLineWidth (1);
              //g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (!plotAsSystematic) g->Draw ("P");
            else {
              ((TGAE*)g->Clone ())->Draw ("5P");
              g->Draw ("2P");
            }
          } // end loop over iPhi
        }

        if (!canvasExists) {
          if (iCent == 0) myText (0.25, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
          else            myText (0.25, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

          if (iCent == 2 && iSpc == 0)      myText (0.72, 0.87, kBlack, "Z #rightarrow #it{ee}", 0.06);
          else if (iCent == 2 && iSpc == 1) myText (0.72, 0.87, kBlack, "Z #rightarrow #it{#mu#mu}", 0.06);

          if (iCent == 0) myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
          else if (iCent == 1) {
            if (iPtZ == nPtZBins-1) myText (0.485, 0.87, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.06);
            else                    myText (0.485, 0.87, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.06);
          }
          else if (iCent == numCentBins-1) {
            myText (0.43, 0.90, kBlack, "MB", 0.06);
            myText (0.52, 0.90, kBlack, "Z-tag", 0.06);
          }
          if (iCent == numCentBins-1) {
            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              myLineText (0.52, 0.85-0.06*(iPhi-1), colors[iPhi], iPhi+1, "", 2.0, 0.054) ;
              myMarkerTextNoLine (0.59, 0.85-0.06*(iPhi-1), colors[iPhi], kFullCircle, "", 1.5, 0.054); // for plotting data vs bkg.

              const char* lo = GetPiString (phiLowBins[iPhi]);
              const char* hi = GetPiString (phiHighBins[iPhi]);

              myText (0.61, 0.84-0.06*(iPhi-1), kBlack, Form ("(%s, %s)", lo, hi), 0.054);
            } // end loop over iPhi
          }
        }


        bottomPad->cd ();
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

        if (!canvasExists) {
          TH1D* h = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);
          useTrkPt ? h->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->SetLineWidth (0);
          h->SetMarkerSize (0);

          h->GetXaxis ()->SetMoreLogLabels ();

          h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
          h->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ} d#Delta#phi");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (1.5 * h->GetYaxis ()->GetTitleOffset ());

          h->DrawCopy ("hist ][");
          delete h;
        }

        if (!plotSignal)
          continue;

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* h = (useTrkPt ? h_trk_pt_dphi_sub : h_trk_xhz_dphi_sub)[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            h->DrawCopy ("bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          } // end loop over iPhi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));

            TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_sub : h_trk_xhz_dphi_sub)[iSpc][iPtZ][iPhi][iCent]);
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
              //g->SetLineWidth (1);
              //g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            if (!plotAsSystematic) g->Draw ("P");
            else {
              ((TGAE*)g->Clone ())->Draw ("5P");
              g->Draw ("2P");
            }
          } // end loop over iPhi
        }
      } // end loop over iCent
      
      c->SaveAs (Form ("%s/TrkYields/allYields_%s_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllYields_dPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_AllYields_%s_dPtZ_%s", useTrkPt ? "pTch" : "xhZ", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 395*numCentBins, 685);
      gDirectory->Add (c);
    }

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      c->cd ();

      const char* topPadName = Form ("p_AllYields_top_%s_dPtZ_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iCent, spc);
      const char* bottomPadName = Form ("p_AllYields_bottom_%s_dPtZ_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iCent, spc);

      TPad* topPad = nullptr, *bottomPad = nullptr;
      if (!canvasExists) {
        topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
        bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

        gDirectory->Add (topPad);
        gDirectory->Add (bottomPad);

        topPad->SetTopMargin (0.04);
        topPad->SetBottomMargin (0);
        topPad->SetLeftMargin (0.20);
        topPad->SetRightMargin (0.06);
        bottomPad->SetTopMargin (0);
        bottomPad->SetBottomMargin (0.20);
        bottomPad->SetLeftMargin (0.20);
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
      gPad->SetLogx ();
      gPad->SetLogy ();

      double min = 1e30, max = 0;
      if (useTrkPt) {
        if (iCent == 0)       { min = 1e-5; max = 4e1; }
        else if (iCent == 1)  { min = 1e-5; max = 4e2; }
        else if (iCent == 2)  { min = 5e-5; max = 2e3; }
        else if (iCent == 3)  { min = 8e-5; max = 6e3; }
      }
      else {
        if (iCent == 0)       { min = 2e-5; max = 2e3; }
        else if (iCent == 1)  { min = 1e-3; max = 5e4; }
        else if (iCent == 2)  { min = 1e-3; max = 9e4; }
        else if (iCent == 3)  { min = 1e-3; max = 9e5; }
      }

      if (!canvasExists) {
        TH1D* h = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);
        useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
        h->GetYaxis ()->SetRangeUser (min, max);

        h->SetLineWidth (0);
        h->SetMarkerSize (0);

        h->GetXaxis ()->SetMoreLogLabels ();

        h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
        h->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

        h->GetXaxis ()->SetTitleFont (43);
        h->GetXaxis ()->SetTitleSize (axisTextSize);
        h->GetXaxis ()->SetLabelFont (43);
        h->GetXaxis ()->SetLabelSize (axisTextSize);

        h->GetYaxis ()->SetTitleFont (43);
        h->GetYaxis ()->SetTitleSize (axisTextSize);
        h->GetYaxis ()->SetLabelFont (43);
        h->GetYaxis ()->SetLabelSize (axisTextSize);

        h->GetYaxis ()->SetTitleOffset (1.5 * h->GetYaxis ()->GetTitleOffset ());
        h->DrawCopy ("hist ][");
        delete h;

        if (iCent == 0) {
          myText (0.25, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
          myText (0.25, 0.14, kBlack, "Before subtraction", 0.060);
        }
        else myText (0.25, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

        if (iCent == 2 && iSpc == 0)      myText (0.72, 0.87, kBlack, "Z#rightarrow #it{ee}", 0.06);
        else if (iCent == 2 && iSpc == 1) myText (0.72, 0.87, kBlack, "Z#rightarrow #it{#mu#mu}", 0.06);

        if (iCent == 0)       myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
        else if (iCent == 1)  myText (0.485, 0.87, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.06);
        else if (iCent == numCentBins-1) {
          myText (0.38, 0.90, kBlack, "MB", 0.06);
          myText (0.47, 0.90, kBlack, "Z-tag", 0.06);
        }

        if (iCent == numCentBins-1) {
          for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
            myLineText (0.47, 0.85-0.06*(iPtZ-2), colors[iPtZ-1], iPtZ, "", 2.0, 0.054) ;
            myMarkerTextNoLine (0.54, 0.85-0.06*(iPtZ-2), colors[iPtZ-1], markerStyles[iPtZ-2], "", 1.5, 0.054); // for plotting data vs bkg.

            if (iPtZ == nPtZBins-1) myText (0.56, 0.84-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.054);
            else                    myText (0.56, 0.84-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
          }
        }
      }

      if (plotFill) {
        for (int iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
          TH1D* h = (useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[iSpc][iPtZ][iCent];

          //h->SetFillColorAlpha (fillColors[iPtZ-2], fillAlpha);
          //h->SetMarkerSize (0);
          h->SetLineColor (colors[iPtZ-1]);
          h->SetLineStyle (iPtZ);
          h->SetLineWidth (1);

          h->Draw ("hist ][ same");
        }
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = markerStyles[iPtZ-2];
          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[iSpc][iPtZ][iCent]);
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
            g->SetLineWidth (0);
            //g->SetLineWidth (1);
            //g->SetLineColor (colors[iPtZ-1]);
            g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
          }

          if (!plotAsSystematic) g->Draw ("P");
          else {
            ((TGAE*)g->Clone ())->Draw ("5P");
            g->Draw ("2P");
          }
        } // end loop over iPtZ
      }


      bottomPad->cd ();
      GetDrawnObjects ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      if (useTrkPt) { min = 1e-3; max = 2e1; }
      else          { min = 1e-2; max = 8e2; }

      if (!canvasExists) {
        TH1D* h = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);
        useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
        h->GetYaxis ()->SetRangeUser (min, max);

        h->SetLineWidth (0);
        h->SetMarkerSize (0);

        h->GetXaxis ()->SetMoreLogLabels ();

        h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
        h->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

        h->GetXaxis ()->SetTitleFont (43);
        h->GetXaxis ()->SetTitleSize (axisTextSize);
        h->GetXaxis ()->SetLabelFont (43);
        h->GetXaxis ()->SetLabelSize (axisTextSize);

        h->GetYaxis ()->SetTitleFont (43);
        h->GetYaxis ()->SetTitleSize (axisTextSize);
        h->GetYaxis ()->SetLabelFont (43);
        h->GetYaxis ()->SetLabelSize (axisTextSize);

        h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
        h->GetYaxis ()->SetTitleOffset (1.5 * h->GetYaxis ()->GetTitleOffset ());
        h->DrawCopy ("hist ][");
        delete h;

        if (iCent == 0) myText (0.25, 0.24, kBlack, "After subtraction", 0.060/padRatio);
      }

      if (!plotSignal)
        continue;

      if (plotFill) {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][iCent];

          if (!h) continue;

          h->SetFillColorAlpha (fillColors[iPtZ-1], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          h->SetMarkerStyle (kFullCircle);

          h->DrawCopy ("bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        } // end loop over iPtZ
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          const Style_t markerStyle = markerStyles[iPtZ-2];

          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][iCent]);
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
            g->SetLineWidth (0);
            //g->SetLineWidth (1);
            //g->SetLineColor (colors[iPtZ-1]);
            g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
          }

          if (!plotAsSystematic) g->Draw ("P");
          else {
            ((TGAE*)g->Clone ())->Draw ("5P");
            g->Draw ("2P");
          }
        } // end loop over iPtZ
      }
    } // end loop over iCent
    
    c->SaveAs (Form ("%s/TrkYields/allYields_%s_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", spc));
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllYields_dPtZ_SpcComp (const bool useTrkPt, const bool plotAsSystematic) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  const char* canvasName = Form ("c_AllYields_%s_dPtZ", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 395*numCentBins, 685);
    gDirectory->Add (c);
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* topPadName = Form ("p_AllYields_top_%s_dPtZ_iCent%i", useTrkPt ? "pTch" : "xhZ", iCent);
    const char* bottomPadName = Form ("p_AllYields_bottom_%s_dPtZ_iCent%i", useTrkPt ? "pTch" : "xhZ", iCent);

    TPad* topPad = nullptr, *bottomPad = nullptr;
    if (!canvasExists) {
      topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
      bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

      gDirectory->Add (topPad);
      gDirectory->Add (bottomPad);

      topPad->SetTopMargin (0.04);
      topPad->SetBottomMargin (0);
      topPad->SetLeftMargin (0.20);
      topPad->SetRightMargin (0.06);
      bottomPad->SetTopMargin (0);
      bottomPad->SetBottomMargin (0.20);
      bottomPad->SetLeftMargin (0.20);
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
    gPad->SetLogx ();
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    if (useTrkPt) {
      if (iCent == 0)       { min = 1e-5; max = 4e1; }
      else if (iCent == 1)  { min = 1e-5; max = 4e2; }
      else if (iCent == 2)  { min = 5e-5; max = 2e3; }
      else if (iCent == 3)  { min = 8e-5; max = 6e3; }
    }
    else {
      if (iCent == 0)       { min = 2e-5; max = 2e3; }
      else if (iCent == 1)  { min = 1e-3; max = 5e4; }
      else if (iCent == 2)  { min = 1e-3; max = 9e4; }
      else if (iCent == 3)  { min = 1e-3; max = 9e5; }
    }

    if (!canvasExists) {
      TH1D* h = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);
      useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
      h->GetYaxis ()->SetRangeUser (min, max);

      h->SetLineWidth (0);
      h->SetMarkerSize (0);

      h->GetXaxis ()->SetMoreLogLabels ();

      h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      h->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetTitleSize (axisTextSize);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetXaxis ()->SetLabelSize (axisTextSize);

      h->GetYaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleSize (axisTextSize);
      h->GetYaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelSize (axisTextSize);

      h->GetYaxis ()->SetTitleOffset (1.5 * h->GetYaxis ()->GetTitleOffset ());
      h->DrawCopy ("hist ][");
      delete h;

      if (iCent == 0) {
        myText (0.25, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
        myText (0.25, 0.14, kBlack, "Before subtraction", 0.060);
      }
      else myText (0.25, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

      if (iCent == 0)       myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      else if (iCent == 1)  myText (0.485, 0.87, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.06);
      else if (iCent == numCentBins - 1) {
        for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          if (iPtZ == 2) {
            myText (0.37, 0.87, kBlack, "#it{ee}", 0.06);
            myText (0.47, 0.87, kBlack, "#it{#mu#mu}", 0.06);
          }
          myMarkerTextNoLine (0.42, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kFullCircle, "", 1.5, 0.054); // for plotting MC reco vs truth
          myMarkerTextNoLine (0.52, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kOpenCircle, "", 1.5, 0.054);

          if (iPtZ == nPtZBins-1) myText (0.56, 0.85-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.054);
          else myText (0.56, 0.85-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
        }
      }
    }

    for (short iSpc = 0; iSpc < 2; iSpc++) {
      for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        const Style_t markerStyle = (iSpc == 1 ? kOpenCircle : kFullCircle);

        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[iSpc][iPtZ][iCent]);
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
          g->SetLineWidth (0);
          //g->SetLineWidth (1);
          //g->SetLineColor (colors[iPtZ-1]);
          g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
        }

        if (!plotAsSystematic) g->Draw ("P");
        else {
          ((TGAE*)g->Clone ())->Draw ("5P");
          g->Draw ("2P");
        }
      } // end loop over iPtZ
    } // end loop over iSpc

    bottomPad->cd ();
    GetDrawnObjects ();
    //plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
    gPad->SetLogx ();

    if (!canvasExists) {
      TH1D* h = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);
      useTrkPt ? h->GetXaxis ()->SetLimits (trk_min_pt, trk_max_pt) : h->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
      h->GetYaxis ()->SetRangeUser (0, 2);

      h->SetLineWidth (0);
      h->SetMarkerSize (0);

      h->GetXaxis ()->SetMoreLogLabels ();

      h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      h->GetYaxis ()->SetTitle ("Electrons / Muons");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetTitleSize (axisTextSize);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetXaxis ()->SetLabelSize (axisTextSize);

      h->GetYaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleSize (axisTextSize);
      h->GetYaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelSize (axisTextSize);

      h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
      h->GetYaxis ()->SetTitleOffset (1.5 * h->GetYaxis ()->GetTitleOffset ());
      h->DrawCopy ("hist ][");
      delete h;
    }

    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      TH1D* h = (useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[0][iPtZ][iCent];
      h->Divide ((useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[1][iPtZ][iCent]);

      TGAE* g = GetTGAE (h);
      delete h;
      RecenterGraph (g);

      if (!plotAsSystematic) {
        ResetXErrors (g);
        g->SetMarkerColor (colors[iPtZ-1]);
        g->SetLineColor (colors[iPtZ-1]);
        g->SetMarkerSize (1);
        g->SetLineWidth (2);
      } else {
        g->SetMarkerSize (0); 
        g->SetLineWidth (0);
        //g->SetLineWidth (1);
        //g->SetLineColor (colors[iPtZ-1]);
        g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
      }

      if (!plotAsSystematic) g->Draw ("P");
      else {
        ((TGAE*)g->Clone ())->Draw ("5P");
        g->Draw ("2P");
      }
    } // end loop over iPtZ
  } // end loop over iCent
  
  c->SaveAs (Form ("%s/TrkYields/allYields_%s_SpcComp.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots all track yield covariance matrices.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllCovMatrices () {
  for (short iSpc : {0, 1})
    for (bool iPt : {true, false})
      for (short iCent : {0, 1, 2, 3})
        for (short iPtZ : {2, 3, 4})
          PlotCovMatrix (iPt, iSpc, iPtZ, iCent);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the track yield covariance matrix.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotCovMatrix (const bool useTrkPt, const short pSpc, const short pPtZ, const short pCent) {
  const char* canvasName = Form ("c_CovMatrix_%s_iSpc%i_iPtZ%i_iCent%i", useTrkPt ? "pTch" : "xhZ", pSpc, pPtZ, pCent);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
  }

  c->cd ();

  gPad->SetLogx ();
  gPad->SetLogy ();
  gPad->SetLogz ();

  TH2D* h2 = (useTrkPt ? h2_trk_pt_ptz_cov : h2_trk_xhz_ptz_cov)[pSpc][pPtZ][pCent];
  //TH1D* h1 = (useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[pSpc][pPtZ][pCent];

  //for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
  //  for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
  //    h2->SetBinContent (iX, iY, sqrt (h2->GetBinContent (iX, iY) / ((h1->GetBinContent (iX)) * (h1->GetBinContent (iY)))));
  //    //h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / ((h1->GetBinContent (iX)) * (h1->GetBinContent (iY))));
  //  }
  //}

  h2->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
  h2->GetYaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
  //h2->GetZaxis ()->SetTitle (useTrkPt ? "[cov (Y(x), Y(y)) / N_{evt} #LTY(x)#GT#LTY(y)#GT]^{1/2}" : "[cov (Y(x), Y(y)) / N_{evt} #LTY(x)#GT#LTY(y)#GT]^{1/2}");
  //h2->GetZaxis ()->SetTitle (useTrkPt ? "cov (Y(x), Y(y)) / N_{evt} #LTY(x)#GT#LTY(y)#GT" : "cov (Y(x), Y(y)) / N_{evt} #LTY(x)#GT#LTY(y)#GT");
  h2->GetZaxis ()->SetTitle (useTrkPt ? "cov (Y(x), Y(y)) / N_{evt} [GeV^{-2}]" : "cov (Y(x), Y(y)) / N_{evt}");

  h2->GetZaxis ()->SetTitleOffset (1.2 * h2->GetZaxis ()->GetTitleOffset ());

  h2->GetZaxis ()->SetRangeUser (0.5 * h2->GetMinimum (0), 2 * h2->GetMaximum ()); // log z scale
  //h2->GetZaxis ()->SetRangeUser (1.2 * h2->GetMinimum (), 1.2 * h2->GetMaximum ()); // linear z scale
  h2->Draw ("colz");

  myText (0.23, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  if (pCent == 0) myText (0.23, 0.80, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  else            myText (0.23, 0.80, kBlack, Form ("Pb+Pb, 5.02 TeV, %i-%i%%", (int)centCuts[pCent], (int)centCuts[pCent-1]), 0.04);
  if (pPtZ == nPtZBins-1) myText (0.23, 0.75, kBlack, Form ("Z #rightarrow #it{%s}, #it{p}_{T}^{Z} > %g GeV", (pSpc == 0 ? "ee" : (pSpc == 1 ? "#mu#mu" : "ll")), zPtBins[pPtZ]), 0.04);
  else                    myText (0.23, 0.75, kBlack, Form ("Z #rightarrow #it{%s}, %g < #it{p}_{T}^{Z} < %g GeV", (pSpc == 0 ? "ee" : (pSpc == 1 ? "#mu#mu" : "ll")), zPtBins[pPtZ], zPtBins[pPtZ+1]), 0.04);

  c->SaveAs (Form ("%s/CovMatrices/covMatrix_%s_iSpc%i_iPtZ%i_iCent%i.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", pSpc, pPtZ, pCent));
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
      } // end loop over iCent

      ppHist = h_trk_xhz_ptz_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        PbPbHist = (TH1D*)(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())));
        PbPbHist->Divide (ppHist);
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      } // end loop over iCent

      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        TH1D* ppHist = h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* PbPbHist = (TH1D*)(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())));
          PbPbHist->Divide (ppHist);
          h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
        } // end loop over iCent
        ppHist = h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* PbPbHist = (TH1D*)(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())));
          PbPbHist->Divide (ppHist);
          h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
        } // end loop over iCent
      } // end loop over iPhi
    } // end loop over iPtZ
  } // end loop over iSpc
  iaaCalculated = true;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each deltaPhi bin in all centralities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dPhi (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_iaa_%s_dPhi_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-1), 502);
        c->Divide (numCentBins-1, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        gPad->SetLogx ();

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_iaa : h_trk_xhz_dphi_iaa)[iSpc][iPtZ][iPhi][iCent]);
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
            //g->SetLineWidth (1);
            //g->SetLineColor (colors[iPhi]);
            g->SetFillColorAlpha (fillColors[iPhi], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          g->GetYaxis ()->SetRangeUser (0, max_iaa);

          g->GetXaxis ()->SetMoreLogLabels ();

          g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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

          if (!canvasExists) LabelIAA_dPhi (iCent, iPhi, iPtZ);
        } // end loop over iPhi
      } // end loop over iCent
      c->cd (1);

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over iCent

      c->SaveAs (Form ("%s/IAA/iaa_%s_dPhi_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for IAA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAA_dPhi (const short iCent, const short iPhi, const short iPtZ) {

  if (iCent == 1)
    myText (0.50, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
  else if (iCent == 2) {
    if (iPtZ == nPtZBins-1) myText (0.50, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    else myText (0.50, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
  }
  else if (iCent == numCentBins-1) {
    if (iPhi == 1) {
      myText (0.655, 0.85, kBlack, "#Delta#phi", 0.05);
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
    myMarkerTextNoLine (0.63, 0.925-0.06*(iPhi+1), colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myText (0.65, 0.91-0.06*(iPhi+1), kBlack, Form ("(%s, %s)", lo, hi), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each centrality in all deltaPhi bins
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dCent (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_iaa_%s_dCent_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numPhiBins-1), 500);
        c->Divide (numPhiBins-1, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        gPad->SetLogx ();

        for (int iCent = 1; iCent < numCentBins; iCent++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_iaa : h_trk_xhz_dphi_iaa)[iSpc][iPtZ][iPhi][iCent]);
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
            g->SetLineWidth (0);
            //g->SetLineWidth (1);
            //g->SetLineColor (colors[iCent]);
            g->SetFillColorAlpha (fillColors[iCent], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          g->GetYaxis ()->SetRangeUser (0, max_iaa);
          xmin = g->GetXaxis ()->GetXmin ();
          xmax = g->GetXaxis ()->GetXmax ();

          g->GetXaxis ()->SetMoreLogLabels ();

          g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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

          if (!canvasExists) LabelIAA_dCent (iCent, iPhi, iPtZ);
        } // end loop over iPhi
      } // end loop over iCent
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
      } // end loop over iCent

      c->SaveAs (Form ("%s/IAA/iaa_%s_dCent_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for IAA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAA_dCent (const short iCent, const short iPhi, const short iPtZ) {

  if (iPhi == 1 && iCent == 1) {
    myText (0.50, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
    if (iPtZ == nPtZBins-1) myText (0.50, 0.78, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    else myText (0.50, 0.78, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
  }
  else if (iPhi == numPhiBins-1) {
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerTextNoLine (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    //myMarkerTextNoLine (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myMarkerTextNoLine (0.63, 0.925-0.06*iCent, colors[iCent], kFullCircle, "", 1.4, 0.05);
    myText (0.65, 0.91-0.06*iCent, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each pT^Z bin in all centralities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_iaa_%s_dPtZ_%s", useTrkPt ? "pTch" : "xhZ", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
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

        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[iSpc][iPtZ][iCent]);
        RecenterGraph (g);

        if (!plotAsSystematic) {
          ResetXErrors (g);
          deltaize (g, 1+(iPtZ-3)*0.04, true); // 2.5 = 0.5*(numPhiBins-1)
          g->SetLineColor (colors[iPtZ-1]);
          g->SetMarkerColor (colors[iPtZ-1]);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize ((markerStyle == kOpenDiamond || markerStyle == kFullDiamond) ? 2.0 : 1.2);
          g->SetLineWidth (3);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (0);
          //g->SetLineWidth (1);
          //g->SetLineColor (colors[iPtZ-1]);
          g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
        g->GetYaxis ()->SetRangeUser (0, max_iaa);

        g->GetXaxis ()->SetMoreLogLabels ();

        g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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

        if (!canvasExists) LabelIAA_dPtZ (iCent, iPtZ);
      } // end loop over iPtZ
    } // end loop over iCent

    for (short iCent = 1; iCent < numCentBins; iCent++) {
      c->cd (iCent);
      myText (0.22, 0.22, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    } // end loop over iCent

    c->SaveAs (Form ("%s/IAA/iaa_%s_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", spc));
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for IAA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAA_dPtZ (const short iCent, const short iPtZ) {

  if (iCent == 1 && iPtZ == 2) {
    myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.065);
    myText (0.22, 0.78, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  }
  else if (iCent == 2 && iPtZ == 2) {
    const char* lo = GetPiString (phiLowBins[1]);
    const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
    myText (0.50, 0.88, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.05);
  }
  else if (iCent == numCentBins-1) {
    myMarkerTextNoLine (0.55, 0.900-0.06*(iPtZ-2), colors[iPtZ-1], markerStyles[iPtZ-2], "", (markerStyles[iPtZ-2] == kOpenDiamond || markerStyles[iPtZ-2] == kFullDiamond) ? 2.0 : 1.2, 0.05);
    if (iPtZ == nPtZBins-1) myText (0.565, 0.89-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    else myText (0.565, 0.89-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each pT^Z bin in a single centrality
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSingleIAA_dPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pPtZ, const short iCent, const short pSpc) {
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
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_iaa_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", pPtZ, iCent, spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 800, 800);
      gDirectory->Add (c);
    }

    double xmin = 0, xmax = 0;
    gPad->SetLogx ();
    gPad->SetLogy ();

    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      const Style_t markerStyle = markerStyles[iPtZ-3];

      TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[iSpc][iPtZ][iCent]);
      RecenterGraph (g);

      if (!plotAsSystematic) {
        //deltaize (g, 1+((nPtZBins-2)*((int)useAltMarker)-iPtZ)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
        ResetXErrors (g);
        //deltaize (g, 1+(iPtZ-3)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
        g->SetLineColor (colors[iPtZ-iPtZLo+1]);
        g->SetMarkerColor (colors[iPtZ-iPtZLo+1]);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize ((markerStyle == kOpenDiamond || markerStyle == kFullDiamond) ? 2.0 : 1.2);
        g->SetLineWidth (2);
      } else {
        g->SetMarkerSize (0); 
        g->SetLineWidth (0);
        //g->SetLineWidth (1);
        //g->SetLineColor (colors[iPtZ-iPtZLo+1]);
        g->SetFillColorAlpha (fillColors[iPtZ-iPtZLo+1], 0.3);
      }

      useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
      //g->GetYaxis ()->SetRangeUser (0.1, max_iaa);
      g->GetYaxis ()->SetRangeUser (0.15, max_iaa);

      g->GetXaxis ()->SetMoreLogLabels ();
      g->GetYaxis ()->SetMoreLogLabels ();

      g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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
    } // end loop over iPtZ

    if (!canvasExists) {
      myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
      myText (0.20, 0.83, kBlack, "Pb+Pb, 5.02 TeV", 0.045);
      const char* lo = GetPiString (phiLowBins[1]);
      const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
      myText (0.20, 0.22, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
      myText (0.20, 0.77, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.04);
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        myMarkerTextNoLine (0.674, 0.868-0.05*(iPtZ-iPtZLo), colors[iPtZ-iPtZLo+1], markerStyles[iPtZ-3], "", (markerStyles[iPtZ-3] == kOpenDiamond || markerStyles[iPtZ-3] == kFullDiamond) ? 3.2 : 2.0, 0.04);
        if (iPtZ == nPtZBins-1) myText (0.690, 0.86-0.05*(iPtZ-iPtZLo), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.032);
        else myText (0.690, 0.86-0.05*(iPtZ-iPtZLo), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.032);
      }

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    }

    c->SaveAs (Form ("%s/IAA/iaa_%s_iCent%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iCent, spc));
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between some Pb+Pb centrality and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dPtZ_SpcComp (const bool useTrkPt, const bool plotAsSystematic) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  short iPtZLo = 2;
  short iPtZHi = nPtZBins;
  short iCentLo = 1;
  short iCentHi = numCentBins;

  const int axisTextSize = 22;

  const char* canvasName = Form ("c_iaa_%s_SpcComp", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 400*(numCentBins-1), 1200);
    c->Divide (numCentBins-1, iPtZHi - iPtZLo);
    gDirectory->Add (c);
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {

    for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        const Style_t markerStyle = kFullCircle;
        c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);

        gPad->SetLogx ();

        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[iSpc][iPtZ][iCent]);
        RecenterGraph (g);

        if (!plotAsSystematic) {
          ResetXErrors (g);
          deltaize (g, 1+iSpc*0.05-0.025, true); // 2.5 = 0.5*(numPhiBins-1)
          g->SetLineColor (colors[iSpc+1]);
          g->SetMarkerColor (colors[iSpc+1]);
          //g->SetLineColor (colors[2*(iPtZ-iPtZLo)+iSpc+1]);
          //g->SetMarkerColor (colors[2*(iPtZ-iPtZLo)+iSpc+1]);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (3);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (0);
          //g->SetLineWidth (1);
          //g->SetLineColor (colors[2*(iPtZ-iPtZLo)+iSpc+1]);
          g->SetFillColorAlpha (fillColors[iSpc+1], 0.3);
          //g->SetFillColorAlpha (fillColors[2*(iPtZ-iPtZLo)+iSpc+1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
        g->GetYaxis ()->SetRangeUser (0, max_iaa);

        g->GetXaxis ()->SetMoreLogLabels ();
        g->GetYaxis ()->SetMoreLogLabels ();

        g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
        g->GetYaxis ()->SetTitle ("I_{AA}");

        g->GetXaxis ()->SetTitleFont (43);
        g->GetXaxis ()->SetTitleSize (axisTextSize);
        g->GetXaxis ()->SetLabelFont (43);
        g->GetXaxis ()->SetLabelSize (axisTextSize);

        g->GetYaxis ()->SetTitleFont (43);
        g->GetYaxis ()->SetTitleSize (axisTextSize);
        g->GetYaxis ()->SetLabelFont (43);
        g->GetYaxis ()->SetLabelSize (axisTextSize);

        g->GetXaxis ()->SetTitleOffset (3.0 * g->GetXaxis ()->GetTitleOffset ());
        g->GetYaxis ()->SetTitleOffset (3.0 * g->GetYaxis ()->GetTitleOffset ());

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
            myText (0.65, 0.77, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
          if (iCent == iCentLo && iPtZ == iPtZLo) {
            myText (0.25, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.08);
            myText (0.25, 0.75, kBlack, "Pb+Pb, 5.02 TeV", 0.07);
          }
          //if (iCent == iCentLo+1 && iPtZ == iPtZLo) {
          //  const char* lo = "3#pi/4";
          //  const char* hi = "#pi";
          //  myText (0.25, 0.85, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.06);
          //}
        }
        for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
          if (iCent == iCentLo+1) {
            c->cd (iCent-iCentLo + 1);
            if (iPtZ == iPtZLo) {
              myText (0.656, 0.875, kBlack, "#it{ee}", 0.06);
              myText (0.736, 0.875, kBlack, "#it{#mu#mu}", 0.06);
              myMarkerTextNoLine (0.790, 0.810-0.10*(iPtZ-iPtZLo), colors[2*(iPtZ-iPtZLo)+1], kFullCircle, "", 1.3, 0.06);
              myMarkerTextNoLine (0.710, 0.810-0.10*(iPtZ-iPtZLo), colors[2*(iPtZ-iPtZLo)+2], kFullCircle, "", 1.3, 0.06);
            }
          }
          else if (iCent == iCentLo+2) {
            c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo)+iCent);
            if (iPtZ == nPtZBins-1) myText (0.510, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.07);
            else                    myText (0.430, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.07);
          }
        } // end loop over iPtZ
      } // end loop over iCent
    }
  } // end loop over iSpc
  for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);
      const float xmin = useTrkPt ? trk_min_pt : allXhZBins[0];
      const float xmax = useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : allXhZBins[maxNXhZBins];
      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    }
  } // end loop over iCent

  //if (!plotAsSystematic) {
  //  for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
  //    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
  //      c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);

  //      for (short iSpc = 0; iSpc < 2; iSpc++) {
  //        TGAE* g1 = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[0][iPtZ][iCent]);
  //        TGAE* g2 = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[1][iPtZ][iCent]);

  //        double chisq = 0;
  //        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
  //        for (int ix = 0; ix < g1->GetN (); ix++) {
  //          g1->GetPoint (ix, x1, y1);
  //          g2->GetPoint (ix, x2, y2);

  //          if (y1 > y2) {
  //            chisq += pow (y2-y1, 2)/ (pow (g2->GetErrorYhigh (ix), 2) + pow (g1->GetErrorYlow (ix), 2));
  //          
  //          }
  //          else if (y1 < y2) {
  //            chisq += pow (y2-y1, 2)/ (pow (g1->GetErrorYhigh (ix), 2) + pow (g2->GetErrorYlow (ix), 2));
  //          }
  //        }
  //        if (iPtZ == 4)
  //          myText (0.50, 0.61, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chisq, 6), 0.06);
  //        else
  //          myText (0.50, 0.61, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chisq, 5), 0.06);
  //      } // end loop over iSpc
  //    } // end loop over iPtZ
  //  } // end loop over iCent
  //}

  c->SaveAs (Form ("%s/IAA/iaa_%s_dPtZ_SpcComp.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
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
        } // end loop over iCent

        periphHist = h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent]) {
            centHist = (TH1D*)(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            h_trk_xhz_dphi_icp[iSpc][iPtZ][iPhi][iCent] = centHist;
          }
        } // end loop over iCent
      } // end loop over iPhi
    } // end loop over iPtZ
  } // end loop over iSpc
  icpCalculated = true;
  return;
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots ICP for each deltaPhi bin in all centralities
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: PlotICP_dPhi (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_icp_%s_dPhi_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-2), 500);
        c->Divide (numCentBins-2, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        gPad->SetLogx ();

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
          
          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_icp : h_trk_xhz_dphi_icp)[iSpc][iPtZ][iPhi][iCent]);
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
            g->SetLineWidth (0);
            //g->SetLineWidth (1);
            //g->SetLineColor (colors[iPhi]);
            g->SetFillColorAlpha (fillColors[iPhi], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          g->GetYaxis ()->SetRangeUser (0, max_icp);
          xmin = g->GetXaxis ()->GetXmin ();
          xmax = g->GetXaxis ()->GetXmax ();

          g->GetXaxis ()->SetMoreLogLabels ();

          g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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

          if (!canvasExists) LabelICP_dPhi (iCent, iPhi, iPtZ);
        } // end loop over iPhi
      } // end loop over iCent

      c->cd (1);
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over iCent

      c->SaveAs (Form ("%s/ICP/icp_%s_dPhi_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for ICP measurements
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: LabelICP_dPhi (const short iCent, const short iPhi, const short iPtZ) {
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
// Plots ICP for each centrality in all deltaPhi bins
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: PlotICP_dCent (const bool useTrkPt, const bool plotAsSystematic, const short pSpc, const short pPtZ) {
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

      const char* canvasName = Form ("c_icp_%s_dCent_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numPhiBins-1), 500);
        c->Divide (numPhiBins-1, 1);
        gDirectory->Add (c);
      }

      double xmin = 0, xmax = 0;
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        gPad->SetLogx ();

        for (int iCent = 2; iCent < numCentBins; iCent++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
          
          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_icp : h_trk_xhz_dphi_icp)[iSpc][iPtZ][iPhi][iCent]);
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
            g->SetLineWidth (0);
            //g->SetLineWidth (1);
            //g->SetLineColor (colors[iCent-1]);
            g->SetFillColorAlpha (fillColors[iCent-1], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          g->GetYaxis ()->SetRangeUser (0, max_icp);
          xmin = g->GetXaxis ()->GetXmin ();
          xmax = g->GetXaxis ()->GetXmax ();

          g->GetXaxis ()->SetMoreLogLabels ();

          g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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

          if (!canvasExists) LabelICP_dCent (iCent, iPhi, iPtZ);
        } // end loop over iPhi
      } // end loop over iCent

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
      } // end loop over iCent

      c->SaveAs (Form ("%s/ICP/icp_%s_dCent_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for ICP measurements
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: LabelICP_dCent (const short iCent, const short iPhi, const short iPtZ) {
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
// Plots IAA for each pT^Z bin in all centralities
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: PlotICP_dPtZ (const bool useTrkPt, const bool plotAsSystematic, const short pSpc) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!icpCalculated)
    CalculateICP ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_icp_%s_dPtZ_%s", useTrkPt ? "pTch" : "xhZ", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 500*(numCentBins-2), 500);
      c->Divide (numCentBins-2, 1);
      gDirectory->Add (c);
    }

    double xmin = 0, xmax = 0;
    for (short iCent = 2; iCent < numCentBins; iCent++) {
      c->cd (iCent-1);
      gPad->SetLogx ();

      for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
        
        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_icp : h_trk_xhz_ptz_icp)[iSpc][iPtZ][iCent]);
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
          g->SetLineWidth (0);
          //g->SetLineWidth (1);
          //g->SetLineColor (colors[iPtZ-1]);
          g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
        g->GetYaxis ()->SetRangeUser (0, max_icp);
        xmin = g->GetXaxis ()->GetXmin ();
        xmax = g->GetXaxis ()->GetXmax ();

        g->GetXaxis ()->SetMoreLogLabels ();

        g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
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

        if (!canvasExists) LabelICP_dPtZ (iCent, iPtZ);
      } // end loop over iPhi
    } // end loop over iCent

    for (short iCent = 2; iCent < numCentBins; iCent++) {
      c->cd (iCent-1);
      myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    } // end loop over iCent

    c->SaveAs (Form ("%s/ICP/icp_%s_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", spc));
  } // end loop over iSpc
}
*/




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for ICP measurements
////////////////////////////////////////////////////////////////////////////////////////////////
/*
void PhysicsAnalysis :: LabelICP_dPtZ (const short iCent, const short iPtZ) {
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

  TF1**** f_trk_pt_ptz_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);
  TF1**** f_trk_xhz_ptz_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Fit (f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Fit (f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
      }
    }
  }

  outFile->Close ();

  Delete3DArray (f_trk_pt_ptz_iaa, 3, nPtZBins, numCentBins);
  Delete3DArray (f_trk_xhz_ptz_iaa, 3, nPtZBins, numCentBins);

  _gDirectory->cd ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots signal-to-background panels for data
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSignalToBkg (const bool useTrkPt, const short iSpc) {
  TCanvas* c = new TCanvas ("c_stb", "", 1200, 400);
  c->Divide (3, 1);

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    c->cd (iPtZ-2+1);
    gPad->SetLogx ();
    gPad->SetLogy ();

    TH1D* htemp = new TH1D (Form ("htemp_iPtZ%i", iPtZ), "", (useTrkPt ? nPtchBins : nXhZBins)[iPtZ], (useTrkPt ? pTchBins : xhZBins)[iPtZ]);
    htemp->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    htemp->GetYaxis ()->SetTitle ("Y / Y_{bkg}");

    htemp->SetMarkerSize (0);
    htemp->SetLineWidth (0);

    htemp->GetXaxis ()->SetMoreLogLabels ();

    htemp->GetYaxis ()->SetRangeUser (8e-4, 2e4);
    htemp->Draw ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_sig_to_bkg : h_trk_xhz_ptz_sig_to_bkg)[iSpc][iPtZ][iCent]);
      RecenterGraph (g);
      ResetXErrors (g);
      deltaize (g, 1 + 0.01*(iCent - (numCentBins-1)), true);
      
      g->SetMarkerStyle (markerStyles[iCent]);
      g->SetMarkerSize ((markerStyles[iCent] == kOpenDiamond || markerStyles[iCent] == kFullDiamond) ? 2.0 : 1.2);
      g->SetMarkerColor (colors[iCent]);
      g->SetLineColor (colors[iCent]);
      g->SetLineWidth (2);
      g->Draw ("P");
    }

    if (iPtZ == 2) {
      myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.060);
      myText (0.22, 0.81, kBlack, "Pb+Pb, 5.02 TeV, 1.7 nb^{-1}", 0.050);
      myText (0.22, 0.74, kBlack, "#it{pp}, 5.02 TeV, 260 pb^{-1}", 0.050);
      myText (0.56, 0.24, kBlack, "15 < #it{p}_{T}^{Z} < 30 GeV", 0.050);
    }
    else if (iPtZ == 3) {
      myText (0.56, 0.24, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.050);
    }
    else if (iPtZ == 4) {
      myMarkerTextNoLine (0.25, 0.88, colors[0], markerStyles[0], "#it{pp}", 1.2, 0.05);
      myMarkerTextNoLine (0.25, 0.815, colors[1], markerStyles[1], "30-80%", 1.2, 0.05);
      myMarkerTextNoLine (0.25, 0.75, colors[2], markerStyles[2], "10-30%", 2.0, 0.05);
      myMarkerTextNoLine (0.25, 0.685, colors[3], markerStyles[3], "0-10%", 1.2, 0.05);
      myText (0.62, 0.24, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.05);
    }
  }

  c->SaveAs (Form ("%s/TrkYields/sigToBkg_%s_dPtZ.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints IAA values for the bin specified
////////////////////////////////////////////////////////////////////////////////////////////////
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

