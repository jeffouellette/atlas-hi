#ifndef __FullAnalysis_h__
#define __FullAnalysis_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

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

#include <iostream>
#include <string>

using namespace atlashi;
using namespace std;

class FullAnalysis : public PhysicsAnalysis {

  public:

  // Analysis checks
  TH1D*   h_fcal_et               = nullptr;
  TH1D*   h_fcal_et_reweighted    = nullptr;

  TH1D**  h_q2                  = Get1DArray <TH1D*> (numFinerCentBins);
  TH1D**  h_q2_reweighted       = Get1DArray <TH1D*> (numFinerCentBins);
  TH1D**  h_psi2                = Get1DArray <TH1D*> (numFinerCentBins);
  TH1D**  h_psi2_reweighted     = Get1DArray <TH1D*> (numFinerCentBins);
  TH1D*   h_PbPb_vz             = nullptr;
  TH1D*   h_PbPb_vz_reweighted  = nullptr;
  TH1D*   h_pp_vz               = nullptr;
  TH1D*   h_pp_vz_reweighted    = nullptr;
  TH1D*   h_pp_nch              = nullptr;
  TH1D*   h_pp_nch_reweighted   = nullptr;

  TH1D*** h_z_phi         = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_z_pt          = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_z_pt_ratio    = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH2D**** h_z_y_phi      = Get3DArray <TH2D*> (numCentBins, 3, nPtZBins);   // iCent, iSpc, iPtZ
  TH1D**** h_z_eta        = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins);   // iCent, iSpc, iPtZ
  TH1D**** h_z_eta_ratio  = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins);   // iCent, iSpc, iPtZ
  TH1D**** h_z_y          = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins);   // iCent, iSpc, iPtZ
  TH1D**** h_z_y_ratio    = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins);   // iCent, iSpc, iPtZ
  TH1D**** h_z_m          = Get3DArray <TH1D*> (numCentBins, 3, 3);          // iCent, iSpc, iReg
  TH1D**** h_z_m_ratio    = Get3DArray <TH1D*> (numCentBins, 3, 3);          // iCent, iSpc, iReg
  TH1D*** h_z_lepton_dphi = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_pt     = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_eta    = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_trk_pt = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_trk_pt        = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH2D*** h_lepton_trk_dr = Get2DArray <TH2D*> (numCentBins, 3);             // iCent, iSpc


  FullAnalysis () : PhysicsAnalysis () { }

  FullAnalysis (const char* _name) {//, const char* subDir) {
    FullAnalysis ();
    name = _name;
    //directory = Form ("%s/", subDir);
    plotFill = false;
  }

  virtual ~FullAnalysis () {

    Delete1DArray (h_q2,            numFinerCentBins);
    Delete1DArray (h_q2_reweighted, numFinerCentBins);

    Delete2DArray (h_z_phi,         numCentBins, 3);
    Delete2DArray (h_z_pt,          numCentBins, 3);
    Delete2DArray (h_z_pt_ratio,    numCentBins, 3);
    Delete3DArray (h_z_y_phi,       numCentBins, 3, nPtZBins);
    Delete3DArray (h_z_eta,         numCentBins, 3, nPtZBins);
    Delete3DArray (h_z_eta_ratio,   numCentBins, 3, nPtZBins);
    Delete3DArray (h_z_y,           numCentBins, 3, nPtZBins);
    Delete3DArray (h_z_y_ratio,     numCentBins, 3, nPtZBins);
    Delete3DArray (h_z_m,           numCentBins, 3, 3);
    Delete3DArray (h_z_m_ratio,     numCentBins, 3, 3);
    Delete2DArray (h_z_lepton_dphi, numCentBins, 3);
    Delete2DArray (h_lepton_pt,     numCentBins, 3);
    Delete2DArray (h_lepton_eta,    numCentBins, 3);
    Delete2DArray (h_lepton_trk_pt, numCentBins, 3);
    Delete2DArray (h_trk_pt,        numCentBins, 3);
    Delete2DArray (h_lepton_trk_dr, numCentBins, 3);
  }


  protected:
  void LabelZMassSpectra (const short iSpc, const short iCent, const short iReg);

  public:

  virtual void CreateHists () override;
  virtual void CopyAnalysis (FullAnalysis* a, const bool copyBkgs = false);
  virtual void CombineHists () override;
  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  virtual void SaveHists (const char* histFileName = "savedHists.root") override;
  virtual void ScaleHists () override;
  virtual void Execute (const char* inFileName = "outFile.root", const char* outFileName = "savedHists.root") override;

  //void PrintZYields ();

  void PlotFCalDists (const bool _treatAsData = true);
  void PlotQ2Dists (const bool _treatAsData = true);
  void PlotQ2Weights ();
  void PlotPsi2Dists (const bool _treatAsData = true);
  void PlotPsi2Weights ();
  void PlotVZDists (const bool _treatAsData = true);
  void PlotNchDists (const bool _treatAsData = true);

  void PlotLeptonPtSpectra ();
  void PlotLeptonEtaSpcComp (FullAnalysis* a = nullptr);
  void PlotLeptonTrackPtSpectra ();
  void PlotLeptonTrackDR ();
  void PlotLeptonTrackDRProjX ();
  void PlotZPtSpectra ();
  void PlotZYPhiMap (const short pPtZ);
  void PlotZEtaMap (const short pPtZ);
  void PlotZYMap (const short pPtZ);
  void PlotZYMapSpcComp (const short pPtZ, FullAnalysis* a = nullptr);
  void PlotZMassSpectra ();
  void PlotZPhiYield ();
  void PlotZLeptonDPhi ();

  void CalculateZPtDistRatio (FullAnalysis* a);
  void CalculateZEtaDistRatio ();
  void CalculateZYDistRatio (FullAnalysis* a);
  void CalculateZMassSpectraRatio (FullAnalysis* a);

};


////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CreateHists () {
  PhysicsAnalysis :: CreateHists ();

  h_fcal_et = new TH1D (Form ("h_fcal_et_%s", name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins); 
  h_fcal_et->Sumw2 ();
  h_fcal_et_reweighted = new TH1D (Form ("h_fcal_et_reweighted_%s", name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
  h_fcal_et_reweighted->Sumw2 ();

  for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
    h_q2[iCent]               = new TH1D (Form ("h_q2_iCent%i_%s", iCent, name.c_str ()), "", 50, 0, 0.3);
    h_q2[iCent]->Sumw2 ();
    h_q2_reweighted[iCent]    = new TH1D (Form ("h_q2_reweighted_iCent%i_%s", iCent, name.c_str ()), "", 50, 0, 0.3);
    h_q2_reweighted[iCent]->Sumw2 ();
    h_psi2[iCent]             = new TH1D (Form ("h_psi2_iCent%i_%s", iCent, name.c_str ()), "", 50, -pi/2, pi/2);
    h_psi2[iCent]->Sumw2 ();
    h_psi2_reweighted[iCent]  = new TH1D (Form ("h_psi2_reweighted_iCent%i_%s", iCent, name.c_str ()), "", 50, -pi/2, pi/2);
    h_psi2_reweighted[iCent]->Sumw2 ();
  }
  h_PbPb_vz = new TH1D (Form ("h_PbPb_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_PbPb_vz_reweighted = new TH1D (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz = new TH1D (Form ("h_pp_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz_reweighted = new TH1D (Form ("h_pp_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_nch = new TH1D (Form ("h_pp_nch_%s", name.c_str ()), "", 80, -0.5, 160.5);
  h_pp_nch_reweighted = new TH1D (Form ("h_pp_nch_reweighted_%s", name.c_str ()), "", 80, -0.5, 160.5);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = new TH1D (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 80, 0, pi);
      h_z_pt[iCent][iSpc]           = new TH1D (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 300, 0, 300);
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]  = new TH2D (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 20, -2.5, 2.5, 20, 0, 2*pi);
        h_z_eta[iCent][iSpc][iPtZ]    = new TH1D (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 16, -5.0, 5.0);
        h_z_y[iCent][iSpc][iPtZ]      = new TH1D (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 16, -2.5, 2.5);
      }
      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg]    = new TH1D (Form ("h_z_m_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()), "", 40, 76, 106);
        h_z_m[iCent][iSpc][iReg]->Sumw2 ();
      }
      h_z_lepton_dphi[iCent][iSpc]  = new TH1D (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 45, 0, pi);
      h_lepton_pt[iCent][iSpc]      = new TH1D (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 250, 0, 250);
      h_lepton_eta[iCent][iSpc]     = new TH1D (Form ("h_lepton_eta_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 50, -2.5, 2.5);
      h_lepton_trk_pt[iCent][iSpc]  = new TH1D (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 250, 0, 250);
      h_trk_pt[iCent][iSpc]         = new TH1D (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 100, 0, 100);
      h_lepton_trk_dr[iCent][iSpc]  = new TH2D (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 100, 0., 1., 40, 0, 2);
      
      h_z_phi[iCent][iSpc]->Sumw2 ();
      h_z_pt[iCent][iSpc]->Sumw2 ();
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]->Sumw2 ();
        h_z_eta[iCent][iSpc][iPtZ]->Sumw2 ();
        h_z_y[iCent][iSpc][iPtZ]->Sumw2 ();
      }
      h_z_lepton_dphi[iCent][iSpc]->Sumw2 ();
      h_lepton_pt[iCent][iSpc]->Sumw2 ();
      h_lepton_eta[iCent][iSpc]->Sumw2 ();
      h_lepton_trk_pt[iCent][iSpc]->Sumw2 ();
      h_trk_pt[iCent][iSpc]->Sumw2 ();
      h_lepton_trk_dr[iCent][iSpc]->Sumw2 ();

    }
  }

  histsLoaded = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms as clones from another analysis, where appropriate
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CopyAnalysis (FullAnalysis* a, const bool copyBkgs) {
  if (name == "")
    cout << "Warning in FullAnalysis :: CopyAnalysis: name of analysis not set!" << endl;

  PhysicsAnalysis :: CopyAnalysis ((PhysicsAnalysis*)a, copyBkgs);

  // Don't need to clone these histograms

  // Should clone these histograms
  h_fcal_et               = (TH1D*) a->h_fcal_et->Clone (Form ("h_fcal_et_%s", name.c_str ()));
  h_fcal_et_reweighted    = (TH1D*) a->h_fcal_et_reweighted->Clone (Form ("h_fcal_et_reweighted_%s", name.c_str ()));

  for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
    h_q2[iCent]               = (TH1D*) a->h_q2[iCent]->Clone (Form ("h_q2_iCent%i_%s", iCent, name.c_str ()));
    h_q2_reweighted[iCent]    = (TH1D*) a->h_q2_reweighted[iCent]->Clone (Form ("h_q2_reweighted_iCent%i_%s", iCent, name.c_str ()));
    h_psi2[iCent]             = (TH1D*) a->h_psi2[iCent]->Clone (Form ("h_psi2_iCent%i_%s", iCent, name.c_str ()));
    h_psi2_reweighted[iCent]  = (TH1D*) a->h_psi2_reweighted[iCent]->Clone (Form ("h_psi2_reweighted_iCent%i_%s", iCent, name.c_str ()));
  }
  h_PbPb_vz             = (TH1D*) a->h_PbPb_vz->Clone (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted  = (TH1D*) a->h_PbPb_vz_reweighted->Clone (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz               = (TH1D*) a->h_pp_vz->Clone (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted    = (TH1D*) a->h_pp_vz_reweighted->Clone (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  h_pp_nch              = (TH1D*) a->h_pp_nch->Clone (Form ("h_pp_nch_%s", name.c_str ()));
  h_pp_nch_reweighted   = (TH1D*) a->h_pp_nch_reweighted->Clone (Form ("h_pp_nch_reweighted_%s", name.c_str ()));

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = (TH1D*) a->h_z_phi[iCent][iSpc]->Clone (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt[iCent][iSpc]           = (TH1D*) a->h_z_pt[iCent][iSpc]->Clone (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_z_pt_ratio[iCent][iSpc])
        h_z_pt_ratio[iCent][iSpc]   = (TH1D*) a->h_z_pt_ratio[iCent][iSpc]->Clone (Form ("h_z_pt_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]        = (TH2D*) a->h_z_y_phi[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_eta[iCent][iSpc][iPtZ]          = (TH1D*) a->h_z_eta[iCent][iSpc][iPtZ]->Clone (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        if (a->h_z_eta_ratio[iCent][iSpc][iPtZ])
          h_z_eta_ratio[iCent][iSpc][iPtZ]  = (TH1D*) a->h_z_eta_ratio[iCent][iSpc][iPtZ]->Clone (Form ("h_z_eta_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_y[iCent][iSpc][iPtZ]            = (TH1D*) a->h_z_y[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        if (a->h_z_y_ratio[iCent][iSpc][iPtZ])
          h_z_y_ratio[iCent][iSpc][iPtZ]    = (TH1D*) a->h_z_y_ratio[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
      }

      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg]    = (TH1D*) a->h_z_m[iCent][iSpc][iReg]->Clone (Form ("h_z_m_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
        if (a->h_z_m_ratio[iCent][iSpc][iReg])
          h_z_m_ratio[iCent][iSpc][iReg] = (TH1D*) a->h_z_m_ratio[iCent][iSpc][iReg]->Clone (Form ("h_z_m_ratio_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
      }

      if (a->h_z_lepton_dphi[iCent][iSpc])  h_z_lepton_dphi[iCent][iSpc]  = (TH1D*) a->h_z_lepton_dphi[iCent][iSpc]->Clone (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_pt[iCent][iSpc])      h_lepton_pt[iCent][iSpc]      = (TH1D*) a->h_lepton_pt[iCent][iSpc]->Clone (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_eta[iCent][iSpc])     h_lepton_eta[iCent][iSpc]     = (TH1D*) a->h_lepton_eta[iCent][iSpc]->Clone (Form ("h_lepton_eta_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_trk_pt[iCent][iSpc])  h_lepton_trk_pt[iCent][iSpc]  = (TH1D*) a->h_lepton_trk_pt[iCent][iSpc]->Clone (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_trk_pt[iCent][iSpc])         h_trk_pt[iCent][iSpc]         = (TH1D*) a->h_trk_pt[iCent][iSpc]->Clone (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_trk_dr[iCent][iSpc])  h_lepton_trk_dr[iCent][iSpc]  = (TH2D*) a->h_lepton_trk_dr[iCent][iSpc]->Clone (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()));

    } // end loop over iSpc
  } // end loop over iCent

  histsLoaded = true;
  histsScaled = true;
  return;    
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  SetupDirectories ("", "ZTrackAnalysis/");
  if (histsLoaded)
    return;

  PhysicsAnalysis :: LoadHists (histFileName, false);

  TDirectory* _gDirectory = gDirectory;
  if (!histFile->IsOpen ()) {
    cout << "Error in FullAnalysis :: LoadHists: histFile not open after calling parent function, exiting." << endl;
    return;
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = (TH1D*) histFile->Get (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt[iCent][iSpc]           = (TH1D*) histFile->Get (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt_ratio[iCent][iSpc]     = (TH1D*) histFile->Get (Form ("h_z_pt_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]        = (TH2D*) histFile->Get (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_eta[iCent][iSpc][iPtZ]          = (TH1D*) histFile->Get (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_eta_ratio[iCent][iSpc][iPtZ]    = (TH1D*) histFile->Get (Form ("h_z_eta_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_y[iCent][iSpc][iPtZ]            = (TH1D*) histFile->Get (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_y_ratio[iCent][iSpc][iPtZ]      = (TH1D*) histFile->Get (Form ("h_z_y_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
      }

      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg]        = (TH1D*) histFile->Get (Form ("h_z_m_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
        h_z_m_ratio[iCent][iSpc][iReg]  = (TH1D*) histFile->Get (Form ("h_z_m_ratio_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
      }
      h_z_lepton_dphi[iCent][iSpc]  = (TH1D*) histFile->Get (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_pt[iCent][iSpc]      = (TH1D*) histFile->Get (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_eta[iCent][iSpc]     = (TH1D*) histFile->Get (Form ("h_lepton_eta_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_pt[iCent][iSpc]  = (TH1D*) histFile->Get (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_trk_pt[iCent][iSpc]         = (TH1D*) histFile->Get (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_dr[iCent][iSpc]  = (TH2D*) histFile->Get (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()));

    }
  }
  h_fcal_et               = (TH1D*) histFile->Get (Form ("h_fcal_et_%s", name.c_str ()));
  h_fcal_et_reweighted    = (TH1D*) histFile->Get (Form ("h_fcal_et_reweighted_%s", name.c_str ()));

  for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
    h_q2[iCent]               = (TH1D*) histFile->Get (Form ("h_q2_iCent%i_%s", iCent, name.c_str ()));
    h_q2_reweighted[iCent]    = (TH1D*) histFile->Get (Form ("h_q2_reweighted_iCent%i_%s", iCent, name.c_str ()));
    h_psi2[iCent]             = (TH1D*) histFile->Get (Form ("h_psi2_iCent%i_%s", iCent, name.c_str ()));
    h_psi2_reweighted[iCent]  = (TH1D*) histFile->Get (Form ("h_psi2_reweighted_iCent%i_%s", iCent, name.c_str ()));
  }
  h_PbPb_vz             = (TH1D*) histFile->Get (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted  = (TH1D*) histFile->Get (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz               = (TH1D*) histFile->Get (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted    = (TH1D*) histFile->Get (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  h_pp_nch               = (TH1D*) histFile->Get (Form ("h_pp_nch_%s", name.c_str ()));
  h_pp_nch_reweighted    = (TH1D*) histFile->Get (Form ("h_pp_nch_reweighted_%s", name.c_str ()));
  
  histsLoaded = true;

  if (_finishHists) {
    CombineHists ();
    ScaleHists ();
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: SaveHists (const char* histFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");
  if (!histsLoaded)
    return;

  PhysicsAnalysis :: SaveHists (histFileName);

  TDirectory* _gDirectory = gDirectory;
  if (!histFile) {
    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "update");
    histFile->cd ();
  }
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      SafeWrite (h_z_phi[iCent][iSpc]);
      SafeWrite (h_z_pt[iCent][iSpc]);
      SafeWrite (h_z_pt_ratio[iCent][iSpc]);
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_z_y_phi[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_eta[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_eta_ratio[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_y[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_y_ratio[iCent][iSpc][iPtZ]);
      }
      for (short iReg = 0; iReg < 3; iReg++) {
        SafeWrite (h_z_m[iCent][iSpc][iReg]);
        SafeWrite (h_z_m_ratio[iCent][iSpc][iReg]);
      }
      SafeWrite (h_z_lepton_dphi[iCent][iSpc]);
      SafeWrite (h_lepton_pt[iCent][iSpc]);
      SafeWrite (h_lepton_eta[iCent][iSpc]);
      SafeWrite (h_lepton_trk_pt[iCent][iSpc]);
      SafeWrite (h_trk_pt[iCent][iSpc]);
      SafeWrite (h_lepton_trk_dr[iCent][iSpc]);
    }
  }
  SafeWrite (h_fcal_et);
  SafeWrite (h_fcal_et_reweighted);

  for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
    SafeWrite (h_q2[iCent]);
    SafeWrite (h_q2_reweighted[iCent]);
    SafeWrite (h_psi2[iCent]);
    SafeWrite (h_psi2_reweighted[iCent]);
  }
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
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CombineHists () {

  PhysicsAnalysis :: CombineHists ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      if (h_z_phi[iCent][iSpc]) h_z_phi[iCent][2]->Add (h_z_phi[iCent][iSpc]);
      if (h_z_pt[iCent][iSpc]) h_z_pt[iCent][2]->Add (h_z_pt[iCent][iSpc]);
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        if (h_z_y_phi[iCent][iSpc][iPtZ]) h_z_y_phi[iCent][2][iPtZ]->Add (h_z_y_phi[iCent][iSpc][iPtZ]);
        if (h_z_eta[iCent][iSpc][iPtZ]) h_z_eta[iCent][2][iPtZ]->Add (h_z_eta[iCent][iSpc][iPtZ]);
        if (h_z_y[iCent][iSpc][iPtZ]) h_z_y[iCent][2][iPtZ]->Add (h_z_y[iCent][iSpc][iPtZ]);
      }
      for (short iReg = 0; iReg < 2; iReg++) {
        if (h_z_m[iCent][iSpc][iReg]) h_z_m[iCent][2][iReg]->Add (h_z_m[iCent][iSpc][iReg]);
        if (h_z_m[iCent][iSpc][iReg]) h_z_m[iCent][iSpc][2]->Add (h_z_m[iCent][iSpc][iReg]);
        if (h_z_m[iCent][iSpc][iReg]) h_z_m[iCent][2][2]->Add (h_z_m[iCent][iSpc][iReg]);
      }
      if (h_lepton_pt[iCent][iSpc]) h_lepton_pt[iCent][2]->Add (h_lepton_pt[iCent][iSpc]);
      if (h_lepton_eta[iCent][iSpc]) h_lepton_eta[iCent][2]->Add (h_lepton_eta[iCent][iSpc]);
      if (h_lepton_pt[iCent][iSpc]) h_lepton_trk_pt[iCent][2]->Add (h_lepton_pt[iCent][iSpc]);
      if (h_trk_pt[iCent][iSpc]) h_trk_pt[iCent][2]->Add (h_trk_pt[iCent][iSpc]);
      if (h_lepton_trk_dr[iCent][iSpc]) h_lepton_trk_dr[iCent][2]->Add (h_lepton_trk_dr[iCent][iSpc]);

    } // end loop over species
  } // end loop over centralities
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  PhysicsAnalysis :: ScaleHists ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {

      double normFactor = 0;
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
        normFactor += h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2);

      if (h_lepton_pt[iCent][iSpc]) {
        h_lepton_pt[iCent][iSpc]->Rebin (5);
        if (normFactor > 0)
          h_lepton_pt[iCent][iSpc]->Scale (1. / normFactor, "width");
      }

      if (h_lepton_eta[iCent][iSpc]) {
        h_lepton_eta[iCent][iSpc]->Rebin (2);
        if (normFactor > 0)
          h_lepton_eta[iCent][iSpc]->Scale (1. / normFactor, "width");
      }

      if (h_lepton_trk_pt[iCent][iSpc]) {
        h_lepton_trk_pt[iCent][iSpc]->Rebin (5);
        if (normFactor > 0)
          h_lepton_trk_pt[iCent][iSpc]->Scale (1. / normFactor, "width");
      }

      if (h_trk_pt[iCent][iSpc]) {
        h_trk_pt[iCent][iSpc]->Rebin (5);
        if (normFactor > 0)
          h_trk_pt[iCent][iSpc]->Scale (1. / normFactor, "width");
      }

      if ( h_z_pt[iCent][iSpc]) {
        h_z_pt[iCent][iSpc]->Rebin (10);
        h_z_pt[iCent][iSpc]->Scale (0.1);
        if (h_z_pt[iCent][iSpc]->Integral () > 0)
          h_z_pt[iCent][iSpc]->Scale (1. / h_z_pt[iCent][iSpc]->Integral (), "width");
      }

      for (short iReg = 0; iReg < 3; iReg++) {
        //if (h_z_m[iCent][iSpc]->GetMaximum () > 0)
        //  h_z_m[iCent][iSpc]->Scale (1. / (h_z_m[iCent][iSpc]->GetMaximum ()));
        if (h_z_m[iCent][iSpc][iReg]) {
          if (h_z_m[iCent][iSpc][iReg]->Integral () > 0)
            h_z_m[iCent][iSpc][iReg]->Scale (1. / h_z_m[iCent][iSpc][iReg]->Integral ());
        }
      }

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        if (h_z_eta[iCent][iSpc][iPtZ]) {
          if (h_z_eta[iCent][iSpc][iPtZ]->Integral () > 0)
            h_z_eta[iCent][iSpc][iPtZ]->Scale (1. / h_z_eta[iCent][iSpc][iPtZ]->Integral (), "width");
        }

        if (h_z_y[iCent][iSpc][iPtZ]) {
          h_z_y[iCent][iSpc][iPtZ]->Rebin (2);
          if (h_z_y[iCent][iSpc][iPtZ]->Integral () > 0)
            h_z_y[iCent][iSpc][iPtZ]->Scale (1. / h_z_y[iCent][iSpc][iPtZ]->Integral (), "width");
        }
      }

      if (h_z_phi[iCent][iSpc]) {
        h_z_phi[iCent][iSpc]->Rebin (8);
        if (h_z_phi[iCent][iSpc]->Integral () > 0)
          h_z_phi[iCent][iSpc]->Scale (1. / h_z_phi[iCent][iSpc]->Integral (), "width");
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
void FullAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_eta = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;


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

      if (event_weight == 0)
        continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      const short iFinerCent = GetFinerCentBin (fcal_et);
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFinerCent]->Fill (q2);
      h_q2_reweighted[iFinerCent]->Fill (q2, event_weight);
      h_psi2[iFinerCent]->Fill (psi2);
      h_psi2_reweighted[iFinerCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt, event_weight);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt || trk_max_pt < trkpt)
          continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta->at (iTrk), l1_trk_phi, trk_phi->at (iTrk)) < 0.03 || DeltaR (l2_trk_eta, trk_eta->at (iTrk), l2_trk_phi, trk_phi->at (iTrk)) < 0.03))
          continue;

        {
          float mindr = pi;
          float ptdiff = 0;
          float dr = DeltaR (trk_eta->at (iTrk), l1_trk_eta, trk_phi->at (iTrk), l1_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l1_trk_pt) / (trkpt + l1_trk_pt);
          }
          dr = DeltaR (trk_eta->at (iTrk), l2_trk_eta, trk_phi->at (iTrk), l2_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l2_trk_pt) / (trkpt + l2_trk_pt);
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
        }

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), true);
        const double trkPur = doTrackPurVar ? GetTrackingPurity (fcal_et, trkpt, trk_eta->at (iTrk), true) : 1.;
        if (trkEff == 0 || trkPur == 0)
          continue;
        const double trkWeight = event_weight * trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, trkWeight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
        }
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
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

      if (event_weight == 0)
        continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp
      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt || trk_max_pt < trkpt)
          continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta->at (iTrk), l1_trk_phi, trk_phi->at (iTrk)) < 0.03 || DeltaR (l2_trk_eta, trk_eta->at (iTrk), l2_trk_phi, trk_phi->at (iTrk)) < 0.03))
          continue;

        {
          float mindr = pi;
          float ptdiff = 0;
          float dr = DeltaR (trk_eta->at (iTrk), l1_trk_eta, trk_phi->at (iTrk), l1_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l1_trk_pt) / (trkpt + l1_trk_pt);
          }
          dr = DeltaR (trk_eta->at (iTrk), l2_trk_eta, trk_phi->at (iTrk), l2_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l2_trk_pt) / (trkpt + l2_trk_pt);
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
        }

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), false);
        const double trkPur = doTrackPurVar ? GetTrackingPurity (fcal_et, trkpt, trk_eta->at (iTrk), false) : 1.;
        if (trkEff == 0 || trkPur == 0)
          continue;
        const double trkWeight = event_weight * trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, trkWeight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
        }
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
        }
      } // end loop over tracks

    } // end loop over pp tree
    cout << "Done primary pp loop." << endl;
  }

  //CombineHists ();
  //ScaleHists ();
  
  SaveHists (outFileName);

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotFCalDists (const bool _treatAsData) {
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

  if (!_treatAsData) {
    //h_fcal_et->Scale (1./h_fcal_et->Integral ());
    h_fcal_et->SetLineColor (kGray+1);

    h_fcal_et->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
    h_fcal_et->GetYaxis ()->SetTitle ("Counts");

    h_fcal_et->Draw (canvasExists ? "same hist" : "hist");

    h_fcal_et_reweighted->SetLineColor (kBlue);
    //h_fcal_et_reweighted->Scale (1./h_fcal_et_reweighted->Integral ());

    h_fcal_et_reweighted->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
    h_fcal_et_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_fcal_et_reweighted->Draw ("same hist");

    myText (0.72, 0.81, kGray+1, "Unweighted", 0.04);
    myText (0.72, 0.74, kBlue, "Reweighted", 0.04);

    myText (0.22, 0.21, kBlack, "Minimum bias", 0.04);
  }
  else {
    h_fcal_et->SetLineColor (kBlack);
    //h_fcal_et->Scale (1./h_fcal_et->Integral ());

    h_fcal_et->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
    h_fcal_et->GetYaxis ()->SetTitle ("Counts");

    h_fcal_et->Draw (canvasExists ? "same hist" : "hist");

    myText (0.72, 0.88, kBlack, "Z-tagged data", 0.04);
  }

  myText (0.22, 0.28, kBlack, "Pb+Pb, 5.02 TeV", 0.04);

  c->SaveAs (Form ("%s/FCalDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Q2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotQ2Dists (const bool _treatAsData) {
  const char* canvasName = "c_q2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (3, 3);
  }

  c->cd ();

  for (short iCent = 1; iCent < numFinerCentBins; iCent++) {
    c->cd (numFinerCentBins-iCent);
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetDrawnObjects ();
    GetMinAndMax (min, max, true);
    SetMinAndMax (min, max);

    if (!_treatAsData) {

      float max = std::fmax (h_q2[iCent]->GetMaximum (), h_q2_reweighted[iCent]->GetMaximum ());
      float min = std::fmin (h_q2[iCent]->GetMinimum (0), h_q2_reweighted[iCent]->GetMinimum (0));
      TH1D* h = h_q2[iCent];

      h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

      //h->SetLineColor (colors[iCent]);
      h->SetLineColor (kGray+1);

      h->GetXaxis ()->SetTitle ("#it{q}_{2}");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw (!canvasExists ? "hist" : "same hist");

      myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)finerCentCuts[iCent], (int)finerCentCuts[iCent-1]), 0.06);

      h = h_q2_reweighted[iCent];

      h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

      h->SetLineColor (colors[iCent]);
      //h->SetLineStyle (2);

      h->GetXaxis ()->SetTitle ("#it{q}_{2}");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw ("same hist");

      myText (0.61, 0.74, colors[iCent], "Reweighted", 0.05);
      myText (0.61, 0.68, kGray+1, "Unweighted", 0.05);
    }

    else {
      TH1D* h = h_q2[iCent];

      if (h->GetMinimum (0) < min)
        min = h->GetMinimum (0);
      if (h->GetMaximum () > max)
        max = h->GetMaximum ();

      h->SetLineColor (kBlack);

      h->GetXaxis ()->SetTitle ("#it{q}_{2}");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw (!canvasExists ? "hist" : "same hist");

      myText (0.61, 0.80, kBlack, "Z-tagged Data", 0.05);
    }

    SetMinAndMax (0.5*min, 2*max);

  }

  myText (0.52, 0.31, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.52, 0.23, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/Q2Dists.pdf", plotPath.Data ()));
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Plot Q2 weights
//////////////////////////////////////////////////////////////////////////////////////////////////
//void FullAnalysis :: PlotQ2Weights () {
//  const char* canvasName = "c_q2_weights";
//  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
//  TCanvas* c = nullptr;
//  if (canvasExists)
//    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
//  else {
//    c = new TCanvas (canvasName, "", 1200, 1200);
//    gDirectory->Add (c);
//    c->Divide (3,3);
//  }
//
//  for (short iCent = 1; iCent < numFinerCentBins; iCent++) {
//    c->cd (numFinerCentBins-iCent);
//
//    TH1D* h = h_PbPbQ2_weights[iCent][5];
//    h->GetYaxis ()->SetRangeUser (0.5, 1.5);
//
//    h->SetMarkerColor (colors[iCent]);
//    h->SetLineColor (colors[iCent]);
//
//    h->GetXaxis ()->SetTitle ("#left|#it{q}_{2}#right|");
//    h->GetYaxis ()->SetTitle ("Event Weight");
//
//    h->Draw (!canvasExists ? "e1" : "same e1");
//
//    myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)finerCentCuts[iCent], (int)finerCentCuts[iCent-1]), 0.06);
//  }
//
//  c->cd (1);
//
//  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
//  myText (0.22, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
//  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
//  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);
//
//  c->SaveAs (Form ("%s/Q2Weights.pdf", plotPath.Data ()));
//
//}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Psi2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotPsi2Dists (const bool _treatAsData) {
  const char* canvasName = "c_psi2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (3, 3);
  }

  c->cd ();

  for (short iCent = 1; iCent < numFinerCentBins; iCent++) {
    c->cd (numFinerCentBins-iCent);
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetDrawnObjects ();
    GetMinAndMax (min, max, true);
    SetMinAndMax (min, max);

    if (!_treatAsData) {

      float max = std::fmax (h_psi2[iCent]->GetMaximum (), h_psi2_reweighted[iCent]->GetMaximum ());
      float min = std::fmin (h_psi2[iCent]->GetMinimum (0), h_psi2_reweighted[iCent]->GetMinimum (0));
      TH1D* h = h_psi2[iCent];

      h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

      //h->SetLineColor (colors[iCent]);
      h->SetLineColor (kGray+1);

      h->GetXaxis ()->SetTitle ("#Psi_{2}");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw (!canvasExists ? "hist" : "same hist");

      myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)finerCentCuts[iCent], (int)finerCentCuts[iCent-1]), 0.06);

      h = h_psi2_reweighted[iCent];

      h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

      h->SetLineColor (colors[iCent]);
      //h->SetLineStyle (2);

      h->GetXaxis ()->SetTitle ("#Psi_{2}");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw ("same hist");

      myText (0.61, 0.74, colors[iCent], "Reweighted", 0.05);
      myText (0.61, 0.68, kGray+1, "Unweighted", 0.05);
    }

    else {
      TH1D* h = h_psi2[iCent];

      if (h->GetMinimum (0) < min)
        min = h->GetMinimum (0);
      if (h->GetMaximum () > max)
        max = h->GetMaximum ();

      h->SetLineColor (kBlack);

      h->GetXaxis ()->SetTitle ("#Psi_{2}");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw (!canvasExists ? "hist" : "same hist");

      myText (0.61, 0.80, kBlack, "Z-tagged Data", 0.05);
    }

    SetMinAndMax (0.5*min, 2*max);

  }

  myText (0.52, 0.31, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.52, 0.23, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/Psi2Dists.pdf", plotPath.Data ()));
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Plot Psi2 weights
//////////////////////////////////////////////////////////////////////////////////////////////////
//void FullAnalysis :: PlotPsi2Weights () {
//  const char* canvasName = "c_psi2_weights";
//  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
//  TCanvas* c = nullptr;
//  if (canvasExists)
//    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
//  else {
//    c = new TCanvas (canvasName, "", 1200, 1200);
//    gDirectory->Add (c);
//    c->Divide (3,3);
//  }
//
//  for (short iCent = 1; iCent < numFinerCentBins; iCent++) {
//    c->cd (numFinerCentBins-iCent);
//
//    TH1D* h = h_PbPbPsi2_weights[iCent][5];
//    h->GetYaxis ()->SetRangeUser (0.5, 1.5);
//
//    h->SetMarkerColor (colors[iCent]);
//    h->SetLineColor (colors[iCent]);
//
//    h->GetXaxis ()->SetTitle ("#Psi_{2}");
//    h->GetYaxis ()->SetTitle ("Event Weight");
//
//    h->Draw (!canvasExists ? "e1" : "same e1");
//
//    myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)finerCentCuts[iCent], (int)finerCentCuts[iCent-1]), 0.06);
//  }
//
//  c->cd (1);
//
//  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
//  myText (0.22, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
//  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
//  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);
//
//  c->SaveAs (Form ("%s/Psi2Weights.pdf", plotPath.Data ()));
//
//}





////////////////////////////////////////////////////////////////////////////////////////////////
// Plot VZ distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotVZDists (const bool _treatAsData) {
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
void FullAnalysis :: PlotNchDists (const bool _treatAsData) {
  const char* canvasName = "c_nch";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
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

  h_pp_nch->Scale (1./h_pp_nch->Integral ());
  h_pp_nch->GetXaxis ()->SetTitle ("N_{ch}");
  h_pp_nch->GetYaxis ()->SetTitle ("dN_{evt} / N_{evt}");

  h_pp_nch->Draw (canvasExists ? "same hist" : "hist");

  if (!_treatAsData) {
    h_pp_nch_reweighted->Scale (1./h_pp_nch_reweighted->Integral ());
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
// Plot lepton Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonPtSpectra () {
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
    c->SaveAs (Form ("%s/LeptonPtSpectra/%sPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : "Muon"));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot lepton eta spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonEtaSpcComp (FullAnalysis* a) {
  const char* canvasName = "c_lepton_eta_SpcComp";
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
    c->cd ();
  }

  const short iCent = 0;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    TH1D* h = (TH1D*) h_lepton_eta[iCent][iSpc]->Clone ();

    uPad->cd ();
    if (plotFill) {
      h->SetFillColorAlpha (fillColors[iSpc+1], 0.2);
      h->SetLineColor (kBlack);
      h->SetMarkerSize (0);
      h->SetLineWidth (0);
      //h->GetYaxis ()->SetRangeUser (0, 1.3);
      h->GetYaxis ()->SetRangeUser (0, 1.0);

      h->GetXaxis ()->SetTitle ("Lepton #eta");
      //h->GetYaxis ()->SetTitle ("Arb. Units");
      h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{l}/d#eta");
      h->GetXaxis ()->SetTitleSize (0.04/0.6);
      h->GetYaxis ()->SetTitleSize (0.04/0.6);
      h->GetXaxis ()->SetLabelSize (0.04/0.6);
      h->GetYaxis ()->SetLabelSize (0.04/0.6);
      h->GetXaxis ()->SetTitleOffset (1.5*0.6);
      h->GetYaxis ()->SetTitleOffset (1.5*0.6);

      h->DrawCopy (iSpc == 0 && !canvasExists ? "bar" : "bar same");
      h->SetLineWidth (1);
      h->Draw ("hist same");

      gPad->RedrawAxis ();
    }
    else {
      TGraphAsymmErrors* g = GetTGAE (h);
      ResetXErrors (g);
      //deltaize (g, 0.1*(-1.5+iCent));

      const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (1);
      g->SetLineWidth (1);
      g->SetLineColor (colors[iSpc+1]);
      g->SetMarkerColor (colors[iSpc+1]);
      g->GetYaxis ()->SetRangeUser (0, 1.0);

      g->GetXaxis ()->SetTitle ("Lepton #eta");
      g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{l}/d#eta");
      g->GetXaxis ()->SetTitleSize (0.04/0.6);
      g->GetYaxis ()->SetTitleSize (0.04/0.6);
      g->GetXaxis ()->SetLabelSize (0.04/0.6);
      g->GetYaxis ()->SetLabelSize (0.04/0.6);
      g->GetXaxis ()->SetTitleOffset (1.5*0.6);
      g->GetYaxis ()->SetTitleOffset (1.5*0.6);
      g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");
    }

    if (a) {
      dPad->cd ();
      //while (a->h_lepton_eta[iCent][iSpc]->GetNbinsX () >= 2*h->GetNbinsX ()) {
      //  a->h_lepton_eta[iCent][iSpc]->Rebin (2);
      //  a->h_lepton_eta[iCent][iSpc]->Scale (0.5);
      //}
      h->Divide (a->h_lepton_eta[iCent][iSpc]);
      if (h) {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (colors[iSpc+1]);
        g->SetMarkerColor (colors[iSpc+1]);
        //g->GetYaxis ()->SetRangeUser (0.76, 1.24);
        g->GetYaxis ()->SetRangeUser (0.9, 1.1);

        g->GetXaxis ()->SetTitle ("Lepton #eta");
        g->GetYaxis ()->SetTitle ("Data / MC Reco");
        g->GetXaxis ()->SetTitleSize (0.04/0.4);
        g->GetYaxis ()->SetTitleSize (0.04/0.4);
        g->GetXaxis ()->SetLabelSize (0.04/0.4);
        g->GetYaxis ()->SetLabelSize (0.04/0.4);
        g->GetXaxis ()->SetTitleOffset (2.5*0.4);
        g->GetYaxis ()->SetTitleOffset (1.5*0.4);

        g->GetYaxis ()->CenterTitle ();
        g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");

        if (iSpc == 0 && !canvasExists) {
          TLine* l = new TLine (-2.5, 1, 2.5, 1);
          l->SetLineColor (46);
          l->SetLineWidth (2);
          l->SetLineStyle (5);
          l->Draw ("same");
        }
      }
    }
  }
    
  uPad->cd ();
  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
  if (iCent == 0)
    myText (0.22, 0.80, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
  else
    myText (0.22, 0.83-iCent*0.06, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

  myText (0.58, 0.90, kBlack, "Data", 0.032/0.6);
  myText (0.66, 0.90, kBlack, "MC Reco.", 0.032/0.6);
  for (int iSpc = 0; iSpc < 2; iSpc++) {
    const char* spcLabel = iSpc == 0 ? "Z #rightarrow e^{+}e^{-}" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-}" : "Z #rightarrow l^{+}l^{-}");
    myMarkerTextNoLine (0.64, 0.823-iSpc*0.06, colors[iSpc+1], kFullCircle, "", 1.8, 0.04/0.6);
    //myBoxText (0.72, 0.823-iSpc*0.06, colors[iSpc+1], kOpenCircle, "", 1.5, 0.04/0.6);
    myOnlyBoxText (0.75, 0.823-iSpc*0.06, 1.2, fillColors[iSpc+1], kBlack, 1, "", 0.06, 1001);
    myText (0.79, 0.82-iSpc*0.06, colors[iSpc+1], spcLabel, 0.04/0.6);
  }

  c->SaveAs (Form ("%s/LeptonEtaSpcComp.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot track Pt spectra for each lepton species
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonTrackPtSpectra () {
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
    c->SaveAs (Form ("%s/LeptonTrackPtSpectra/%sTrackPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Comb")));
  }

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dR between leptons and tracks
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonTrackDR () {
  
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
      //h2->RebinX (2);
      //h2->RebinY (2);
      h2->GetXaxis ()->SetTitle (Form ("min (#DeltaR (track, %s))", iSpc == 0 ? "electrons" : (iSpc == 1 ? "muons" : "leptons")));
      //h2->GetYaxis ()->SetTitle ("#it{p}_{T}^{h} [GeV]");
      h2->GetYaxis ()->SetTitle ("|#Delta#it{p}_{T}| / <#it{p}_{T}>");
      //h2->GetYaxis ()->SetTitle ("#Delta#phi");
      h2->GetZaxis ()->SetTitle ("Counts");

      h2->GetYaxis ()->SetTitleOffset (1.1);

      h2->Draw ("colz");

      myText (0.56, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);

      if (iCent != 0) {
        myText (0.56, 0.82, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
        myText (0.56, 0.76, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
      }
      else
        myText (0.56, 0.82, kBlack, "#it{pp}, 5.02 TeV", 0.04);

      c->SaveAs (Form ("%s/LeptonTrackDists/%sTrackDist_iCent%i.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Lepton"), iCent));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dR between leptons and tracks
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonTrackDRProjX () {
  
  const char* canvasName = Form ("c_trk_dr");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1400, 600);
    gDirectory->Add (c);
    c->Divide (2, 1);
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    c->cd (iSpc+1);
    //gPad->SetLogy ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {

      TH2D* h2 = h_lepton_trk_dr[iCent][iSpc];

      TH1D* h = h2->ProjectionX ();

      int nEvt = 0;
      for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        nEvt += h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2);
      }
      h->Scale (1./nEvt);
      

      //TH1D* h = h2->ProjectionX ("temp0", h2->GetYaxis ()->FindBin (2), h2->GetYaxis ()->FindBin (3)-1);
      //h->Rebin (2);
      h->GetXaxis ()->SetTitle (Form ("Minimum #DeltaR (track, %s)", iSpc == 0 ? "electrons" : (iSpc == 1 ? "muons" : "leptons")));
      h->GetYaxis ()->SetTitle ("Tracks / Event");
      h->SetLineColor (colors[iCent]);
      h->SetMarkerColor (colors[iCent]);
      h->GetXaxis ()->SetRangeUser (0, 0.16);
      //h->GetYaxis ()->SetRangeUser (0.5, 0.5e4);
      h->GetYaxis ()->SetRangeUser (0, 1);
      h->GetYaxis ()->SetTitleOffset (1.1);

      h->DrawCopy (iCent == 0 ? "e1" : "e1 same");
      delete h;

      myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);

      if (iCent != 0) {
        myText (0.22, 0.80-iCent*0.06, colors[iCent], Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
      }
      else
        myText (0.22, 0.80, kBlack, "#it{pp}, 5.02 TeV", 0.04);

    }
  }
  c->SaveAs (Form ("%s/LeptonTrackDists/LeptonTrackDist_ProjX.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates ratio of Z pT distributions between data and MC
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CalculateZPtDistRatio (FullAnalysis* a) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      h_z_pt_ratio[iCent][iSpc] = (TH1D*) h_z_pt[iCent][iSpc]->Clone (Form ("h_z_pt_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt_ratio[iCent][iSpc]->Divide (a->h_z_pt[iCent][iSpc]);
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZPtSpectra () {
  //const char* canvasName = "c_z_pt";
  //const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  //TCanvas* c = nullptr;
  ////TPad* uPad = nullptr, *dPad = nullptr;
  //if (canvasExists) {
  //  c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  //  //uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
  //  //dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
  //}
  //else {
  //  c = new TCanvas (canvasName, "", 800, 800);
  //  c->cd ();
  //  //uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
  //  //dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
  //  //uPad->SetBottomMargin (0);
  //  //dPad->SetTopMargin (0);
  //  //dPad->SetBottomMargin (0.25);
  //  //uPad->Draw ();
  //  //dPad->Draw ();
  //  gDirectory->Add (c);
  //  //gDirectory->Add (uPad);
  //  //gDirectory->Add (dPad);
  //}
  //c->cd ();
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    const char* canvasName = Form ("c_z_pt_iCent%i", iCent);
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
    c->cd ();

    uPad->cd ();
    gPad->SetLogy ();

    TH1D* h = h_z_pt[iCent][2];

    if (plotFill) {
      h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
      h->SetLineColor (kBlack);
      h->SetMarkerSize (0);
      h->SetLineWidth (0);
      h->GetYaxis ()->SetRangeUser (1e-6, 0.06);

      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");
      h->GetYaxis ()->SetTitle ("1/N_{Z} dN/d#it{p}_{T} [GeV^{-1}]");
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
    } else {
      TGraphAsymmErrors* g = GetTGAE (h);
      ResetXErrors (g);

      const int markerStyle = kFullCircle;
      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (1.25);
      g->SetLineWidth (2);
      g->SetLineColor (kBlack);
      g->SetMarkerColor (kBlack);
      //g->SetLineColor (colors[iCent]);
      //g->SetMarkerColor (colors[iCent]);
      g->GetYaxis ()->SetRangeUser (1e-6, 0.06);

      g->GetXaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
      g->GetYaxis ()->SetTitle ("1/N_{Z} dN/d#it{p}_{T} [GeV^{-1}]");
      g->GetXaxis ()->SetTitleSize (0.04/0.6);
      g->GetYaxis ()->SetTitleSize (0.04/0.6);
      g->GetXaxis ()->SetLabelSize (0.04/0.6);
      g->GetYaxis ()->SetLabelSize (0.04/0.6);
      g->GetXaxis ()->SetTitleOffset (1.5*0.6);
      g->GetYaxis ()->SetTitleOffset (1.5*0.6);
      g->Draw (!canvasExists/* && iCent == 0*/ ? "AP" : "P");
    }

    dPad->cd ();
    h = h_z_pt_ratio[iCent][2];
    if (h) {
      TGraphAsymmErrors* g = GetTGAE (h);
      ResetXErrors (g);

      const int markerStyle = kFullCircle;
      g->SetMarkerStyle (markerStyle);
      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (1);
      g->SetLineWidth (1);
      g->SetLineColor (kBlack);
      g->SetMarkerColor (kBlack);
      g->GetYaxis ()->SetRangeUser (0.5, 1.5);

      g->GetXaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
      g->GetYaxis ()->SetTitle ("Data / MC");
      g->GetXaxis ()->SetTitleSize (0.04/0.4);
      g->GetYaxis ()->SetTitleSize (0.04/0.4);
      g->GetXaxis ()->SetLabelSize (0.04/0.4);
      g->GetYaxis ()->SetLabelSize (0.04/0.4);
      g->GetXaxis ()->SetTitleOffset (2.5*0.4);
      g->GetYaxis ()->SetTitleOffset (1.5*0.4);
      g->GetYaxis ()->CenterTitle ();
      g->Draw (!canvasExists/* && iCent == 0*/ ? "AP" : "P");
    }
    else {
      cout << "Warning in FullAnalysis :: PlotZPtSpectra: Z pT spectra ratio not stored, needs to be calculated!" << endl;
    }

    uPad->cd ();

    myText (0.66, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/0.6);
    myText (0.26, 0.85, kBlack, "Z #rightarrow l^{+}l^{-} Events", 0.04/0.6);
    myMarkerText (0.753, 0.65, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
    myOnlyBoxText (0.76, 0.55, 1.2, fillColors[iCent], kBlack, 1, "MC", 0.04/0.6, 1001);

    if (iCent == 0)
      myText (0.66, 0.75, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
    else
      myText (0.66, 0.75, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

    c->SaveAs (Form ("%s/ZPtSpectra/z_pt_spectrum_iCent%i.pdf", plotPath.Data (), iCent));
  }
  //c->SaveAs (Form ("%s/ZPtSpectra/z_pt_spectrum.pdf", plotPath.Data ()));
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Prints yield of Z's that meet the event selection criteria in each centrality
//////////////////////////////////////////////////////////////////////////////////////////////////
//void FullAnalysis :: PrintZYields () {
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    const char* spc = (iSpc == 0 ? "Z->ee" : (iSpc == 1 ? "Z->mumu" : "Z->ll"));
//    for (short iCent = 0; iCent < numCentBins; iCent++) {
//      float yield = h_z_counts[iSpc][2][iCent]->GetBinContent (2);
//      if (iCent == 0) 
//        cout << "pp " << spc << " # Z's > 25 GeV  =  " << yield << endl;
//      else
//        cout << Form ("Pb+Pb %i-%i%% ", (int)centCuts[iCent], (int)centCuts[iCent-1]) << spc << " # Z's > 25 GeV  =  " << yield << endl;
//    }
//  }
//}





////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the RAPIDITY-phi distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZYPhiMap (const short pPtZ) {
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    if (pPtZ != -1 && pPtZ != iPtZ)
      continue;

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

        const char* canvasName = Form ("c_z_y_phi_%s_iCent%i_iPtZ%i", spc, iCent, iPtZ);
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

        TH2D* h = h_z_y_phi[iCent][iSpc][iPtZ];

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
          myText (0.18, 0.84, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04);
        }
        else
          myText (0.18, 0.84, kBlack, Form ("Pb+Pb %i-%i%%, 5.02 TeV", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);

        c->SaveAs (Form ("%s/ZYPhiDists/z%s_y_phi_iCent%i_iPtZ%i.pdf", plotPath.Data (), spc, iCent, iPtZ));
      }
    }
  }

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates ratio of Z eta distributions between central Pb+Pb to pp
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CalculateZEtaDistRatio () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_eta_ratio[iCent][iSpc][iPtZ] = (TH1D*) h_z_eta[iCent][iSpc][iPtZ]->Clone (Form ("h_z_eta_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_eta_ratio[iCent][iSpc][iPtZ]->Divide (h_z_eta[0][iSpc][iPtZ]);
      }
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the PSEUDORAPIDITY distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZEtaMap (const short pPtZ) {
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    if (pPtZ != -1 && pPtZ != iPtZ)
      continue;

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      const char* canvasName = Form ("c_z_eta_%s_iPtZ%i", spc, iPtZ);
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
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        TH1D* h = h_z_eta[iCent][iSpc][iPtZ];

        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (colors[iCent]);
        g->SetMarkerColor (colors[iCent]);
        g->GetYaxis ()->SetRangeUser (0, 0.22);

        g->GetXaxis ()->SetTitle ("#eta");
        g->GetYaxis ()->SetTitle ("1/N_{Z} dN/d#eta");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (!canvasExists && iCent == 0 ? "AP" : "P");

        myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/0.6);

        const char* spcLabel = iSpc == 0 ? "Z #rightarrow e^{+}e^{-} Events" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-} Events" : "Z #rightarrow l^{+}l^{-} Events");
        myText (0.66, 0.85, kBlack, spcLabel, 0.04/0.6);
        if (iPtZ == nPtZBins-1)
          myText (0.66, 0.75, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04/0.6);
        else
          myText (0.66, 0.75, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04/0.6);
        if (iCent == 0)
          myText (0.28, 0.26, kBlack, Form ("#it{pp}"), 0.03/0.6);
        else
          myText (0.28, 0.26-iCent*0.06, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.03/0.6);
      }

      dPad->cd ();
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        TH1D* h = h_z_eta_ratio[iCent][iSpc][iPtZ];
        if (h) {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (colors[iCent]);
          g->SetMarkerColor (colors[iCent]);
          g->GetYaxis ()->SetRangeUser (0.76, 1.24);

          g->GetXaxis ()->SetTitle ("#eta");
          g->GetYaxis ()->SetTitle ("Pb+Pb / #it{pp}");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);

          g->GetYaxis ()->CenterTitle ();
          g->Draw (!canvasExists && iCent == 1 ? "AP" : "P");
        }
        else {
          cout << "Warning in FullAnalysis :: PlotZEtaDst: Z Eta spectra ratio not stored, needs to be calculated!" << endl;
        }
      }

      c->SaveAs (Form ("%s/ZEtaDists/z%s_eta_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates ratio of Z y distributions between central Pb+Pb to pp
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CalculateZYDistRatio (FullAnalysis* a) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_ratio[iCent][iSpc][iPtZ] = (TH1D*) h_z_y[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_y_ratio[iCent][iSpc][iPtZ]->Rebin (2);
        h_z_y_ratio[iCent][iSpc][iPtZ]->Scale (0.5);
        TH1D* temp = (TH1D*)a->h_z_y[iCent][iSpc][iPtZ]->Clone ("temp");
        temp->Rebin (2);
        temp->Scale (1./temp->Integral (), "width");
        h_z_y_ratio[iCent][iSpc][iPtZ]->Divide (temp);
        delete temp;
      }
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the PSEUDORAPIDITY distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZYMap (const short pPtZ) {
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    if (pPtZ != -1 && pPtZ != iPtZ)
      continue;

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      const char* canvasName = Form ("c_z_y_%s_iPtZ%i", spc, iPtZ);
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
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        TH1D* h = h_z_y[iCent][iSpc][iPtZ];

        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (colors[iCent]);
        g->SetMarkerColor (colors[iCent]);
        g->GetYaxis ()->SetRangeUser (0, 0.4);

        g->GetXaxis ()->SetTitle ("y_{Z}");
        g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/dy_{Z}");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (iCent == 0 ? "AP" : "P");

        myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/0.6);

        const char* spcLabel = iSpc == 0 ? "Z #rightarrow e^{+}e^{-} Events" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-} Events" : "Z #rightarrow l^{+}l^{-} Events");
        myText (0.66, 0.85, kBlack, spcLabel, 0.04/0.6);
        if (iPtZ == nPtZBins-1)
          myText (0.66, 0.75, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04/0.6);
        else
          myText (0.66, 0.75, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04/0.6);
        if (iCent == 0)
          myText (0.28, 0.26, kBlack, Form ("#it{pp}"), 0.03/0.6);
        else
          myText (0.28, 0.26-iCent*0.06, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.03/0.6);
      }

      dPad->cd ();
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        TH1D* h = h_z_y_ratio[iCent][iSpc][iPtZ];
        if (h) {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (colors[iCent]);
          g->SetMarkerColor (colors[iCent]);
          g->GetYaxis ()->SetRangeUser (0.76, 1.24);
          //g->GetYaxis ()->SetRangeUser (0, 2);

          g->GetXaxis ()->SetTitle ("y_{Z}");
          g->GetYaxis ()->SetTitle ("Pb+Pb / #it{pp}");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);

          g->GetYaxis ()->CenterTitle ();
          g->Draw (iCent == 1 ? "AP" : "P");
        }
        else {
          cout << "Warning in FullAnalysis :: PlotZEtaDst: Z Eta spectra ratio not stored, needs to be calculated!" << endl;
        }
      }

      c->SaveAs (Form ("%s/ZYDists/z%s_y_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the rapidity distribution of reconstructed Z->ee vs Z->mumu bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZYMapSpcComp (const short pPtZ, FullAnalysis* a) {
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    if (pPtZ != -1 && pPtZ != iPtZ)
      continue;

    const char* canvasName = Form ("c_z_y_SpcComp_iPtZ%i", iPtZ);
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
      c->cd ();
    }

    const short iCent = 0;

    for (short iSpc = 0; iSpc < 2; iSpc++) {
      //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      TH1D* h = (TH1D*) h_z_y[iCent][iSpc][iPtZ]->Clone ();
      //h->Rebin (2);
      //h->Scale (0.5);

      uPad->cd ();
      if (plotFill) {
        h->SetFillColorAlpha (fillColors[iSpc+1], 0.2);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);
        //h->GetYaxis ()->SetRangeUser (0, 1.3);
        h->GetYaxis ()->SetRangeUser (0, 0.5);

        h->GetXaxis ()->SetTitle ("y_{Z}");
        h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/dy_{Z}");
        h->GetXaxis ()->SetTitleSize (0.04/0.6);
        h->GetYaxis ()->SetTitleSize (0.04/0.6);
        h->GetXaxis ()->SetLabelSize (0.04/0.6);
        h->GetYaxis ()->SetLabelSize (0.04/0.6);
        h->GetXaxis ()->SetTitleOffset (1.5*0.6);
        h->GetYaxis ()->SetTitleOffset (1.5*0.6);

        if (!useAltMarker) {
          h->DrawCopy (iSpc == 0 && !canvasExists ? "bar" : "bar same");
          h->SetLineWidth (1);
        }
        else {
          h->SetLineStyle (2);
          h->SetLineColor (colors[iSpc+1]);
          h->SetFillColorAlpha (fillColors[iSpc+1], 0);
          h->SetLineWidth (2);
        }
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      }
      else {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (colors[iSpc+1]);
        g->SetMarkerColor (colors[iSpc+1]);
        g->GetYaxis ()->SetRangeUser (0, 0.5);

        g->GetXaxis ()->SetTitle ("y_{Z}");
        g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/dy_{Z}");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");
      }

      if (a) {
        dPad->cd ();
        //while (a->h_z_y[iCent][iSpc][iPtZ]->GetNbinsX () >= 2*h->GetNbinsX ()) {
        //  a->h_z_y[iCent][iSpc][iPtZ]->Rebin (2);
        //  a->h_z_y[iCent][iSpc][iPtZ]->Scale (0.5);
        //}
        h->Divide (a->h_z_y[iCent][iSpc][iPtZ]);

        dPad->cd ();
        if (h) {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);

          const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (colors[iSpc+1]);
          g->SetMarkerColor (colors[iSpc+1]);
          g->GetYaxis ()->SetRangeUser (0.76, 1.24);
          //g->GetYaxis ()->SetRangeUser (0, 2);

          g->GetXaxis ()->SetTitle ("y_{Z}");
          g->GetYaxis ()->SetTitle ("Data / MC Truth");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);

          g->GetYaxis ()->CenterTitle ();
          g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");

          if (iSpc == 0 && !canvasExists) {
            TLine* l = new TLine (-2.5, 1, 2.5, 1);
            l->SetLineColor (46);
            l->SetLineWidth (2);
            l->SetLineStyle (5);
            l->Draw ("same");
          }
        }
      }
      else {
        cout << "Warning in FullAnalysis :: PlotZYMapSpcComp: Z y spectra ratio not stored, needs to be calculated!" << endl;
      }
    }
      
    uPad->cd ();
    myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
    if (iCent == 0)
      myText (0.22, 0.82, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
    else
      myText (0.22, 0.82, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);
    if (iPtZ == nPtZBins-1)
      myText (0.22, 0.74, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04/0.6);
    else
      myText (0.22, 0.74, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04/0.6);

    myText (0.58, 0.90, kBlack, "Data", 0.032/0.6);
    myText (0.66, 0.90, kBlack, "MC Reco.", 0.032/0.6);
    for (int iSpc = 0; iSpc < 2; iSpc++) {
      const char* spcLabel = iSpc == 0 ? "Z #rightarrow e^{+}e^{-}" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-}" : "Z #rightarrow l^{+}l^{-}");
      myMarkerTextNoLine (0.64, 0.823-iSpc*0.06, colors[iSpc+1], kFullCircle, "", 1.8, 0.04/0.6);
      //myBoxText (0.72, 0.823-iSpc*0.06, colors[iSpc+1], kOpenCircle, "", 1.5, 0.04/0.6);
      myOnlyBoxText (0.75, 0.823-iSpc*0.06, 1.2, fillColors[iSpc+1], kBlack, 1, "", 0.06, 1001);
      myText (0.79, 0.82-iSpc*0.06, colors[iSpc+1], spcLabel, 0.04/0.6);
    }

    c->SaveAs (Form ("%s/ZYDists/z_y_iPtZ%i.pdf", plotPath.Data (), iPtZ));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates ratio of Z mass spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CalculateZMassSpectraRatio (FullAnalysis* a) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      //if (iCent > 0) {
      //  h_z_m[iCent][iSpc]->Rebin (2);
      //  h_z_m[iCent][iSpc]->Scale (0.5);
      //  a->h_z_m[iCent][iSpc]->Rebin (2);
      //  a->h_z_m[iCent][iSpc]->Scale (0.5);
      //}

      for (short iReg = 0; iReg < 3; iReg++) {

        h_z_m_ratio[iCent][iSpc][iReg] = (TH1D*) h_z_m[iCent][iSpc][iReg]->Clone (Form ("h_z_m_ratio_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
        h_z_m_ratio[iCent][iSpc][iReg]->Divide (a->h_z_m[iCent][iSpc][iReg]);
      }
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z mass spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZMassSpectra () {
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    for (short iReg = 0; iReg < 2; iReg++) {
      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        const char* canvasName = Form ("c_z_m_%s_iCent%i_iReg%i", spc, iCent, iReg);
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
        TH1D* h = h_z_m[iCent][iSpc][iReg];
        if (plotFill) {
          h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          //h->GetYaxis ()->SetRangeUser (0, 1.3);
          h->GetYaxis ()->SetRangeUser (0, 0.12);

          h->GetXaxis ()->SetTitle (Form ("m_{%s} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          //h->GetYaxis ()->SetTitle ("Arb. Units");
          h->GetYaxis ()->SetTitle ("Counts / Total");
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
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);
          //deltaize (g, 0.1*(-1.5+iCent));

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (kBlack);
          g->SetMarkerColor (kBlack);
          //g->GetYaxis ()->SetRangeUser (0, 1.3);
          g->GetYaxis ()->SetRangeUser (0, 0.12);

          g->GetXaxis ()->SetTitle (Form ("m_{%s} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          //g->GetYaxis ()->SetTitle ("Arb. Units");
          g->GetYaxis ()->SetTitle ("Counts / Total");
          g->GetXaxis ()->SetTitleSize (0.04/0.6);
          g->GetYaxis ()->SetTitleSize (0.04/0.6);
          g->GetXaxis ()->SetLabelSize (0.04/0.6);
          g->GetYaxis ()->SetLabelSize (0.04/0.6);
          g->GetXaxis ()->SetTitleOffset (1.5*0.6);
          g->GetYaxis ()->SetTitleOffset (1.5*0.6);
          g->Draw (!canvasExists ? "AP" : "P");
        }
        LabelZMassSpectra (iSpc, iCent, iReg);
        

        dPad->cd ();
        h = h_z_m_ratio[iCent][iSpc][iReg];
        if (h) {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (kBlack);
          g->SetMarkerColor (kBlack);
          g->GetYaxis ()->SetRangeUser (0.1, 2.4);

          g->GetXaxis ()->SetTitle (Form ("m_{%s} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          g->GetYaxis ()->SetTitle ("Data / MC");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);
          g->GetYaxis ()->CenterTitle ();
          g->Draw (!canvasExists ? "AP" : "P");

          if (!canvasExists) {
            TLine* l = new TLine (76, 1, 106, 1);
            l->SetLineColor (46);
            l->SetLineWidth (2);
            l->SetLineStyle (5);
            l->Draw ("same");
          }
        }
        else {
          cout << "Warning in FullAnalysis :: PlotZMassSpectra: Z mass spectra ratio not stored, needs to be calculated!" << endl;
        }

        if (iReg == 2) 
          c->SaveAs (Form ("%s/ZMassSpectra/z%s_mass_spectrum_iCent%i.pdf", plotPath.Data (), spc, iCent));
        else
          c->SaveAs (Form ("%s/ZMassSpectra/z%s_mass_spectrum_iCent%i_iReg%i.pdf", plotPath.Data (), spc, iCent, iReg));
      }
    }
  }
}



void FullAnalysis :: LabelZMassSpectra (const short iSpc, const short iCent, const short iReg) {
  myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/0.6);
  const char* spc = iSpc == 0 ? "Z #rightarrow e^{+}e^{-}" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-}" : "Z #rightarrow l^{+}l^{-}");
  myText (0.71, 0.85, kBlack, spc, 0.04/0.6);
  if (iCent == 0) {
    myText (0.22, 0.76, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
  }
  else
    myText (0.22, 0.76, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

  myOnlyBoxText (0.76, 0.67, 1.2, fillColors[iCent], kBlack, 1, "MC", 0.04/0.6, 1001);

  //TVirtualPad* cPad = gPad; // store current pad
  //TBox* b = TBoxNDC (0.4+0.6*(0.598-0.025), 0.67-0.06*numPhiBins-0.018, 0.4+0.6*(0.598+0.025), 0.67-0.06*numPhiBins+0.018);
  //b->SetFillColorAlpha (fillColors[iCent], fillAlpha);
  //b->Draw ("l");
  //cPad->cd ();
  //myText (0.753, 0.67, kBlack, "MC", 0.04/0.6);
  myMarkerText (0.753, 0.76, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);

  if (iReg == 0)
    myText (0.22, 0.67, kBlack, Form ("#left|y^{%s}#right| < 1", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))), 0.04/0.6);
  else if (iReg == 1)
    myText (0.22, 0.67, kBlack, Form ("#left|y^{%s}#right| > 1", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))), 0.04/0.6);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZPhiYield () {
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
    float v2 = 0, v2err = 0;
    float c = 0, cerr = 0;

    for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
      c += h->GetBinContent (ix) * h->GetBinWidth (ix);
      cerr += pow (h->GetBinError (ix) * h->GetBinWidth (ix), 2);

      v2 += h->GetBinContent (ix) * cos (h->GetBinCenter (ix)) * h->GetBinWidth (ix);
      v2err += pow (h->GetBinError (ix) * cos (h->GetBinCenter (ix)) * h->GetBinWidth (ix), 2);
    }
    cerr = sqrt (cerr);
    v2err = sqrt (v2err);
    v2 = v2 / (c);
    v2err = sqrt (pow (v2err / c, 2) + pow (v2 * cerr / (c*c), 2));
    c = c / pi;
    cerr = cerr / pi;

    TGraphAsymmErrors* g = GetTGAE (h);
    deltaize (g, (1.5-iCent)*0.02, false);

    g->GetXaxis ()->SetTitle ("2#left|#phi_{Z} - #Psi_{2}#right|");
    g->GetYaxis ()->SetTitle ("1/N_{Z} dN/d#Delta#phi");

    g->GetYaxis ()->SetRangeUser (0.16, 0.5);

    g->SetLineColor (colors[iCent]);
    g->SetMarkerColor (colors[iCent]);
    g->SetMarkerSize (1);
    if (iCent == 1)
      g->Draw ("AP");
    else
      g->Draw ("P");

    TF1* fit = new TF1 (Form ("fit_iCent%i", iCent), "[0]*(1+2*[1]*cos(x))", 0, pi);
    fit->FixParameter (0, c);
    fit->FixParameter (1, v2);
    fit->SetLineColor (colors[iCent]);
    fit->SetLineStyle (2);
    fit->SetLineWidth (2);
    fit->Draw ("same");

    myText (0.22, 0.38-0.055*iCent, colors[iCent], Form ("v_{2} = %s", FormatMeasurement (v2, v2err, 2)), 0.04);
    myText (0.66, 0.88-0.055*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);

  } // end loop over cents

  //myText (0.66, 0.88, colors[0], "#it{pp}", 0.04);
  myText (0.25, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.056);
  myText (0.25, 0.80, kBlack, "#it{p}_{T}^{Z} > 5 GeV", 0.05);
  myText (0.66, 0.88, kBlack, "Pb+Pb, 5.02 TeV", 0.04);

  c->SaveAs (Form ("%s/ZPhiYields.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZLeptonDPhi () {
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
      //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

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


#endif
