#include "Intercalibration.h"
#include "Params.h"
#include "CalibUtils.h"
#include "TreeVariables.h"

#include <Utilities.h>

#include <TH2D.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

int nTruthPhotons = 0;

/**
 * Returns true iff the object at eta, phi is truth matched to a photon within dR_thres.
 */
bool isTruthPhoton (const int n_t_p, // number of truth photons
                    const vector<float>* t_p_pt, // truth photon pts
                    const vector<float>* t_p_eta, // truth photon etas
                    const vector<float>* t_p_phi, // truth photon phis
                    const float eta,
                    const float phi,
                    const float dR_thres = 0.3)
{
  int l_t_p = -1;
  for (int tp = 0; tp < n_t_p; tp++) {
   if (l_t_p == -1 || t_p_pt->at (l_t_p) < t_p_pt->at (tp))
    l_t_p = tp;
  }
  if (l_t_p == -1)
   return false; // return false iff there are no truth photons

  const float dR = DeltaR (t_p_eta->at (l_t_p), eta, t_p_phi->at (l_t_p), phi);

  if (dR < dR_thres)
   nTruthPhotons++;

  return dR < dR_thres;
}


/**
 * Returns true iff the object at eta, phi is truth matched to an electron within dR_thres.
 */
bool isTruthElectron (const int n_t_e, // number of truth electrons
                    const vector<float>* t_e_pt, // truth electron pts
                    const vector<float>* t_e_eta, // truth electron etas
                    const vector<float>* t_e_phi, // truth electron phis
                    const float eta,
                    const float phi,
                    const float dR_thres = 0.3)
{
  int l_t_e = -1;
  for (int te = 0; te < n_t_e; te++) {
   if (l_t_e == -1 || t_e_pt->at (l_t_e) < t_e_pt->at (te))
    l_t_e = te;
  }
  if (l_t_e == -1)
   return false; // return false iff there are no truth electrons

  const float dR = DeltaR (t_e_eta->at (l_t_e), eta, t_e_phi->at (l_t_e), phi);

  return dR < dR_thres;
}


/**
 * Main analysis routine.
 */
void Intercalibration (const char* directory,
                       const int dataSet,
                       const bool isPeriodA,
                       const double pt_low,
                       const double pt_high,
                       const char* inFileName,
                       const bool cutTruthEgammas,
                       const double crossSection_microbarns,
                       const double filterEfficiency,
                       const int numberEvents)
{

  SetupDirectories ("", "JetCalibration/");

  const bool isSignalOnlySample = TString (inFileName).Contains ("signalonly");
  const TString identifier = GetIdentifier (dataSet, inFileName, true, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, true, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In Intercalibration.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, true);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets (false); // make sure we are getting the full jet collections
  t->SetGetHIJets ();
  t->SetGetTruthJets ();
  t->SetGetTruthElectrons ();
  t->SetGetTruthPhotons ();
  t->SetBranchAddresses ();

  int jet_n;
  vector<float> *calib_jet_pt, *reco_jet_pt; // calib is at the EtaJES scale in MC, reco is at the EM scale
  vector<float> *calib_jet_eta, *reco_jet_eta;
  vector<float> *calib_jet_phi, *reco_jet_phi;
  vector<float> *calib_jet_e, *reco_jet_e;

  // initialize histograms
  TH3D* jetRecoRespHist = new TH3D (Form ("jetRecoResp_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
  TH3D* jetCalibRespHist = new TH3D (Form ("jetCalibResp_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
  TH2D* jetRecoRespCounts = new TH2D (Form ("jetRecoCounts_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins);
  TH2D* jetCalibRespCounts = new TH2D (Form ("jetCalibCounts_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins);

  jetRecoRespHist->Sumw2 ();
  jetCalibRespHist->Sumw2 ();
  jetRecoRespCounts->Sumw2 ();
  jetCalibRespCounts->Sumw2 ();


  TH1D* truthPhotonsPerEvent = new TH1D (Form ("truthPhotonsPerEvent_%s", identifier.Data ()), "", 9, -0.5, 8.5);
  TH1D* truthRatePerEvent = new TH1D (Form ("truthRatePerEvent_%s", identifier.Data ()), "", 50, 0, 1);
  truthPhotonsPerEvent->Sumw2 ();
  truthRatePerEvent->Sumw2();
  
  const int numEntries = tree->GetEntries ();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (int entry = 0; entry < numEntries; entry++) {
   tree->GetEntry (entry);


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;


   jet_n = *(& (t->akt4hi_jet_n));
   calib_jet_pt = t->akt4hi_em_xcalib_jet_pt;
   reco_jet_pt = t->akt4hi_em_jet_pt;
   calib_jet_eta = t->akt4hi_em_xcalib_jet_eta;
   reco_jet_eta = t->akt4hi_em_jet_eta;
   calib_jet_phi = t->akt4hi_em_xcalib_jet_phi;
   reco_jet_phi = t->akt4hi_em_jet_phi;
   calib_jet_e = t->akt4hi_em_xcalib_jet_e;
   reco_jet_e = t->akt4hi_em_jet_e;

   const float weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;


   /////////////////////////////////////////////////////////////////////////////
   // main calibrated jet loop
   /////////////////////////////////////////////////////////////////////////////
   for (int j = 0; j < jet_n; j++) {
    const float jpt = calib_jet_pt->at (j);
    const float jeta = calib_jet_eta->at (j);
    const float jphi = calib_jet_phi->at (j);
    const float je = calib_jet_e->at (j);

    if (!InHadCal (jeta, 0.4))
     continue; // require jets to be completely inside the hadronic calorimeter
    if (InDisabledHEC (jeta, jphi))
     continue; // Reject event on additional HEC cuts
    if (cutTruthEgammas && isTruthPhoton (t->truth_photon_n, t->truth_photon_pt, t->truth_photon_eta, t->truth_photon_phi, jeta, jphi))
     continue; // If rejecting truth egammas, calculate dR isolation from truth egammas, rejecting if threshold is met.

    // truth match this jet
    double minDeltaR = 1000;
    int truth_jet = -1;
    for (int tj = 0; tj < t->truth_jet_n; tj++) {
     double dR = DeltaR (jeta, t->truth_jet_eta->at (tj), jphi, t->truth_jet_phi->at (tj));
     if (dR < minDeltaR) {
      minDeltaR = dR;
      truth_jet = tj;
     }
    }

    if (0 <= truth_jet && truth_jet < t->truth_jet_n && minDeltaR < 0.2) {
     if (t->truth_jet_pt->at (truth_jet) < pt_low || pt_high < t->truth_jet_pt->at (truth_jet))
      continue; // Only look at jets which are within the pT bounds of this MC sample 

     double ratio = 0;
     double xval = 0;
     if (calcPtClosure) {
      ratio = jpt / t->truth_jet_pt->at (truth_jet);
      xval = t->truth_jet_pt->at (truth_jet);
     } 
     else {
      ratio = je / t->truth_jet_e->at (truth_jet);
      xval = t->truth_jet_e->at (truth_jet);
     }
     
     jetCalibRespHist->Fill (xval, jeta, ratio, weight);
     jetCalibRespCounts->Fill (xval, jeta);
    }
    
   } // end calibrated jet loop

   nTruthPhotons = 0;

   /////////////////////////////////////////////////////////////////////////////
   // main reconstructed jet loop
   /////////////////////////////////////////////////////////////////////////////
   for (int j = 0; j < jet_n; j++) {
    const float jpt = reco_jet_pt->at (j);
    const float jeta = reco_jet_eta->at (j);
    const float jphi = reco_jet_phi->at (j);
    const float je = reco_jet_e->at (j);

    if (!InHadCal (jeta, 0.4))
     continue; // require jets to be completely inside the hadronic calorimeter
    if (InDisabledHEC (jeta, jphi))
     continue; // Reject event on additional HEC cuts
    if (cutTruthEgammas && isTruthPhoton (t->truth_photon_n, t->truth_photon_pt, t->truth_photon_eta, t->truth_photon_phi, jeta, jphi))
     continue; // If rejecting truth egammas, calculate dR isolation from truth egammas, rejecting if threshold is met.

    // truth match this jet
    double minDeltaR = 1000;
    int truth_jet = -1;
    for (int tj = 0; tj < t->truth_jet_n; tj++) {
     double dR = DeltaR (jeta, t->truth_jet_eta->at (tj), jphi, t->truth_jet_phi->at (tj));
     if (dR < minDeltaR) {
      minDeltaR = dR;
      truth_jet = tj;
     }
    }

    if (0 <= truth_jet && truth_jet < t->truth_jet_n && minDeltaR < 0.2) {
     if (t->truth_jet_pt->at (truth_jet) < pt_low || pt_high < t->truth_jet_pt->at (truth_jet))
      continue; // Only look at jets which are within the pT bounds of this MC sample 

     double ratio = 0;
     double xval = 0;
     if (calcPtClosure) {
      ratio = jpt / t->truth_jet_pt->at (truth_jet);
      xval = t->truth_jet_pt->at (truth_jet);
     } 
     else {
      ratio = je / t->truth_jet_e->at (truth_jet);
      xval = t->truth_jet_e->at (truth_jet);
     }

     jetRecoRespHist->Fill (xval, jeta, ratio, weight);
     jetRecoRespCounts->Fill (xval, jeta);
    }
   } // end reconstructed jet loop

   truthPhotonsPerEvent->Fill (nTruthPhotons);
   truthRatePerEvent->Fill ((float)nTruthPhotons / (float)jet_n);
   nTruthPhotons = 0;

  } // end loop over events

  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////


  const char* outFileName = Form ("%s/Intercalibration/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  // Write histograms to output and clean memory
  jetCalibRespHist->Write ();
  jetCalibRespCounts->Write ();
  jetRecoRespHist->Write ();
  jetRecoRespCounts->Write ();

  truthPhotonsPerEvent->Write ();
  truthRatePerEvent->Write ();

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
