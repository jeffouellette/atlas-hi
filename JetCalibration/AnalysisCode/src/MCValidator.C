#include "MCValidator.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <Trigger.h>
#include <TreeVariables.h>
#include <ArrayTemplates.h>

#include <TH2D.h>
#include <TH3D.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {


/**
 * Main analysis routine.
 */
void MCValidator (const char* directory,
                  const int dataSet,
                  const bool isPeriodA,
                  const double pt_low,
                  const double pt_high,
                  const char* inFileName,
                  const double crossSection_microbarns,
                  const double filterEfficiency,
                  const int numberEvents)
{
  const int numetabins_mc = 26;
  const double etabins_mc[27] = {-4.5, -4.0, -3.6, -3.2, -2.8, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5};

  SetupDirectories ("", "JetCalibration/");

  const bool isSignalOnlySample = TString (inFileName).Contains ("signalonly");
  const TString identifier = GetIdentifier (dataSet, inFileName, true, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, true, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In MCValidator.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  int jet_n;
  vector<float> *calib_jet_pt = NULL, *reco_jet_pt = NULL; // calib is at the EtaJES scale in MC, reco is at the EM scale
  vector<float> *calib_jet_eta = NULL, *reco_jet_eta = NULL;
  vector<float> *calib_jet_phi = NULL, *reco_jet_phi = NULL;
  vector<float> *calib_jet_e = NULL, *reco_jet_e = NULL;
  vector<vector<double>> *jet_sampling = NULL;

  TreeVariables* t = new TreeVariables (tree, true);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets (false); // make sure we are getting the full jet collections
  t->SetGetHIJets ();
  t->SetGetTruthJets ();
  t->SetGetTruthElectrons ();
  t->SetGetTruthPhotons ();
  tree->SetBranchAddress ("akt4hi_sampling", &jet_sampling);
  t->SetBranchAddresses ();

  // initialize histograms
  TH3D* jetRecoRespHist = new TH3D (Form ("jetRecoResp_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
  TH3D* jetCalibRespHist = new TH3D (Form ("jetCalibResp_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
  TH2D* jetRecoRespCounts = new TH2D (Form ("jetRecoCounts_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins);
  TH2D* jetCalibRespCounts = new TH2D (Form ("jetCalibCounts_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins);

  jetRecoRespHist->Sumw2 ();
  jetCalibRespHist->Sumw2 ();
  jetRecoRespCounts->Sumw2 ();
  jetCalibRespCounts->Sumw2 ();

  TH3D** jetSamplingHist = Get1DArray <TH3D*> (28);
  for (short iJS = 0; iJS < 28; iJS++) {
   jetSamplingHist[iJS] = new TH3D (Form ("jetSamplingHist_iJS%i_dataSet%s", iJS, identifier.Data ()), "", numetabins_mc, etabins_mc, numphibins, phibins, numpbins, pbins);
   jetSamplingHist[iJS]->Sumw2 ();
  }

  TH2D* jetEtaPhiCorr = new TH2D (Form ("jetEtaPhiCorr_dataSet%s", identifier.Data ()), "", numetabins_mc, etabins_mc, 48, -pi, pi);
  jetEtaPhiCorr->Sumw2 ();

  TH2D* jetPtSpectrum = new TH2D (Form ("jetPtSpectrum_dataSet%s", identifier.Data ()), "", numpbins, pbins, numetabins, etabins);
  jetPtSpectrum->Sumw2 ();

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

    jetEtaPhiCorr->Fill (jeta, jphi, weight);
    jetPtSpectrum->Fill (jpt, jeta, weight);

    for (short iJS = 0; iJS < 28; iJS++) {
     jetSamplingHist[iJS]->Fill (jeta, jphi, 1e-3 * (jet_sampling->at (j).at(iJS)));
    }

    if (!InHadCal (jeta, 0.4))
     continue; // require jets to be completely inside the hadronic calorimeter
    if (InDisabledHEC (jeta, jphi))
     continue; // Reject event on additional HEC cuts

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
      xval = t->truth_jet_pt->at (truth_jet);
      ratio = jpt / xval;
     } 
     else {
      xval = t->truth_jet_e->at (truth_jet);
      ratio = je / xval;
     }
     jetCalibRespHist->Fill (xval, jeta, ratio, weight);
     jetCalibRespCounts->Fill (xval, jeta);
    }
    
   } // end calibrated jet loop

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
      xval = t->truth_jet_pt->at (truth_jet);
      ratio = jpt / xval;
     } 
     else {
      xval = t->truth_jet_e->at (truth_jet);
      ratio = je / xval;
     }

     jetRecoRespHist->Fill (xval, jeta, ratio, weight);
     jetRecoRespCounts->Fill (xval, jeta);
    }
   } // end reconstructed jet loop

  } // end loop over events

  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////


  const char* outFileName = Form ("%s/MCValidator/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");
  outFile->cd ();

  // Write histograms to output and clean memory
  jetCalibRespHist->Write ();
  jetCalibRespCounts->Write ();
  jetRecoRespHist->Write ();
  jetRecoRespCounts->Write ();

  for (short iJS = 0; iJS < 28; iJS++) {
   jetSamplingHist[iJS]->Write ();
  }

  jetEtaPhiCorr->Write ();
  jetPtSpectrum->Write ();

  outFile->Close ();
  return;
}

} // end namespace
