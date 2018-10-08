#include "GammaJets.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Utils.h"

#include <Trigger.h>
#include <ArrayTemplates.h>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

vector<Trigger*> photonTriggers = {};

void GammaJets (const char* directory,
                const int dataSet,
                const double luminosity,
                const bool isMC, 
                const bool isPeriodA,
                const char* inFileName,
                const double crossSection_microbarns,
                const double filterEfficiency,
                const int numberEvents)
{

  SetupDirectories ("", "pPb_8TeV_2016_jet_calibration/");

  const bool isSignalOnlySample = isMC && (TString (inFileName).Contains ("signalonly"));
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In GammaJets.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets ();
  t->SetGetHIJets ();
  t->SetGetPhotons ();
  t->SetBranchAddresses ();

  if (!isMC) {
   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* temp = new Trigger (photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    temp->minPt = photonTriggerMinPtCuts[photonTriggerN];
    temp->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
    photonTriggers.push_back (temp);
    tree->SetBranchAddress (photonTriggerNames[photonTriggerN], & (temp->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", photonTriggerNames[photonTriggerN]), & (temp->trigPrescale));
   }
  } // end branch triggers

  // initialize histograms
  TH3D** gJetHists = Get1DArray <TH3D*> (3);
  TH2D* gJetCounts;
  //TH2D* gJetHistsSys[3][numetabins+1];

  {
   TString data = "data";
   if (isMC && !isSignalOnlySample) data = "mc_overlay";
   else if (isMC && isSignalOnlySample) data = "mc_signal";

   gJetCounts = new TH2D (Form ("gJetCounts_dataSet%s_%s", identifier.Data (), data.Data ()), "", numpbins, pbins, numetabins, etabins);
   gJetCounts->Sumw2();

   for (short iErr = 0; iErr < 3; iErr++) {
    TString error = "sys_lo";
    if (iErr == 1) error = "stat";
    else if (iErr == 2) error = "sys_hi";

    gJetHists[iErr] = new TH3D (Form ("gJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
    gJetHists[iErr]->Sumw2 ();

    //if (iErr == 1) {
    // gJetHistsSys[iErr][iEta] = new TH2D (Form ("gJetPtRatioSys_dataSet%s_iEta%i_%s_%s", identifier.Data (), iEta, data.Data (), error.Data ()), "", numpbins, pbins, numSigmaBins, -maxSigma, maxSigma);
    // gJetHistsSys[iErr][iEta]->Sumw2 ();
    //}
   }
  }

  //int** nGammaJet = Get2DArray <int> (numetabins+1, numpbins+1);

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");

  const long long numEntries = tree->GetEntries ();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (long long entry = 0; entry < numEntries; entry++) {
   tree->GetEntry (entry);


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if (t->nvert <= 0 || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;


   /////////////////////////////////////////////////////////////////////////////
   // photon + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    TLorentzVector photon;
    const double photon_pt = t->photon_pt->at (p);
    const double photon_eta = t->photon_eta->at (p);
    const double photon_phi = t->photon_phi->at (p);
    photon.SetPtEtaPhiM (photon_pt, photon_eta, photon_phi, 0);

    // photon cuts
    if (photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at (p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*photon_pt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (photon_eta) || InDisabledHEC (photon_eta, photon_phi))
     continue; // require photon to be in EMCal

    // triggering and event weighting
    double weight = 1;
    if (!isMC) {
     Trigger* photonTrigger = NULL;
     for (Trigger* trig : photonTriggers) {
      if (trig->trigPrescale > 0 &&
          (photonTrigger == NULL ||
           (trig->trigPrescale < photonTrigger->trigPrescale &&
            trig->minPt <= photon_pt &&
            photon_pt <= trig->maxPt)))
       photonTrigger = trig;
     }
     if (photonTrigger == NULL ||
         !photonTrigger->trigBool)
      continue;
     weight = photonTrigger->trigPrescale;
    }
    else weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

    // Put the photon in the right pt bin
    short iP = 0;
    if (pbins[0] < photon_pt &&
        photon_pt < pbins[numpbins]) {
     while (pbins[iP] < photon_pt) iP++;
    }
    iP--;

    // jet finding
    int lj = -1; // allowed to be in HEC when finding
    for (int j = 0; j < t->jet_n; j++) {
     // cuts on leading jet
     if (t->jet_pt->at (j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (t->jet_eta->at (j), 0.4))
      continue; // require jets inside hadronic calorimeter
     double minDeltaR = 1000;
     for (int p2 = 0; p2 < t->photon_n; p2++) { // find minimum dR to a reco photon
      if (!t->photon_tight->at (p2))
       continue; // only compare with tight photons
      const double dR = DeltaR (t->photon_eta->at (p2), t->jet_eta->at (j), t->photon_phi->at (p2), t->jet_phi->at (j));
      if (dR < minDeltaR)
       minDeltaR = dR;
     }
     if (minDeltaR < 0.4)
      continue; // require jet is not reconstructed as a photon
     if (DeltaPhi (photon_phi, t->jet_phi->at (j)) < 7*pi/8)
      continue; // cut on gamma+jet samples not back-to-back in the transverse plane

     // compare to the leading jet
     else if (lj == -1 ||
              t->jet_pt->at (lj) < t->jet_pt->at (j)) {
      lj = j;
     }
    } // end jet finding loop
    if (lj == -1) // true iff there are no jets opposite photon
     continue; // reject on no jets

    // store relevant jet kinematics
    const double ljet_pt = t->jet_pt->at (lj);
    const double ljet_eta = t->jet_eta->at (lj);
    const double ljet_phi = t->jet_phi->at (lj);

    // Put the jet in the right eta bin
    short iEta = 0;
    if (etabins[0] < ljet_eta &&
        ljet_eta < etabins[numetabins]) {
     while (etabins[iEta] < ljet_eta) iEta++;
    }
    iEta--;

    // jet cuts
    if (iEta == -1)
     continue; // Reject jets outside eta bounds
    if (InDisabledHEC (ljet_eta, ljet_phi))
     continue; // require jet to be outside of disabled HEC
    bool hasOtherJet = false;
    for (int j = 0; j < t->jet_n; j++) {
     if (j == lj)
      continue; // don't look at the leading jet, its our candidate :)
     if (DeltaR (t->jet_eta->at (j), photon_eta, t->jet_phi->at (j), photon_phi) < 0.4)
      continue; // reject jets that are just this photon
     if (t->jet_pt->at (j) < 12 || InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
      continue; // basic jet pT cut, also reject on the disabled HEC
     //const double s_dphi = DeltaPhi (t->jet_phi->at (j), photon_phi);
     const double s_xjref = t->jet_pt->at (j) / (photon_pt);// * TMath::Cos (pi - s_dphi));
     if (0.1 < s_xjref) {
      hasOtherJet = true;
      break;
     }
    }
    if (hasOtherJet)
     continue; // cut on other jets that look back-to-back with gamma

    // Calculate opening angle in the transverse plane
    const double dPhi = DeltaPhi (ljet_phi, photon_phi);

    // Calculate systematics
    const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError (ljet_pt, ljet_eta));

    // Calculate ptref and xjrefs
    const double ptref = photon_pt * TMath::Cos (pi - dPhi);
    const double ptratio[3] = { (ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

    // Fill xjref histograms
    for (short iErr = 0; iErr < 3; iErr++) {
     gJetHists[iErr]->Fill (photon_pt, ljet_eta, ptratio[iErr], weight);
    }
    gJetCounts->Fill (photon_pt, ljet_eta);
   } // end loop over photons
    
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/GammaJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  // Write histograms to output and clean memory
  for (short iErr = 0; iErr < 3; iErr++) {
   gJetHists[iErr]->Write ();
  }

  Delete1DArray (gJetHists, 3);

  gJetCounts->Write();
  if (gJetCounts) { delete gJetCounts; gJetCounts = NULL; }

  //gJetHistsSys[1][iEta]->Write ();
  //if (gJetHistsSys[1][iEta]) delete gJetHistsSys[1][iEta];

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
