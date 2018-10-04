#include "ZmumuJets.h"
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

vector<Trigger*> muonTriggers = {};

void ZmumuJets (const char* directory,
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
   cout << "Error: In ZmumuJets.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets ();
  t->SetGetHIJets ();
  t->SetGetMuons ();
  t->SetBranchAddresses ();

  if (!isMC) {
   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* temp = new Trigger (muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    temp->minPt = muonTriggerMinPtCuts[muonTriggerN];
    temp->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
    muonTriggers.push_back (temp);
    tree->SetBranchAddress (muonTriggerNames[muonTriggerN], & (temp->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", muonTriggerNames[muonTriggerN]), & (temp->trigPrescale));
   }
  } // end branch triggers

  // initialize histograms
  TH3D** zmumuJetHists = Get1DArray <TH3D*> (3);
  //TH2D* zmumuJetHistsSys[3][numetabins];

  {
   TString data = "data";
   if (isMC && !isSignalOnlySample) data = "mc_overlay";
   else if (isMC && isSignalOnlySample) data = "mc_signal";

   for (short iErr = 0; iErr < 3; iErr++) {
    TString error = "sys_lo";
    if (iErr == 1) error = "stat";
    else if (iErr == 2) error = "sys_hi";

    zmumuJetHists[iErr] = new TH3D (Form ("zmumuJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpzbins, pzbins, numetabins, etabins, numxjrefbins, xjrefbins);
    zmumuJetHists[iErr]->Sumw2 ();
    //zmumuJetHistsSys[iErr][iEta] = new TH2D (Form ("zmumuJetPtRatioSys_dataSet%s_iEta%i_%s_%s", identifier.Data (), iEta, data.Data (), error.Data ()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zmumuJetHistsSys[iErr][iEta]->Sumw2 ();
   }
  }

  int** nZmumuJet = Get2DArray <int> (numetabins+1, numpzbins+1);

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");
  dataOverMCFile = new TFile (rootPath + "cc_difference.root", "READ");

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
   // now reject events with less than 2 muons
   /////////////////////////////////////////////////////////////////////////////
   if (t->muon_n < 2) continue;


   /////////////////////////////////////////////////////////////////////////////
   // Z->mumu + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int m1 = 0; m1 < t->muon_n; m1++) { // loop over primary muon
    // primary muon cuts
    if (t->muon_pt->at (m1) < muon_pt_cut)
     continue; // basic muon pT cuts
    if (!t->muon_loose->at (m1))
     continue; // require loose muons
    if (2.4 < abs (t->muon_eta->at (m1)))
     continue; // reject muons reconstructed outside muon spectrometer

    for (int m2 = 0; m2 < m1; m2++) { // loop over secondary muon
     // secondary muon cuts
     if (t->muon_pt->at (m2) < muon_pt_cut)
      continue; // basic muon pT cuts
     if (!t->muon_loose->at (m2))
      continue; // require loose muons
     if (2.4 < abs (t->muon_eta->at (m2))) 
      continue; // reject muons reconstructed outside muon spectrometer

     // relevant muon kinematic data
     TLorentzVector muon1, muon2;
     muon1.SetPtEtaPhiM (t->muon_pt->at (m1), t->muon_eta->at (m1), t->muon_phi->at (m1), muon_mass);
     muon2.SetPtEtaPhiM (t->muon_pt->at (m2), t->muon_eta->at (m2), t->muon_phi->at (m2), muon_mass);
     const int lm = (t->muon_pt->at (m1) > t->muon_pt->at (m2) ? m1 : m2);
     const double leading_muon_pt = t->muon_pt->at (lm);
     //const double leading_muon_eta = t->muon_eta->at (lm);

     // triggering and event weighting
     double weight = 1;
     if (!isMC) {
      Trigger* muonTrigger = NULL;
      for (Trigger* trig : muonTriggers) {
       if (trig->trigPrescale > 0 &&
           (muonTrigger == NULL ||
            (trig->trigPrescale < muonTrigger->trigPrescale &&
             trig->minPt <= leading_muon_pt &&
             leading_muon_pt <= trig->maxPt)))
        muonTrigger = trig;
      }
      if (muonTrigger == NULL ||
          !muonTrigger->trigBool ||
          muonTrigger->trigPrescale <= 0.)
       continue;
      weight = muonTrigger->trigPrescale;
     }
     else weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

     // Reco mumu invariant 4-momentum (best guess for Z boson)
     const TLorentzVector Z = muon1 + muon2;
     const double Z_pt = Z.Pt (); 
     const double Z_phi = Z.Phi ();
     const double Z_m = Z.M ();

     // Z boson, dimuon cuts
     if (t->muon_charge->at (m1) == t->muon_charge->at (m2))
      continue; // require oppositely charged muons
     if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m)
      continue; // cut on our sample Z boson mass
     if (Z_pt < Z_pt_cut)
      continue; // pt cut on Z bosons

     // Put the Z in the right pt bin
     short iP = 0;
     if (pzbins[0] < Z_pt &&
         Z_pt < pzbins[numpzbins]) {
      while (pzbins[iP] < Z_pt) iP++;
     }
     iP--;

     // jet finding
     int lj = -1;
     int sj = -1;
     for (int j = 0; j < t->jet_n; j++) {
      if (t->jet_pt->at (j) < jet_pt_cut)
       continue; // basic jet pT cut
      if (!InHadCal (t->jet_eta->at (j), 0.4))
       continue; // require jets inside hadronic calorimeter
      if (DeltaR (t->muon_eta->at (m1), t->jet_eta->at (j), t->muon_phi->at (m1), t->jet_phi->at (j)) < 0.2 ||
          DeltaR (t->muon_eta->at (m2), t->jet_eta->at (j), t->muon_phi->at (m2), t->jet_phi->at (j)) < 0.2)
       continue; // require jets to be isolated from both muons
      if (DeltaPhi (t->jet_phi->at (j), Z_phi) < 3*pi/4)
       continue; // require jet to be back-to-back with Z in transverse plane

      // compare to leading jet
      else if (lj == -1 ||
               t->jet_pt->at (lj) < t->jet_pt->at (j)) {
       if (lj == -1 || !InDisabledHEC (t->jet_eta->at (lj), t->jet_phi->at (lj)))
        sj = lj;
       lj = j;
      }

      // compare to subleading jet
      else if ( (sj == -1 ||
                t->jet_pt->at (sj) < t->jet_pt->at (j)) &&
               !InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j))) {
       sj = j;
      }
     } // end jet finding loop
     if (lj == -1) // true iff no candidate jet is found
      continue; // reject on no candidate jet

     // relevant jet kinematic data
     const double ljet_pt = t->jet_pt->at (lj);
     const double ljet_eta = t->jet_eta->at (lj);
     const double ljet_phi = t->jet_phi->at (lj);
     const double sjet_pt = ( (0 <= sj && sj < t->jet_n) ? t->jet_pt->at (sj) : 0);
     const double sjet_phi = ( (0 <= sj && sj < t->jet_n) ? t->jet_phi->at (sj) : 0);

     // Put the jet in the right eta, pt bins
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
      continue; // Reject event on additional HEC cuts
     if (sjet_pt > 12) {
      const double subleading_dPhi = DeltaPhi (sjet_phi, Z_phi);
      //if (sjet_pt / Z_pt > 0.3)
      if (sjet_pt / (Z_pt * TMath::Cos (pi - subleading_dPhi)) > 0.2)
      //if (sjet_pt * TMath::Cos (pi - subleading_dPhi) / Z_pt > 0.2)
       continue; // suppress dijets by requiring leading jet to dominate ptref
     }

     // Calculate opening angle in the transverse plane
     const double dPhi = DeltaPhi (ljet_phi, Z_phi);

     // Calculate systematics on jet pT
     const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError (ljet_pt, ljet_eta));

     // Calculate ptref and xjrefs
     const double ptref = Z_pt * TMath::Cos (pi - dPhi);
     const double ptratio[3] = { (ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

     // Fill dimuon xjref histograms
     for (short iErr = 0; iErr < 3; iErr++) {
      zmumuJetHists[iErr]->Fill (Z_pt, ljet_eta, ptratio[iErr], weight);
      //zmumuJetHists[iErr]->Fill (ljet_pt, ljet_eta, ptratio[iErr], weight);
     }

     // Increment Z+jet counters
     if (iEta != -1 && iP != -1) nZmumuJet[iEta][iP]++;
     if (iEta != -1) nZmumuJet[iEta][numpzbins]++;
     if (iP != -1) nZmumuJet[numetabins][iP]++;
     nZmumuJet[numetabins][numpzbins]++;

     // Fill additional systematics histograms
     //for (short iErr = 0; iErr < 3; iErr++)
     // zmumuJetHistsSys[iErr][iEta]->Fill (ljet_pt, ptratio[iErr]/ptref, weight);
    }
   } // end loop over muon pairs
   // end Z->mumu type events

  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  dataOverMCFile->Close ();
  if (dataOverMCFile) delete dataOverMCFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/ZmumuJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");


  // Write histograms to output and clean memory
  for (short iErr = 0; iErr < 3; iErr++) {
   zmumuJetHists[iErr]->Write ();
  }

  Delete1DArray (zmumuJetHists, 3);

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  TVectorD nZmumuJetVec ( (numetabins+1)* (numpzbins+1));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (int iP = 0; iP <= numpzbins; iP++) {
    nZmumuJetVec[iEta+ (numetabins+1)*iP] = (double)nZmumuJet[iEta][iP];
   }
  }
  nZmumuJetVec.Write (Form ("nZmumuJetVec_%s", identifier.Data ()));

  Delete2DArray (nZmumuJet, numetabins+1, numpzbins+1);

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
