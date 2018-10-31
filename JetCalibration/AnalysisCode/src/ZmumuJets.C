#include "ZmumuJets.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <Trigger.h>
#include <ArrayTemplates.h>
#include <TreeVariables.h>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

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

  SetupDirectories ("", "JetCalibration/");

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
  TH2D* zmumuJetCounts;

  {
   TString data = "data";
   if (isMC && !isSignalOnlySample) data = "mc_overlay";
   else if (isMC && isSignalOnlySample) data = "mc_signal";

   zmumuJetCounts = new TH2D (Form ("zmumuJetCounts_dataSet%s_%s", identifier.Data (), data.Data ()), "", numpbins, pbins, numetabins, etabins);
   zmumuJetCounts->Sumw2();

   for (short iErr = 0; iErr < 3; iErr++) {
    TString error = "sys_lo";
    if (iErr == 1) error = "stat";
    else if (iErr == 2) error = "sys_hi";

    zmumuJetHists[iErr] = new TH3D (Form ("zmumuJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
    zmumuJetHists[iErr]->Sumw2 ();
   }
  }

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

     // Z boson, dimuon cuts
     if (t->muon_charge->at (m1) == t->muon_charge->at (m2))
      continue; // require oppositely charged muons
     if (Z.M () < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z.M ())
      continue; // cut on our sample Z boson mass
     if (Z.Pt () < Z_pt_cut)
      continue; // pt cut on Z bosons

     // jet finding
     int lj = -1;
     for (int j = 0; j < t->jet_n; j++) {
      if (t->jet_pt->at (j) < jet_pt_cut)
       continue; // basic jet pT cut
      if (!InHadCal (t->jet_eta->at (j), 0.4))
       continue; // require jets inside hadronic calorimeter
      if (InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
       continue; // Reject event on additional HEC cuts
      if (DeltaR (t->muon_eta->at (m1), t->jet_eta->at (j), t->muon_phi->at (m1), t->jet_phi->at (j)) < 0.2 ||
          DeltaR (t->muon_eta->at (m2), t->jet_eta->at (j), t->muon_phi->at (m2), t->jet_phi->at (j)) < 0.2)
       continue; // require jets to be isolated from both muons
      if (DeltaPhi (t->jet_phi->at (j), Z.Phi ()) < 3*pi/4)
       continue; // require jet to be back-to-back with Z in transverse plane

      // compare to leading jet
      else if (lj == -1 || t->jet_pt->at (lj) < t->jet_pt->at (j)) {
       lj = j;
      }
     } // end jet finding loop
     if (lj == -1) // true iff no candidate jet is found
      continue; // reject on no candidate jet

     // relevant jet kinematic data
     const double ljet_pt = t->jet_pt->at (lj);
     const double ljet_eta = t->jet_eta->at (lj);
     const double ljet_phi = t->jet_phi->at (lj);

     // jet cuts
     //if (InDisabledHEC (ljet_eta, ljet_phi))
     // continue; // Reject event on additional HEC cuts
     bool hasOtherJet = false;
     for (int j = 0; j < t->jet_n; j++) {
      if (j == lj)
       continue; // don't look at the leading jet, its our candidate :)
      if (DeltaR (t->muon_eta->at (m1), t->jet_eta->at (j), t->muon_phi->at (m1), t->jet_phi->at (j)) < 0.2 ||
          DeltaR (t->muon_eta->at (m2), t->jet_eta->at (j), t->muon_phi->at (m2), t->jet_phi->at (j)) < 0.2)
       continue; // require jets to be isolated from both muons
      if (t->jet_pt->at (j) < 12 || InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
       continue; // basic jet pT cut, also reject on the disabled HEC
      const double s_dphi = DeltaPhi (t->jet_phi->at (j), Z.Phi ());
      if (0.1 < t->jet_pt->at (j) / (Z.Pt () * cos (pi - s_dphi))) {
       hasOtherJet = true;
       break;
      }
     }
     if (hasOtherJet)
      continue; // cut on other jets that look back-to-back with gamma

     // Calculate opening angle in the transverse plane
     const double dPhi = DeltaPhi (ljet_phi, Z.Phi ());

     // Calculate systematics on jet pT
     const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError (ljet_pt, ljet_eta));

     // Calculate ptref and xjrefs
     const double ptref = Z.Pt () * TMath::Cos (pi - dPhi);
     const double ptratio[3] = { (ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

     // Fill dimuon xjref histograms
     for (short iErr = 0; iErr < 3; iErr++) {
      zmumuJetHists[iErr]->Fill (Z.Pt (), ljet_eta, ptratio[iErr], weight);
     }
     zmumuJetCounts->Fill (Z.Pt (), ljet_eta);
    }
   } // end loop over muon pairs
   // end Z->mumu type events

  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
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

  zmumuJetCounts->Write ();
  if (zmumuJetCounts) { delete zmumuJetCounts; zmumuJetCounts = NULL; }

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
