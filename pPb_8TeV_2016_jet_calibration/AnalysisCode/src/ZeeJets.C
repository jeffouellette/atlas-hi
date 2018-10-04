#include "ZeeJets.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Utils.h"

#include <Trigger.h>
#include <ArrayTemplates.h>

#include <TSystemDirectory.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

vector<Trigger*> electronTriggers = {};

void ZeeJets (const char* directory,
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
   cout << "Error: In ZeeJets.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets ();
  t->SetGetHIJets ();
  t->SetGetElectrons ();
  t->SetBranchAddresses ();

  if (!isMC) {
   for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* temp = new Trigger (electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
    temp->minPt = electronTriggerMinPtCuts[electronTriggerN];
    temp->maxPt = electronTriggerMaxPtCuts[electronTriggerN];
    electronTriggers.push_back (temp);
    tree->SetBranchAddress (electronTriggerNames[electronTriggerN], & (temp->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", electronTriggerNames[electronTriggerN]), & (temp->trigPrescale));
   }
  } // end branch triggers

  // initialize histograms
  TH3D** zeeJetHists = Get1DArray <TH3D*> (3);
  //TH2D* zeeJetHistsSys[3][numetabins];

  {
   TString data = "data";
   if (isMC && !isSignalOnlySample) data = "mc_overlay";
   else if (isMC && isSignalOnlySample) data = "mc_signal";

   for (short iErr = 0; iErr < 3; iErr++) {
    TString error = "sys_lo";
    if (iErr == 1) error = "stat";
    else if (iErr == 2) error = "sys_hi";

    zeeJetHists[iErr] = new TH3D (Form ("zeeJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpzbins, pzbins, numetabins, etabins, numxjrefbins, xjrefbins);
    zeeJetHists[iErr]->Sumw2 ();
    //zeeJetHistsSys[iErr][iEta] = new TH2D (Form ("zeeJetPtRatioSys_dataSet%s_iEta%i_%s_%s", identifier.Data (), iEta, data.Data (), error.Data ()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zeeJetHistsSys[iErr][iEta]->Sumw2 ();
   }
  }

  int** nZeeJet = Get2DArray <int> (numetabins+1, numpzbins+1);

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
   // now reject events with less than 2 electrons
   /////////////////////////////////////////////////////////////////////////////
   if (t->electron_n < 2) continue;


   /////////////////////////////////////////////////////////////////////////////
   // Z (ee) + jet events
   /////////////////////////////////////////////////////////////////////////////
   for (int e1 = 0; e1 < t->electron_n; e1++) { // loop over primary electron

    // electron cuts
    if (t->electron_pt->at (e1) < electron_pt_cut)
     continue; // basic electron pT cuts
    if (!InEMCal (t->electron_eta->at (e1)))
     continue; // reject electrons reconstructed outside EMCal
    if (!t->electron_loose->at (e1))
     continue; // reject non-loose electrons
    if (t->electron_d0sig->at (e1) > 5)
     continue; // d0 (transverse impact parameter) significance cut
    if (t->electron_delta_z0_sin_theta->at (e1) > 0.5)
     continue; // z0 (longitudinal impact parameter) vertex compatibility cut

    for (int e2 = 0; e2 < e1; e2++) { // loop over secondary electron
     // electron cuts
     if (t->electron_pt->at (e2) < electron_pt_cut)
      continue; // basic electron pT cut
     if (!InEMCal (t->electron_eta->at (e2)))
      continue; // reject electrons reconstructed outside the EMCal
     if (!t->electron_loose->at (e2))
      continue; // reject non-loose electrons
     if (t->electron_d0sig->at (e2) > 5)
      continue; // d0 (transverse impact parameter) significance cut
     if (t->electron_delta_z0_sin_theta->at (e2) > 0.5)
      continue; // z0 (longitudinal impact parameter) vertex compatibility cut

     // relevant electron kinematic data
     TLorentzVector electron1, electron2;
     electron1.SetPtEtaPhiM (t->electron_pt->at (e1), t->electron_eta->at (e1), t->electron_phi->at (e1), 0);//electron_mass);
     electron2.SetPtEtaPhiM (t->electron_pt->at (e2), t->electron_eta->at (e2), t->electron_phi->at (e2), 0);//electron_mass);
     const int le = (t->electron_pt->at (e1) > t->electron_pt->at (e2) ? e1 : e2);
     const double leading_electron_pt = t->electron_pt->at (le);
     //const double leading_electron_eta = t->electron_eta->at (le);

     // triggering and event weighting
     double weight = 1;
     if (!isMC) {
      Trigger* electronTrigger = NULL;
      for (Trigger* trig : electronTriggers) {
       if (trig->trigPrescale > 0 &&
           (electronTrigger == NULL ||
            (trig->trigPrescale < electronTrigger->trigPrescale &&
             trig->minPt <= leading_electron_pt &&
             leading_electron_pt <= trig->maxPt)))
        electronTrigger = trig;
      }
      if (electronTrigger == NULL ||
          !electronTrigger->trigBool ||
          electronTrigger->trigPrescale <= 0.)
       continue;
      //bool nonElectronTrigger = false;
      //for (Trigger* trig : photonTriggers) {
      // if (trig->trigBool) nonElectronTrigger = true;
      //}
      //for (Trigger* trig : muonTriggers) {
      // if (trig->trigBool) nonElectronTrigger = true;
      //}
      //if (nonElectronTrigger)
      // continue; // reject events where a non-electron trigger also fired.
      weight = electronTrigger->trigPrescale;
     }
     else weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

     // Reco ee invariant 4-momentum (best guess for Z boson)
     const TLorentzVector Z = electron1 + electron2;
     const double Z_pt = Z.Pt (); 
     const double Z_phi = Z.Phi ();
     const double Z_m = Z.M ();

     // Z boson, dielectron cuts
     if (t->electron_charge->at (e1) == t->electron_charge->at (e2))
      continue; // opposite charge requirement
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

     // Jet finding
     int lj = -1; // allowed to be in disabled HEC while finding
     int sj = -1; // not allowed to be in disabled HEC while finding
     for (int j = 0; j < t->jet_n; j++) {
      if (t->jet_pt->at (j) < jet_pt_cut)
       continue; // basic jet pT cut
      if (!InHadCal (t->jet_eta->at (j), 0.4))
       continue; // require jets inside hadronic calorimeter
      // Require jets to be isolated from both electrons
      if (DeltaR (t->electron_eta->at (e1), t->jet_eta->at (j), t->electron_phi->at (e1), t->jet_phi->at (j)) < 0.2 ||
          DeltaR (t->electron_eta->at (e2), t->jet_eta->at (j), t->electron_phi->at (e2), t->jet_phi->at (j)) < 0.2)
       continue; // require jets to be isolated from both electrons
      if (DeltaPhi (t->jet_phi->at (j), Z_phi) < 3*pi/4)
       continue; // require jet to be back-to-back with Z in transverse plane

      // compare to leading jet
      else if (lj == -1 ||
               t->jet_pt->at (lj) < t->jet_pt->at (j)) {
       if (lj != -1 && !InDisabledHEC (t->jet_eta->at (lj), t->jet_phi->at (lj)))
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

     // exclude uncorrected jets with pt < 20 GeV
     //if (precalib_jet_pt->at (lj) < 20) continue;

     // relevant jet kinematic data
     const double ljet_pt = t->jet_pt->at (lj);
     const double ljet_eta = t->jet_eta->at (lj);
     const double ljet_phi = t->jet_phi->at (lj);
     const double sjet_pt = ( (0 <= sj && sj < t->jet_n) ? t->jet_pt->at (sj) : 0);
     const double sjet_phi = ( (0 <= sj && sj < t->jet_n) ? t->jet_phi->at (sj) : 0);

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

     // Fill xjref histograms 
     for (short iErr = 0; iErr < 3; iErr++) {
      zeeJetHists[iErr]->Fill (Z_pt, ljet_eta, ptratio[iErr], weight);
      //zeeJetHists[iErr]->Fill (ljet_pt, ljet_eta, ptratio[iErr], weight);
     }

     // Increment Z+jet counters
     if (iEta != -1 && iP != -1) nZeeJet[iEta][iP]++;
     if (iEta != -1) nZeeJet[iEta][numpzbins]++;
     if (iP != -1) nZeeJet[numetabins][iP]++;
     nZeeJet[numetabins][numpzbins]++;

     // Fill additional systematics histograms
     //for (short iErr = 0; iErr < 3; iErr++) zeeJetHistsSys[iErr][iEta]->Fill (ljet_pt, ptratio[iErr]/ptref, weight);
    }
   } // end loop over electron pairs
   // end Z->ee type events
    
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  dataOverMCFile->Close ();
  if (dataOverMCFile) delete dataOverMCFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/ZeeJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");


  // Write histograms to output and clean memory
  for (short iErr = 0; iErr < 3; iErr++) {
   zeeJetHists[iErr]->Write ();
  }

  Delete1DArray (zeeJetHists, 3);

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  TVectorD nZeeJetVec ( (numetabins+1)* (numpzbins+1));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (int iP = 0; iP <= numpzbins; iP++) {
    nZeeJetVec[iEta+ (numetabins+1)*iP] = (double)nZeeJet[iEta][iP];
   }
  }
  nZeeJetVec.Write (Form ("nZeeJetVec_%s", identifier.Data ()));

  Delete2DArray (nZeeJet, numetabins+1, numpzbins+1);

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
