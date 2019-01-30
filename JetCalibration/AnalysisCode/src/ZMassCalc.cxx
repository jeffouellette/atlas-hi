#include "ZMassCalc.h"
#include "Params.h"
#include "CalibUtils.h"
#include "TreeVariables.h"

#include <Utilities.h>
#include <ArrayTemplates.h>
#include <Trigger.h>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

vector<Trigger*> electronTriggers = {};
vector<Trigger*> muonTriggers = {};

void ZMassCalc (const char* directory,
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
  if (isSignalOnlySample) {
    cout << "Not setup for signal MC samples! Quitting." << endl;
    return;
  }

  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Find the relevant TTree for this run
  //////////////////////////////////////////////////////////////////////////////
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In ZMassCalc.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup local branches for tree 
  //////////////////////////////////////////////////////////////////////////////
  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
    t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetMuons ();
  t->SetBranchAddresses ();

  //////////////////////////////////////////////////////////////////////////////
  // Setup triggers 
  //////////////////////////////////////////////////////////////////////////////
  if (!isMC) {
    for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
      Trigger* temp = new Trigger (electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
      temp->minPt = electronTriggerMinPtCuts[electronTriggerN];
      temp->maxPt = electronTriggerMaxPtCuts[electronTriggerN];
      electronTriggers.push_back (temp);
      tree->SetBranchAddress (electronTriggerNames[electronTriggerN], & (temp->trigBool));
      tree->SetBranchAddress (Form ("%s_prescale", electronTriggerNames[electronTriggerN]), & (temp->trigPrescale));
    }

    for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
      Trigger* temp = new Trigger (muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
      temp->minPt = muonTriggerMinPtCuts[muonTriggerN];
      temp->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
      muonTriggers.push_back (temp);
      tree->SetBranchAddress (muonTriggerNames[muonTriggerN], & (temp->trigBool));
      tree->SetBranchAddress (Form ("%s_prescale", muonTriggerNames[muonTriggerN]), & (temp->trigPrescale));
    }
  } // end branch triggers

  //////////////////////////////////////////////////////////////////////////////
  // Setup output tree
  //////////////////////////////////////////////////////////////////////////////
  const char* outFileName = Form ("%s/ZMassCalc/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  TTree* outTree = new TTree ("jeffsztree", "jeffsztree");
  outTree->SetDirectory (outFile);

  float zpt = 0, zeta = 0, zphi = 0, zm = 0, l1pt = 0, l2pt = 0, l1eta = 0, l2eta = 0, l1phi = 0, l2phi = 0;
  double evtWeight = 0;
  bool _isMC = isMC, _isPeriodA = isPeriodA, isZee = false;

  outTree->Branch ("evt_weight", &evtWeight, "evt_weight/D");
  outTree->Branch ("isMC", &_isMC, "isMC/O");
  outTree->Branch ("isPeriodA", &_isPeriodA, "isPeriodA/O");
  outTree->Branch ("isZee", &isZee, "isZee/O");
  outTree->Branch ("Z_pt", &zpt, "Z_pt/F");
  outTree->Branch ("Z_eta", &zeta, "Z_eta/F");
  outTree->Branch ("Z_phi", &zphi, "Z_phi/F");
  outTree->Branch ("Z_m", &zm, "Z_m/F");
  outTree->Branch ("l1_pt", &l1pt, "l1_pt/F");
  outTree->Branch ("l2_pt", &l2pt, "l2_pt/F");
  outTree->Branch ("l1_eta", &l1eta, "l1_eta/F");
  outTree->Branch ("l2_eta", &l2eta, "l2_eta/F");
  outTree->Branch ("l1_phi", &l1phi, "l1_phi/F");
  outTree->Branch ("l2_phi", &l2phi, "l2_phi/F");

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  const long long numEntries = tree->GetEntries ();
  for (long long entry = 0; entry < numEntries; entry++) {
    tree->GetEntry (entry);

    /////////////////////////////////////////////////////////////////////////////
    // basic event selection: e.g., require a primary vertex
    /////////////////////////////////////////////////////////////////////////////
    if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1))
      continue;

    /////////////////////////////////////////////////////////////////////////////
    // now reject events with less than 2 muons or 2 electrons
    /////////////////////////////////////////////////////////////////////////////
    if (t->muon_n < 2 && t->electron_n < 2)
      continue;

    /////////////////////////////////////////////////////////////////////////////
    // Z (ee) events
    /////////////////////////////////////////////////////////////////////////////
    for (int e1 = 0; e1 < t->electron_n; e1++) { // loop over primary electron
      // electron cuts
      if (t->electron_pt->at (e1) < electron_pt_cut)
        continue; // basic electron pT cuts
      if (!InEMCal (t->electron_eta->at (e1)))
        continue; // reject electrons reconstructed outside EMCal
      if (!t->electron_tight->at (e1))
        continue; // reject non-tight electrons
      if (fabs (t->electron_d0sig->at (e1)) > 5)
        continue; // d0 (transverse impact parameter) significance cut
      if (fabs (t->electron_delta_z0_sin_theta->at (e1)) > 0.5)
        continue; // z0 (longitudinal impact parameter) vertex compatibility cut

      for (int e2 = 0; e2 < e1; e2++) { // loop over secondary electron
        // electron cuts
        if (t->electron_pt->at (e2) < electron_pt_cut)
          continue; // basic electron pT cut
        if (!InEMCal (t->electron_eta->at (e2)))
          continue; // reject electrons reconstructed outside the EMCal
        if (!t->electron_tight->at (e2))
          continue; // reject non-tight electrons
        if (fabs (t->electron_d0sig->at (e2)) > 5)
          continue; // d0 (transverse impact parameter) significance cut
        if (fabs (t->electron_delta_z0_sin_theta->at (e2)) > 0.5)
          continue; // z0 (longitudinal impact parameter) vertex compatibility cut

        // relevant electron kinematic data
        l1pt = t->electron_pt->at (e1);
        l2pt = t->electron_pt->at (e2);
        l1eta = t->electron_eta->at (e1);
        l2eta = t->electron_eta->at (e2);
        l1phi = t->electron_phi->at (e1);
        l2phi = t->electron_phi->at (e2);

        if (!isMC) {
          if (fabs (l1eta) < 0.8) l1pt *= 0.9941;
          else if (fabs (l1eta) < 1.475) l1eta *= 0.9933;
          if (fabs (l2eta) < 0.8) l2pt *= 0.9941;
          else if (fabs (l2eta) < 1.475) l2eta *= 0.9933;
        }

        TLorentzVector electron1, electron2;
        electron1.SetPtEtaPhiM (l1pt, l1eta, l1phi, electron_mass);
        electron2.SetPtEtaPhiM (l2pt, l2eta, l2phi, electron_mass);
        const int le = (l1pt > l2pt ? e1 : e2);
        const double leading_electron_pt = t->electron_pt->at (le);

        // triggering and event weighting
        evtWeight = -1;
        if (!isMC) {
          Trigger* electronTrigger = NULL;
          for (Trigger* trig : electronTriggers) {
            if (trig->trigPrescale > 0 && (electronTrigger == NULL || (trig->trigPrescale < electronTrigger->trigPrescale && trig->minPt <= leading_electron_pt && leading_electron_pt <= trig->maxPt)))
              electronTrigger = trig;
          }
          if (electronTrigger == NULL || !electronTrigger->trigBool || electronTrigger->trigPrescale <= 0.)
            continue;
          evtWeight = electronTrigger->trigPrescale;
        }
        else
          evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;
        if (evtWeight <= 0)
          continue;

        // Reco ee invariant 4-momentum (best guess for Z boson)
        const TLorentzVector Z = electron1 + electron2;
        zpt = Z.Pt (); 
        zeta = Z.Eta ();
        zphi = Z.Phi ();
        zm = Z.M ();

        // Z boson, dielectron cuts
        if (zpt < Z_pt_cut)
          continue; // pt cut on Z bosons
        if (t->electron_charge->at (e1) == t->electron_charge->at (e2))
          continue; // opposite charge requirement

        isZee = true;
        outTree->Fill ();
      }
    } // end loop over electron pairs
    // end Z->ee type events

    /////////////////////////////////////////////////////////////////////////////
    // Z (mumu) events
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
        l1pt = t->muon_pt->at (m1);
        l2pt = t->muon_pt->at (m2);
        l1eta = t->muon_eta->at (m1);
        l2eta = t->muon_eta->at (m2);
        l1phi = t->muon_phi->at (m1);
        l2phi = t->muon_phi->at (m2);
        TLorentzVector muon1, muon2;
        muon1.SetPtEtaPhiM (l1pt, l1eta, l1phi, muon_mass);
        muon2.SetPtEtaPhiM (l2pt, l2eta, l2phi, muon_mass);
        const int lm = (l1pt > l2pt ? m1 : m2);
        const double leading_muon_pt = t->muon_pt->at (lm);

        // triggering and event weighting
        evtWeight = -1;
        if (!isMC) {
          Trigger* muonTrigger = NULL;
          for (Trigger* trig : muonTriggers) {
            if (trig->trigPrescale > 0 && (muonTrigger == NULL || (trig->trigPrescale < muonTrigger->trigPrescale && trig->minPt <= leading_muon_pt && leading_muon_pt <= trig->maxPt)))
              muonTrigger = trig;
          }
          if (muonTrigger == NULL || !muonTrigger->trigBool || muonTrigger->trigPrescale <= 0.)
            continue;
          bool nonMuonTrigger = false;
          for (Trigger* trig : electronTriggers) {
            if (trig->trigBool) nonMuonTrigger = true;
          }
          if (nonMuonTrigger)
            continue; // reject events where a non-electron trigger also fired.
          evtWeight = muonTrigger->trigPrescale;
        }
        else
          evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;
        if (evtWeight <= 0)
          continue;
      
        // Reco mumu invariant 4-momentum (best guess for Z boson)
        const TLorentzVector Z = muon1 + muon2;
        zpt = Z.Pt (); 
        zeta = Z.Eta ();
        zphi = Z.Phi ();
        zm = Z.M ();

        // Z boson, dimuon cuts
        if (zpt < Z_pt_cut)
          continue; // pt cut on Z bosons
        if (t->muon_charge->at (m1) == t->muon_charge->at (m2))
          continue; // require oppositely charged muons

        isZee = false;
        outTree->Fill ();
      }
    } // end loop over muon pairs
    // end Z->mumu type events

  } // end loop over events


  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  outFile->Write ();

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
