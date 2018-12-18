#include "ZeeJets.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <Trigger.h>
#include <ArrayTemplates.h>
#include <TreeVariables.h>

#include <TSystemDirectory.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

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
    cout << "Error: In ZeeJets.C: TTree not obtained for given data set. Quitting." << endl;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup local branches for tree 
  //////////////////////////////////////////////////////////////////////////////
  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
    t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets ();
  t->SetGetHIJets ();
  t->SetGetElectrons ();
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
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup output tree
  //////////////////////////////////////////////////////////////////////////////
  const char* outFileName = Form ("%s/ZeeJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  TTree* outTree = new TTree ("jeffsztree", "jeffsztree");
  outTree->SetDirectory (outFile);

  float zpt = 0, zeta = 0, zphi = 0, zm = 0, jpt = 0, jeta = 0, jphi = 0, je = 0, dPhi = 0, jpterr = 0;
  double evtWeight = 0;
  bool _isMC = isMC, _isPeriodA = isPeriodA;

  outTree->Branch ("evt_weight", &evtWeight, "evt_weight/D");
  outTree->Branch ("isMC", &_isMC, "isMC/O");
  outTree->Branch ("isPeriodA", &_isPeriodA, "isPeriodA/O");
  outTree->Branch ("Z_pt", &zpt, "Z_pt/F");
  outTree->Branch ("Z_eta", &zeta, "Z_eta/F");
  outTree->Branch ("Z_phi", &zphi, "Z_phi/F");
  outTree->Branch ("Z_m", &zm, "Z_m/F");
  outTree->Branch ("jet_pt", &jpt, "jet_pt/F");
  outTree->Branch ("jet_eta", &jeta, "jet_eta/F");
  outTree->Branch ("jet_phi", &jphi, "jet_phi/F");
  outTree->Branch ("jet_e", &je, "jet_e/F");
  outTree->Branch ("delta_phi", &dPhi, "delta_phi/F");
  outTree->Branch ("jet_pt_sys", &jpterr, "jet_pt_sys/F");

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  const long long numEntries = tree->GetEntries ();
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

      /////////////////////////////////////////////////////////////////////////////
      // electron cuts
      /////////////////////////////////////////////////////////////////////////////
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
        /////////////////////////////////////////////////////////////////////////////
        // electron cuts
        /////////////////////////////////////////////////////////////////////////////
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

        /////////////////////////////////////////////////////////////////////////////
        // relevant electron kinematic data
        /////////////////////////////////////////////////////////////////////////////
        TLorentzVector electron1, electron2;
        electron1.SetPtEtaPhiM (t->electron_pt->at (e1), t->electron_eta->at (e1), t->electron_phi->at (e1), electron_mass);
        electron2.SetPtEtaPhiM (t->electron_pt->at (e2), t->electron_eta->at (e2), t->electron_phi->at (e2), electron_mass);
        const int le = (t->electron_pt->at (e1) > t->electron_pt->at (e2) ? e1 : e2);
        const double leading_electron_pt = t->electron_pt->at (le);

        /////////////////////////////////////////////////////////////////////////////
        // triggering and event weighting
        /////////////////////////////////////////////////////////////////////////////
        evtWeight = -1;
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
          evtWeight = electronTrigger->trigPrescale;
        }
        else
          evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;
        if (evtWeight <= 0)
          continue;

        /////////////////////////////////////////////////////////////////////////////
        // Reco ee invariant 4-momentum (best guess for Z boson)
        /////////////////////////////////////////////////////////////////////////////
        const TLorentzVector Z = electron1 + electron2;
        zpt = Z.Pt ();
        zeta = Z.Eta ();
        zphi = Z.Phi ();
        zm = Z.M ();

        /////////////////////////////////////////////////////////////////////////////
        // Z boson, dielectron cuts
        /////////////////////////////////////////////////////////////////////////////
        if (t->electron_charge->at (e1) == t->electron_charge->at (e2))
          continue; // opposite charge requirement
        if (Z.M () < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z.M ())
          continue; // cut on our sample Z boson mass
        if (zpt < Z_pt_cut)
          continue; // pt cut on Z bosons

        /////////////////////////////////////////////////////////////////////////////
        // Jet finding
        /////////////////////////////////////////////////////////////////////////////
        int lj = -1;
        for (int j = 0; j < t->jet_n; j++) {
          if (t->jet_pt->at (j) < jet_pt_cut)
            continue; // basic jet pT cut
          if (!InHadCal (t->jet_eta->at (j), 0.4))
            continue; // require jets inside hadronic calorimeter
          if (InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
            continue; // Reject event on additional HEC cuts
          if (DeltaR (t->electron_eta->at (e1), t->jet_eta->at (j), t->electron_phi->at (e1), t->jet_phi->at (j)) < 0.2 ||
              DeltaR (t->electron_eta->at (e2), t->jet_eta->at (j), t->electron_phi->at (e2), t->jet_phi->at (j)) < 0.2)
            continue; // require jets to be isolated from both electrons
          if (DeltaPhi (t->jet_phi->at (j), zphi) < 3*pi/4)
            continue; // require jet to be back-to-back with Z in transverse plane

          // compare to leading jet
          else if (lj == -1 || t->jet_pt->at (lj) < t->jet_pt->at (j)) {
            lj = j;
          }
        } // end jet finding loop
        if (lj == -1) // true iff no candidate jet is found
          continue; // reject on no candidate jet

        /////////////////////////////////////////////////////////////////////////////
        // relevant jet kinematic data
        /////////////////////////////////////////////////////////////////////////////
        jpt = t->jet_pt->at (lj);
        jeta = t->jet_eta->at (lj);
        jphi = t->jet_phi->at (lj);
        je = t->jet_e->at (lj);

        /////////////////////////////////////////////////////////////////////////////
        // jet cuts
        /////////////////////////////////////////////////////////////////////////////
        bool hasOtherJet = false;
        for (int j = 0; j < t->jet_n; j++) {
         if (j == lj)
           continue; // don't look at the leading jet, its our candidate :)
         if (DeltaR (t->electron_eta->at (e1), t->jet_eta->at (j), t->electron_phi->at (e1), t->jet_phi->at (j)) < 0.2 ||
             DeltaR (t->electron_eta->at (e2), t->jet_eta->at (j), t->electron_phi->at (e2), t->jet_phi->at (j)) < 0.2)
           continue; // require jets to be isolated from both electrons
         if (t->jet_pt->at (j) < 12 || InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
           continue; // basic jet pT cut, also reject on the disabled HEC
         const double s_dphi = DeltaPhi (t->jet_phi->at (j), zphi);
         if (0.1 < t->jet_pt->at (j) / (zpt * cos (pi - s_dphi))) {
           hasOtherJet = true;
           break;
         }
        }
        if (hasOtherJet)
          continue; // cut on other jets that look back-to-back with gamma

        /////////////////////////////////////////////////////////////////////////////
        // Calculate opening angle in the transverse plane
        /////////////////////////////////////////////////////////////////////////////
        dPhi = DeltaPhi (jphi, zphi);

        /////////////////////////////////////////////////////////////////////////////
        // Calculate systematics on jet pT
        /////////////////////////////////////////////////////////////////////////////
        jpterr = (isMC ? 0:GetXCalibSystematicError (jpt, jeta));

        outTree->Fill ();
      }
    } // end loop over electron pairs
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  outFile->Write ();

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
