#include "SubtractionStudy.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <Trigger.h>
#include <ArrayTemplates.h>
#include <TreeVariables.h>

#include <TH3D.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

vector<Trigger*> photonTriggers = {};

void SubtractionStudy (const char* directory,
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
    cout << "Error: In SubtractionStudy.C: TTree not obtained for given data set. Quitting." << endl;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup local branches for tree 
  //////////////////////////////////////////////////////////////////////////////
  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
    t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetHIJets ();
  t->SetGetSimpleJets (false);
  t->SetGetPhotons ();
  t->SetBranchAddresses ();

  //////////////////////////////////////////////////////////////////////////////
  // Setup triggers 
  //////////////////////////////////////////////////////////////////////////////
  if (!isMC) {
    for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
      Trigger* temp = new Trigger (photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
      temp->minPt = photonTriggerMinPtCuts[photonTriggerN];
      temp->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
      photonTriggers.push_back (temp);
      tree->SetBranchAddress (photonTriggerNames[photonTriggerN], & (temp->trigBool));
      tree->SetBranchAddress (Form ("%s_prescale", photonTriggerNames[photonTriggerN]), & (temp->trigPrescale));
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Setup output tree
  //////////////////////////////////////////////////////////////////////////////
  const char* outFileName = Form ("%s/SubtractionStudy/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  TTree* outTree = new TTree ("SubtractionTree", "SubtractionTree");
  outTree->SetDirectory (outFile);

  float unsubjpt = 0, unsubjeta = 0, unsubjphi = 0, unsubje = 0, subjpt = 0, subjeta = 0, subjphi = 0, subje = 0, fcalScaleFactor = 0, ppt = 0, peta = 0, pphi = 0;
  double evtWeight = 0;
  bool _isMC = isMC, _isPeriodA = isPeriodA;

  outTree->Branch ("evt_weight", &evtWeight, "evt_weight/D");
  outTree->Branch ("isMC", &_isMC, "isMC/O");
  outTree->Branch ("isPeriodA", &_isPeriodA, "isPeriodA/O");
  outTree->Branch ("unsub_jet_pt", &unsubjpt, "unsub_jet_pt/F");
  outTree->Branch ("unsub_jet_eta", &unsubjeta, "unsub_jet_eta/F");
  outTree->Branch ("unsub_jet_phi", &unsubjphi, "unsub_jet_phi/F");
  outTree->Branch ("unsub_jet_e", &unsubje, "unsub_jet_e/F");
  outTree->Branch ("sub_jet_pt", &subjpt, "sub_jet_pt/F");
  outTree->Branch ("sub_jet_eta", &subjeta, "sub_jet_eta/F");
  outTree->Branch ("sub_jet_phi", &subjphi, "sub_jet_phi/F");
  outTree->Branch ("sub_jet_e", &subje, "sub_jet_e/F");
  outTree->Branch ("photon_pt", &ppt, "photon_pt/F");
  outTree->Branch ("photon_eta", &peta, "photon_eta/F");
  outTree->Branch ("photon_phi", &pphi, "photon_phi/F");
  outTree->Branch ("fcalScaleFactor", &fcalScaleFactor, "fcalScaleFactor/F");

  fcalFile = new TFile (rootPath + "fcalHistograms.root", "READ");

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  const long long numEntries = tree->GetEntries ();
  for (long long entry = 0; entry < numEntries; entry++) {
    tree->GetEntry (entry);

    /////////////////////////////////////////////////////////////////////////////
    // basic event selection: e.g., require a primary vertex
    /////////////////////////////////////////////////////////////////////////////
    if (t->nvert <= 0 || (t->nvert >= 1 && t->vert_type->at (0) != 1))
      continue;

    /////////////////////////////////////////////////////////////////////////////
    // photon + jet type events
    /////////////////////////////////////////////////////////////////////////////
    for (short iP = 0; iP < t->photon_n; iP++) { // loop over all photons
      // relevant photon kinematic data
      ppt = t->photon_pt->at (iP);
      peta = t->photon_eta->at (iP);
      pphi = t->photon_phi->at (iP);

      // photon cuts
      if (ppt < photon_pt_cut)
        continue; // basic pT cut on photons
      if (!t->photon_tight->at (iP))
        continue; // require tight photons
      if (t->photon_topoetcone40->at (iP) > isolationEnergyIntercept + isolationEnergySlope*ppt)
        continue; // require maximum isolation energy on gammas
      if (!InEMCal (peta) || InDisabledHEC (peta, pphi))
        continue; // require photon to be in EMCal

      // triggering and event weighting
      evtWeight = -1;
      if (!isMC) {
        Trigger* photonTrigger = NULL;
        for (Trigger* trig : photonTriggers) {
          if (trig->trigPrescale > 0 &&
              (photonTrigger == NULL ||
               (trig->trigPrescale < photonTrigger->trigPrescale &&
                trig->minPt <= ppt &&
                ppt <= trig->maxPt)))
            photonTrigger = trig;
        }
        if (photonTrigger == NULL || !photonTrigger->trigBool)
          continue;
        evtWeight = photonTrigger->trigPrescale;
      }
      else
        evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;
      if (evtWeight <= 0)
        continue;

      // jet finding
      short lJ = -1;
      for (short iJ = 0; iJ < t->akt4hi_jet_n; iJ++) {
        // cuts on leading jet
        if (t->akt4hi_em_xcalib_jet_pt->at (iJ) < jet_pt_cut)
          continue; // basic jet pT cut
        if (!InHadCal (t->akt4hi_em_xcalib_jet_eta->at (iJ), 0.4))
          continue; // require jets inside hadronic calorimeter
        if (InDisabledHEC (t->akt4hi_em_xcalib_jet_eta->at (iJ), t->akt4hi_em_xcalib_jet_phi->at(iJ)))
          continue; // require jet to be outside of disabled HEC
        if (DeltaPhi (pphi, t->akt4hi_em_xcalib_jet_phi->at (iJ)) < 7*pi/8)
          continue; // cut on gamma+jet samples not back-to-back in the transverse plane

        // compare to the leading jet
        else if (lJ == -1 || t->akt4hi_em_xcalib_jet_pt->at (lJ) < t->akt4hi_em_xcalib_jet_pt->at (iJ)) {
          lJ = iJ;
        }
      } // end jet finding loop
      if (lJ == -1) // true iff there are no jets opposite photon
        continue; // reject on no jets

      // store relevant jet kinematics
      subjpt = t->akt4hi_em_jet_pt->at (lJ);
      subjeta = t->akt4hi_em_jet_eta->at (lJ);
      subjphi = t->akt4hi_em_jet_phi->at (lJ);
      subje = t->akt4hi_em_jet_e->at (lJ);
      unsubjpt = t->akt4hi_constit_jet_pt->at (lJ);
      unsubjeta = t->akt4hi_constit_jet_eta->at (lJ);
      unsubjphi = t->akt4hi_constit_jet_phi->at (lJ);
      unsubje = t->akt4hi_constit_jet_e->at (lJ);

      // jet cuts
      bool hasOtherJet = false;
      for (short iJ = 0; iJ < t->akt4hi_jet_n; iJ++) {
        if (iJ == lJ)
          continue; // don't look at the leading jet, its our candidate :)
        if (DeltaR (t->akt4hi_em_xcalib_jet_eta->at (iJ), peta, t->akt4hi_em_xcalib_jet_phi->at (iJ), pphi) < 0.4)
          continue; // reject jets that are just this photon
        if (t->akt4hi_em_xcalib_jet_pt->at (iJ) < 12 || InDisabledHEC (t->akt4hi_em_xcalib_jet_eta->at (iJ), t->akt4hi_em_xcalib_jet_phi->at (iJ)))
          continue; // basic jet pT cut, also reject on the disabled HEC
        const double s_dphi = DeltaPhi (t->akt4hi_em_xcalib_jet_phi->at (iJ), pphi);
        const double s_xjref = t->akt4hi_em_xcalib_jet_pt->at (iJ) / (ppt * TMath::Cos (pi - s_dphi));
        if (0.1 < s_xjref) {
          hasOtherJet = true;
          break;
        }
      }
      if (hasOtherJet)
        continue; // cut on other jets that look back-to-back with gamma

      if (isMC) fcalScaleFactor = GetFCalScaleFactor (isPeriodA ? t->fcalA_et : t->fcalC_et);
      else fcalScaleFactor = 1.;

      outTree->Fill ();
    } // end loop over photons
    
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
