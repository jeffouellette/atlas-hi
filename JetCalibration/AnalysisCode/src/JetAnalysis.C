#include "JetAnalysis.h"
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

vector<Trigger*> photonTriggers = {};

void JetAnalysis (const char* directory,
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
    cout << "Error: In JetAnalysis.C: TTree not obtained for given data set. Quitting." << endl;
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
  t->SetGetSimpleJets ();
  t->SetGetTruthJets ();
  t->SetGetPhotons ();
  t->SetGetTruthPhotons ();
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
  const char* outFileName = Form ("%s/JetAnalysis/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  TTree* outTree = new TTree ("jeffsjets", "jeffsjets");
  outTree->SetDirectory (outFile);

  float jpt = 0, jeta = 0, jphi = 0, je = 0, ppt = 0, peta = 0, pphi = 0;
  double evtWeight = 0;
  bool _isPeriodA = isPeriodA, _isMC = isMC;
  outTree->Branch ("evt_weight", &evtWeight, "evt_weight/D");

  outTree->Branch ("jet_pt", &jpt, "jet_pt/F");
  outTree->Branch ("jet_eta", &jeta, "jet_eta/F");
  outTree->Branch ("jet_phi", &jphi, "jet_phi/F");
  outTree->Branch ("jet_e", &je, "jet_e/F");
  outTree->Branch ("photon_pt", &ppt, "photon_pt/F");
  outTree->Branch ("photon_eta", &peta, "photon_eta/F");
  outTree->Branch ("photon_phi", &pphi, "photon_phi/F");
  outTree->Branch ("isPeriodA", &_isPeriodA, "isPeriodA/O");
  outTree->Branch ("isMC", &_isMC, "isMC/O");

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
    // find leading photon
    /////////////////////////////////////////////////////////////////////////////
    short lP = -1;
    for (short iP = 0; iP < t->photon_n; iP++) {

      /////////////////////////////////////////////////////////////////////////////
      // photon kinematic info
      /////////////////////////////////////////////////////////////////////////////
      ppt = t->photon_pt->at (iP);
      peta = t->photon_eta->at (iP);
      pphi = t->photon_phi->at (iP);

      /////////////////////////////////////////////////////////////////////////////
      // photon cuts
      /////////////////////////////////////////////////////////////////////////////
      if (ppt < photon_pt_cut)
        continue; // basic pT cut on photons
      if (!t->photon_tight->at (iP))
        continue; // require tight photons
      if (t->photon_topoetcone40->at (iP) > isolationEnergyIntercept + isolationEnergySlope*ppt)
        continue; // require maximum isolation energy on gammas
      if (!InEMCal (peta) || InDisabledHEC (peta, pphi))
        continue; // require photon to be in EMCal

      if (isMC && 1 <= dataSet && dataSet <= numdpbins) {
        short tp = -1;
        float minDeltaR = 1000;
        for (short iTP = 0; iTP < t->truth_photon_n; iTP++) {
          const float deltaR = DeltaR (t->truth_photon_eta->at (iTP), peta, t->truth_photon_phi->at (iTP), pphi);
          if (deltaR < minDeltaR) {
            tp = iTP;
            minDeltaR = deltaR;
          }
        }
        if (minDeltaR > 0.2)
          continue; // require photons to be truth-matched in MC
        if (t->truth_photon_pt->at (tp) < dpbins[dataSet-1] || dpbins[dataSet] < t->truth_photon_pt->at (tp))
          continue; // require matched truth photons to be in the DP slice
      }

      if (lP == -1 || t->photon_pt->at (lP) < ppt)
        lP = iP;
    }
    if (lP == -1)
      continue; // require a leading photon
    ppt = t->photon_pt->at (lP);
    peta = t->photon_eta->at (lP);
    pphi = t->photon_phi->at (lP);

    /////////////////////////////////////////////////////////////////////////////
    // triggering and event weighting
    /////////////////////////////////////////////////////////////////////////////
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
      if (photonTrigger == NULL ||
          !photonTrigger->trigBool)
        continue;
      evtWeight = photonTrigger->trigPrescale;
    }
    else
      evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;
    if (evtWeight <= 0)
      continue;


    /////////////////////////////////////////////////////////////////////////////
    // loop over jets
    /////////////////////////////////////////////////////////////////////////////
    for (short iJ = 0; iJ < t->jet_n; iJ++) { // loop over all jets

      /////////////////////////////////////////////////////////////////////////////
      // relevant jet kinematic data
      /////////////////////////////////////////////////////////////////////////////
      jpt = t->jet_pt->at (iJ);
      jeta = t->jet_eta->at (iJ);
      jphi = t->jet_phi->at (iJ);
      je = t->jet_e->at (iJ);

      /////////////////////////////////////////////////////////////////////////////
      // jet cuts
      /////////////////////////////////////////////////////////////////////////////
      if (jpt < jet_pt_cut)
        continue; // basic pT cut on jets
      if (!InHadCal (jeta) || InDisabledHEC (jeta, jphi))
        continue; // require jet to be in EMCal
      if (DeltaPhi (t->photon_phi->at (lP), jphi) < 7*pi/8)
        continue; // require jet to be opposite to photon
      if (isMC && 1 <= dataSet && dataSet <= numdpbins) {
        short tj = -1;
        float minDeltaR = 1000;
        for (short iTJ = 0; iTJ < t->truth_jet_n; iTJ++) {
          const float deltaR = DeltaR (t->truth_jet_eta->at (iTJ), jeta, t->truth_jet_phi->at (iTJ), jphi);
          if (deltaR < minDeltaR) {
            tj = iTJ;
            minDeltaR = deltaR;
          }
        }
        if (minDeltaR > 0.2)
          continue; // require jets to be truth-matched in MC
        if (t->truth_jet_pt->at (tj) < dpbins[dataSet-1] || dpbins[dataSet] < t->truth_jet_pt->at (tj))
          continue; // require matched truth jets to be in the DP slice
      }

      outTree->Fill ();
    } // end loop over jets
     
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
