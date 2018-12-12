#include "RtrkComparison.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <Trigger.h>
#include <TreeVariables.h>
#include <ArrayTemplates.h>

#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

vector<Trigger*> triggers = {};

void RtrkComparison (const char* directory,
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

  const bool isSignalOnlySample = isMC && TString (inFileName).Contains ("signalonly");
  //if (isSignalOnlySample) {
  //  cout << "Not setup for signal MC samples! Quitting." << endl;
  //  return;
  //}

  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In RtrkComparison.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetSimpleJets (true);
  t->SetGetHIJets ();
  //t->SetGetTracks ();
  //t->SetGetElectrons ();
  t->SetGetPhotons ();
  if (isMC) t->SetGetTruthPhotons ();
  t->SetBranchAddresses ();

  vector<float>* jtrkpt500_vec = 0, *jtrkpt1000_vec = 0;
  tree->SetBranchAddress ("akt4hi_jet_SumPtTrkPt500", &jtrkpt500_vec);
  tree->SetBranchAddress ("akt4hi_jet_SumPtTrkPt1000", &jtrkpt1000_vec);

  if (!isMC) {
   //for (int jetTriggerN = 0; jetTriggerN < jetTrigLength; jetTriggerN++) {
   // Trigger* trig = new Trigger (jetTriggerNames[jetTriggerN], jetTriggerMinPtCuts[jetTriggerN], -2.47, 2.47);
   // trig->minPt = jetTriggerMinPtCuts[jetTriggerN];
   // trig->maxPt = jetTriggerMaxPtCuts[jetTriggerN];
   // triggers.push_back (trig);
   // tree->SetBranchAddress (jetTriggerNames[jetTriggerN], & (trig->trigBool));
   // tree->SetBranchAddress (Form ("%s_prescale", jetTriggerNames[jetTriggerN]), & (trig->trigPrescale));
   //}
   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* trig = new Trigger (photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    trig->minPt = photonTriggerMinPtCuts[photonTriggerN];
    trig->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
    triggers.push_back (trig);
    tree->SetBranchAddress (photonTriggerNames[photonTriggerN], & (trig->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", photonTriggerNames[photonTriggerN]), & (trig->trigPrescale));
   }
  } // end branch triggers

  const char* outFileName = Form ("%s/RtrkComparison/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  TTree* outTree = new TTree ("RtrkTree", "RtrkTree");
  outTree->SetDirectory (outFile);

  float jpt = 0, jeta = 0, jphi = 0, je = 0, jpterr = 0, jtrkpt500 = 0, jtrkpt1000 = 0;
  double evtWeight = 0;
  bool _isMC = isMC, _isPeriodA = isPeriodA;

  outTree->Branch ("evt_weight", &evtWeight, "evt_weight/D");
  outTree->Branch ("isMC", &_isMC, "isMC/O");
  outTree->Branch ("isPeriodA", &_isPeriodA, "isPeriodA/O");
  outTree->Branch ("jet_pt", &jpt, "jet_pt/F");
  outTree->Branch ("jet_eta", &jeta, "jet_eta/F");
  outTree->Branch ("jet_phi", &jphi, "jet_phi/F");
  outTree->Branch ("jet_e", &je, "jet_e/F");
  outTree->Branch ("jet_SumPtTrkPt500", &jtrkpt500, "jet_SumPtTrkPt500/F");
  outTree->Branch ("jet_SumPtTrkPt1000", &jtrkpt1000, "jet_SumPtTrkPt1000/F");
  outTree->Branch ("jet_pt_sys", &jpterr, "jet_pt_sys/F");

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");

  const int numEntries = tree->GetEntries ();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (int entry = 0; entry < numEntries; entry++) {

    tree->GetEntry (entry);

    /////////////////////////////////////////////////////////////////////////////
    // basic event selection: e.g., require a primary vertex
    /////////////////////////////////////////////////////////////////////////////
    if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1))
      continue;


    /////////////////////////////////////////////////////////////////////////////
    // find the leading photon and require it to trigger 
    /////////////////////////////////////////////////////////////////////////////
    int lP = -1;
    for (int iP = 0; iP < t->photon_n; iP++) {
      if (!InEMCal (t->photon_eta->at (iP)))
        continue;
      if (!t->photon_tight->at (iP))
        continue;
      if (t->photon_topoetcone40->at (iP) > isolationEnergyIntercept + isolationEnergySlope*t->photon_pt->at (iP))
        continue;
      if (isMC && 1 <= dataSet && dataSet <= numdpbins && (t->photon_pt->at (iP) < dpbins[dataSet-1] || dpbins[dataSet] < t->photon_pt->at (iP)))
        continue; // require matched truth photons to be in the DP slice

      if (lP == -1 || t->photon_pt->at (lP) < t->photon_pt->at (iP))
        lP = iP;
    }
    if (lP == -1)
      continue;


    /////////////////////////////////////////////////////////////////////////////
    // main event selection - check if a trigger fired (data only)
    /////////////////////////////////////////////////////////////////////////////
    evtWeight = -1; 
    if (!isMC) {
      for (Trigger* trig : triggers) {
        if (t->photon_pt->at (lP) < trig->minPt || trig->maxPt < t->photon_pt->at (lP))
          continue; // if leading jet not in pT bounds of trigger
        if (trig->trigBool && trig->trigPrescale > 0)
          evtWeight = trig->trigPrescale; // store prescale from trigger if fired
          break; // only 1 trigger will correspond to each pT range
      }
    }
    else { // MC weight
      evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;
    }
    if (evtWeight <= 0)
      continue; // reject events which are weighted to 0 or disable


    ///////////////////////////////////////////////////////////////////////////////
    //// find the leading jet 
    ///////////////////////////////////////////////////////////////////////////////
    //int lj = -1;
    //for (int j = 0; j < jet_n; j++) {
    //  if (lj == -1 || jet_pt->at (lj) < jet_pt->at (j))
    //    lj = j;
    //}
    //if (lj == -1)
    //  continue; // reject when can't find leading jet


    /////////////////////////////////////////////////////////////////////////////
    // main jet loop
    /////////////////////////////////////////////////////////////////////////////
    for (int j = 0; j < t->jet_n; j++) {
      jpt = t->jet_pt->at (j);
      jeta = t->jet_eta->at (j);
      jphi = t->jet_phi->at (j);
      je = t->jet_e->at (j);
      jtrkpt500 = jtrkpt500_vec->at (j);
      jtrkpt1000 = jtrkpt1000_vec->at (j);

      if (jpt < jet_pt_cut)
        continue; // basic jet pT cut
      if (!InHadCal (jeta, 2.8)) // 2.8 since jets must be within 2.1 to have tracks
        continue; // require jets inside hadronic calorimeter
      if (InDisabledHEC (jeta, jphi, 0.4))
        continue; // Reject event on additional HEC cuts
      //if (isMC && 1 <= dataSet && dataSet <= numdpbins && (jpt < dpbins[dataSet-1] || dpbins[dataSet] < jpt))
      //  continue;
      if (DeltaPhi (jphi, t->photon_phi->at (lP)) < 7*pi/8)
        continue;

      bool isPhoton = false;
      for (int iP = 0; iP < t->photon_n; iP++) {
        if (!InEMCal (t->photon_eta->at (iP)))
          continue;
        if (!t->photon_tight->at (iP))
          continue;
        if (t->photon_topoetcone40->at (iP) > isolationEnergyIntercept + isolationEnergySlope*t->photon_pt->at (iP))
          continue;
        if (isMC && 1 <= dataSet && dataSet <= numdpbins && (t->photon_pt->at (iP) < dpbins[dataSet-1] || dpbins[dataSet] < t->photon_pt->at (iP)))
          continue;

        if (DeltaR (jeta, t->photon_eta->at (iP), jphi, t->photon_phi->at (iP)) < 0.4) {
          isPhoton = true;
          break;
        }
      }
      if (isPhoton)
        continue;

      ///////////////////////////////////////////////////////////////////////////////
      //// loop over tracks to calculate Rtrk
      ///////////////////////////////////////////////////////////////////////////////
      //jtrkpt500 = 0;
      //jtrkpt1000 = 0;
      //for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      //  //if (1e-3 * (t->trk_pt->at (iTrk)) < trk_pt_cut)
      //  //  continue; // reject tracks below pT threshold for consistency in data & MC
      //  if (t->trk_quality_4->at (iTrk))
      //    continue; // cut on track quality
      //  if (DeltaR (jeta, t->trk_eta->at (iTrk), jphi, t->trk_phi->at (iTrk)) > 0.4) { // if track is within jet cone
      //    continue;
      //  }
      //  if (t->trk_pt->at (iTrk) < 0.5)
      //    continue;
      //  jtrkpt500 += t->trk_pt->at (iTrk);
      //  if (t->trk_pt->at (iTrk) < 1.0)
      //    continue;
      //  jtrkpt1000 += t->trk_pt->at (iTrk);
      //}

      jpterr = GetXCalibSystematicError (jpt, jeta);

      outTree->Fill ();
     
    } // end jet loop

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
