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
  t->SetGetSimpleJets (false);
  t->SetGetHIJets ();
  t->SetGetEMTopoJets ();
  t->SetGetTracks ();
  //t->SetGetElectrons ();
  //t->SetGetPhotons ();
  t->SetBranchAddresses ();

  if (!isMC) {
   for (int jetTriggerN = 0; jetTriggerN < jetTrigLength; jetTriggerN++) {
    Trigger* trig = new Trigger (jetTriggerNames[jetTriggerN], jetTriggerMinPtCuts[jetTriggerN], -2.47, 2.47);
    trig->minPt = jetTriggerMinPtCuts[jetTriggerN];
    trig->maxPt = jetTriggerMaxPtCuts[jetTriggerN];
    triggers.push_back (trig);
    tree->SetBranchAddress (jetTriggerNames[jetTriggerN], & (trig->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", jetTriggerNames[jetTriggerN]), & (trig->trigPrescale));
   }
  } // end branch triggers

  int jet_n;
  vector<float>* jet_pt;
  vector<float>* jet_eta;
  vector<float>* jet_phi;
  vector<float>* jet_e;

  // initialize histograms
  TH3D*** jetRtrkHists = Get2DArray <TH3D*> (2, 3); // iAlgo, iErr
  TH2D** jetRtrkCounts = Get1DArray <TH2D*> (2);

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
 
   for (short iErr = 0; iErr < 3; iErr++) {
    TString error = "sys_lo";
    if (iErr == 1) error = "stat";
    else if (iErr == 2) error = "sys_hi";

    jetRtrkHists[iAlgo][iErr] = new TH3D (Form ("jetRtrkDist_dataSet%s_%s_%s_%s", identifier.Data (), algo.Data (), (isMC ? "mc":"data"), error.Data ()), "", numpbins, pbins, numetabins, etabins, numrtrkbins, rtrkbins);
    jetRtrkHists[iAlgo][iErr]->Sumw2 ();
   }

   jetRtrkCounts[iAlgo] = new TH2D (Form ("jetRtrkCounts_dataSet%s_%s_%s", identifier.Data (), algo.Data (), (isMC ? "mc":"data")), "", numpbins, pbins, numetabins, etabins);
   jetRtrkCounts[iAlgo]->Sumw2 ();
   
  }

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");

  const int numEntries = tree->GetEntries ();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (int entry = 0; entry < numEntries; entry++) {
   if (isMC && entry >= 2900652) continue;

   tree->GetEntry (entry);

   // basically just do the full analysis twice, once for each algorithm
   for (short iAlgo = 0; iAlgo < 2; iAlgo++) {

    /////////////////////////////////////////////////////////////////////////////
    // basic event selection: e.g., require a primary vertex
    /////////////////////////////////////////////////////////////////////////////
    if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;


    if (iAlgo == 0) { // HI algorithm at the EM scale
     jet_n = *(&(t->akt4hi_jet_n));
     jet_pt = t->akt4hi_em_xcalib_jet_pt;
     jet_eta = t->akt4hi_em_xcalib_jet_eta;
     jet_phi = t->akt4hi_em_xcalib_jet_phi;
     jet_e = t->akt4hi_em_xcalib_jet_e;
    }
    else { // EMTopo algorithm at the EM scale
     jet_n = *(&(t->akt4emtopo_jet_n));
     jet_pt = t->akt4emtopo_calib_jet_pt;
     jet_eta = t->akt4emtopo_calib_jet_eta;
     jet_phi = t->akt4emtopo_calib_jet_phi;
     jet_e = t->akt4emtopo_calib_jet_e;
    }


    /////////////////////////////////////////////////////////////////////////////
    // find the leading jet 
    /////////////////////////////////////////////////////////////////////////////
    int lj = -1;
    for (int j = 0; j < jet_n; j++) {
     if (lj == -1 || jet_pt->at (lj) < jet_pt->at (j))
      lj = j;
    }
    if (lj == -1)
     continue; // reject when can't find leading jet

    /////////////////////////////////////////////////////////////////////////////
    // main event selection - check if a trigger fired (data only)
    /////////////////////////////////////////////////////////////////////////////
    float weight = 0; // weight is 0 unless a trigger fired
    if (!isMC) {
     for (Trigger* trig : triggers) {
      if (jet_pt->at (lj) < trig->minPt || trig->maxPt < jet_pt->at (lj))
       continue; // if leading jet not in pT bounds of trigger
      if (trig->trigBool && trig->trigPrescale > 0)
       weight = trig->trigPrescale; // store prescale from trigger if fired
      break; // only 1 trigger will correspond to each pT range
     }
    }
    else { // MC weight
     weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;
    }
    if (weight == 0)
     continue; // reject events which are weighted to 0


    /////////////////////////////////////////////////////////////////////////////
    // main jet loop
    /////////////////////////////////////////////////////////////////////////////
    for (int j = 0; j < jet_n; j++) {
     const float jpt = jet_pt->at (j);
     const float jeta = jet_eta->at (j);
     const float jphi = jet_phi->at (j);

     if (jpt < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (jeta, 0.4))
      continue; // require jets inside hadronic calorimeter
     if (InDisabledHEC (jeta, jphi))
      continue; // Reject event on additional HEC cuts

     /////////////////////////////////////////////////////////////////////////////
     // loop over tracks to calculate Rtrk
     /////////////////////////////////////////////////////////////////////////////
     TLorentzVector sumTracks, newTrack;
     sumTracks.SetPtEtaPhiM (0, 0, 0, 0);
     for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      if (1e-3 * (t->trk_pt->at (iTrk)) < trk_pt_cut)
       continue; // reject tracks below pT threshold for consistency in data & MC
      if (t->trk_quality_4->at (iTrk))
       continue; // cut on track quality
      if (DeltaR (jeta, t->trk_eta->at (iTrk), jphi, t->trk_phi->at (iTrk)) < 0.4) { // if track is within jet cone
       newTrack.SetPtEtaPhiM (1e-3 * (t->trk_pt->at (iTrk)), t->trk_eta->at (iTrk), t->trk_phi->at (iTrk), 0);
       sumTracks = sumTracks + newTrack;
      }
     }
     float sum_trk_pt = sumTracks.Pt ();

     const float jptsys = GetXCalibSystematicError (jpt, jeta);
     const float jpts[3] = {jpt-jptsys, jpt, jpt+jptsys};

     for (short iErr = 0; iErr < 3; iErr++) {
      jetRtrkHists[iAlgo][iErr]->Fill (jpt, jeta, sum_trk_pt / jpts[iErr], weight);
     }
     jetRtrkCounts[iAlgo]->Fill (jpt, jeta);
     
    } // end jet loop

   } // end loop over jet algorithms
    
  } // end loop over events

  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/RtrkComparison/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "recreate");

  // Write histograms to output and clean memory
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   for (short iErr = 0; iErr < 3; iErr++) {
    jetRtrkHists[iAlgo][iErr]->Write ();
   }
   jetRtrkCounts[iAlgo]->Write ();
  }

  Delete2DArray (jetRtrkHists, 2, 3);
  Delete1DArray (jetRtrkCounts, 2);

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
