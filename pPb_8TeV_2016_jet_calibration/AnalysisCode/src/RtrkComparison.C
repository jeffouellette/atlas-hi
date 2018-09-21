#include "RtrkComparison.h"
#include "Params.h"

#include <TFile.h>
#include <TSystemDirectory.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

TFile* xCalibSystematicsFile = NULL;

vector<Trigger*> triggers = {};


double GetXCalibSystematicError (const double jpt, const double jeta) {
  TFile* file = xCalibSystematicsFile;
  if (!file || !file->IsOpen ())
   return 0;

  if (TMath::Abs (jeta) < xcalibEtabins[0] ||
      xcalibEtabins[sizeof (xcalibEtabins)/sizeof (xcalibEtabins[0]) -1] < TMath::Abs (jeta)) {
   return 0;
  }

  short iEta = 0;
  while (xcalibEtabins[iEta] < TMath::Abs (jeta)) iEta++;
  iEta--;

  const TString hname = TString ("fsys_rel_") + Form ("%i", iEta);
  TH1D* fsys_rel = (TH1D*)file->Get (hname.Data ());

  return TMath::Abs (fsys_rel->GetBinContent (fsys_rel->FindBin (jpt)) - 1) * jpt;
}


TString GetIdentifier (const int dataSet, const TString inFileName, const bool isMC, const bool isSignalOnlySample, const bool periodA) {
  if (!isMC) return to_string (dataSet);

  TString id = (periodA ? "pPb_" : "Pbp_");

  id = id + (isSignalOnlySample ? "Signal_" : "Overlay_");

  if (inFileName.Contains ("jetjet")) { // dijet
   if (dataSet <= 0) return "";
   id = id + "Dijet_Slice" + to_string (dataSet);
  }
  else if (inFileName.Contains ("42310") && inFileName.Contains ("Slice")) { // gamma+jet
   if (dataSet <= 0) return "";
   id = id + "GammaJet_Slice" + to_string (dataSet);
  }
  else if (inFileName.Contains ("ZeeJet")) { // Zee+jet
   if (dataSet < 0) return "";
   id = id + "ZeeJet" + (dataSet == 0 ? "" : "_Slice" + to_string (dataSet));
  }
  else if (inFileName.Contains ("ZmumuJet")) { // Zmumu+jet
   if (dataSet != 0) return "";
   id = id + "ZmumuJet";
  }

  return id;
}


void RtrkComparison (const int dataSet,
                     const double luminosity,
                     const bool isMC, 
                     const bool isPeriodA,
                     const TString inFileName,
                     const double crossSection_microbarns,
                     const double filterEfficiency,
                     const int numberEvents)
{

  SetupDirectories ("", "pPb_8TeV_2016_jet_calibration/");

  const bool isSignalOnlySample = isMC && TString (inFileName).Contains ("signalonly");
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   if (inFileName == "") {
    if (!isMC) fileIdentifier = to_string (dataSet);
    else fileIdentifier = TString (dataSet > 0 ? ("Slice" + to_string (dataSet)) : (dataSet==0 ? "ZmumuJet" : ("ZeeJet"+to_string (-dataSet)))) + (isPeriodA ? ".pPb":".Pbp");
   }
   else fileIdentifier = inFileName;

   const TString path = dataPath + "/rtrk_data/";
   TSystemDirectory dir (path.Data (), path.Data ());
   TList* sysfiles = dir.GetListOfFiles ();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile* sysfile;
   TString fname;
   TIter next (sysfiles);

   while ( (sysfile = (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In RtrkComparison.C (breakpoint B): Found " << fname.Data () << endl;
     
     if (fname.Contains (fileIdentifier)) {
      file = new TFile (path+fname, "READ");
      tree = (TTree*)file->Get ("tree");
      break;
     }
    }
   }
  }
  if (tree == NULL || file == NULL) {
   cout << "Error: In RtrkComparison.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetHIJets ();
  t->SetGetEMTopoJets ();
  t->SetGetTracks ();
  //t->SetGetElectrons ();
  //t->SetGetPhotons ();
  t->SetBranchAddresses ();

  if (!isMC) {
   //for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
   // Trigger* trig = new Trigger (electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
   // trig->minPt = electronTriggerMinPtCuts[electronTriggerN];
   // trig->maxPt = electronTriggerMaxPtCuts[electronTriggerN];
   // triggers.push_back (trig);
   // tree->SetBranchAddress (electronTriggerNames[electronTriggerN], & (trig->trigBool));
   // tree->SetBranchAddress (Form ("%s_prescale", electronTriggerNames[electronTriggerN]), & (trig->trigPrescale));
   //}

   //for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
   // Trigger* trig = new Trigger (muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
   // trig->minPt = muonTriggerMinPtCuts[muonTriggerN];
   // trig->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
   // triggers.push_back (trig);
   // tree->SetBranchAddress (muonTriggerNames[muonTriggerN], & (trig->trigBool));
   // tree->SetBranchAddress (Form ("%s_prescale", muonTriggerNames[muonTriggerN]), & (trig->trigPrescale));
   //}

   //for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
   // Trigger* trig = new Trigger (photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
   // trig->minPt = photonTriggerMinPtCuts[photonTriggerN];
   // trig->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
   // triggers.push_back (trig);
   // tree->SetBranchAddress (photonTriggerNames[photonTriggerN], & (trig->trigBool));
   // tree->SetBranchAddress (Form ("%s_prescale", photonTriggerNames[photonTriggerN]), & (trig->trigPrescale));
   //}

   for (int jetTriggerN = 0; jetTriggerN < jetTrigLength; jetTriggerN++) {
    Trigger* trig = new Trigger (jetTriggerNames[jetTriggerN], jetTriggerMinPtCuts[jetTriggerN], -2.47, 2.47);
    trig->minPt = jetTriggerMinPtCuts[jetTriggerN];
    trig->maxPt = jetTriggerMaxPtCuts[jetTriggerN];
    triggers.push_back (trig);
    tree->SetBranchAddress (jetTriggerNames[jetTriggerN], & (trig->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", jetTriggerNames[jetTriggerN]), & (trig->trigPrescale));
   }

  } // end branch triggers

  int* jet_n;
  vector<float>* jet_pt;
  vector<float>* jet_eta;
  vector<float>* jet_phi;
  vector<float>* jet_e;

  // initialize histograms
  TH2D* jetRtrkHists[2][3][numetabins+1];

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
 
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iErr = 0; iErr < 3; iErr++) {
     TString error = "sys_lo";
     if (iErr == 1) error = "stat";
     else if (iErr == 2) error = "sys_hi";

     jetRtrkHists[iAlgo][iErr][iEta] = new TH2D (Form ("jetRtrkDist_dataSet%s_%s_iEta%i_%s_%s", identifier.Data (), algo.Data (), iEta, (isMC ? "mc":"data"), error.Data ()), "", numpzbins, pzbins, numrtrkbins, rtrkbins);
     jetRtrkHists[iAlgo][iErr][iEta]->Sumw2 ();
    }
   }
  }

  int** nJet = new int*[2]; //[numetabins+1] = {{}, {}};
  for (int i = 0; i < 2; i++) {
   nJet[i] = new int[numetabins+1];
   for (int iEta = 0; iEta <= numetabins; iEta++) {
    nJet[i][iEta] = 0;
   }
  }

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");

  const int numEntries = tree->GetEntries ();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (int entry = 0; entry < numEntries; entry++) {
   tree->GetEntry (entry);

   // basically just do the full analysis twice, once for each algorithm
   for (short iAlgo = 0; iAlgo < 2; iAlgo++) {

    if (iAlgo == 0) { // HI algorithm at the EM scale
     jet_n = & (t->akt4hi_jet_n);
     jet_pt = t->akt4hi_em_xcalib_jet_pt;
     jet_eta = t->akt4hi_em_xcalib_jet_eta;
     jet_phi = t->akt4hi_em_xcalib_jet_phi;
     jet_e = t->akt4hi_em_xcalib_jet_e;
    }
    else { // EMTopo algorithm at the EM scale
     jet_n = & (t->akt4emtopo_jet_n);
     jet_pt = t->akt4emtopo_calib_jet_pt;
     jet_eta = t->akt4emtopo_calib_jet_eta;
     jet_phi = t->akt4emtopo_calib_jet_phi;
     jet_e = t->akt4emtopo_calib_jet_e;
    }

    /////////////////////////////////////////////////////////////////////////////
    // basic event selection: e.g., require a primary vertex
    /////////////////////////////////////////////////////////////////////////////
    if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;


    /////////////////////////////////////////////////////////////////////////////
    // main event selection - check if a trigger fired (data only)
    /////////////////////////////////////////////////////////////////////////////
    float weight = 0; // weight is 0 unless a trigger fired
    if (!isMC) {
     for (Trigger* trig : triggers) {
      if (trig->trigBool &&
          trig->trigPrescale > 0 &&
          (weight == 0 || trig->trigPrescale < weight))
       weight = trig->trigPrescale;
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
    for (int j = 0; j < *jet_n; j++) {
     const float jpt = jet_pt->at (j);
     const float jeta = jet_eta->at (j);
     const float jphi = jet_phi->at (j);

     if (jpt < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (jeta, 0.4))
      continue; // require jets inside hadronic calorimeter
     if (InDisabledHEC (jeta, jphi))
      continue; // Reject event on additional HEC cuts

     //// reject jets near selected photons & electrons
     //double minDeltaR = 1000;
     //for (int p = 0; p < t->photon_n; p++) {
     // // photon cuts
     // if (t->photon_pt->at (p) < photon_pt_cut)
     //  continue; // basic pT cut on photons
     // if (!t->photon_tight->at (p))
     //  continue; // require tight cuts on photons
     // if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*t->photon_pt->at (p))
     //  continue; // require maximum isolation energy on gammas
     // if (!InEMCal (t->photon_eta->at (p)) || InDisabledHEC (t->photon_eta->at (p), t->photon_phi->at (p)))
     //  continue; // require photon to be in EMCal
     // const double deltaR = DeltaR (jeta, t->photon_eta->at (p), jphi, t->photon_phi->at (p));
     // if (deltaR < minDeltaR) {
     //  minDeltaR = deltaR;
     // }
     //}
     //for (int e = 0; e < t->electron_n; e++) {
     // // electron cuts
     // if (t->electron_pt->at (e) < electron_pt_cut)
     //  continue; // basic electron pT cuts
     // if (!InEMCal (t->electron_eta->at (e)))
     //  continue; // reject electrons reconstructed outside EMCal
     // if (!t->electron_loose->at (e))
     //  continue; // reject non-loose electrons
     // if (t->electron_d0sig->at (e) > 5)
     //  continue; // d0 (transverse impact parameter) significance cut
     // if (t->electron_delta_z0_sin_theta->at (e) > 0.5)
     //  continue; // z0 (longitudinal impact parameter) vertex compatibility cut
     // const double deltaR = DeltaR (jeta, t->electron_eta->at (e), jphi, t->electron_phi->at (e));
     // if (deltaR < minDeltaR) {
     //  minDeltaR = deltaR;
     // }
     //}

     // Put the jet in the right eta bin
     short iEta = 0;
     if (etabins[0] < jeta &&
         jeta < etabins[numetabins]) {
      while (etabins[iEta] < jeta) iEta++;
     }
     iEta--;

     if (iEta == -1)
      continue; // Reject jets outside eta bounds

     /////////////////////////////////////////////////////////////////////////////
     // loop over tracks to calculate Rtrk
     /////////////////////////////////////////////////////////////////////////////
     TLorentzVector tlv;
     tlv.SetPtEtaPhiM (0, 0, 0, 0);
     for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      if (1e-3 * (t->trk_pt->at (iTrk)) < trk_pt_cut)
       continue; // reject tracks below pT threshold for consistency in data & MC
      if (t->trk_quality_4->at (iTrk))
       continue; // cut on track quality
      if (DeltaR (jeta, t->trk_eta->at (iTrk), jphi, t->trk_phi->at (iTrk)) < 0.4) { // if track is within jet cone
       TLorentzVector newTrack;
       newTrack.SetPtEtaPhiM (1e-3 * (t->trk_pt->at (iTrk)), t->trk_eta->at (iTrk), t->trk_phi->at (iTrk), 0);
       tlv = tlv + newTrack;
      }
     }
     float sum_trk_pt = tlv.Pt ();
//     cout << "sum_trk_pt / calo_pt = " << sum_trk_pt << ", " << jpt << " = " << sum_trk_pt / jpt << endl;

     const float jptsys = GetXCalibSystematicError (jpt, jeta);
     const float jpts[3] = {jpt-jptsys, jpt, jpt+jptsys};

     for (short iErr = 0; iErr < 3; iErr++) {
      if (jpts[iErr] > 0) 
       jetRtrkHists[iAlgo][iErr][iEta]->Fill (jpt, sum_trk_pt / jpts[iErr], weight);
     }
     nJet[iAlgo][iEta]++;

     //if (sum_trk_pt == 0 && abs (jeta) < 2.1) cout << entry << endl;
     
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
  TFile* outFile = new TFile (outFileName, "RECREATE");

  // Write histograms to output and clean memory
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iErr = 0; iErr < 3; iErr++) {
     jetRtrkHists[iAlgo][iErr][iEta]->Write ();
     if (jetRtrkHists[iAlgo][iErr][iEta]) delete jetRtrkHists[iAlgo][iErr][iEta];
    }
   }
  }

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  TVectorD nJetVec (2* (numetabins+1));
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    nJetVec[iEta + iAlgo* (numetabins+1)] = (double)nJet[iAlgo][iEta];
   }
  }
  nJetVec.Write (Form ("nJetVec_%s", identifier.Data ()));

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
