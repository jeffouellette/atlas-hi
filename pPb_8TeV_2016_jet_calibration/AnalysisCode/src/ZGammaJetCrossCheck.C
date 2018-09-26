#include "ZGammaJetCrossCheck.h"
#include "Params.h"

#include <ArrayTemplates.h>

#include <TFile.h>
#include <TSystemDirectory.h>
//#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

TFile* xCalibSystematicsFile = NULL;
TFile* dataOverMCFile = NULL;

vector<Trigger*> electronTriggers = {};
vector<Trigger*> muonTriggers = {};
vector<Trigger*> photonTriggers = {};

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


double GetNewXCalibSystematicError (const double jeta, const double refpt, const bool periodA) {
  return 0;
  TFile* file = dataOverMCFile;
  if (!file || !file->IsOpen ())
   return 0;

  if (jeta < etabins[0] ||
      etabins[numetabins] < jeta) {
   return 0;
  }

  short bin = 0;
  while (bin <= numetabins && etabins[bin] < jeta) bin++;
  bin--;

  const char* period = (periodA ?  "periodA" : "periodB");
  const TString hname = TString (Form ("gJetPtRatio_diff%i_stat_%s", bin, period));
  TH1D* hist = (TH1D*)file->Get (hname.Data ());

  return hist->GetBinContent (hist->FindBin (refpt)) * refpt;
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


void ZGammaJetCrossCheck (const int dataSet,
                          const double luminosity,
                          const bool isMC, 
                          const bool isPeriodA,
                          const TString inFileName,
                          const double crossSection_microbarns,
                          const double filterEfficiency,
                          const int numberEvents)
{

  SetupDirectories ("", "pPb_8TeV_2016_jet_calibration/");

  const bool isSignalOnlySample = isMC && (TString (inFileName).Contains ("signalonly"));
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   if (inFileName == "") {
    if (!isMC) fileIdentifier = to_string (dataSet);
    else {
     cout << "Cannot identify this MC file! Exiting." << endl;
     return;
    }
   } else fileIdentifier = inFileName;

   const TString dataPathTemp = dataPath;// + "/data810/";
   TSystemDirectory dir (dataPathTemp.Data (), dataPathTemp.Data ());
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
     if (debugStatements) cout << "Status: In ZGammaJetCrossCheck.C (breakpoint B): Found " << fname.Data () << endl;
     
     if (fname.Contains (fileIdentifier)) {
      file = new TFile (dataPathTemp+fname, "READ");
      tree = (TTree*)file->Get ("tree");
      break;
     }
    }
   }
  }
  if (tree == NULL || file == NULL) {
   cout << "Error: In ZGammaJetCrossCheck.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
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
  t->SetGetPhotons ();
  t->SetGetMuons ();
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

   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* temp = new Trigger (muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    temp->minPt = muonTriggerMinPtCuts[muonTriggerN];
    temp->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
    muonTriggers.push_back (temp);
    tree->SetBranchAddress (muonTriggerNames[muonTriggerN], & (temp->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", muonTriggerNames[muonTriggerN]), & (temp->trigPrescale));
   }

   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* temp = new Trigger (photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    temp->minPt = photonTriggerMinPtCuts[photonTriggerN];
    temp->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
    photonTriggers.push_back (temp);
    tree->SetBranchAddress (photonTriggerNames[photonTriggerN], & (temp->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", photonTriggerNames[photonTriggerN]), & (temp->trigPrescale));
   }
  } // end branch triggers

  // initialize histograms
  TH3D** zeeJetHists = Get1DArray <TH3D*> (3);
  //TH2D* zeeJetHistsSys[3][numetabins];
  TH3D** zmumuJetHists = Get1DArray <TH3D*> (3);
  //TH2D* zmumuJetHistsSys[3][numetabins];
  TH3D** gJetHists = Get1DArray <TH3D*> (3);
  //TH2D* gJetHistsSys[3][numetabins+1];

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

    zmumuJetHists[iErr] = new TH3D (Form ("zmumuJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpzbins, pzbins, numetabins, etabins, numxjrefbins, xjrefbins);
    zmumuJetHists[iErr]->Sumw2 ();
    //zmumuJetHistsSys[iErr][iEta] = new TH2D (Form ("zmumuJetPtRatioSys_dataSet%s_iEta%i_%s_%s", identifier.Data (), iEta, data.Data (), error.Data ()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zmumuJetHistsSys[iErr][iEta]->Sumw2 ();

    gJetHists[iErr] = new TH3D (Form ("gJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
    gJetHists[iErr]->Sumw2 ();
    //if (iErr == 1) {
    // gJetHistsSys[iErr][iEta] = new TH2D (Form ("gJetPtRatioSys_dataSet%s_iEta%i_%s_%s", identifier.Data (), iEta, data.Data (), error.Data ()), "", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
    // gJetHistsSys[iErr][iEta]->Sumw2 ();
    //}
   }
  }

  int** nZeeJet = Get2DArray <int> (numetabins+1, numpzbins+1);
  int** nZmumuJet = Get2DArray <int> (numetabins+1, numpzbins+1);
  int** nGammaJet = Get2DArray <int> (numetabins+1, numpbins+1);

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
   if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;


   /////////////////////////////////////////////////////////////////////////////
   // now reject events with less than 2 muons, 2 electrons AND 1 photon
   /////////////////////////////////////////////////////////////////////////////
   if (t->photon_n < 1 && t->muon_n < 2 && t->electron_n < 2) continue;


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
      bool nonMuonTrigger = false;
      for (Trigger* trig : electronTriggers) {
       if (trig->trigBool) nonMuonTrigger = true;
      }
      for (Trigger* trig : photonTriggers) {
       if (trig->trigBool) nonMuonTrigger = true;
      }
      if (nonMuonTrigger)
       continue; // reject events where a non-electron trigger also fired.
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


   /////////////////////////////////////////////////////////////////////////////
   // photon + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    TLorentzVector photon;
    const double photon_pt = t->photon_pt->at (p);
    const double photon_eta = t->photon_eta->at (p);
    const double photon_phi = t->photon_phi->at (p);
    photon.SetPtEtaPhiM (photon_pt, photon_eta, photon_phi, 0);

    // photon cuts
    if (photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at (p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*photon_pt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (photon_eta) || InDisabledHEC (photon_eta, photon_phi))
     continue; // require photon to be in EMCal

    // triggering and event weighting
    double weight = 1;
    if (!isMC) {
     Trigger* photonTrigger = NULL;
     for (Trigger* trig : photonTriggers) {
      if (trig->trigPrescale > 0 &&
          (photonTrigger == NULL ||
           (trig->trigPrescale < photonTrigger->trigPrescale &&
            trig->minPt <= photon_pt &&
            photon_pt <= trig->maxPt)))
       photonTrigger = trig;
      //if (trig->trigPrescale > 0 &&
      //    trig->minPt <= photon_pt &&
      //    photon_pt <= trig->maxPt)
      // photonTrigger = trig;
     }
     if (photonTrigger == NULL ||
         !photonTrigger->trigBool)
      continue;
     bool nonPhotonTrigger = false;
     for (Trigger* trig : electronTriggers) {
      if (trig->trigBool) nonPhotonTrigger = true;
     }
     for (Trigger* trig : muonTriggers) {
      if (trig->trigBool) nonPhotonTrigger = true;
     }
     if (nonPhotonTrigger)
      continue; // reject events where a non-electron trigger also fired.
     weight = photonTrigger->trigPrescale;
    }
    else weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

    // Put the photon in the right pt bin
    short iP = 0;
    if (pbins[0] < photon_pt &&
        photon_pt < pbins[numpzbins]) {
     while (pbins[iP] < photon_pt) iP++;
    }
    iP--;

    // jet finding
    int lj = -1; // allowed to be in HEC when finding
    int sj = -1; // not allowed to be in HEC when finding
    for (int j = 0; j < t->jet_n; j++) {
     // cuts on leading jet
     if (t->jet_pt->at (j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (t->jet_eta->at (j), 0.4))
      continue; // require jets inside hadronic calorimeter
     double minDeltaR = 1000;
     for (int p2 = 0; p2 < t->photon_n; p2++) { // find minimum dR to a reco photon
      if (!t->photon_tight->at (p2))
       continue; // only compare with tight photons
      const double dR = DeltaR (t->photon_eta->at (p2), t->jet_eta->at (j), t->photon_phi->at (p2), t->jet_phi->at (j));
      if (dR < minDeltaR)
       minDeltaR = dR;
     }
     if (minDeltaR < 0.4)
      continue; // require jet is not reconstructed as a photon
     if (DeltaPhi (photon_phi, t->jet_phi->at (j)) < 3*pi/4)
      continue; // cut on gamma+jet samples not back-to-back in the transverse plane

     // compare to the leading jet
     else if (lj == -1 ||
              t->jet_pt->at (lj) < t->jet_pt->at (j)) {
      if (lj != -1 &&
          !InDisabledHEC (t->jet_eta->at (lj), t->jet_phi->at (lj)))
       sj = lj;
      lj = j;
     }

     // compare to the subleading jet
     else if ( (sj == -1 ||
               t->jet_pt->at (sj) < t->jet_pt->at (j)) &&
              !InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j))) {
      sj = j;
     }
    } // end jet finding loop
    if (lj == -1) // true iff there are no jets opposite photon
     continue; // reject on no jets

    // store relevant jet kinematics
    const double ljet_pt = t->jet_pt->at (lj);
    const double ljet_eta = t->jet_eta->at (lj);
    const double ljet_phi = t->jet_phi->at (lj);
    //const double sjet_pt = ( (0 <= sj && sj < t->jet_n) ? t->jet_pt->at (sj) : 0);
    //const double sjet_phi = ( (0 <= sj && sj < t->jet_n) ? t->jet_phi->at (sj) : 0);

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
     continue; // require jet to be outside of disabled HEC
    bool hasOtherJet = false;
    for (int j = 0; j < t->jet_n; j++) {
     if (j == lj)
      continue; // don't look at the leading jet, its our candidate :)
     if (DeltaR (t->jet_eta->at (j), photon_eta, t->jet_phi->at (j), photon_phi) < 0.4)
      continue; // reject jets that are just this photon
     if (t->jet_pt->at (j) < 12 || InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
      continue; // basic jet pT cut, also reject on the disabled HEC
     //const double s_dphi = DeltaPhi (t->jet_phi->at (j), photon_phi);
     const double s_xjref = t->jet_pt->at (j) / (photon_pt);// * TMath::Cos (pi - s_dphi));
     if (0.1 < s_xjref) {
      hasOtherJet = true;
      break;
     }
    }
    if (hasOtherJet)
     continue; // cut on other jets that look back-to-back with gamma
    //if (sjet_pt > 12) {
    // const double subleading_dPhi = DeltaPhi (sjet_phi, photon_phi);
    // //if (sjet_pt / photon_pt > 0.1)
    // if (sjet_pt / (photon_pt * TMath::Cos (pi - subleading_dPhi)) > 0.2)
    // //if (sjet_pt * TMath::Cos (pi - subleading_dPhi) / photon_pt > 0.2)
    // //if (sjet_pt / photon_pt > 0.02)
    //  continue; // suppress dijets by requiring leading jet to dominate ptref
    //}

    // Calculate opening angle in the transverse plane
    const double dPhi = DeltaPhi (ljet_phi, photon_phi);

    // Calculate systematics
    const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError (ljet_pt, ljet_eta));

    // Calculate ptref and xjrefs
    const double ptref = photon_pt * TMath::Cos (pi - dPhi);
    const double ptratio[3] = { (ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

    // Fill xjref histograms
    for (short iErr = 0; iErr < 3; iErr++) {
     gJetHists[iErr]->Fill (photon_pt, ljet_eta, ptratio[iErr], weight);
     //gJetHists[iErr]->Fill (ljet_pt, ljet_eta, ptratio[iErr], weight);
    }

    //// if data, calculate additional systematics for this jet
    //if (!isMC) {
    // const double newJetPtSys = GetNewXCalibSystematicError (ljet_eta, photon_pt, isPeriodA) / ljet_pt;
    // if (iEta != -1) gJetHistsSys[1][iEta]->Fill (ljet_pt, newJetPtSys, weight);
    // gJetHistsSys[1][numetabins]->Fill (ljet_pt, newJetPtSys, weight);
    //}

    // Increment gamma+jet counters
    if (iEta != -1 && iP != -1) nGammaJet[iEta][iP]++;
    if (iEta != -1) nGammaJet[iEta][numpbins]++;
    if (iP != -1) nGammaJet[numetabins][iP]++;
    nGammaJet[numetabins][numpbins]++;

    //if (90 < ptref && ptref < 110 && 0.6 < ptratio[1] && ptratio[1] < 0.8) {
    // cout << "Found a weird event! Dataset: " << identifier.Data () << ", entry: " << entry << endl;
    //}
   }
    
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  dataOverMCFile->Close ();
  if (dataOverMCFile) delete dataOverMCFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%sdataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");


  // Write histograms to output and clean memory
  for (short iErr = 0; iErr < 3; iErr++) {
   zmumuJetHists[iErr]->Write ();
   zeeJetHists[iErr]->Write ();
   gJetHists[iErr]->Write ();
  }

  Delete1DArray (zmumuJetHists, 3);
  Delete1DArray (zeeJetHists, 3);
  Delete1DArray (gJetHists, 3);

  //gJetHistsSys[1][iEta]->Write ();
  //if (gJetHistsSys[1][iEta]) delete gJetHistsSys[1][iEta];

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  TVectorD nZeeJetVec ( (numetabins+1)* (numpzbins+1));
  TVectorD nZmumuJetVec ( (numetabins+1)* (numpzbins+1));
  TVectorD nGammaJetVec ( (numetabins+1)* (numpbins+1));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (int iP = 0; iP <= numpzbins; iP++) {
    nZeeJetVec[iEta+ (numetabins+1)*iP] = (double)nZeeJet[iEta][iP];
    nZmumuJetVec[iEta+ (numetabins+1)*iP] = (double)nZmumuJet[iEta][iP];
   }
   for (int iP = 0; iP <= numpbins; iP++) {
    nGammaJetVec[iEta+ (numetabins+1)*iP] = (double)nGammaJet[iEta][iP];
   }
  }
  nZeeJetVec.Write (Form ("nZeeJetVec_%s", identifier.Data ()));
  nZmumuJetVec.Write (Form ("nZmumuJetVec_%s", identifier.Data ()));
  nGammaJetVec.Write (Form ("nGammaJetVec_%s", identifier.Data ()));

  Delete2DArray (nZeeJet, numetabins+1, numpzbins+1);
  Delete2DArray (nZmumuJet, numetabins+1, numpzbins+1);
  Delete2DArray (nGammaJet, numetabins+1, numpbins+1);

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
