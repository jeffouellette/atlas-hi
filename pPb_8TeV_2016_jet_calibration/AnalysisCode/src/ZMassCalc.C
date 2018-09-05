#include "ZGammaJetCrossCheck.h"

#include <ArrayTemplates.h>

#include <TFile.h>
#include <TSystemDirectory.h>
//#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

vector<Trigger*> electronTriggers = {};
vector<Trigger*> muonTriggers = {};

TString GetIdentifier (const int dataSet, const bool isMC, const bool isValidationSample, const bool periodA) {
  if (!isMC) return to_string(dataSet);
  TString id = "";
  if (periodA) id = "pPb_";
  else id = "Pbp_";
  if (dataSet > 0) { // true for GammaJet samples
   if (isValidationSample) id = id + "Valid_";
   else id = id + "Overlay_";
   id = id + "GammaJet_Slice" + to_string(dataSet);
  }
  else {
   if (dataSet == 0) { // true for Zmumu samples
    id = id + "ZmumuJet";
   }
   else if (dataSet == -6) { // true for Zee overlay samples
    id = id + "ZeeJet_Overlay";
   }
   else { // true for Zee signal-only samples
    id = id + "ZeeJet_Slice" + to_string(-dataSet);
   }
  }
  return id;
}


void ZGammaJetCrossCheck (const int dataSet,
                          const double luminosity,
                          const bool isMC, 
                          const bool isMCperiodAflag,
                          const TString inFileName)
{

  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  bool isPeriodA;
  if (!isMC) isPeriodA = dataSet < 313500;
  else isPeriodA = isMCperiodAflag;

  const bool isValidationSample = isMC && (TString(inFileName).Contains("valid") || TString(inFileName).Contains("Zee"));
  const TString identifier = GetIdentifier(dataSet, isMC, isValidationSample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   if (inFileName == "") {
    if (!isMC) fileIdentifier = to_string(dataSet);
    else fileIdentifier = TString(dataSet > 0 ? ("Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet" : ("ZeeJet"+to_string(-dataSet)))) + (isMCperiodAflag ? ".pPb":".Pbp");
   } else fileIdentifier = inFileName;

   const TString dataPathTemp = dataPath;// + "/data810/";
   TSystemDirectory dir(dataPathTemp.Data(), dataPathTemp.Data());
   TList* sysfiles = dir.GetListOfFiles();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile* sysfile;
   TString fname;
   TIter next(sysfiles);

   while ((sysfile = (TSystemFile*)next())) {
    fname = sysfile->GetName();
    if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
     if (debugStatements) cout << "Status: In ZGammaJetCrossCheck.C (breakpoint B): Found " << fname.Data() << endl;
     
     if (fname.Contains(fileIdentifier)) {
      file = new TFile(dataPathTemp+fname, "READ");
      tree = (TTree*)file->Get("tree");
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

  TreeVariables* t = new TreeVariables(tree, isMC);
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetMuons ();
  t->SetBranchAddresses();

  if (!isMC) {
   for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* temp = new Trigger(electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
    temp->minPt = electronTriggerMinPtCuts[electronTriggerN];
    temp->maxPt = electronTriggerMaxPtCuts[electronTriggerN];
    electronTriggers.push_back(temp);
    tree->SetBranchAddress(electronTriggerNames[electronTriggerN], &(temp->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", electronTriggerNames[electronTriggerN]), &(temp->trigPrescale));
   }

   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* temp = new Trigger(muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    temp->minPt = muonTriggerMinPtCuts[muonTriggerN];
    temp->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
    muonTriggers.push_back(temp);
    tree->SetBranchAddress(muonTriggerNames[muonTriggerN], &(temp->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", muonTriggerNames[muonTriggerN]), &(temp->trigPrescale));
   }
  } // end branch triggers

  // initialize histograms
  TH1D*** zMassSpectra = Get2DArray <TH1D*> (2, numetabins+1);

  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (short spcType = 0; spcType < 2; spcType++) {
    const TString species = (spcType == 0 ? "mumu":"ee");

    zMassSpectra[spcType][iEta] = new TH1D(Form("z%sMassSpectrum_dataSet%s_%s_iEta%i", species.Data(), identifier.Data(), (isMC ? "mc":"data"), iEta), "", 50, 60, 110);
    zMassSpectra[spcType][iEta]->Sumw2();
   }
  }

  int* nZeeMass = Get1DArray <int> (numetabins+1);
  int* nZmumuMass = Get1DArray <int> (numetabins+1);

  const long long numEntries = tree->GetEntries();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (long long entry = 0; entry < numEntries; entry++) {
   tree->GetEntry(entry);


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if ((t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at(0) != 1)) continue;


   /////////////////////////////////////////////////////////////////////////////
   // now reject events with less than 2 muons, 2 electrons AND 1 photon
   /////////////////////////////////////////////////////////////////////////////
   if (t->photon_n < 1 && t->muon_n < 2 && t->electron_n < 2) continue;


   /////////////////////////////////////////////////////////////////////////////
   // Z (ee) + jet events
   /////////////////////////////////////////////////////////////////////////////
   for (int e1 = 0; e1 < t->electron_n; e1++) { // loop over primary electron

    // electron cuts
    if (t->electron_pt->at(e1) < electron_pt_cut)
     continue; // basic electron pT cuts
    if (!InEMCal (t->electron_eta->at(e1)))
     continue; // reject electrons reconstructed outside EMCal
    if (!t->electron_loose->at(e1))
     continue; // reject non-loose electrons
    if (t->electron_d0sig->at(e1) > 5)
     continue; // d0 (transverse impact parameter) significance cut
    if (t->electron_delta_z0_sin_theta->at(e1) > 0.5)
     continue; // z0 (longitudinal impact parameter) vertex compatibility cut

    for (int e2 = 0; e2 < e1; e2++) { // loop over secondary electron
     // electron cuts
     if (t->electron_pt->at(e2) < electron_pt_cut)
      continue; // basic electron pT cut
     if (!InEMCal (t->electron_eta->at(e2)))
      continue; // reject electrons reconstructed outside the EMCal
     if (!t->electron_loose->at(e2))
      continue; // reject non-loose electrons
     if (t->electron_d0sig->at(e2) > 5)
      continue; // d0 (transverse impact parameter) significance cut
     if (t->electron_delta_z0_sin_theta->at(e2) > 0.5)
      continue; // z0 (longitudinal impact parameter) vertex compatibility cut

     // relevant electron kinematic data
     TLorentzVector electron1, electron2;
     electron1.SetPtEtaPhiM (t->electron_pt->at(e1), t->electron_eta->at(e1), t->electron_phi->at(e1), 0);//electron_mass);
     electron2.SetPtEtaPhiM (t->electron_pt->at(e2), t->electron_eta->at(e2), t->electron_phi->at(e2), 0);//electron_mass);
     const int le = (t->electron_pt->at(e1) > t->electron_pt->at(e2) ? e1 : e2);
     const double leading_electron_pt = t->electron_pt->at(le);
     //const double leading_electron_eta = t->electron_eta->at(le);

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
     const double Z_pt = Z.Pt(); 
     const double Z_eta = Z.Eta();
     const double Z_m = Z.M();

     // Z boson, dielectron cuts
     if (Z_pt < Z_pt_cut)
      continue; // pt cut on Z bosons
     if (t->electron_charge->at(e1) == t->electron_charge->at(e2))
      continue; // opposite charge requirement

     // Fill dielectron mass spectra
     zMassSpectra[1][numetabins]->Fill (Z_m, weight);

     // Increment total Z count
     nZeeMass[numetabins]++;

     // Put the Z in the right eta bin
     short iEta = 0;
     if (etabins[0] < Z_eta &&
         Z_eta < etabins[numetabins]) {
      while (etabins[iEta] < Z_eta) iEta++;
     }
     iEta--;

     if (iEta != -1) {
      // Fill additional dielectron mass spectrum
      zMassSpectra[1][iEta]->Fill (Z_m, weight);

      // Fill appropriate Z counter
      nZeeMass[iEta]++;
     }

    }
   } // end loop over electron pairs
   // end Z->ee type events


   /////////////////////////////////////////////////////////////////////////////
   // Z->mumu + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int m1 = 0; m1 < t->muon_n; m1++) { // loop over primary muon
    // primary muon cuts
    if (t->muon_pt->at(m1) < muon_pt_cut)
     continue; // basic muon pT cuts
    if (!t->muon_loose->at(m1))
     continue; // require loose muons
    if (2.4 < abs(t->muon_eta->at(m1)))
     continue; // reject muons reconstructed outside muon spectrometer

    for (int m2 = 0; m2 < m1; m2++) { // loop over secondary muon
     // secondary muon cuts
     if (t->muon_pt->at(m2) < muon_pt_cut)
      continue; // basic muon pT cuts
     if (!t->muon_loose->at(m2))
      continue; // require loose muons
     if (2.4 < abs(t->muon_eta->at(m2))) 
      continue; // reject muons reconstructed outside muon spectrometer

     // relevant muon kinematic data
     TLorentzVector muon1, muon2;
     muon1.SetPtEtaPhiM (t->muon_pt->at(m1), t->muon_eta->at(m1), t->muon_phi->at(m1), muon_mass);
     muon2.SetPtEtaPhiM (t->muon_pt->at(m2), t->muon_eta->at(m2), t->muon_phi->at(m2), muon_mass);
     const int lm = (t->muon_pt->at(m1) > t->muon_pt->at(m2) ? m1 : m2);
     const double leading_muon_pt = t->muon_pt->at(lm);
     //const double leading_muon_eta = t->muon_eta->at(lm);

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
      if (nonMuonTrigger)
       continue; // reject events where a non-electron trigger also fired.
      weight = muonTrigger->trigPrescale;
     }
     else weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

     // Reco mumu invariant 4-momentum (best guess for Z boson)
     const TLorentzVector Z = muon1 + muon2;
     const double Z_pt = Z.Pt(); 
     const double Z_eta = Z.Eta();
     const double Z_m = Z.M();

     // Z boson, dimuon cuts
     if (Z_pt < Z_pt_cut)
      continue; // pt cut on Z bosons
     if (t->muon_charge->at(m1) == t->muon_charge->at(m2))
      continue; // require oppositely charged muons

     // Fill dimuon mass spectrum
     zMassSpectra[0][numetabins]->Fill (Z_m, weight);

     // Increment total Z counter
     nZmumuMass[numetabins]++;

     // Put the Z boson in the right eta bin
     short iEta = 0;
     if (etabins[0] < Z_eta &&
         Z_eta < etabins[numetabins]) {
      while (etabins[iEta] < Z_eta) iEta++;
     }
     iEta--;

     if (iEta != -1) {
      // Fill additional dimuon mass spectrum specific to opposite jet
      zMassSpectra[0][iEta]->Fill (Z_m, weight);

      // Increment appropriate Z counter
      nZmumuMass[iEta]++;
     }

    }
   } // end loop over muon pairs
   // end Z->mumu type events

  } // end loop over events


  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form("%sdataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");

  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (short spcType = 0; spcType < 2; spcType++) {
    zMassSpectra[spcType][iEta]->Write();
   }
  }
  Delete2DArray (zMassSpectra, 2, numetabins+1);

  TVectorD infoVec(2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write(Form("infoVec_%s", identifier.Data()));

  TVectorD nZeeMassVec (numetabins+1);
  TVectorD nZmumuMassVec (numetabins+1);
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   nZeeMassVec[iEta] = (double)nZeeMass[iEta];
   nZmumuMassVec[iEta] = (double)nZmumuMass[iEta];
  }
  nZeeMassVec.Write(Form("nZeeMassVec_%s", identifier.Data()));
  nZmumuMassVec.Write(Form("nZmumuMassVec_%s", identifier.Data()));

  Delete1DArray (nZeeMass, numetabins+1);
  Delete1DArray (nZmumuMass, numetabins+1);

  outFile->Close();
  if (outFile) delete outFile;
  return;
}

} // end namespace
