#include "GammaJets.h"
#include "Params.h"
#include "TreeVariables.h"

#include <Trigger.h>
#include <ArrayTemplates.h>

#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

TFile* xCalibSystematicsFile = NULL;
TFile* dataOverMCFile = NULL;

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


void GammaJets (const int dataSet,
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

   const TString dataPathTemp = dataPath;// + "/hion5_300/";
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
  t->SetGetPhotons ();
  t->SetBranchAddresses ();

  if (!isMC) {
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

    gJetHists[iErr] = new TH3D (Form ("gJetPtRatio_dataSet%s_%s_%s", identifier.Data (), data.Data (), error.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
    gJetHists[iErr]->Sumw2 ();
    //if (iErr == 1) {
    // gJetHistsSys[iErr][iEta] = new TH2D (Form ("gJetPtRatioSys_dataSet%s_iEta%i_%s_%s", identifier.Data (), iEta, data.Data (), error.Data ()), "", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
    // gJetHistsSys[iErr][iEta]->Sumw2 ();
    //}
   }
  }

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
   if (t->nvert <= 0 || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;


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

  const char* outFileName = Form ("%s/GammaJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");


  // Write histograms to output and clean memory
  for (short iErr = 0; iErr < 3; iErr++) {
   gJetHists[iErr]->Write ();
  }

  Delete1DArray (gJetHists, 3);

  //gJetHistsSys[1][iEta]->Write ();
  //if (gJetHistsSys[1][iEta]) delete gJetHistsSys[1][iEta];

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  TVectorD nGammaJetVec ( (numetabins+1)* (numpbins+1));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (int iP = 0; iP <= numpbins; iP++) {
    nGammaJetVec[iEta+ (numetabins+1)*iP] = (double)nGammaJet[iEta][iP];
   }
  }
  nGammaJetVec.Write (Form ("nGammaJetVec_%s", identifier.Data ()));

  Delete2DArray (nGammaJet, numetabins+1, numpbins+1);

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
