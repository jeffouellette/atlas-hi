#include "FCalDistribution.h"

#include "TreeVariables.h"
#include "Params.h"

#include <GlobalParams.h>
#include <Trigger.h>

#include <TFile.h>
#include <TSystemDirectory.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

vector<Trigger*> triggers = {};

TString GetIdentifier (const int dataSet, const TString inFileName, const bool isMC, const bool isSignalOnlySample, const bool periodA) {
  if (!isMC) return to_string(dataSet);

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


void FCalDistribution (const int dataSet,
                       const double luminosity,
                       const bool isMC, 
                       const bool isPeriodA,
                       const TString inFileName)
{

  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  const bool isSignalOnlySample = isMC && (TString(inFileName).Contains("signalonly"));
  const TString identifier = GetIdentifier(dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   if (inFileName == "") {
    if (!isMC) fileIdentifier = to_string(dataSet);
    else fileIdentifier = TString(dataSet > 0 ? ("Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet" : ("ZeeJet"+to_string(-dataSet)))) + (isPeriodA ? ".pPb":".Pbp");
   } else fileIdentifier = inFileName;

   TSystemDirectory dir(dataPath.Data(), dataPath.Data());
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
     if (debugStatements) cout << "Status: In FCalDistributions.C (breakpoint B): Found " << fname.Data() << endl;
     
     if (fname.Contains(fileIdentifier)) {
      file = new TFile(dataPath+fname, "READ");
      tree = (TTree*)file->Get("tree");
      break;
     }
    }
   }
  }
  if (tree == NULL || file == NULL) {
   cout << "Error: In FCalDistributions.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetVertices ();
  t->SetGetFCals ();
  t->SetBranchAddresses();

  if (!isMC) {
   for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* trig = new Trigger(electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
    trig->minPt = electronTriggerMinPtCuts[electronTriggerN];
    trig->maxPt = electronTriggerMaxPtCuts[electronTriggerN];
    triggers.push_back(trig);
    tree->SetBranchAddress(electronTriggerNames[electronTriggerN], &(trig->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", electronTriggerNames[electronTriggerN]), &(trig->trigPrescale));
   }

   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* trig = new Trigger(muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    trig->minPt = muonTriggerMinPtCuts[muonTriggerN];
    trig->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
    triggers.push_back(trig);
    tree->SetBranchAddress(muonTriggerNames[muonTriggerN], &(trig->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", muonTriggerNames[muonTriggerN]), &(trig->trigPrescale));
   }

   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* trig = new Trigger(photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    trig->minPt = photonTriggerMinPtCuts[photonTriggerN];
    trig->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
    triggers.push_back(trig);
    tree->SetBranchAddress(photonTriggerNames[photonTriggerN], &(trig->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", photonTriggerNames[photonTriggerN]), &(trig->trigPrescale));
   }
  } // end branch triggers

  TH1D* fCal_p_et = new TH1D (Form ("fCal_p_et_dataSet%s", identifier.Data()), "", 125, -50, 200);
  fCal_p_et->Sumw2();
  TH1D* fCal_Pb_et = new TH1D (Form ("fCal_Pb_et_dataSet%s", identifier.Data()), "", 125, -50, 200);
  fCal_Pb_et->Sumw2();

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

   double weight = 0;
   if (isMC) // get weight from MC generation properties if MC
    weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;
   else { // get weight from trigger prescale if data
    Trigger* bestTrigger = NULL;
    for (Trigger* trig : triggers) {
     if (trig->trigPrescale > 0 &&
         (bestTrigger == NULL ||
          trig->trigPrescale < bestTrigger->trigPrescale))
      bestTrigger = trig;
    }
    if (bestTrigger == NULL ||
        !bestTrigger->trigBool ||
        bestTrigger->trigPrescale <= 0.)
     continue;
    weight = bestTrigger->trigPrescale;
   }

   if (isPeriodA) {
     fCal_p_et->Fill (t->fcalC_et, weight);
     fCal_Pb_et->Fill (t->fcalA_et, weight);
   }
   else {
     fCal_p_et->Fill (t->fcalA_et, weight);
     fCal_Pb_et->Fill (t->fcalC_et, weight);
   }

  } // end event loop

  const char* outFileName = Form("%s/FCalDistribution/dataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");

  fCal_p_et->Write();
  if (fCal_p_et) delete fCal_p_et;
  fCal_Pb_et->Write();
  if (fCal_Pb_et) delete fCal_Pb_et;

  outFile->Close();
  if (outFile) delete outFile;

}

} // end namespace
