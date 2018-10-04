#include "FCalDistribution.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Utils.h"

#include <GlobalParams.h>
#include <Trigger.h>

#include <TH2D.h>
#include <TLorentzVector.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

vector<Trigger*> triggers = {};

void FCalDistribution (const char* directory,
                       const int dataSet,
                       const double luminosity,
                       const bool isMC, 
                       const bool isPeriodA,
                       const char* inFileName,
                       const double crossSection_microbarns,
                       const double filterEfficiency,
                       const int numberEvents)
{

  SetupDirectories ("", "pPb_8TeV_2016_jet_calibration/");

  const bool isSignalOnlySample = isMC && (TString (inFileName).Contains ("signalonly"));
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In FCalDistributions.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetFCals ();
  t->SetBranchAddresses ();

  if (!isMC) {
   for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* trig = new Trigger (electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
    trig->minPt = electronTriggerMinPtCuts[electronTriggerN];
    trig->maxPt = electronTriggerMaxPtCuts[electronTriggerN];
    triggers.push_back (trig);
    tree->SetBranchAddress (electronTriggerNames[electronTriggerN], & (trig->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", electronTriggerNames[electronTriggerN]), & (trig->trigPrescale));
   }

   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* trig = new Trigger (muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    trig->minPt = muonTriggerMinPtCuts[muonTriggerN];
    trig->maxPt = muonTriggerMaxPtCuts[muonTriggerN];
    triggers.push_back (trig);
    tree->SetBranchAddress (muonTriggerNames[muonTriggerN], & (trig->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", muonTriggerNames[muonTriggerN]), & (trig->trigPrescale));
   }

   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* trig = new Trigger (photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    trig->minPt = photonTriggerMinPtCuts[photonTriggerN];
    trig->maxPt = photonTriggerMaxPtCuts[photonTriggerN];
    triggers.push_back (trig);
    tree->SetBranchAddress (photonTriggerNames[photonTriggerN], & (trig->trigBool));
    tree->SetBranchAddress (Form ("%s_prescale", photonTriggerNames[photonTriggerN]), & (trig->trigPrescale));
   }
  } // end branch triggers

  TH1D* fCal_p_et = new TH1D (Form ("fCal_p_et_dataSet%s", identifier.Data ()), "", 125, -50, 200);
  fCal_p_et->Sumw2 ();
  TH1D* fCal_Pb_et = new TH1D (Form ("fCal_Pb_et_dataSet%s", identifier.Data ()), "", 125, -50, 200);
  fCal_Pb_et->Sumw2 ();

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

  const char* outFileName = Form ("%s/FCalDistribution/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  fCal_p_et->Write ();
  if (fCal_p_et) delete fCal_p_et;
  fCal_Pb_et->Write ();
  if (fCal_Pb_et) delete fCal_Pb_et;

  outFile->Close ();
  if (outFile) delete outFile;

}

} // end namespace
