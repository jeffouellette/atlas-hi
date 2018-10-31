#include "HadronicContamination.h"
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

void HadronicContamination (const char* directory,
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
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In HadronicContamination.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
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
  TH2D** contHists = Get1DArray <TH2D*> (2); // 0=tight, 1=loose non tight
  TH2D** contCounts = Get1DArray <TH2D*> (2);
  //TH2D* gJetHistsSys[3][numetabins+1];

  for (short iCont = 0; iCont < 2; iCont++) {
   TString data = "data";
   if (isMC && !isSignalOnlySample) data = "mc_overlay";
   else if (isMC && isSignalOnlySample) data = "mc_signal";

   const TString contType = (iCont == 0 ? "tight":"loosenontight");

   contHists[iCont] = new TH2D (Form ("contCounts_dataSet%s_%s_%s", identifier.Data (), data.Data (), contType.Data ()), "", numpbins, pbins, numetabins, etabins);
   contHists[iCont]->Sumw2();

   contCounts[iCont] = new TH2D (Form ("contCounts_dataSet%s_%s_%s", identifier.Data (), data.Data (), contType.Data ()), "", numpbins, pbins, numetabins, etabins);
   contCounts[iCont]->Sumw2();

  }

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
    const double photon_pt = t->photon_pt->at (p);
    const double photon_eta = t->photon_eta->at (p);
    const double photon_phi = t->photon_phi->at (p);

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
     }
     if (photonTrigger == NULL ||
         !photonTrigger->trigBool)
      continue;
     weight = photonTrigger->trigPrescale;
    }
    else weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

   } // end loop over photons
    
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/HadronicContamination/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  // Write histograms to output and clean memory
  for (short iErr = 0; iErr < 3; iErr++) {
   gJetHists[iErr]->Write ();
  }

  Delete1DArray (gJetHists, 3);

  gJetCounts->Write();
  if (gJetCounts) { delete gJetCounts; gJetCounts = NULL; }

  //gJetHistsSys[1][iEta]->Write ();
  //if (gJetHistsSys[1][iEta]) delete gJetHistsSys[1][iEta];

  TVectorD infoVec (2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write (Form ("infoVec_%s", identifier.Data ()));

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
