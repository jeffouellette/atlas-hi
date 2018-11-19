#include "PhotonAnalysis.h"
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

void PhotonAnalysis (const char* directory,
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
   cout << "Error: In PhotonAnalysis.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, isMC);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetHIJets ();
  t->SetGetSimpleJets ();
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
  TH2D*** photonSpectrum = Get2DArray <TH2D*> (2, 2); // tight/lnt, weighted/non-weighted
  TH2D*** photonEtaPhi = Get2DArray <TH2D*> (2, 2); // tight/lnt, weighted/non-weighted

  {
   TString data = "data";
   if (isMC && !isSignalOnlySample) data = "mc_overlay";
   else if (isMC && isSignalOnlySample) data = "mc_signal";

   for (short iPCut = 0; iPCut < 2; iPCut++) {
    const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

    for (short iWgt = 0; iWgt < 2; iWgt++) {
     const TString weight = (iWgt == 0 ? "unweighted":"weighted");

     photonSpectrum[iPCut][iWgt] = new TH2D (Form ("photonSpectrum_dataSet%s_%s_%s", identifier.Data (), data.Data (), pCut.Data (), weight.Data ()), "", numpbins, pbins, numetabins, etabins);
     photonSpectrum[iPCut][iWgt]->Sumw2 ();

     photonEtaPhi[iPCut][iWgt] = new TH2D (Form ("photonEtaPhi_dataSet%s_%s_%s", identifier.Data (), data.Data (), pCut.Data (), weight.Data ()), "", 48, -2.4, 2.4, numphibins, phibins);
     photonEtaPhi[iPCut][iWgt]->Sumw2 ();
    }
   }
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
    TLorentzVector photon;
    const double photon_pt = t->photon_pt->at (p);
    const double photon_eta = t->photon_eta->at (p);
    const double photon_phi = t->photon_phi->at (p);
    photon.SetPtEtaPhiM (photon_pt, photon_eta, photon_phi, 0);

    // photon cuts
    if (photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
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

    // find what contribution this photon makes to the signal
    short iPCut = 0;
    if (t->photon_tight->at (p))
     iPCut = 0;
    else if (t->photon_loose->at (p))
     iPCut = 1;
    else
     continue;

    // find leading jet
    int lj = -1;
    for (int j = 0; j < t->jet_n; j++) {
     // cuts on leading jet
     if (t->jet_pt->at (j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (t->jet_eta->at (j), 0.4))
      continue; // require jets inside hadronic calorimeter
     if (InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at(j)))
      continue; // require jet to be outside of disabled HEC
     if (DeltaPhi (photon_phi, t->jet_phi->at (j)) < 7*pi/8)
      continue; // cut on gamma+jet samples not back-to-back in the transverse plane

     // compare to the leading jet
     else if (lj == -1 || t->jet_pt->at (lj) < t->jet_pt->at (j)) {
      lj = j;
     }
    } // end jet finding loop
    if (lj == -1) // true iff there are no jets opposite photon
     continue; // reject on no jets

    // store relevant jet kinematics
    const double ljet_pt = t->jet_pt->at (lj);
    const double ljet_eta = t->jet_eta->at (lj);
    const double ljet_phi = t->jet_phi->at (lj);

    // jet cuts
    bool hasOtherJet = false;
    for (int j = 0; j < t->jet_n; j++) {
     if (j == lj)
      continue; // don't look at the leading jet, its our candidate :)
     if (DeltaR (t->jet_eta->at (j), photon_eta, t->jet_phi->at (j), photon_phi) < 0.4)
      continue; // reject jets that are just this photon
     if (t->jet_pt->at (j) < 12)
      continue; // basic jet pT cut
     if (InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
      continue; // reject on the disabled HEC
     const double s_dphi = DeltaPhi (t->jet_phi->at (j), photon_phi);
     const double s_xjref = t->jet_pt->at (j) / (photon_pt * TMath::Cos (pi - s_dphi));
     if (0.1 < s_xjref) {
      hasOtherJet = true;
      break;
     }
    }
    if (hasOtherJet)
     continue; // cut on other jets that look back-to-back with gamma

    // Fill xjref histograms
    photonSpectrum[iPCut][0]->Fill (photon_pt, photon_eta, weight);
    photonSpectrum[iPCut][1]->Fill (photon_pt, photon_eta);
    photonEtaPhi[iPCut][0]->Fill (photon_eta, photon_phi, weight);
    photonEtaPhi[iPCut][1]->Fill (photon_eta, photon_phi);
   } // end loop over photons
    
  } // end loop over events


  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/PhotonAnalysis/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  // Write histograms to output and clean memory
  for (short iPCut = 0; iPCut < 2; iPCut++) {
   for (short iWgt = 0; iWgt < 2; iWgt++) {
    photonSpectrum[iPCut][iWgt]->Write ();
    photonEtaPhi[iPCut][iWgt]->Write ();
   }
  }

  Delete2DArray (photonSpectrum, 2, 2);
  Delete2DArray (photonEtaPhi, 2, 2);

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
