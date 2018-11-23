#include "PhotonJets.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <Trigger.h>
#include <ArrayTemplates.h>
#include <TreeVariables.h>

#include <TH3D.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

vector<Trigger*> photonTriggers = {};

void PhotonJets (const char* directory,
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
  if (isSignalOnlySample) {
    cout << "Not setup for signal MC samples! Quitting." << endl;
    return;
  }

  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
   cout << "Error: In PhotonJets.C: TTree not obtained for given data set. Quitting." << endl;
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

  //// initialize histograms
  //TH3D*** gJetHists = Get2DArray <TH3D*> (3, 2); // error type, photonic/non-photonic
  //TH2D*** gJetCounts = Get2DArray <TH2D*> (2, 2); // photonic/non-photonic, weight

  //{
  // TString data = "data";
  // if (isMC && !isSignalOnlySample) data = "mc_overlay";
  // else if (isMC && isSignalOnlySample) data = "mc_signal";

  // for (short iPCut = 0; iPCut < 2; iPCut++) {
  //  const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

  //  for (short iWgt = 0; iWgt < 2; iWgt++) {
  //   const TString weight = (iWgt == 0 ? "unweighted":"weighted");

  //   gJetCounts[iPCut][iWgt] = new TH2D (Form ("gJetCounts_dataSet%s_%s_%s_%s", identifier.Data (), data.Data (), pCut.Data (), weight.Data ()), "", numpbins, pbins, numetabins, etabins);
  //   gJetCounts[iPCut][iWgt]->Sumw2();
  //  }

  //  for (short iErr = 0; iErr < 3; iErr++) {
  //   TString error = "sys_lo";
  //   if (iErr == 1) error = "stat";
  //   else if (iErr == 2) error = "sys_hi";

  //   gJetHists[iErr][iPCut] = new TH3D (Form ("gJetPtRatio_dataSet%s_%s_%s_%s", identifier.Data (), data.Data (), error.Data (), pCut.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
  //   gJetHists[iErr][iPCut]->Sumw2 ();
  //  }
  // }
  //}

  const char* outFileName = Form ("%s/PhotonJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  TTree* outTree = new TTree ("PhotonJetTree", "PhotonJetTree");
  outTree->SetDirectory (outFile);

  float evtWeight = 0, purityFactor = 0, ppt = 0, peta = 0, pphi = 0, jpt = 0, jeta = 0, jphi = 0, je = 0, dPhi = 0, jpterr = 0;
  bool _isMC = isMC, _isPeriodA = isPeriodA, ptight = false;

  outTree->Branch ("evt_weight", &evtWeight, "evt_weight/F");
  outTree->Branch ("purity_factor", &purityFactor, "purity_factor/F");
  outTree->Branch ("isMC", &_isMC, "isMC/O");
  outTree->Branch ("isPeriodA", &_isPeriodA, "isPeriodA/O");
  outTree->Branch ("photon_pt", &ppt, "photon_pt/F");
  outTree->Branch ("photon_eta", &peta, "photon_eta/F");
  outTree->Branch ("photon_phi", &pphi, "photon_phi/F");
  outTree->Branch ("tight_photon", &ptight, "tight_photon/O");
  outTree->Branch ("jet_pt", &jpt, "jet_pt/F");
  outTree->Branch ("jet_eta", &jeta, "jet_eta/F");
  outTree->Branch ("jet_phi", &jphi, "jet_phi/F");
  outTree->Branch ("jet_e", &je, "jet_e/F");
  outTree->Branch ("delta_phi", &dPhi, "delta_phi/F");
  outTree->Branch ("jet_pt_sys", &jpterr, "jet_pt_sys/F");

  xCalibSystematicsFile = new TFile (rootPath + "cc_sys_090816.root", "READ");
  purityFile = new TFile (rootPath + "PhotonPurities.root", "READ");

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
    ppt = t->photon_pt->at (p);
    peta = t->photon_eta->at (p);
    pphi = t->photon_phi->at (p);

    // photon cuts
    if (ppt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*ppt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (peta) || InDisabledHEC (peta, pphi))
     continue; // require photon to be in EMCal

    // triggering and event weighting
    evtWeight = 1;
    if (!isMC) {
     Trigger* photonTrigger = NULL;
     for (Trigger* trig : photonTriggers) {
      if (trig->trigPrescale > 0 &&
          (photonTrigger == NULL ||
           (trig->trigPrescale < photonTrigger->trigPrescale &&
            trig->minPt <= ppt &&
            ppt <= trig->maxPt)))
       photonTrigger = trig;
     }
     if (photonTrigger == NULL || !photonTrigger->trigBool)
      continue;
     evtWeight = photonTrigger->trigPrescale;
    }
    else evtWeight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

    // find what contribution this photon makes to the signal
    if (t->photon_tight->at (p))
     ptight = true;
    else if (t->photon_loose->at (p))
     ptight = false;
    else
     continue;

    // jet finding
    int lj = -1;
    for (int j = 0; j < t->jet_n; j++) {
     // cuts on leading jet
     if (t->jet_pt->at (j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (t->jet_eta->at (j), 0.4))
      continue; // require jets inside hadronic calorimeter
     if (InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at(j)))
      continue; // require jet to be outside of disabled HEC
     if (DeltaPhi (pphi, t->jet_phi->at (j)) < 7*pi/8)
      continue; // cut on gamma+jet samples not back-to-back in the transverse plane

     // compare to the leading jet
     else if (lj == -1 || t->jet_pt->at (lj) < t->jet_pt->at (j)) {
      lj = j;
     }
    } // end jet finding loop
    if (lj == -1) // true iff there are no jets opposite photon
     continue; // reject on no jets

    // store relevant jet kinematics
    jpt = t->jet_pt->at (lj);
    jeta = t->jet_eta->at (lj);
    jphi = t->jet_phi->at (lj);
    je = t->jet_e->at (lj);

    // jet cuts
    //if (InDisabledHEC (jeta, jphi))
    // continue; // require jet to be outside of disabled HEC
    bool hasOtherJet = false;
    for (int j = 0; j < t->jet_n; j++) {
     if (j == lj)
      continue; // don't look at the leading jet, its our candidate :)
     if (DeltaR (t->jet_eta->at (j), peta, t->jet_phi->at (j), pphi) < 0.4)
      continue; // reject jets that are just this photon
     if (t->jet_pt->at (j) < 12 || InDisabledHEC (t->jet_eta->at (j), t->jet_phi->at (j)))
      continue; // basic jet pT cut, also reject on the disabled HEC
     const double s_dphi = DeltaPhi (t->jet_phi->at (j), pphi);
     const double s_xjref = t->jet_pt->at (j) / (ppt * TMath::Cos (pi - s_dphi));
     if (0.1 < s_xjref) {
      hasOtherJet = true;
      break;
     }
    }
    if (hasOtherJet)
     continue; // cut on other jets that look back-to-back with gamma

    // Calculate opening angle in the transverse plane
    dPhi = DeltaPhi (jphi, pphi);

    // Calculate systematics
    jpterr = (isMC ? 0:GetXCalibSystematicError (jpt, jeta));

    //// Calculate ptref and xjrefs
    //const double ptref = ppt * TMath::Cos (pi - dPhi);
    //const double ptratio[3] = { (jpt-jpterr)/ptref, jpt/ptref, (jpt+jpterr)/ptref};

    purityFactor = (ptight ? 1. : 1-GetPurity (ppt, peta));

    //// Fill xjref histograms
    //for (short iErr = 0; iErr < 3; iErr++) {
    // gJetHists[iErr][iPCut]->Fill (ppt, jeta, ptratio[iErr], evtWeight * (1.-purity));
    //}
    //gJetCounts[iPCut][0]->Fill (ppt, jeta);
    //gJetCounts[iPCut][1]->Fill (ppt, jeta, evtWeight);

    outTree->Fill ();
   } // end loop over photons
    
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close ();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;

  purityFile->Close ();
  if (purityFile) delete purityFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  //const char* outFileName = Form ("%s/PhotonJets/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  //TFile* outFile = new TFile (outFileName, "RECREATE");

  //// Write histograms to output and clean memory
  //for (short iPCut = 0; iPCut < 2; iPCut++) {
  // for (short iErr = 0; iErr < 3; iErr++) {
  //  gJetHists[iErr][iPCut]->Write ();
  // }
  // for (short iWgt = 0; iWgt < 2; iWgt++) {
  //  gJetCounts[iPCut][iWgt]->Write();
  // }
  //}

  //Delete2DArray (gJetHists, 3, 2);
  //Delete2DArray (gJetCounts, 2, 2);

  outFile->Write ();
  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
