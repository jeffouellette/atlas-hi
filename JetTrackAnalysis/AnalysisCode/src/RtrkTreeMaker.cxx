#include "RtrkTreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"

#include <Utilities.h>

#include <iostream>

using namespace std;

namespace JetTrackAnalysis {

void RtrkTreeMaker (const char* directory,
                    const int dataSet,
                    const bool isMC,
                    const bool isPeriodA,
                    const char* inFileName,
                    const double crossSection_microbarns,
                    const double filterEfficiency,
                    const int numberEvents) {
 
  SetupDirectories ("", "JetTrackAnalysis/");

  const bool isSignalOnlySample = isMC && TString (inFileName).Contains ("signalonly");
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = nullptr;
  if (file) tree = (TTree*)file->Get ("bush");
  if (file == nullptr || tree == nullptr) {
    cout << "Error: In RtrkTreeMaker.cxx: TTree not obtained for given data set. Quitting." << endl;
    return;
  }

  //First sort jets & tracks into many, smaller TTrees.
  //This is where the sorting based on event information (e.g. centrality, Ntrk, jet pT) will go.
  //Event mixing will take place based on these categories so that total memory usage at any point in time is minimized.

  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetJets ();
  t->SetGetTracks ();
  if (crossSection_microbarns != 0)
    t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetBranchAddresses ();

  TFile* outFile = new TFile (Form ("%s/RtrkStudy/%s.root", rootPath.Data (), identifier.Data ()), "update");

  float jet_pt = 0, jet_eta = 0, jet_phi = 0, jet_e = 0, sum_trk_pt = 0, evtWeight = 0;
  TTree* outTree = new TTree (isSignalOnlySample ? "signalJets" : "overlayJets", isSignalOnlySample ? "signalJets" : "overlayJets");
  outTree->Branch ("jet_pt",      &jet_pt,      "jet_pt/F");
  outTree->Branch ("jet_eta",     &jet_eta,     "jet_eta/F");
  outTree->Branch ("jet_phi",     &jet_phi,     "jet_phi/F");
  outTree->Branch ("jet_e",       &jet_e,       "jet_e/F");
  outTree->Branch ("sum_trk_pt",  &sum_trk_pt,  "sum_trk_pt/F");
  outTree->Branch ("evtWeight",   &evtWeight,   "evtWeight/F");
  
  const int nEvts = tree->GetEntries ();

  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    tree->GetEntry (iEvt);

    // Find leading jet
    int lJ = -1;
    for (int iJ = 0; iJ < t->jet_n; iJ++) {
      if (lJ == -1 || t->jet_pt->at (lJ) < t->jet_pt->at (iJ))
        lJ = iJ;
    }
    if (lJ == -1)
      continue; 

    // Find subleading jet
    int sJ = -1;
    for (int iJ = 0; iJ < t->jet_n; iJ++) {
      if (iJ == lJ)
        continue;
      if (sJ == -1 || t->jet_pt->at (sJ) < t->jet_pt->at (sJ))
        sJ = iJ;
    }
    if (sJ == -1)
      continue;

    evtWeight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;

    // Cuts on jets 
    if (InDisabledHEC (t->jet_eta->at (lJ), t->jet_phi->at (lJ)))
      continue; // leading jet not in disabled HEC
    if (InDisabledHEC (t->jet_eta->at (sJ), t->jet_phi->at (sJ)))
      continue; // subleading jet not in disabled HEC
    
    if (DeltaPhi (t->jet_phi->at (lJ), t->jet_phi->at (sJ)) < 3.*pi / 4.)
      continue; // leading and subleading at least delta phi apart

    if (isPeriodA && t->jet_eta->at (lJ) > 3.2)
      continue; // reject leading jets outside barrel to avoid a centrality bias
    else if (!isPeriodA && t->jet_eta->at (lJ) < -3.2)
      continue;
 
    if (isPeriodA && t->jet_eta->at (sJ) > 3.2)
      continue; // reject subleading jets outside barrel to avoid a centrality bias
    else if (!isPeriodA && t->jet_eta->at (sJ) < -3.2)
      continue;

    // Loop over tracks
    sum_trk_pt = 0;
    for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      // Cuts on tracks
      if (1e-3 * t->trk_pt->at (iTrk) < trk_pt_cut)
        continue; // track pT cut
      if (!t->trk_quality_4->at (iTrk))
        continue; // cut on track quality

      if (DeltaR (t->jet_eta->at (lJ), t->trk_eta->at (iTrk), t->jet_phi->at (lJ), t->trk_phi->at (iTrk)) < 0.4)
        sum_trk_pt += t->trk_pt->at (iTrk);
    }

    jet_pt = t->jet_pt->at (lJ);
    jet_eta = t->jet_eta->at (lJ);
    jet_phi = t->jet_phi->at (lJ);
    jet_e = t->jet_e->at (lJ);

    outTree->Fill ();

  }    

  outFile->Write ();
  outFile->Close ();

}

} // end namespace
