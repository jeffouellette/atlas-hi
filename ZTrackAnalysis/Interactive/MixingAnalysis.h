#ifndef __MixingAnalysis_h__
#define __MixingAnalysis_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace atlashi;

class MixingAnalysis : public PhysicsAnalysis {

  public:

  bool takeNonTruthTracks = false;
  bool doCentMixing = true;
  bool doQ2Mixing = false;
  bool doPsi2Mixing = true;
  bool doPsi3Mixing = false;
  bool doPPTransMinMixing = true; // by default analyses are not performing trans-min mixing. Only really applies to pp bkg.
  bool doPPTransMaxMixing = false; // variation for analyses to use trans-max mixing. Only really applies to pp bkg.

  int nQ2MixBins = 1;
  double* q2MixBins = nullptr;
  int nPsi2MixBins = 16;
  double* psi2MixBins = nullptr;
  int nPsi3MixBins = 1;
  double* psi3MixBins = nullptr;

  short GetQ2MixBin (const float q2) {
    if (!q2MixBins)
      return -1;
    short i = 0;
    while (i < nQ2MixBins) {
      if (q2 < q2MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetPsi2MixBin (const float psi2) {
    if (!psi2MixBins)
      return -1;
    short i = 0;
    while (i < nPsi2MixBins) {
      if (psi2 < psi2MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetPsi3MixBin (const float psi3) {
    if (!psi3MixBins)
      return -1;
    short i = 0;
    while (i < nPsi3MixBins) {
      if (psi3 < psi3MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  MixingAnalysis (const char* _name = "bkg") : PhysicsAnalysis () {
    name = _name;
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    isBkg = true;
    hasBkg = false;
    backgroundSubtracted = true;
    histsUnfolded = true;
    iaaCalculated = true;
    //icpCalculated = true;
    doPPMBMixing = false;
  }

  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  void Execute (const char* inFileName, const char* outFileName) override {
    cout << "Error: In MixingAnalysis::Execute: Called invalid function! 4 arguments are required to specify the collision system and mixing file name. Exiting." << endl;
  }
  void Execute (const bool isPbPb, const char* inFileName, const char* mbInFileName, const char* outFileName) {
    cout << "Warning: In MixingAnalysis::Execute: No mixed event file name specified. Proceeding without saving mixed events." << endl;
    Execute (isPbPb, inFileName, mbInFileName, outFileName, "");
  }
  void Execute (const bool isPbPb, const char* inFileName, const char* mbInFileName, const char* outFileName, const char* mixedFileName);

  void PlotCentralityDists ();
};




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load histograms into memory, then combine channels.
//////////////////////////////////////////////////////////////////////////////////////////////////
void MixingAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);

  PhysicsAnalysis :: CombineHists ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MixingAnalysis :: Execute (const bool isPbPb, const char* inFileName, const char* mbInFileName, const char* outFileName, const char* mixedFileName) {

  LoadEventWeights ();
  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));

  SetupDirectories ("", "ZTrackAnalysis/");

  cout << "Arguments provided: " << endl;
  cout << "isPbPb = " << isPbPb << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "mbInFileName = " << mbInFileName << endl;
  cout << "outFileName = " << outFileName << endl;
  cout << "mixedFileName = " << mixedFileName << endl;

  const bool saveMixedEvents = (string (mixedFileName) != "");
  if (saveMixedEvents) cout << "Note: will save mixed event files!" << endl;

  cout << "Mixing fraction provided is " << mixingFraction << endl;
  const bool doSameFileMixing = (string (inFileName) == string (mbInFileName));
  if (doSameFileMixing) cout << "Attempting to mix events within the same file!" << endl;

  if (doQ2Mixing != (nQ2MixBins > 1)) {
    cout << "|q2| mixing strategy inconsistent, please check nQ2MixBins is consistent with doQ2Mixing." << endl;
    nQ2MixBins = 1;
    doQ2Mixing = false;
  }
  if (doPsi2Mixing != (nPsi2MixBins > 1)) {
    cout << "Psi2 angle mixing strategy inconsistent, please check nPsi2MixBins is consistent with doPsi2Mixing." << endl;
    nPsi2MixBins = 1;
    doPsi2Mixing = false;
  }
  if (doPsi3Mixing != (nPsi3MixBins > 1)) {
    cout << "Psi3 angle mixing strategy inconsistent, please check nPsi3MixBins is consistent with doPsi3Mixing." << endl;
    nPsi3MixBins = 1;
    doPsi3Mixing = false;
  }

  if (doQ2Mixing) {
    cout << "Attempting to mix events in |q2| with " << nQ2MixBins << " bins" << endl;
    q2MixBins = linspace (0, 0.2, nQ2MixBins);
  }
  if (doPsi2Mixing) {
    cout << "Attempting to mix events in psi2 angle with " << nPsi2MixBins << " bins" << endl;
    psi2MixBins = linspace (-pi/2, pi/2, nPsi2MixBins);
  }
  if (doPsi3Mixing) {
    cout << "Attempting to mix events in psi3 angle with " << nPsi3MixBins << " bins" << endl;
    psi3MixBins = linspace (-pi/3, pi/3, nPsi3MixBins);
  }

  CreateHists ();


  // Get TTree with mixed event pool
  TFile* mbInFile = new TFile (Form ("%s/%s", rootPath.Data (), mbInFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), mbInFileName) << endl;
  TTree* mbTree = (TTree*) mbInFile->Get (doSameFileMixing ? (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree") : (isPbPb ? "PbPbTrackTree" : "ppZTrackTree"));
  const int nMBEvts = mbTree->GetEntries ();
  cout << "N events in mixing pool = " << nMBEvts << endl;

  int iMBEvt = 0;
  mbTree->GetEntry (iMBEvt);

  std::vector <int> mbEventOrder = {};
  std::vector <bool> mbEventsUsed = {};
  for (int i = 0; i < nMBEvts; i++) {
    mbEventOrder.push_back (i);
    mbEventsUsed.push_back (false);
  }
  if (!doSameFileMixing) {
    std::srand (1);
    std::random_shuffle (mbEventOrder.begin (), mbEventOrder.end ());

    cout << "Shuffled minimum bias events randomly; first ten events are ";
    for (int i = 0; i < 10; i++) {
      cout << mbEventOrder[i];
      if (i < 9)
        cout << ", ";
    }
    cout << endl;
  }


  // Get TTree with Z bosons
  cout << Form ("Reading Z bosons from %s/%s", rootPath.Data (), inFileName) << endl;
  TFile* inFile = nullptr;
  TTree* zTree = nullptr;
  if (doSameFileMixing)
    zTree = mbTree;
  else {
    inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
    zTree = (TTree*) inFile->Get (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree");
  }
  if (!zTree) cout << "Got a null Z boson tree!" << endl;

  int nZEvts = zTree->GetEntries ();

  if (nZEvts == 0)
    cout << "Warning! No Z's to mix with in this run!" << endl;
  cout << "N Z-tagged events = " << nZEvts << endl;
  cout << "For this Z tree, maximum mixing fraction = " << nMBEvts / nZEvts << endl;

  bool doShuffle = false;
  //if (mixingFraction * nZEvts > nMBEvts && !doSameFileMixing) {
  //  cout << "Warning! Mixing fraction too high, will use " << (float)(nMBEvts / mixingFraction) / (float)(nZEvts) * 100. << "% of Z events" << endl;
  //  nZEvts = nMBEvts / mixingFraction;
  //  doShuffle = true;
  //}

  std::vector <int> zEventOrder = {};
  std::vector <int> zEventsUsed = {};
  for (int i = 0; i < nZEvts; i++) {
    zEventOrder.push_back (i);
    zEventsUsed.push_back (0);
  }

  if (doShuffle) {
    std::srand (1);
    std::random_shuffle (zEventOrder.begin (), zEventOrder.end ());

    cout << "Shuffled Z boson events randomly; first ten events are ";
    for (int i = 0; i < 10; i++) {
      cout << mbEventOrder[i];
      if (i < 9)
        cout << ", ";
    }
    cout << endl;
  }


  //// output TTree to store mixed event
  TFile* mixedEventsFile = (saveMixedEvents ? new TFile (Form ("%s/%s", rootPath.Data (), mixedFileName), "recreate") : nullptr);
  TTree* mixedEventsTree = (saveMixedEvents ? new TTree (isPbPb ? "PbPbMixedTree" : "ppMixedTree", isPbPb ? "PbPbMixedTree" : "ppMixedTree") : nullptr);


  // variables for branches
  int event_number = 0, z_event_number = 0;
  int run_number = 0, z_run_number = 0;
  bool isEE = false;//, passes_toroid = false;
  bool _isEE = false;
  float event_weight = 1., z_event_weight = 1., fcal_weight = 1., q2_weight = 1., psi2_weight = 1.;
  float fcal_et = 0, vz = 0, zdcEnergy = 0, z_fcal_et = 0, z_vz = 0, z_zdcEnergy = 0;
  float phi_transmin = 0, phi_transmax = 0;
  //float q2x_a = 0, q2y_a = 0, q2x_c = 0, q2y_c = 0, z_q2x_a = 0, z_q2y_a = 0, z_q2x_c = 0, z_q2y_c = 0;
  float q2 = 0, q3 = 0, q4 = 0, z_q2 = 0, z_q3 = 0, z_q4 = 0;
  float psi2 = 0, psi3 = 0, psi4 = 0, z_psi2 = 0, z_psi3 = 0, z_psi4 = 0;
  float _z_pt = 0, _z_phi = 0, _z_y = 0, _z_m = 0;
  float z_pt = 0, z_phi = 0, z_y = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  float _l1_pt = 0, _l1_eta = 0, _l1_phi = 0, _l1_trk_pt = 0, _l1_trk_eta = 0, _l1_trk_phi = 0, _l2_pt = 0, _l2_eta = 0, _l2_phi = 0, _l2_trk_pt = 0, _l2_trk_eta = 0, _l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  int _l1_charge = 0, _l2_charge = 0;
  int ntrk = 0, z_ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  bool trk_truth_matched[10000];
  float z_trk_pt[10000], z_trk_eta[10000], z_trk_phi[10000];
  bool HLT_mb_sptrk_L1ZDC_A_C_VTE50 = false;
  bool HLT_noalg_pc_L1TE50_VTE600_0ETA49 = false;
  bool HLT_noalg_cc_L1TE600_0ETA49 = false;

  int***    trks_counts   = Get3DArray <int> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (maxNPtchBins, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (maxNPtchBins, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (maxNPtchBins, 40);

  int nMixed = 0, nNotMixed = 0;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is PbPb
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (isPbPb) {
    mbTree->LoadBaskets (4000000000); //2000000000 = 2GB
    mbTree->SetBranchAddress ("run_number",   &run_number);
    mbTree->SetBranchAddress ("event_number", &event_number);
    mbTree->SetBranchAddress ("event_weight", &event_weight);
    mbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    mbTree->SetBranchAddress ("zdcEnergy",    &zdcEnergy);
    //mbTree->SetBranchAddress ("q2x_a",         &q2x_a);
    //mbTree->SetBranchAddress ("q2y_a",         &q2y_a);
    //mbTree->SetBranchAddress ("q2x_c",         &q2x_c);
    //mbTree->SetBranchAddress ("q2y_c",         &q2y_c);
    mbTree->SetBranchAddress ("q2",           &q2);
    mbTree->SetBranchAddress ("psi2",         &psi2);
    mbTree->SetBranchAddress ("q3",           &q3);
    mbTree->SetBranchAddress ("psi3",         &psi3);
    mbTree->SetBranchAddress ("q4",           &q4);
    mbTree->SetBranchAddress ("psi4",         &psi4);
    mbTree->SetBranchAddress ("vz",           &vz);
    mbTree->SetBranchAddress ("ntrk",         &ntrk);
    mbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    mbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    mbTree->SetBranchAddress ("trk_phi",      &trk_phi);
    mbTree->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",       &HLT_mb_sptrk_L1ZDC_A_C_VTE50);
    mbTree->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49",  &HLT_noalg_pc_L1TE50_VTE600_0ETA49);
    mbTree->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49",        &HLT_noalg_cc_L1TE600_0ETA49);

    if (takeNonTruthTracks)
      mbTree->SetBranchAddress ("trk_truth_matched", &trk_truth_matched);


    if (!doSameFileMixing) {
      zTree->LoadBaskets (2000000000);
      zTree->SetBranchAddress ("run_number",    &z_run_number);
      zTree->SetBranchAddress ("event_number",  &z_event_number);
      zTree->SetBranchAddress ("event_weight",  &z_event_weight);
      zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
      zTree->SetBranchAddress ("zdcEnergy",     &z_zdcEnergy);
      //zTree->SetBranchAddress ("q2x_a",          &z_q2x_a);
      //zTree->SetBranchAddress ("q2y_a",          &z_q2y_a);
      //zTree->SetBranchAddress ("q2x_c",          &z_q2x_c);
      //zTree->SetBranchAddress ("q2y_c",          &z_q2y_c);
      zTree->SetBranchAddress ("q2",            &z_q2);
      zTree->SetBranchAddress ("psi2",          &z_psi2);
      zTree->SetBranchAddress ("q3",            &z_q3);
      zTree->SetBranchAddress ("psi3",          &z_psi3);
      zTree->SetBranchAddress ("q4",            &z_q4);
      zTree->SetBranchAddress ("psi4",          &z_psi4);
      zTree->SetBranchAddress ("vz",            &z_vz);
      zTree->SetBranchAddress ("ntrk",          &z_ntrk);
      zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
      zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
      zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    }
    zTree->SetBranchAddress ("isEE",          &_isEE);
    zTree->SetBranchAddress ("z_pt",          &_z_pt);
    zTree->SetBranchAddress ("z_phi",         &_z_phi);
    zTree->SetBranchAddress ("z_y",           &_z_y);
    zTree->SetBranchAddress ("z_m",           &_z_m);
    zTree->SetBranchAddress ("l1_pt",         &_l1_pt);
    zTree->SetBranchAddress ("l1_eta",        &_l1_eta);
    zTree->SetBranchAddress ("l1_phi",        &_l1_phi);
    zTree->SetBranchAddress ("l1_trk_pt",     &_l1_trk_pt);
    zTree->SetBranchAddress ("l1_trk_eta",    &_l1_trk_eta);
    zTree->SetBranchAddress ("l1_trk_phi",    &_l1_trk_phi);
    zTree->SetBranchAddress ("l1_charge",     &_l1_charge);
    zTree->SetBranchAddress ("l2_pt",         &_l2_pt);
    zTree->SetBranchAddress ("l2_eta",        &_l2_eta);
    zTree->SetBranchAddress ("l2_phi",        &_l2_phi);
    zTree->SetBranchAddress ("l2_trk_pt",     &_l2_trk_pt);
    zTree->SetBranchAddress ("l2_trk_eta",    &_l2_trk_eta);
    zTree->SetBranchAddress ("l2_trk_phi",    &_l2_trk_phi);
    zTree->SetBranchAddress ("l2_charge",     &_l2_charge);

    if (saveMixedEvents) {
      mixedEventsTree->Branch ("run_number",      &run_number,      "run_number/I");
      mixedEventsTree->Branch ("event_number",    &event_number,    "event_number/I");
      mixedEventsTree->Branch ("fcal_et",         &fcal_et,         "fcal_et/F");
      mixedEventsTree->Branch ("zdcEnergy",       &zdcEnergy,       "zdcEnergy/F");
      mixedEventsTree->Branch ("q2",              &q2,              "q2/F");
      mixedEventsTree->Branch ("psi2",            &psi2,            "psi2/F");
      mixedEventsTree->Branch ("vz",              &vz,              "vz/F");
      mixedEventsTree->Branch ("ntrk",            &ntrk,            "ntrk/I");
      //mixedEventsTree->Branch ("trk_pt",          &trk_pt,          "trk_pt[ntrk]/F");
      //mixedEventsTree->Branch ("trk_eta",         &trk_eta,         "trk_eta[ntrk]/F");
      //mixedEventsTree->Branch ("trk_phi",         &trk_phi,         "trk_phi[ntrk]/F");
      mixedEventsTree->Branch ("isEE",            &isEE,            "isEE/O");
      mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/I");
      mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/I");
      mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
      mixedEventsTree->Branch ("z_fcal_et",       &z_fcal_et,       "z_fcal_et/F");
      mixedEventsTree->Branch ("z_zdcEnergy",     &z_zdcEnergy,     "z_zdcEnergy/F");
      mixedEventsTree->Branch ("z_q2",            &z_q2,            "z_q2/F");
      mixedEventsTree->Branch ("z_psi2",          &z_psi2,          "z_psi2/F");
      mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
      mixedEventsTree->Branch ("z_pt",            &z_pt,            "z_pt/F");
      mixedEventsTree->Branch ("z_phi",           &z_phi,           "z_phi/F");
      mixedEventsTree->Branch ("z_y",             &z_y,             "z_y/F");
      mixedEventsTree->Branch ("z_m",             &z_m,             "z_m/F");
      mixedEventsTree->Branch ("l1_pt",           &l1_pt,           "l1_pt/F");
      mixedEventsTree->Branch ("l1_eta",          &l1_eta,          "l1_eta/F");
      mixedEventsTree->Branch ("l1_phi",          &l1_phi,          "l1_phi/F");
      mixedEventsTree->Branch ("l1_trk_pt",       &l1_trk_pt,       "l1_trk_pt/F");
      mixedEventsTree->Branch ("l1_trk_eta",      &l1_trk_eta,      "l1_trk_eta/F");
      mixedEventsTree->Branch ("l1_trk_phi",      &l1_trk_phi,      "l1_trk_phi/F");
      mixedEventsTree->Branch ("l1_charge",       &l1_charge,       "l1_charge/I");
      mixedEventsTree->Branch ("l2_pt",           &l2_pt,           "l2_pt/F");
      mixedEventsTree->Branch ("l2_eta",          &l2_eta,          "l2_eta/F");
      mixedEventsTree->Branch ("l2_phi",          &l2_phi,          "l2_phi/F");
      mixedEventsTree->Branch ("l2_trk_pt",       &l2_trk_pt,       "l2_trk_pt/F");
      mixedEventsTree->Branch ("l2_trk_eta",      &l2_trk_eta,      "l2_trk_eta/F");
      mixedEventsTree->Branch ("l2_trk_phi",      &l2_trk_phi,      "l2_trk_phi/F");
      mixedEventsTree->Branch ("l2_charge",       &l2_charge,       "l2_charge/I");
      if (!doSameFileMixing) {
        mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
        mixedEventsTree->Branch ("z_trk_pt",        &z_trk_pt,        "z_trk_pt[z_ntrk]/F");
        mixedEventsTree->Branch ("z_trk_eta",       &z_trk_eta,       "z_trk_eta[z_ntrk]/F");
        mixedEventsTree->Branch ("z_trk_phi",       &z_trk_phi,       "z_trk_phi[z_ntrk]/F");
      }
    }


    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (zEventOrder[iZEvt % nZEvts]);
      zEventsUsed[iZEvt % nZEvts]++;

      isEE = _isEE;
      z_pt = _z_pt;
      z_phi = _z_phi;
      z_y = _z_y;
      z_m = _z_m;
      if (saveMixedEvents) {
        l1_pt = _l1_pt;
        l1_eta = _l1_eta;
        l1_phi = _l1_phi;
        l1_trk_pt = _l1_trk_pt;
        l1_trk_eta = _l1_trk_eta;
        l1_trk_phi = _l1_trk_phi;
        l1_charge = _l1_charge;
        l2_pt = _l2_pt;
        l2_eta = _l2_eta;
        l2_phi = _l2_phi;
        l2_trk_pt = _l2_trk_pt;
        l2_trk_eta = _l2_trk_eta;
        l2_trk_phi = _l2_trk_phi;
        l2_charge = _l2_charge;
      }
      if (doSameFileMixing) {
        z_run_number = run_number;
        z_event_number = event_number;
        z_event_weight = event_weight;
        z_fcal_et = fcal_et;
        z_zdcEnergy = zdcEnergy;
        //z_q2x_a = q2x_a;
        //z_q2x_c = q2x_c;
        //z_q2y_a = q2y_a;
        //z_q2y_c = q2y_c;
        z_q2 = q2;
        z_psi2 = psi2;
        z_q3 = q3;
        z_psi3 = psi3;
        z_q4 = q4;
        z_psi4 = psi4;
        z_vz = vz;
        z_ntrk = 0; // placeholder
      }
     

      if (fabs (z_vz) > 150) continue;

      if (z_event_weight == 0) continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      if (iPtZ < 2) continue; // skip unneeded events


      // Find the next unused minimum bias event
      {
        const short iFCalEt = GetSuperFineCentBin (z_fcal_et);
        if (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1) continue;
        const short iQ2 = GetQ2MixBin (z_q2);
        if (doQ2Mixing && (iQ2 < 0 || iQ2 > nQ2MixBins-1)) {
          cout << "Out-of-bounds q2, skipping this Z!" << endl;
          continue;
        }

        // disable big branches during event mixing
        mbTree->SetBranchStatus ("trk_pt", 0);
        mbTree->SetBranchStatus ("trk_eta", 0);
        mbTree->SetBranchStatus ("trk_phi", 0);
        
        bool goodMixEvent = false;
        const int _iMBEvt = iMBEvt;
        do {
          iMBEvt = (iMBEvt+1) % nMBEvts;
          mbTree->GetEntry (mbEventOrder[iMBEvt]);
          goodMixEvent = (fabs (vz) < 150 && event_weight != 0); // always require these conditions
          //goodMixEvent &= (!doSameFileMixing || mixingFraction != 1 || !mbEventsUsed[iMBEvt]); // checks for uniqueness (if applicable)
          goodMixEvent &= (!doSameFileMixing || iMBEvt != iZEvt); // don't mix with the exact same event
          goodMixEvent &= (!doCentMixing || iFCalEt == GetSuperFineCentBin (fcal_et)); // do centrality matching
          goodMixEvent &= (!doQ2Mixing   || iQ2 == GetQ2MixBin (q2)); // do q2 matching
          goodMixEvent &= (!doPsi2Mixing || DeltaPhi (psi2, z_psi2) < (pi / nPsi2MixBins)); // do psi2 matching
          goodMixEvent &= (!doPsi3Mixing || DeltaPhi (psi3, z_psi3) < (2.*pi/3. / nPsi3MixBins)); // do psi3 matching
        } while (!goodMixEvent && iMBEvt != _iMBEvt); // only check each event once
        if (_iMBEvt == iMBEvt) {
          //cout << "No minbias event to mix with!!! Wrapped around on the same Z!!! Sum Et = " << z_fcal_et << ", q2 = " << z_q2 << ", psi2 = " << z_psi2 << endl;
          nNotMixed++;
          continue;
        }
        nMixed++;
        mbEventsUsed[iMBEvt] = true;

        // reenable big branches after finding a suitable event
        mbTree->SetBranchStatus ("trk_pt", 1);
        mbTree->SetBranchStatus ("trk_eta", 1);
        mbTree->SetBranchStatus ("trk_phi", 1);
        mbTree->GetEntry (mbEventOrder[iMBEvt]);
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      if (saveMixedEvents) mixedEventsTree->Fill ();

      const short iSpc = (isEE ? 0 : 1); // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = GetCentBin (z_fcal_et);
      if (iCent < 1 || iCent > numCentBins-1) continue;
      const short iFineCent = GetFineCentBin (z_fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1) continue;

      // do a reweighting procedure
      {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)  dphi = pi - dphi;
        if (useCentWgts)  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (z_fcal_et));
        if (useQ2Wgts)    q2_weight   = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
        if (usePsi2Wgts)  psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

        z_event_weight *= fcal_weight * q2_weight * psi2_weight;
      }

      event_weight = z_event_weight;

      // fill event category histograms
      {
        const int centrality = GetPercentileCentrality (fcal_et);
        h_fcal_et->Fill (fcal_et);
        h_fcal_et_reweighted->Fill (fcal_et, event_weight);
        if (HLT_mb_sptrk_L1ZDC_A_C_VTE50) {
          h_centrality[0]->Fill (centrality);
          h_centrality_reweighted[0]->Fill (centrality, event_weight);
        }
        if (HLT_noalg_pc_L1TE50_VTE600_0ETA49) {
          h_centrality[1]->Fill (centrality);
          h_centrality_reweighted[1]->Fill (centrality, event_weight);

        }
        if (HLT_noalg_cc_L1TE600_0ETA49) {
          h_centrality[2]->Fill (centrality);
          h_centrality_reweighted[2]->Fill (centrality, event_weight);
        }
        h_q2[iFineCent]->Fill (q2);
        h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
        h_psi2[iFineCent]->Fill (psi2);
        h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
        h_PbPb_vz->Fill (vz);
        h_PbPb_vz_reweighted->Fill (vz, event_weight);
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkPur == 0 || trkEff == 0) continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {   
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {   
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks
      
      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (iPtch, iCent); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1
      
      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over Pb+Pb tree

    cout << endl << "Done minbias Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is pp
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (!isPbPb) {
    assert (doSameFileMixing == !doPPMBMixing);

    mbTree->LoadBaskets (3000000000);
    mbTree->SetBranchAddress ("run_number",   &run_number);
    mbTree->SetBranchAddress ("event_number", &event_number);
    mbTree->SetBranchAddress ("event_weight", &event_weight);
    mbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    mbTree->SetBranchAddress ("vz",           &vz);
    mbTree->SetBranchAddress ("ntrk",         &ntrk);
    mbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    mbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    mbTree->SetBranchAddress ("trk_phi",      &trk_phi);

    if (!doSameFileMixing) {
      zTree->LoadBaskets (3000000000);
      zTree->SetBranchAddress ("run_number",    &z_run_number);
      zTree->SetBranchAddress ("event_number",  &z_event_number);
      zTree->SetBranchAddress ("event_weight",  &z_event_weight);
      zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
      //zTree->SetBranchAddress ("zdcEnergy",     &z_zdcEnergy);
      //zTree->SetBranchAddress ("q2",            &z_q2);
      //zTree->SetBranchAddress ("psi2",          &z_psi2);
      zTree->SetBranchAddress ("vz",            &z_vz);
      zTree->SetBranchAddress ("ntrk",          &z_ntrk);
      zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
      zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
      zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    }
    if (doPPTransMinMixing || doPPTransMaxMixing) {
      zTree->SetBranchAddress ("phi_transmin",  &phi_transmin);
      zTree->SetBranchAddress ("phi_transmax",  &phi_transmax);
    }
    zTree->SetBranchAddress ("isEE",          &_isEE);
    zTree->SetBranchAddress ("z_pt",          &_z_pt);
    zTree->SetBranchAddress ("z_phi",         &_z_phi);
    zTree->SetBranchAddress ("z_y",           &_z_y);
    zTree->SetBranchAddress ("z_m",           &_z_m);
    //zTree->SetBranchAddress ("l1_pt",         &_l1_pt);
    //zTree->SetBranchAddress ("l1_eta",        &_l1_eta);
    //zTree->SetBranchAddress ("l1_phi",        &_l1_phi);
    //zTree->SetBranchAddress ("l1_trk_pt",     &_l1_trk_pt);
    //zTree->SetBranchAddress ("l1_trk_eta",    &_l1_trk_eta);
    //zTree->SetBranchAddress ("l1_trk_phi",    &_l1_trk_phi);
    //zTree->SetBranchAddress ("l1_charge",     &_l1_charge);
    //zTree->SetBranchAddress ("l2_pt",         &_l2_pt);
    //zTree->SetBranchAddress ("l2_eta",        &_l2_eta);
    //zTree->SetBranchAddress ("l2_phi",        &_l2_phi);
    //zTree->SetBranchAddress ("l2_trk_pt",     &_l2_trk_pt);
    //zTree->SetBranchAddress ("l2_trk_eta",    &_l2_trk_eta);
    //zTree->SetBranchAddress ("l2_trk_phi",    &_l2_trk_phi);
    //zTree->SetBranchAddress ("l2_charge",     &_l2_charge);

    if (saveMixedEvents) {
      mixedEventsTree->Branch ("run_number",      &run_number,      "run_number/i");
      mixedEventsTree->Branch ("event_number",    &event_number,    "event_number/i");
      mixedEventsTree->Branch ("isEE",            &isEE,            "isEE/O");
      mixedEventsTree->Branch ("fcal_et",         &fcal_et,         "fcal_et/F");
      mixedEventsTree->Branch ("vz",              &vz,              "vz/F");
      mixedEventsTree->Branch ("ntrk",            &ntrk,            "ntrk/I");
      //mixedEventsTree->Branch ("trk_pt",          &trk_pt,          "trk_pt[ntrk]/F");
      //mixedEventsTree->Branch ("trk_eta",         &trk_eta,         "trk_eta[ntrk]/F");
      //mixedEventsTree->Branch ("trk_phi",         &trk_phi,         "trk_phi[ntrk]/F");
      mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/i");
      mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/i");
      mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
      mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
      mixedEventsTree->Branch ("z_pt",            &z_pt,            "z_pt/F");
      mixedEventsTree->Branch ("z_phi",           &z_phi,           "z_phi/F");
      mixedEventsTree->Branch ("z_y",             &z_y,             "z_y/F");
      mixedEventsTree->Branch ("z_m",             &z_m,             "z_m/F");
      //mixedEventsTree->Branch ("l1_pt",           &l1_pt,           "l1_pt/F");
      //mixedEventsTree->Branch ("l1_eta",          &l1_eta,          "l1_eta/F");
      //mixedEventsTree->Branch ("l1_phi",          &l1_phi,          "l1_phi/F");
      //mixedEventsTree->Branch ("l1_trk_pt",       &l1_trk_pt,       "l1_trk_pt/F");
      //mixedEventsTree->Branch ("l1_trk_eta",      &l1_trk_eta,      "l1_trk_eta/F");
      //mixedEventsTree->Branch ("l1_trk_phi",      &l1_trk_phi,      "l1_trk_phi/F");
      //mixedEventsTree->Branch ("l1_charge",       &l1_charge,       "l1_charge/I");
      //mixedEventsTree->Branch ("l2_pt",           &l2_pt,           "l2_pt/F");
      //mixedEventsTree->Branch ("l2_eta",          &l2_eta,          "l2_eta/F");
      //mixedEventsTree->Branch ("l2_phi",          &l2_phi,          "l2_phi/F");
      //mixedEventsTree->Branch ("l2_trk_pt",       &l2_trk_pt,       "l2_trk_pt/F");
      //mixedEventsTree->Branch ("l2_trk_eta",      &l2_trk_eta,      "l2_trk_eta/F");
      //mixedEventsTree->Branch ("l2_trk_phi",      &l2_trk_phi,      "l2_trk_phi/F");
      //mixedEventsTree->Branch ("l2_charge",       &l2_charge,       "l2_charge/I");
      mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
      //mixedEventsTree->Branch ("z_trk_pt",        &z_trk_pt,        "z_trk_pt[z_ntrk]/F");
      //mixedEventsTree->Branch ("z_trk_eta",       &z_trk_eta,       "z_trk_eta[z_ntrk]/F");
      //mixedEventsTree->Branch ("z_trk_phi",       &z_trk_phi,       "z_trk_phi[z_ntrk]/F");
    }


    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (zEventOrder[iZEvt % nZEvts]);
      zEventsUsed[iZEvt % nZEvts]++;

      isEE = _isEE;
      z_pt = _z_pt;
      z_phi = _z_phi;
      z_y = _z_y;
      z_m = _z_m;
      //if (saveMixedEvents) {
      //  l1_pt = _l1_pt;
      //  l1_eta = _l1_eta;
      //  l1_phi = _l1_phi;
      //  l1_trk_pt = _l1_trk_pt;
      //  l1_trk_eta = _l1_trk_eta;
      //  l1_trk_phi = _l1_trk_phi;
      //  l1_charge = _l1_charge;
      //  l2_pt = _l2_pt;
      //  l2_eta = _l2_eta;
      //  l2_phi = _l2_phi;
      //  l2_trk_pt = _l2_trk_pt;
      //  l2_trk_eta = _l2_trk_eta;
      //  l2_trk_phi = _l2_trk_phi;
      //  l2_charge = _l2_charge;
      //}
      if (doSameFileMixing) {
        z_run_number = run_number;
        z_event_number = event_number;
        z_event_weight = event_weight;
        z_fcal_et = fcal_et;
        //z_zdcEnergy = zdcEnergy;
        //z_q2x_a = q2x_a;
        //z_q2x_c = q2x_c;
        //z_q2y_a = q2y_a;
        //z_q2y_c = q2y_c;
        z_vz = vz;
        z_ntrk = 0; // placeholder
      }

      if (fabs (z_vz) > 150) continue;

      if (z_event_weight == 0) continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      if (iPtZ < 2) continue; // skip unneeded events

      // Find the next unused minimum bias event
      {
        // disable big branches during event mixing
        mbTree->SetBranchStatus ("trk_pt", 0);
        mbTree->SetBranchStatus ("trk_eta", 0);
        mbTree->SetBranchStatus ("trk_phi", 0);

        bool goodMixEvent = false;
        const int _iMBEvt = iMBEvt;
        do {
          iMBEvt = (iMBEvt+1) % nMBEvts;
          mbTree->GetEntry (mbEventOrder[iMBEvt]);
          goodMixEvent = (fabs (vz) < 150 && event_weight != 0); // always require these conditions
          goodMixEvent &= (doSameFileMixing || mixingFraction != 1 || !mbEventsUsed[iMBEvt]); // checks for uniqueness (if applicable)
          goodMixEvent &= (!doSameFileMixing || iMBEvt != iZEvt); // don't mix with the exact same event
          goodMixEvent &= (!doPPTransMinMixing || (1 < _z_pt && _z_pt < 12 && DeltaPhi (phi_transmin, z_phi) >= 7.*pi/8.));
          goodMixEvent &= (!doPPTransMaxMixing || (1 < _z_pt && _z_pt < 12 && DeltaPhi (phi_transmax, z_phi) >= 7.*pi/8.)); // alternatively mix with trans-max region
        } while (!goodMixEvent && iMBEvt != _iMBEvt); // only check each event once
        if (_iMBEvt == iMBEvt) {
          //cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          nNotMixed++;
          continue;
        }
        nMixed++;
        mbEventsUsed[iMBEvt] = true;

        // reenable big branches after finding a suitable event
        mbTree->SetBranchStatus ("trk_pt", 1);
        mbTree->SetBranchStatus ("trk_eta", 1);
        mbTree->SetBranchStatus ("trk_phi", 1);
        mbTree->GetEntry (mbEventOrder[iMBEvt]);
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      if (saveMixedEvents) mixedEventsTree->Fill ();

      event_weight = z_event_weight;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkPur == 0 || trkEff == 0) continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {   
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {   
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if ((!doPPMBMixing && 7.*pi/8. <= dPhi) || (doPPMBMixing && 3.*pi/4. <= dPhi)) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks
      
      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (iPtch, iCent); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1
      
      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over pp tree

    cout << endl << "Done minbias pp loop." << endl;
  }

  cout << "Number of events mixed:     " << nMixed << endl;
  cout << "Number of events not mixed: " << nNotMixed << " (" << ((float)nNotMixed) * 100 / ((float)nNotMixed + (float)nMixed) << "\%)" << endl;

  for (int i = 0; i < nZEvts; i++) {
    if (zEventsUsed[i] < mixingFraction) {
      cout << "Warning! Some events used fewer than " << mixingFraction << " times!" << endl;
      break;
    }
  }

  if (saveMixedEvents) {
    mixedEventsTree->SetDirectory (mixedEventsFile);
    mixedEventsTree->Write ("", TObject :: kOverwrite);

    if (mixedEventsFile) mixedEventsFile->Close ();
    SaferDelete (&mixedEventsFile);
  }

  Delete3DArray (&trks_counts, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights1, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights2, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete2DArray (&trks_counts_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights1_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights2_inPhi, maxNPtchBins, 40);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (&inFile);

  if (mbInFile) mbInFile->Close ();
  SaferDelete (&mbInFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot centrality distributions for different MB triggers
////////////////////////////////////////////////////////////////////////////////////////////////
void MixingAnalysis :: PlotCentralityDists () {
  const char* canvasName = "c_centrality";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600, 600);
    gDirectory->Add (c);
    c->cd ();
  }

  c->cd ();

  c->SetLogy ();

  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig]->Scale (1., "width");

    h_centrality[iMBTrig]->GetYaxis ()->SetRangeUser (0.5, 1e5);
    h_centrality[iMBTrig]->SetLineColor (colors[iMBTrig+1]);

    h_centrality[iMBTrig]->GetXaxis ()->SetTitle ("Centrality [%]");
    h_centrality[iMBTrig]->GetYaxis ()->SetTitle ("N_{evt}");

//    h_centrality[iMBTrig]->GetYaxis ()->SetRangeUser (5e-2, 2e3);

    h_centrality[iMBTrig]->Draw (!canvasExists && iMBTrig == 0 ? "hist" : "same hist");
  }

  myMarkerTextNoLine (0.22, 0.35, colors[1], kFullCircle, "HLT_mb_sptrk_L1ZDC_A_C_VTE50 (PC)", 1.25, 0.04);
  myMarkerTextNoLine (0.22, 0.30, colors[2], kFullCircle, "HLT_noalg_pc_L1TE50_VTE600_0ETA49 (PC)", 1.25, 0.04);
  myMarkerTextNoLine (0.22, 0.25, colors[3], kFullCircle, "HLT_noalg_cc_L1TE600_0ETA49 (CC)", 1.25, 0.04);

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.22, 0.81, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);

  c->SaveAs (Form ("%s/CentralityDist.pdf", plotPath.Data ()));
}

#endif
