#ifndef __DoMixing_cxx__
#define __DoMixing_cxx__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_set>

using namespace std;
using namespace atlashi;

bool doCentMixing = true;
bool doQ2Mixing = false;
bool doPsi2Mixing = true;
bool doPsi3Mixing = false;
bool doPPTransMinMixing = true; // by default analyses are not performing trans-min mixing. Only really applies to pp bkg.
bool doPPTransMaxMixing = false; // variation for analyses to use trans-max mixing. Only really applies to pp bkg.
bool doVZMixing = false;
bool useImpactParameter = true;

int mixingFraction = 10;

int nQ2MixBins = 1;
double* q2MixBins = nullptr;
int nPsi2MixBins = 16;
double* psi2MixBins = nullptr;
int nPsi3MixBins = 1;
double* psi3MixBins = nullptr;
const int nIPBins = 156;//2640;//508;//298;//478;//298;//156;
double* ipBins = nullptr;

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

short GetIPBin (const float ip) {
  if (!useImpactParameter)
    return -1;
  short i = 0;
  while (i < nIPBins-1) {
    if (ip < ipBins[i+1])
      break;
    i++;
  }
  return i;
}

void GenerateIPMixBins () {
  // impact parameter centralities
  // 0%     1%      10%     20%     30%     40%     50%     60%     70%     80%
  // 0      1.564   4.953   7.009   8.582   9.911   11.083  12.136  13.111  14.032
  int i = 0;
  ipBins = new double[nIPBins+1];

  for (i = 0; i < 10; i++)
    ipBins[i] = i*(1.564)/10.; // 0-1% central
  for (i = 10; i < 36; i++)
    ipBins[i] = 1.564 + (i-10)*(4.953-1.564)/36.; // 1-10% central
  for (i = 36; i < 126; i++)
    ipBins[i] = 4.953 + (i-36)*(11.083-4.953)/80.; // 10-50% central
  for (i = 126; i < 156; i++)
    ipBins[i] = 11.083 + (i-126)*(14.032-11.083)/30.; // 50-80% central
  ipBins[156] = 14.032;

  //for (i = 0; i < 160; i++)
  //  ipBins[i] = i*(1.564)/160.; // 0-1% central
  //for (i = 160; i < 880; i++)
  //  ipBins[i] = 1.564 + (i-160)*(4.953-1.564)/720.; // 1-10% central
  //for (i = 880; i < 2160; i++)
  //  ipBins[i] = 4.953 + (i-880)*(11.083-4.953)/1280.; // 10-50% central
  //for (i = 2160; i < 2640; i++)
  //  ipBins[i] = 11.083 + (i-2160)*(14.032-11.083)/480.; // 50-80% central
  //ipBins[2640] = 14.032;

  //for (i = 0; i < 40; i++)
  //  ipBins[i] = i*(1.564)/40.;
  //for (i = 40; i < 188; i++)
  //  ipBins[i] = 1.564 + (i-40)*(4.953-1.564)/144.;
  //for (i = 188; i < 268; i++)
  //  ipBins[i] = 4.953 + (i-188)*(11.083-4.953)/80.;
  //for (i = 268; i < 298; i++)
  //  ipBins[i] = 11.083 + (i-268)*(14.032-11.083)/30.;
  //ipBins[298] = 14.032;

  //for (i = 0; i < 80; i++)
  //  ipBins[i] = i*(1.564)/80.;
  //for (i = 80; i < 368; i++)
  //  ipBins[i] = 1.564 + (i-80)*(4.953-1.564)/288.;
  //for (i = 368; i < 448; i++)
  //  ipBins[i] = 4.953 + (i-368)*(11.083-4.953)/80.;
  //for (i = 448; i < 508; i++)
  //  ipBins[i] = 11.083 + (i-448)*(14.032-11.083)/60.;
  //ipBins[508] = 14.032;
}




void DoMixing (const bool isPbPb, const char* inFileName, const char* bkgInFileName, const char* mixedFileName) {

  SetupDirectories ("", "ZTrackAnalysis/");

  cout << "Arguments provided: " << endl;
  cout << "isPbPb = " << isPbPb << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "bkgInFileName = " << bkgInFileName << endl;
  cout << "mixedFileName = " << mixedFileName << endl;

  cout << "Mixing fraction provided is " << mixingFraction << endl;
  const bool doSameFileMixing = (string (inFileName) == string (bkgInFileName));
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
  if (useImpactParameter) {
    cout << "Resorting to impact parameter matched mixing with " << nIPBins << " bins in b" << endl;
  //  ipBins = linspace (0, 25, nIPBins);
    GenerateIPMixBins ();
    //for (int i = 0; i <= nIPBins; i++) {
    //  cout << ipBins[i] << endl;
    //}
  }


  // Get TTree with mixed event pool
  TFile* bkgInFile = new TFile (Form ("%s/%s", rootPath.Data (), bkgInFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), bkgInFileName) << endl;
  TTree* bkgTree = (TTree*) bkgInFile->Get (doSameFileMixing ? (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree") : (isPbPb ? "PbPbTrackTree" : "ppZTrackTree"));
  const int nBkgEvts = bkgTree->GetEntries ();
  cout << "N events in mixing pool = " << nBkgEvts << endl;

  int iBkgEvt = 0;
  bkgTree->GetEntry (iBkgEvt);

  vector <int> bkgEventOrder = {};
  vector <bool> bkgEventsUsed = {};
  for (int i = 0; i < nBkgEvts; i++) {
    bkgEventOrder.push_back (i);
    bkgEventsUsed.push_back (false);
  }
  if (!doSameFileMixing) {
    std::srand (1);
    std::random_shuffle (bkgEventOrder.begin (), bkgEventOrder.end ());

    cout << "Shuffled minimum bias events randomly; first ten events are ";
    for (int i = 0; i < 10; i++) {
      cout << bkgEventOrder[i];
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
    zTree = bkgTree;
  else {
    inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
    zTree = (TTree*) inFile->Get (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree");
  }
  if (!zTree) cout << "Got a null Z boson tree!" << endl;

  int nZEvts = zTree->GetEntries ();

  if (nZEvts == 0)
    cout << "Warning! No Z's to mix with in this run!" << endl;
  cout << "N Z-tagged events = " << nZEvts << endl;
  cout << "For this Z tree, maximum mixing fraction = " << nBkgEvts / nZEvts << endl;

  vector <int> zEventsUsed = {};
  vector <vector<int>> mixingMap = vector <vector<int>> (0);
  for (int i = 0; i < nZEvts; i++) {
    zEventsUsed.push_back (0);
    mixingMap.push_back (vector<int> (0));
  }

  // output TTree to store mixed events
  TFile* mixedEventsFile = new TFile (Form ("%s/%s", rootPath.Data (), mixedFileName), "recreate");
  TTree* mixedEventsTree = new TTree (isPbPb ? "PbPbMixedTree" : "ppMixedTree", isPbPb ? "PbPbMixedTree" : "ppMixedTree");


  // variables for branches
  int event_number = 0, z_event_number = 0;
  int run_number = 0, z_run_number = 0;
  bool isEE = false, _isEE = false;//, passes_toroid = false;
  float event_weight = 1., z_event_weight = 1., fcal_weight = 1., q2_weight = 1., psi2_weight = 1.;
  float fcal_et = 0, vz = 0, zdcEnergy = 0, z_fcal_et = 0, z_vz = 0, z_zdcEnergy = 0;
  int ntrk_perp = 0, z_ntrk_perp = 0;
  float ip = 0, eventPlane = 0, z_ip = 0, z_eventPlane = 0;
  float phi_transmin = 0, phi_transmax = 0;
  float q2 = 0, q3 = 0, q4 = 0, z_q2 = 0, z_q3 = 0, z_q4 = 0;
  float psi2 = 0, psi3 = 0, psi4 = 0, z_psi2 = 0, z_psi3 = 0, z_psi4 = 0;
  float z_pt = 0, z_phi = 0, z_y = 0, z_m = 0;
  float _z_pt = 0, _z_phi = 0, _z_y = 0, _z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  //float _l1_pt = 0, _l1_eta = 0, _l1_phi = 0, _l1_trk_pt = 0, _l1_trk_eta = 0, _l1_trk_phi = 0, _l2_pt = 0, _l2_eta = 0, _l2_phi = 0, _l2_trk_pt = 0, _l2_trk_eta = 0, _l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  //int _l1_charge = 0, _l2_charge = 0;
  int ntrk = 0, z_ntrk = 0, z_truth_ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  float z_trk_pt[10000], z_trk_eta[10000], z_trk_phi[10000];
  float z_truth_trk_pt[10000], z_truth_trk_eta[10000], z_truth_trk_phi[10000];

  int nMixed = 0, nNotMixed = 0;
  //std::unordered_set <int> unMixableEvents;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is PbPb
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (isPbPb) {
    bkgTree->LoadBaskets (4000000000); //2000000000 = 2GB
    bkgTree->SetBranchAddress ("event_weight",     &event_weight);
    bkgTree->SetBranchAddress ("fcal_et",          &fcal_et);
    bkgTree->SetBranchAddress ("zdcEnergy",        &zdcEnergy);
    bkgTree->SetBranchAddress ("impactParameter",  &ip);
    bkgTree->SetBranchAddress ("eventPlane",       &eventPlane);
    bkgTree->SetBranchAddress ("q2",               &q2);
    bkgTree->SetBranchAddress ("psi2",             &psi2);
    bkgTree->SetBranchAddress ("q3",               &q3);
    bkgTree->SetBranchAddress ("psi3",             &psi3);
    bkgTree->SetBranchAddress ("q4",               &q4);
    bkgTree->SetBranchAddress ("psi4",             &psi4);
    bkgTree->SetBranchAddress ("vz",               &vz);
    bkgTree->SetBranchAddress ("ntrk_perp",        &ntrk_perp);
    bkgTree->SetBranchAddress ("ntrk",             &ntrk);

    zTree->LoadBaskets (2000000000);
    zTree->SetBranchAddress ("event_weight",    &z_event_weight);
    zTree->SetBranchAddress ("fcal_et",         &z_fcal_et);
    zTree->SetBranchAddress ("zdcEnergy",       &z_zdcEnergy);
    zTree->SetBranchAddress ("impactParameter", &z_ip);
    zTree->SetBranchAddress ("eventPlane",      &z_eventPlane);
    zTree->SetBranchAddress ("q2",              &z_q2);
    zTree->SetBranchAddress ("psi2",            &z_psi2);
    zTree->SetBranchAddress ("q3",              &z_q3);
    zTree->SetBranchAddress ("psi3",            &z_psi3);
    zTree->SetBranchAddress ("q4",              &z_q4);
    zTree->SetBranchAddress ("psi4",            &z_psi4);
    zTree->SetBranchAddress ("vz",              &z_vz);
    zTree->SetBranchAddress ("ntrk_perp",       &z_ntrk_perp);
    zTree->SetBranchAddress ("ntrk",            &z_ntrk);
    zTree->SetBranchAddress ("isEE",            &isEE);
    zTree->SetBranchAddress ("z_pt",            &z_pt);
    zTree->SetBranchAddress ("z_phi",           &z_phi);
    zTree->SetBranchAddress ("z_y",             &z_y);
    zTree->SetBranchAddress ("z_m",             &z_m);
    zTree->SetBranchAddress ("l1_pt",           &l1_pt);
    zTree->SetBranchAddress ("l1_eta",          &l1_eta);
    zTree->SetBranchAddress ("l1_phi",          &l1_phi);
    zTree->SetBranchAddress ("l1_trk_pt",       &l1_trk_pt);
    zTree->SetBranchAddress ("l1_trk_eta",      &l1_trk_eta);
    zTree->SetBranchAddress ("l1_trk_phi",      &l1_trk_phi);
    zTree->SetBranchAddress ("l1_charge",       &l1_charge);
    zTree->SetBranchAddress ("l2_pt",           &l2_pt);
    zTree->SetBranchAddress ("l2_eta",          &l2_eta);
    zTree->SetBranchAddress ("l2_phi",          &l2_phi);
    zTree->SetBranchAddress ("l2_trk_pt",       &l2_trk_pt);
    zTree->SetBranchAddress ("l2_trk_eta",      &l2_trk_eta);
    zTree->SetBranchAddress ("l2_trk_phi",      &l2_trk_phi);
    zTree->SetBranchAddress ("l2_charge",       &l2_charge);
    zTree->SetBranchAddress ("truth_ntrk",      &z_truth_ntrk);

    mixedEventsTree->Branch ("run_number",      &run_number,      "run_number/I");
    mixedEventsTree->Branch ("event_number",    &event_number,    "event_number/I");
    mixedEventsTree->Branch ("fcal_et",         &fcal_et,         "fcal_et/F");
    mixedEventsTree->Branch ("ip",              &ip,              "ip/F");
    mixedEventsTree->Branch ("eventPlane",      &eventPlane,      "eventPlane/F");
    mixedEventsTree->Branch ("zdcEnergy",       &zdcEnergy,       "zdcEnergy/F");
    mixedEventsTree->Branch ("q2",              &q2,              "q2/F");
    mixedEventsTree->Branch ("psi2",            &psi2,            "psi2/F");
    mixedEventsTree->Branch ("q3",              &q3,              "q3/F");
    mixedEventsTree->Branch ("psi3",            &psi3,            "psi3/F");
    mixedEventsTree->Branch ("q4",              &q4,              "q4/F");
    mixedEventsTree->Branch ("psi4",            &psi4,            "psi4/F");
    mixedEventsTree->Branch ("vz",              &vz,              "vz/F");
    mixedEventsTree->Branch ("ntrk_perp",       &ntrk_perp,       "ntrk_perp/I");
    mixedEventsTree->Branch ("ntrk",            &ntrk,            "ntrk/I");
    mixedEventsTree->Branch ("trk_pt",          &trk_pt,          "trk_pt[ntrk]/F");
    mixedEventsTree->Branch ("trk_eta",         &trk_eta,         "trk_eta[ntrk]/F");
    mixedEventsTree->Branch ("trk_phi",         &trk_phi,         "trk_phi[ntrk]/F");
    mixedEventsTree->Branch ("isEE",            &isEE,            "isEE/O");
    mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/I");
    mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/I");
    mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
    mixedEventsTree->Branch ("z_fcal_et",       &z_fcal_et,       "z_fcal_et/F");
    mixedEventsTree->Branch ("z_ip",            &z_ip,            "z_ip/F");
    mixedEventsTree->Branch ("z_eventPlane",    &z_eventPlane,    "z_eventPlane/F");
    mixedEventsTree->Branch ("z_zdcEnergy",     &z_zdcEnergy,     "z_zdcEnergy/F");
    mixedEventsTree->Branch ("z_q2",            &z_q2,            "z_q2/F");
    mixedEventsTree->Branch ("z_psi2",          &z_psi2,          "z_psi2/F");
    mixedEventsTree->Branch ("z_q3",            &z_q3,            "z_q3/F");
    mixedEventsTree->Branch ("z_psi3",          &z_psi3,          "z_psi3/F");
    mixedEventsTree->Branch ("z_q4",            &z_q4,            "z_q4/F");
    mixedEventsTree->Branch ("z_psi4",          &z_psi4,          "z_psi4/F");
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
    mixedEventsTree->Branch ("z_ntrk_perp",     &z_ntrk_perp,     "z_ntrk_perp/I");
    mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
    mixedEventsTree->Branch ("z_trk_pt",        &z_trk_pt,        "z_trk_pt[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_eta",       &z_trk_eta,       "z_trk_eta[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_phi",       &z_trk_phi,       "z_trk_phi[z_ntrk]/F");
    mixedEventsTree->Branch ("z_truth_ntrk",    &z_truth_ntrk,    "z_truth_ntrk/I");
    mixedEventsTree->Branch ("z_truth_trk_pt",  &z_truth_trk_pt,  "z_truth_trk_pt[z_truth_ntrk]/F");
    mixedEventsTree->Branch ("z_truth_trk_eta", &z_truth_trk_eta, "z_truth_trk_eta[z_truth_ntrk]/F");
    mixedEventsTree->Branch ("z_truth_trk_phi", &z_truth_trk_phi, "z_truth_trk_phi[z_truth_ntrk]/F");


    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (iZEvt % nZEvts);
      zEventsUsed[iZEvt % nZEvts]++;

      if (z_event_weight == 0) continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 2 || iPtZ > nPtZBins-1) continue;

      // Find the next unused minimum bias event
      {
        const short iFCalEt = GetSuperFineCentBin (z_fcal_et);
        if (!useImpactParameter && (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1)) continue;
        const short iIP = GetIPBin (z_ip);
        if (useImpactParameter && (iIP < 0 || iIP > nIPBins-1)) {
          cout << "Out-of-bounds impact parameter, skipping this Z!" << endl;
          continue;
        }
        const short iQ2 = GetQ2MixBin (z_q2);
        if (doQ2Mixing && (iQ2 < 0 || iQ2 > nQ2MixBins-1)) {
          cout << "Out-of-bounds q2, skipping this Z!" << endl;
          continue;
        }

        // disable big branches during event mixing
        //bkgTree->SetBranchStatus ("trk_pt", 0);
        //bkgTree->SetBranchStatus ("trk_eta", 0);
        //bkgTree->SetBranchStatus ("trk_phi", 0);
        
        bool goodMixEvent = false;
        int numMixTries = 0;
        do {
          iBkgEvt = (iBkgEvt+1) % nBkgEvts;
          bkgTree->GetEntry (bkgEventOrder[iBkgEvt]);
          goodMixEvent = (fabs (vz) < 150 && event_weight != 0); // always require these conditions
          //goodMixEvent &= (!doSameFileMixing || mixingFraction != 1 || !bkgEventsUsed[iBkgEvt]); // checks for uniqueness (if applicable)
          //goodMixEvent &= (!doSameFileMixing || !bkgEventsUsed[iBkgEvt]); // checks for uniqueness (if applicable)
          goodMixEvent &= (!doSameFileMixing || iBkgEvt != iZEvt); // don't mix with the exact same event
          goodMixEvent &= (!doCentMixing || (!useImpactParameter && iFCalEt == GetSuperFineCentBin (fcal_et)) || (useImpactParameter && iIP == GetIPBin (ip))); // do centrality matching
          //goodMixEvent &= (!doCentMixing || (!useImpactParameter && iFCalEt == GetSuperFineCentBin (fcal_et)) || (useImpactParameter && fabs (z_ntrk_perp - ntrk_perp) <= 3)); // do centrality matching
          //goodMixEvent &= (!doCentMixing || (!useImpactParameter && iFCalEt == GetSuperFineCentBin (fcal_et)) || (useImpactParameter && fabs (z_ip-ip) < 0.001)); // do centrality matching
          goodMixEvent &= (!doQ2Mixing   || iQ2 == GetQ2MixBin (q2)); // do q2 matching
          goodMixEvent &= (!doPsi2Mixing || (!useImpactParameter && DeltaPhi (psi2, z_psi2) < (pi / nPsi2MixBins)) || (useImpactParameter && DeltaPhi (eventPlane, z_eventPlane) < (pi / nPsi2MixBins))); // do psi2 matching
          //goodMixEvent &= (!doPsi2Mixing || DeltaPhi (psi2, z_psi2) < (pi / nPsi2MixBins)); // do psi2 matching
          goodMixEvent &= (!doPsi3Mixing || DeltaPhi (psi3, z_psi3) < (pi / nPsi3MixBins)); // do psi3 matching
          goodMixEvent &= (!doVZMixing || fabs (vz - z_vz) < 15);
          numMixTries++;
        } while (!goodMixEvent && numMixTries <= nBkgEvts); // only check each event once
        if (numMixTries == nBkgEvts+1 || std::find (mixingMap[iZEvt % nZEvts].begin (), mixingMap[iZEvt % nZEvts].end (), iBkgEvt) != mixingMap[iZEvt % nZEvts].end ()) {
          //cout << "No minbias event to mix with!!! Wrapped around on the same Z!!! Sum Et = " << z_fcal_et << ", q2 = " << z_q2 << ", psi2 = " << z_psi2 << endl;
          nNotMixed++;
          continue;
        }
        nMixed++;

        // reenable big branches after finding a suitable event
        //bkgTree->SetBranchStatus ("trk_pt", 1);
        //bkgTree->SetBranchStatus ("trk_eta", 1);
        //bkgTree->SetBranchStatus ("trk_phi", 1);
        //bkgTree->GetEntry (bkgEventOrder[iBkgEvt]);
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      mixingMap[iZEvt % nZEvts].push_back (bkgEventOrder[iBkgEvt]);
    } // end event matching loop


    bkgTree->SetBranchAddress ("run_number",    &run_number);
    bkgTree->SetBranchAddress ("event_number",  &event_number);
    bkgTree->SetBranchAddress ("trk_pt",        &trk_pt);
    bkgTree->SetBranchAddress ("trk_eta",       &trk_eta);
    bkgTree->SetBranchAddress ("trk_phi",       &trk_phi);
    zTree->SetBranchAddress ("run_number",      &z_run_number);
    zTree->SetBranchAddress ("event_number",    &z_event_number);
    zTree->SetBranchAddress ("trk_pt",          &z_trk_pt);
    zTree->SetBranchAddress ("trk_eta",         &z_trk_eta);
    zTree->SetBranchAddress ("trk_phi",         &z_trk_phi);
    zTree->SetBranchAddress ("truth_trk_pt",    &z_truth_trk_pt);
    zTree->SetBranchAddress ("truth_trk_eta",   &z_truth_trk_eta);
    zTree->SetBranchAddress ("truth_trk_phi",   &z_truth_trk_phi);

    for (int iZEvt = 0; iZEvt < nZEvts; iZEvt++) {
      if (nZEvts > 100 && iZEvt % (nZEvts / 100) == 0)
        cout << iZEvt / (nZEvts / 100) << "\% done...\r" << flush;

      if (mixingMap[iZEvt].size () < mixingFraction) continue; // then we couldn't find enough candidate events

      zTree->GetEntry (iZEvt);

      if (z_event_weight == 0) continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 2 || iPtZ > nPtZBins-1) continue;

      for (int iBkgEvt : mixingMap[iZEvt]) {
        bkgTree->GetEntry (iBkgEvt);
        mixedEventsTree->Fill ();
      }
    }

    cout << endl << "Done Pb+Pb loop." << endl;
  } // end tree-filling loop


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is pp
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (!isPbPb) {
    bkgTree->LoadBaskets (3000000000);
    bkgTree->SetBranchAddress ("run_number",   &run_number);
    bkgTree->SetBranchAddress ("event_number", &event_number);
    bkgTree->SetBranchAddress ("event_weight", &event_weight);
    bkgTree->SetBranchAddress ("fcal_et",      &fcal_et);
    bkgTree->SetBranchAddress ("vz",           &vz);
    bkgTree->SetBranchAddress ("ntrk",         &ntrk);
    bkgTree->SetBranchAddress ("trk_pt",       &trk_pt);
    bkgTree->SetBranchAddress ("trk_eta",      &trk_eta);
    bkgTree->SetBranchAddress ("trk_phi",      &trk_phi);

    zTree->LoadBaskets (3000000000);
    zTree->SetBranchAddress ("run_number",    &z_run_number);
    zTree->SetBranchAddress ("event_number",  &z_event_number);
    zTree->SetBranchAddress ("event_weight",  &z_event_weight);
    zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
    zTree->SetBranchAddress ("vz",            &z_vz);
    zTree->SetBranchAddress ("ntrk",          &z_ntrk);
    zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
    zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
    zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    zTree->SetBranchAddress ("phi_transmin",  &phi_transmin);
    zTree->SetBranchAddress ("phi_transmax",  &phi_transmax);
    zTree->SetBranchAddress ("isEE",          &_isEE);
    zTree->SetBranchAddress ("z_pt",          &z_pt);
    zTree->SetBranchAddress ("z_phi",         &z_phi);
    zTree->SetBranchAddress ("z_y",           &z_y);
    zTree->SetBranchAddress ("z_m",           &z_m);

    mixedEventsTree->Branch ("run_number",      &run_number,      "run_number/i");
    mixedEventsTree->Branch ("event_number",    &event_number,    "event_number/i");
    mixedEventsTree->Branch ("isEE",            &isEE,            "isEE/O");
    mixedEventsTree->Branch ("fcal_et",         &fcal_et,         "fcal_et/F");
    mixedEventsTree->Branch ("vz",              &vz,              "vz/F");
    mixedEventsTree->Branch ("ntrk",            &ntrk,            "ntrk/I");
    mixedEventsTree->Branch ("trk_pt",          &trk_pt,          "trk_pt[ntrk]/F");
    mixedEventsTree->Branch ("trk_eta",         &trk_eta,         "trk_eta[ntrk]/F");
    mixedEventsTree->Branch ("trk_phi",         &trk_phi,         "trk_phi[ntrk]/F");
    mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/i");
    mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/i");
    mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
    mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
    mixedEventsTree->Branch ("z_pt",            &_z_pt,            "z_pt/F");
    mixedEventsTree->Branch ("z_phi",           &_z_phi,           "z_phi/F");
    mixedEventsTree->Branch ("z_y",             &_z_y,             "z_y/F");
    mixedEventsTree->Branch ("z_m",             &_z_m,             "z_m/F");
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
    mixedEventsTree->Branch ("z_trk_pt",        &z_trk_pt,        "z_trk_pt[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_eta",       &z_trk_eta,       "z_trk_eta[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_phi",       &z_trk_phi,       "z_trk_phi[z_ntrk]/F");


    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (iZEvt % nZEvts);
      zEventsUsed[iZEvt % nZEvts]++;
      z_m = _z_m;
      isEE = _isEE;
      z_pt = _z_pt;
      z_phi = _z_phi;
      z_y = _z_y;
      z_m = _z_m;
      //l1_pt = _l1_pt;
      //l1_eta = _l1_eta;
      //l1_phi = _l1_phi;
      //l1_trk_pt = _l1_trk_pt;
      //l1_trk_eta = _l1_trk_eta;
      //l1_trk_phi = _l1_trk_phi;
      //l1_charge = _l1_charge;
      //l2_pt = _l2_pt;
      //l2_eta = _l2_eta;
      //l2_phi = _l2_phi;
      //l2_trk_pt = _l2_trk_pt;
      //l2_trk_eta = _l2_trk_eta;
      //l2_trk_phi = _l2_trk_phi;
      //l2_charge = _l2_charge;
      if (doSameFileMixing) {
        z_run_number = run_number;
        z_event_number = event_number;
        z_event_weight = event_weight;
        z_fcal_et = fcal_et;
        //z_zdcEnergy = zdcEnergy;
        z_vz = vz;
        z_ntrk = 0; // placeholder
      }

      if (z_event_weight == 0) continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 2 || iPtZ > nPtZBins-1) continue;

      // Find the next unused minimum bias event
      {
        // disable big branches during event mixing
        bkgTree->SetBranchStatus ("trk_pt", 0);
        bkgTree->SetBranchStatus ("trk_eta", 0);
        bkgTree->SetBranchStatus ("trk_phi", 0);

        bool goodMixEvent = false;
        const int _iBkgEvt = iBkgEvt;
        do {
          iBkgEvt = (iBkgEvt+1) % nBkgEvts;
          bkgTree->GetEntry (bkgEventOrder[iBkgEvt]);
          goodMixEvent = (fabs (vz) < 150 && event_weight != 0); // always require these conditions
          goodMixEvent &= (doSameFileMixing || mixingFraction != 1 || !bkgEventsUsed[iBkgEvt]); // checks for uniqueness (if applicable)
          goodMixEvent &= (!doSameFileMixing || iBkgEvt != iZEvt); // don't mix with the exact same event
          goodMixEvent &= (!doPPTransMinMixing || (1 < _z_pt && _z_pt < 12 && DeltaPhi (phi_transmin, z_phi) >= 7.*pi/8.));
          goodMixEvent &= (!doPPTransMaxMixing || (1 < _z_pt && _z_pt < 12 && DeltaPhi (phi_transmax, z_phi) >= 7.*pi/8.)); // alternatively mix with trans-max region
        } while (!goodMixEvent && iBkgEvt != _iBkgEvt); // only check each event once
        if (_iBkgEvt == iBkgEvt) {
          //cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          nNotMixed++;
          continue;
        }
        nMixed++;
        bkgEventsUsed[iBkgEvt] = true;

        // reenable big branches after finding a suitable event
        bkgTree->SetBranchStatus ("trk_pt", 1);
        bkgTree->SetBranchStatus ("trk_eta", 1);
        bkgTree->SetBranchStatus ("trk_phi", 1);
        bkgTree->GetEntry (bkgEventOrder[iBkgEvt]);
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      mixedEventsTree->Fill ();

    } // end loop over pp tree

    cout << endl << "Done pp loop." << endl;
  }

  cout << "Number of events mixed:     " << nMixed << endl;
  cout << "Number of events not mixed: " << nNotMixed << " (" << ((float)nNotMixed) * 100 / ((float)nNotMixed + (float)nMixed) << "\%)" << endl;

  for (int i = 0; i < nZEvts; i++) {
    if (zEventsUsed[i] < mixingFraction) {
      cout << "Warning! Some events used fewer than " << mixingFraction << " times!" << endl;
      break;
    }
  }

  mixedEventsTree->SetDirectory (mixedEventsFile);
  mixedEventsTree->Write ("", TObject :: kOverwrite);

  if (mixedEventsFile) mixedEventsFile->Close ();
  SaferDelete (&mixedEventsFile);

  if (inFile) inFile->Close ();
  SaferDelete (&inFile);

  if (bkgInFile) bkgInFile->Close ();
  SaferDelete (&bkgInFile);
}





int main (int argc, char** argv) {

  if (argc < 5) {
    cout << "Insufficient arguments!" << endl;
    return 0;
  }

  DoMixing (string (argv[1]) == "true", argv[2], argv[3], argv[4]);

  return 0;
}

#endif
