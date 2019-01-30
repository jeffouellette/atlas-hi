#include "TreeVariables.h"

#include <iostream>
#include <iomanip>

namespace JetTrackAnalysis {

TreeVariables :: TreeVariables (TTree* t, const bool _isMC) {
  tree = t;
  isMC = _isMC;
}


TreeVariables :: ~TreeVariables () {
  if (vert_z) delete vert_z;
  if (vert_ntrk) delete vert_ntrk;
  if (vert_type) delete vert_type;

  if (trk_quality_4) delete trk_quality_4;
  if (trk_d0) delete trk_d0;
  if (trk_z0) delete trk_z0;
  if (trk_theta) delete trk_theta;
  if (trk_charge) delete trk_charge;
  if (trk_pt) delete trk_pt;
  if (trk_eta) delete trk_eta;
  if (trk_phi) delete trk_phi;

  if (jet_pt) delete jet_pt;
  if (jet_eta) delete jet_eta;
  if (jet_phi) delete jet_phi;
  if (jet_e) delete jet_e;

  if (truth_jet_pt) delete truth_jet_pt;
  if (truth_jet_eta) delete truth_jet_eta;
  if (truth_jet_phi) delete truth_jet_phi;
  if (truth_jet_e) delete truth_jet_e;
}


void TreeVariables :: SetBranchAddresses () {

  if (getCollisionRateInfo) {
    tree->SetBranchAddress ("actualInteractionsPerCrossing", &actualInteractionsPerCrossing);
    tree->SetBranchAddress ("averageInteractionsPerCrossing", &averageInteractionsPerCrossing);
  }
  else {
    tree->SetBranchStatus ("actualInteractionsPerCrossing", 0);
    tree->SetBranchStatus ("averageInteractionsPerCrossing", 0);
  }

  if (isMC && getMCInfo) {
    tree->SetBranchAddress ("NumberEvents", &numberEvents);
    tree->SetBranchAddress ("CrossSection_microbarns", &crossSection_microbarns);
    tree->SetBranchAddress ("FilterEfficiency", &filterEfficiency);
  }
  if (!isMC) {
    tree->SetBranchAddress ("eventNumber", &eventNumber);
    tree->SetBranchAddress ("runNumber", &runNumber);
    tree->SetBranchAddress ("lumiBlock", &lumiBlock);
  }

  if (getVertices) {
    tree->SetBranchAddress ("nvert", &nvert);
    //tree->SetBranchAddress ("vert_x", &vert_x);
    tree->SetBranchStatus ("vert_x", 0);
    //tree->SetBranchAddress ("vert_y", &vert_y);
    tree->SetBranchStatus ("vert_y", 0);
    tree->SetBranchAddress ("vert_z", &vert_z);
    //tree->SetBranchAddress ("vert_ntrk", &vert_ntrk);
    tree->SetBranchStatus ("vert_ntrk", 0);
    tree->SetBranchAddress ("vert_type", &vert_type);
  }
  else {
    tree->SetBranchStatus ("nvert", 0);
    tree->SetBranchStatus ("vert_x", 0);
    tree->SetBranchStatus ("vert_y", 0);
    tree->SetBranchStatus ("vert_z", 0);
    tree->SetBranchStatus ("vert_ntrk", 0);
    tree->SetBranchStatus ("vert_type", 0);
  }

  if (getFCals) {
    tree->SetBranchAddress ("fcalA_et", &fcalA_et);
    tree->SetBranchAddress ("fcalC_et", &fcalC_et);
  }
  else {
    tree->SetBranchStatus ("fcalA_et", 0);
    tree->SetBranchStatus ("fcalC_et", 0);
  }

  if (getTracks) {
    tree->SetBranchAddress ("ntrk", &ntrk);
    tree->SetBranchAddress ("trk_quality_4", &trk_quality_4);
    //tree->SetBranchAddress ("trk_d0", &trk_d0);
    tree->SetBranchStatus ("trk_d0", 0);
    //tree->SetBranchAddress ("trk_z0", &trk_z0);
    tree->SetBranchStatus ("trk_z0", 0);
    //tree->SetBranchAddress ("trk_theta", &trk_theta);
    tree->SetBranchStatus ("trk_theta", 0);
    tree->SetBranchAddress ("trk_charge", &trk_charge);
    tree->SetBranchAddress ("trk_pt", &trk_pt);
    tree->SetBranchAddress ("trk_eta", &trk_eta);
    tree->SetBranchAddress ("trk_phi", &trk_phi);
  }
  else {
    tree->SetBranchStatus ("ntrk", 0);
    tree->SetBranchStatus ("trk_quality_4", 0);
    tree->SetBranchStatus ("trk_d0", 0);
    tree->SetBranchStatus ("trk_z0", 0);
    tree->SetBranchStatus ("trk_theta", 0);
    tree->SetBranchStatus ("trk_charge", 0);
    tree->SetBranchStatus ("trk_pt", 0);
    tree->SetBranchStatus ("trk_eta", 0);
    tree->SetBranchStatus ("trk_phi", 0);
  }

  if (getJets) {
    tree->SetBranchAddress ("akt4hi_jet_n", &jet_n);
  }
  else {
    tree->SetBranchStatus ("akt4hi_jet_n", 0);

  }

  if (getJets && jetType == FullCalibrated) {
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_pt", &jet_pt);
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_eta", &jet_eta);
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_phi", &jet_phi);
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_e", &jet_e);
  }
  else {
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_e", 0);
  }

  if (getJets && jetType == EtaJES) {
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_pt", &jet_pt);
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_eta", &jet_eta);
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_phi", &jet_phi);
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_e", &jet_e);
  }
  else {
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_e", 0);
  }

  if (getJets && jetType == EM) {
    tree->SetBranchAddress ("akt4hi_em_jet_pt", &jet_pt);
    tree->SetBranchAddress ("akt4hi_em_jet_eta", &jet_eta);
    tree->SetBranchAddress ("akt4hi_em_jet_phi", &jet_phi);
    tree->SetBranchAddress ("akt4hi_em_jet_e", &jet_e);
  }
  else {
    tree->SetBranchStatus ("akt4hi_em_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_em_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_em_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_em_jet_e", 0);
  }
  

  if (isMC) {

    if (getTruthJets) {
      tree->SetBranchAddress ("akt4_truth_jet_n", &truth_jet_n);
      tree->SetBranchAddress ("akt4_truth_jet_pt", &truth_jet_pt);
      tree->SetBranchAddress ("akt4_truth_jet_eta", &truth_jet_eta);
      tree->SetBranchAddress ("akt4_truth_jet_phi", &truth_jet_phi);
      tree->SetBranchAddress ("akt4_truth_jet_e", &truth_jet_e);
    }
    else {
      tree->SetBranchStatus ("akt4_truth_jet_n", 0);
      tree->SetBranchStatus ("akt4_truth_jet_pt", 0);
      tree->SetBranchStatus ("akt4_truth_jet_eta", 0);
      tree->SetBranchStatus ("akt4_truth_jet_phi", 0);
      tree->SetBranchStatus ("akt4_truth_jet_e", 0);
    }

  }

  return;
}

void TreeVariables :: PrintAll (const long long entry) {
  tree->GetEntry (entry);
  cout << endl << "////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "// Getting event " << entry << "..." << endl;
  cout << "////////////////////////////////////////////////////////////////////////////////" << endl;

  if (getJets) {

    if (jetType == FullCalibrated)
      cout << endl << "Anti-Kt R=0.4 HI Insitu/xCalib + EtaJES Calibratd Jets:" << endl;
    else if (jetType == EtaJES)
      cout << endl << "Anti-Kt R=0.4 HI EtaJES Calibrated Jets:" << endl;
    else if (jetType == EM)
      cout << endl << "Anti-Kt R=0.4 HI Uncalibrated (EM scale) Jets:" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << jet_pt->at (j)
           << setw (10) << jet_eta->at (j)
           << setw (10) << jet_phi->at (j)
           << setw (10) << jet_e->at (j) << endl;
    }
  }

  if (isMC) {
    cout << endl << "Anti-Kt R=0.4 Truth Jets:" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < truth_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << truth_jet_pt->at (j)
           << setw (10) << truth_jet_eta->at (j)
           << setw (10) << truth_jet_phi->at (j)
           << setw (10) << truth_jet_e->at (j) << endl;
    }
  }

  if (getTracks) {
    cout << endl << "Tracks:" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Quality"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi" << endl;

    for (int t = 0; t < ntrk; t++) {
      cout << setw (10) << t
           << setw (10) << trk_quality_4->at (t)
           << setw (10) << trk_pt->at (t)
           << setw (10) << trk_eta->at (t)
           << setw (10) << trk_phi->at (t) << endl;
    }
  }

  return;
}

/** Setter functions **/

void TreeVariables :: SetGetMCInfo (const bool _getMCInfo, const double _crossSection_microbarns, const double _filterEfficiency, const int _numberEvents) {
  getMCInfo = _getMCInfo;
  crossSection_microbarns = _crossSection_microbarns;
  filterEfficiency = _filterEfficiency;
  numberEvents = _numberEvents;
}

void TreeVariables :: SetGetCollisionRateInfo (const bool _getCollisionRateInfo)  { getCollisionRateInfo = _getCollisionRateInfo; }
void TreeVariables :: SetGetVertices (const bool _getVertices)                    { getVertices = _getVertices; }
void TreeVariables :: SetGetFCals (const bool _getFCals)                          { getFCals = _getFCals; }
void TreeVariables :: SetGetTracks (const bool _getTracks)                        { getTracks = _getTracks; }
void TreeVariables :: SetGetJets (const bool _getJets, const JetType _jetType)    { getJets = _getJets; jetType = _jetType; }
void TreeVariables :: SetGetTruthJets (const bool _getTruthJets)                  { getTruthJets = _getTruthJets; }

} // end namespace
