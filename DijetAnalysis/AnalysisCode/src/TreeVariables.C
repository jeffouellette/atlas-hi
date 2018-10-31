#include "TreeVariables.h"
#include <set>

namespace JetAnalysis {

TreeVariables::TreeVariables(TTree* tree, const bool _isMC) {
  t = tree;
  isMC = _isMC;
}


TreeVariables::~TreeVariables() {
  nvert = 0;
  if (vert_x) delete vert_x;
  if (vert_y) delete vert_y;
  if (vert_z) delete vert_z;
  if (vert_ntrk) delete vert_ntrk;
  if (vert_type) delete vert_type;

  ntrk = 0;
  if (trk_quality_4) delete trk_quality_4;
  if (trk_d0) delete trk_d0;
  if (trk_z0) delete trk_z0;
  if (trk_theta) delete trk_theta;
  if (trk_charge) delete trk_charge;
  if (trk_pt) delete trk_pt;
  if (trk_eta) delete trk_eta;
  if (trk_phi) delete trk_phi;

  akt4hi_njet = 0;
  if (akt4hi_em_jet_pt) delete akt4hi_em_jet_pt;
  if (akt4hi_em_jet_eta) delete akt4hi_em_jet_eta;
  if (akt4hi_em_jet_phi) delete akt4hi_em_jet_phi;
  if (akt4hi_em_jet_e) delete akt4hi_em_jet_e;
  if (akt4hi_etajes_jet_pt) delete akt4hi_etajes_jet_pt;
  if (akt4hi_etajes_jet_eta) delete akt4hi_etajes_jet_eta;
  if (akt4hi_etajes_jet_phi) delete akt4hi_etajes_jet_phi;
  if (akt4hi_etajes_jet_e) delete akt4hi_etajes_jet_e;
  if (akt4hi_xcalib_jet_pt) delete akt4hi_xcalib_jet_pt;
  if (akt4hi_xcalib_jet_eta) delete akt4hi_xcalib_jet_eta;
  if (akt4hi_xcalib_jet_phi) delete akt4hi_xcalib_jet_phi;
  if (akt4hi_xcalib_jet_e) delete akt4hi_xcalib_jet_e;

  akt4truth_njet = 0;
  if (akt4truth_jet_pt) delete akt4truth_jet_pt;
  if (akt4truth_jet_eta) delete akt4truth_jet_eta;
  if (akt4truth_jet_phi) delete akt4truth_jet_phi;
  if (akt4truth_jet_e) delete akt4truth_jet_e;

  hlt_njet = 0;
  if (hlt_jet_pt) delete hlt_jet_pt;
  if (hlt_jet_eta) delete hlt_jet_eta;
  if (hlt_jet_phi) delete hlt_jet_phi;
  if (hlt_jet_e) delete hlt_jet_e;
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, int* ptr, const bool get) {
  if (active.count(TString(name)) == 1 && get) t->SetBranchAddress(name, ptr);
  else t->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, float* ptr, const bool get) {
  if (active.count(TString(name)) == 1 && get) t->SetBranchAddress(name, ptr);
  else t->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, double* ptr, const bool get) {
  if (active.count(TString(name)) == 1 && get) t->SetBranchAddress(name, ptr);
  else t->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, vector<int>** ptr, const bool get) {
  if (active.count(TString(name)) == 1 && get) t->SetBranchAddress(name, ptr);
  else t->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, vector<float>** ptr, const bool get) {
  if (active.count(TString(name)) == 1 && get) t->SetBranchAddress(name, ptr);
  else t->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, vector<bool>** ptr, const bool get) {
  if (active.count(TString(name)) == 1 && get) t->SetBranchAddress(name, ptr);
  else t->SetBranchStatus(name, 0);
}


void TreeVariables::SetBranchAddresses (const std::set<TString>& active) {
  BranchHelper (active, "runNumber", &runNumber);
  BranchHelper (active, "eventNumber", &eventNumber);
  BranchHelper (active, "lumiBlock", &lumiBlock);
  //if (isMC) BranchHelper (active, "EventWeight", &eventWeight);

  BranchHelper (active, "actualInteractionsPerCrossing", &actualInteractionsPerCrossing, getCollisionRateInfo);
  BranchHelper (active, "averageInteractionsPerCrossing", &averageInteractionsPerCrossing, getCollisionRateInfo);

  BranchHelper (active, "nvert", &nvert, getVertices);
  //BranchHelper (active, "vert_x", &vert_x, getVertices);
  //BranchHelper (active, "vert_y", &vert_y, getVertices);
  BranchHelper (active, "vert_z", &vert_z, getVertices);
  //BranchHelper (active, "vert_ntrk", &vert_ntrk, getVertices);
  BranchHelper (active, "vert_type", &vert_type, getVertices);

  //BranchHelper (active, "trk_quality_4", &trk_quality_4, getTracks);
  //BranchHelper (active, "trk_d0", &trk_d0, getTracks);
  //BranchHelper (active, "trk_z0", &trk_z0, getTracks);
  //BranchHelper (active, "trk_theta", &trk_theta, getTracks);
  //BranchHelper (active, "trk_charge", &trk_charge, getTracks);
  //BranchHelper (active, "trk_pt", &trk_pt, getTracks);
  //BranchHelper (active, "trk_eta", &trk_eta, getTracks);
  //BranchHelper (active, "trk_phi", &trk_phi, getTracks);

  BranchHelper (active, "fcalA_et", &fcalA_et, getFCals);
  BranchHelper (active, "fcalC_et", &fcalC_et, getFCals);

  BranchHelper (active, "total_jets", &total_jets);
  BranchHelper (active, "clean_jets", &clean_jets);

  if (getSimpleJets) {
    BranchHelper (active, "akt4hi_njet", &njet, getHIJets);
  }
  else {
    BranchHelper (active, "akt4hi_njet", &akt4hi_njet, getHIJets);
  }

  const bool getIntermediateHIJets = getHIJets && !getSimpleJets;
  BranchHelper (active, "akt4hi_em_jet_pt", &akt4hi_em_jet_pt, getIntermediateHIJets);
  BranchHelper (active, "akt4hi_em_jet_eta", &akt4hi_em_jet_eta, getIntermediateHIJets);
  BranchHelper (active, "akt4hi_em_jet_phi", &akt4hi_em_jet_phi, getIntermediateHIJets);
  BranchHelper (active, "akt4hi_em_jet_e", &akt4hi_em_jet_e, getIntermediateHIJets);

  BranchHelper (active, "akt4hi_etajes_jet_pt", &akt4hi_etajes_jet_pt, getIntermediateHIJets);
  BranchHelper (active, "akt4hi_etajes_jet_eta", &akt4hi_etajes_jet_eta, getIntermediateHIJets);
  BranchHelper (active, "akt4hi_etajes_jet_phi", &akt4hi_etajes_jet_phi, getIntermediateHIJets);
  BranchHelper (active, "akt4hi_etajes_jet_e", &akt4hi_etajes_jet_e, getIntermediateHIJets);

  if (getSimpleJets) {
    BranchHelper (active, "akt4hi_xcalib_jet_pt", &akt4hi_xcalib_jet_pt, getHIJets);
    BranchHelper (active, "akt4hi_xcalib_jet_eta", &akt4hi_xcalib_jet_eta, getHIJets);
    BranchHelper (active, "akt4hi_xcalib_jet_phi", &akt4hi_xcalib_jet_phi, getHIJets);
    BranchHelper (active, "akt4hi_xcalib_jet_e", &akt4hi_xcalib_jet_e, getHIJets);
  }
  else {
    BranchHelper (active, "akt4hi_xcalib_jet_pt", &jet_pt, getHIJets);
    BranchHelper (active, "akt4hi_xcalib_jet_eta", &jet_eta, getHIJets);
    BranchHelper (active, "akt4hi_xcalib_jet_phi", &jet_phi, getHIJets);
    BranchHelper (active, "akt4hi_xcalib_jet_e", &jet_e, getHIJets);
  }

  if (isMC) {
   BranchHelper (active, "akt4truth_njet", &akt4truth_njet);
   BranchHelper (active, "akt4truth_jet_pt", &akt4truth_jet_pt);
   BranchHelper (active, "akt4truth_jet_eta", &akt4truth_jet_eta);
   BranchHelper (active, "akt4truth_jet_phi", &akt4truth_jet_phi);
   BranchHelper (active, "akt4truth_jet_e", &akt4truth_jet_e);
  }
  if (!isMC) {
   BranchHelper (active, "hlt_njet", &hlt_njet);
   BranchHelper (active, "hlt_jet_pt", &hlt_jet_pt);
   BranchHelper (active, "hlt_jet_eta", &hlt_jet_eta);
   BranchHelper (active, "hlt_jet_phi", &hlt_jet_phi);
   BranchHelper (active, "hlt_jet_e", &hlt_jet_e);
  }

  return;
}


void TreeVariables :: SetGetCollisionRateInfo (const bool _getCollisionRateInfo) {
  getCollisionRateInfo = _getCollisionRateInfo;
}

void TreeVariables :: SetGetVertices (const bool _getVertices) {
  getVertices = _getVertices;
}

void TreeVariables :: SetGetFCals (const bool _getFCals) {
  getFCals = _getFCals;
}

void TreeVariables :: SetGetTracks (const bool _getTracks) {
  getTracks = _getTracks;
}

void TreeVariables :: SetGetSimpleJets (const bool _getSimpleJets) {
  getSimpleJets = _getSimpleJets;
}

void TreeVariables :: SetGetHIJets (const bool _getHIJets) {
  getHIJets = _getHIJets;
}

void TreeVariables :: SetGetEMTopoJets (const bool _getEMTopoJets) {
  getEMTopoJets = _getEMTopoJets;
}

} // end namespace
