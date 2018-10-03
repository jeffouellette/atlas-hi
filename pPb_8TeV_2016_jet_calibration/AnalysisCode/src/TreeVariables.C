#include "TreeVariables.h"

#include <iostream>
#include <iomanip>

namespace pPb8TeV2016JetCalibration {

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

  if (akt4hi_em_jet_pt) delete akt4hi_em_jet_pt;
  if (akt4hi_em_jet_eta) delete akt4hi_em_jet_eta;
  if (akt4hi_em_jet_phi) delete akt4hi_em_jet_phi;
  if (akt4hi_em_jet_e) delete akt4hi_em_jet_e;

  if (akt4hi_em_etajes_jet_pt) delete akt4hi_em_etajes_jet_pt;
  if (akt4hi_em_etajes_jet_eta) delete akt4hi_em_etajes_jet_eta;
  if (akt4hi_em_etajes_jet_phi) delete akt4hi_em_etajes_jet_phi;
  if (akt4hi_em_etajes_jet_e) delete akt4hi_em_etajes_jet_e;

  if (akt4hi_em_xcalib_jet_pt) delete akt4hi_em_xcalib_jet_pt;
  if (akt4hi_em_xcalib_jet_eta) delete akt4hi_em_xcalib_jet_eta;
  if (akt4hi_em_xcalib_jet_phi) delete akt4hi_em_xcalib_jet_phi;
  if (akt4hi_em_xcalib_jet_e) delete akt4hi_em_xcalib_jet_e;

  if (akt4hi_constit_xcalib_jet_pt) delete akt4hi_constit_xcalib_jet_pt;
  if (akt4hi_constit_xcalib_jet_eta) delete akt4hi_constit_xcalib_jet_eta;
  if (akt4hi_constit_xcalib_jet_phi) delete akt4hi_constit_xcalib_jet_phi;
  if (akt4hi_constit_xcalib_jet_e) delete akt4hi_constit_xcalib_jet_e;

  if (akt4hi_constit_etajes_jet_pt) delete akt4hi_constit_etajes_jet_pt;
  if (akt4hi_constit_etajes_jet_eta) delete akt4hi_constit_etajes_jet_eta;
  if (akt4hi_constit_etajes_jet_phi) delete akt4hi_constit_etajes_jet_phi;
  if (akt4hi_constit_etajes_jet_e) delete akt4hi_constit_etajes_jet_e;

  if (akt4hi_constit_jet_pt) delete akt4hi_constit_jet_pt;
  if (akt4hi_constit_jet_eta) delete akt4hi_constit_jet_eta;
  if (akt4hi_constit_jet_phi) delete akt4hi_constit_jet_phi;
  if (akt4hi_constit_jet_e) delete akt4hi_constit_jet_e;

  if (akt4emtopo_calib_jet_pt) delete akt4emtopo_calib_jet_pt;
  if (akt4emtopo_calib_jet_eta) delete akt4emtopo_calib_jet_eta;
  if (akt4emtopo_calib_jet_phi) delete akt4emtopo_calib_jet_phi;
  if (akt4emtopo_calib_jet_e) delete akt4emtopo_calib_jet_e;

  if (akt4emtopo_em_jet_pt) delete akt4emtopo_em_jet_pt;
  if (akt4emtopo_em_jet_eta) delete akt4emtopo_em_jet_eta;
  if (akt4emtopo_em_jet_phi) delete akt4emtopo_em_jet_phi;
  if (akt4emtopo_em_jet_e) delete akt4emtopo_em_jet_e;

  if (truth_jet_pt) delete truth_jet_pt;
  if (truth_jet_eta) delete truth_jet_eta;
  if (truth_jet_phi) delete truth_jet_phi;
  if (truth_jet_e) delete truth_jet_e;

  if (muon_pt) delete muon_pt;
  if (muon_eta) delete muon_eta;
  if (muon_phi) delete muon_phi;
  if (muon_quality) delete muon_quality;
  if (muon_charge) delete muon_charge;
  if (muon_tight) delete muon_tight;
  if (muon_loose) delete muon_loose;

  if (truth_muon_pt) delete truth_muon_pt;
  if (truth_muon_eta) delete truth_muon_eta;
  if (truth_muon_phi) delete truth_muon_phi;
  if (truth_muon_charge) delete truth_muon_charge;

  if (electron_pt) delete electron_pt;
  if (electron_eta) delete electron_eta;
  if (electron_phi) delete electron_phi;
  if (electron_charge) delete electron_charge;
  if (electron_loose) delete electron_loose;
  if (electron_tight) delete electron_tight;
  if (electron_d0sig) delete electron_d0sig;
  if (electron_delta_z0_sin_theta) delete electron_delta_z0_sin_theta;

  if (truth_electron_pt) delete truth_electron_pt;
  if (truth_electron_eta) delete truth_electron_eta;
  if (truth_electron_phi) delete truth_electron_phi;
  if (truth_electron_charge) delete truth_electron_charge;

  if (photon_pt) delete photon_pt;
  if (photon_eta) delete photon_eta;
  if (photon_phi) delete photon_phi;
  if (photon_tight) delete photon_tight;
  if (photon_isem) delete photon_isem;
  if (photon_convFlag) delete photon_convFlag;
  if (photon_Rconv) delete photon_Rconv;
  if (photon_topoetcone40) delete photon_topoetcone40;

  if (truth_photon_pt) delete truth_photon_pt;
  if (truth_photon_eta) delete truth_photon_eta;
  if (truth_photon_phi) delete truth_photon_phi;
}


void TreeVariables :: SetBranchAddresses () {
  tree->SetBranchStatus ("akt4hi_sampling", 0);
  tree->SetBranchStatus ("akt4emtopo_sampling", 0);

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
   tree->SetBranchAddress ("vert_x", &vert_x);
   tree->SetBranchAddress ("vert_y", &vert_y);
   tree->SetBranchAddress ("vert_z", &vert_z);
   tree->SetBranchAddress ("vert_ntrk", &vert_ntrk);
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
   //tree->SetBranchAddress ("trk_z0", &trk_z0);
   //tree->SetBranchAddress ("trk_theta", &trk_theta);
   //tree->SetBranchAddress ("trk_charge", &trk_charge);
   tree->SetBranchStatus ("trk_d0", 0);
   tree->SetBranchStatus ("trk_z0", 0);
   tree->SetBranchStatus ("trk_theta", 0);
   tree->SetBranchStatus ("trk_charge", 0);
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

  if (!getSimpleJets) {

   tree->SetBranchAddress ("total_jet_n", &total_jet_n);
   tree->SetBranchAddress ("clean_jet_n", &clean_jet_n);

   if (getHIJets) {

    tree->SetBranchAddress ("akt4hi_jet_n", &akt4hi_jet_n);

    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt);
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta);
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi);
    tree->SetBranchAddress ("akt4hi_em_xcalib_jet_e", &akt4hi_em_xcalib_jet_e);

    tree->SetBranchAddress ("akt4hi_em_etajes_jet_pt", &akt4hi_em_etajes_jet_pt);
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_eta", &akt4hi_em_etajes_jet_eta);
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_phi", &akt4hi_em_etajes_jet_phi);
    tree->SetBranchAddress ("akt4hi_em_etajes_jet_e", &akt4hi_em_etajes_jet_e);

    tree->SetBranchAddress ("akt4hi_em_jet_pt", &akt4hi_em_jet_pt);
    tree->SetBranchAddress ("akt4hi_em_jet_eta", &akt4hi_em_jet_eta);
    tree->SetBranchAddress ("akt4hi_em_jet_phi", &akt4hi_em_jet_phi);
    tree->SetBranchAddress ("akt4hi_em_jet_e", &akt4hi_em_jet_e);

    tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_pt", &akt4hi_constit_xcalib_jet_pt);
    tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_eta", &akt4hi_constit_xcalib_jet_eta);
    tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_phi", &akt4hi_constit_xcalib_jet_phi);
    tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_e", &akt4hi_constit_xcalib_jet_e);

    tree->SetBranchAddress ("akt4hi_constit_etajes_jet_pt", &akt4hi_constit_etajes_jet_pt);
    tree->SetBranchAddress ("akt4hi_constit_etajes_jet_eta", &akt4hi_constit_etajes_jet_eta);
    tree->SetBranchAddress ("akt4hi_constit_etajes_jet_phi", &akt4hi_constit_etajes_jet_phi);
    tree->SetBranchAddress ("akt4hi_constit_etajes_jet_e", &akt4hi_constit_etajes_jet_e);

    tree->SetBranchAddress ("akt4hi_constit_jet_pt", &akt4hi_constit_jet_pt);
    tree->SetBranchAddress ("akt4hi_constit_jet_eta", &akt4hi_constit_jet_eta);
    tree->SetBranchAddress ("akt4hi_constit_jet_phi", &akt4hi_constit_jet_phi);
    tree->SetBranchAddress ("akt4hi_constit_jet_e", &akt4hi_constit_jet_e);
   } // end if get HI jets

   else {
    tree->SetBranchStatus ("akt4hi_jet_n", 0);

    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_em_xcalib_jet_e", 0);

    tree->SetBranchStatus ("akt4hi_em_etajes_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_em_etajes_jet_e", 0);

    tree->SetBranchStatus ("akt4hi_em_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_em_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_em_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_em_jet_e", 0);

    tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_e", 0);

    tree->SetBranchStatus ("akt4hi_constit_etajes_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_constit_etajes_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_constit_etajes_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_constit_etajes_jet_e", 0);

    tree->SetBranchStatus ("akt4hi_constit_jet_pt", 0);
    tree->SetBranchStatus ("akt4hi_constit_jet_eta", 0);
    tree->SetBranchStatus ("akt4hi_constit_jet_phi", 0);
    tree->SetBranchStatus ("akt4hi_constit_jet_e", 0);
   }

   if (getEMTopoJets) {
    tree->SetBranchAddress ("akt4emtopo_jet_n", &akt4emtopo_jet_n);

    tree->SetBranchAddress ("akt4emtopo_em_jet_pt", &akt4emtopo_em_jet_pt);
    tree->SetBranchAddress ("akt4emtopo_em_jet_eta", &akt4emtopo_em_jet_eta);
    tree->SetBranchAddress ("akt4emtopo_em_jet_phi", &akt4emtopo_em_jet_phi);
    tree->SetBranchAddress ("akt4emtopo_em_jet_e", &akt4emtopo_em_jet_e);

    tree->SetBranchAddress ("akt4emtopo_calib_jet_pt", &akt4emtopo_calib_jet_pt);
    tree->SetBranchAddress ("akt4emtopo_calib_jet_eta", &akt4emtopo_calib_jet_eta);
    tree->SetBranchAddress ("akt4emtopo_calib_jet_phi", &akt4emtopo_calib_jet_phi);
    tree->SetBranchAddress ("akt4emtopo_calib_jet_e", &akt4emtopo_calib_jet_e);
   } // end if get EMTopo jets
   else {
    tree->SetBranchStatus ("akt4emtopo_jet_n", 0);

    tree->SetBranchStatus ("akt4emtopo_em_jet_pt", 0);
    tree->SetBranchStatus ("akt4emtopo_em_jet_eta", 0);
    tree->SetBranchStatus ("akt4emtopo_em_jet_phi", 0);
    tree->SetBranchStatus ("akt4emtopo_em_jet_e", 0);

    tree->SetBranchStatus ("akt4emtopo_calib_jet_pt", 0);
    tree->SetBranchStatus ("akt4emtopo_calib_jet_eta", 0);
    tree->SetBranchStatus ("akt4emtopo_calib_jet_phi", 0);
    tree->SetBranchStatus ("akt4emtopo_calib_jet_e", 0);
   }
  } // end not getSimpleJets
  else {

   tree->SetBranchStatus ("total_jet_n", 0);
   tree->SetBranchStatus ("clean_jet_n", 0);

   //// temporary solution to running over older data sets
   //tree->SetBranchAddress ("jet_n", &jet_n);

   tree->SetBranchAddress ("akt4hi_jet_n", &jet_n);

   tree->SetBranchAddress ("akt4hi_em_xcalib_jet_pt", &jet_pt);
   tree->SetBranchAddress ("akt4hi_em_xcalib_jet_eta", &jet_eta);
   tree->SetBranchAddress ("akt4hi_em_xcalib_jet_phi", &jet_phi);
   tree->SetBranchAddress ("akt4hi_em_xcalib_jet_e", &jet_e);

   //tree->SetBranchAddress ("akt4hi_em_etajes_jet_pt", &jet_pt);
   //tree->SetBranchAddress ("akt4hi_em_etajes_jet_eta", &jet_eta);
   //tree->SetBranchAddress ("akt4hi_em_etajes_jet_phi", &jet_phi);
   //tree->SetBranchAddress ("akt4hi_em_etajes_jet_e", &jet_e);

   //tree->SetBranchAddress ("akt4hi_em_jet_pt", &jet_pt);
   //tree->SetBranchAddress ("akt4hi_em_jet_eta", &jet_eta);
   //tree->SetBranchAddress ("akt4hi_em_jet_phi", &jet_phi);
   //tree->SetBranchAddress ("akt4hi_em_jet_e", &jet_e);

   //tree->SetBranchAddress ("akt4hi_em_jet_pt", &precalib_jet_pt);
   //tree->SetBranchAddress ("akt4hi_em_jet_eta", &precalib_jet_eta);
   //tree->SetBranchAddress ("akt4hi_em_jet_phi", &precalib_jet_phi);
   //tree->SetBranchAddress ("akt4hi_em_jet_e", &precalib_jet_e);

   //tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_pt", &jet_pt);
   //tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_eta", &jet_eta);
   //tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_phi", &jet_phi);
   //tree->SetBranchAddress ("akt4hi_constit_xcalib_jet_e", &jet_e);

   ////tree->SetBranchAddress ("akt4hi_constit_etajes_jet_pt", &jet_pt);
   ////tree->SetBranchAddress ("akt4hi_constit_etajes_jet_eta", &jet_eta);
   ////tree->SetBranchAddress ("akt4hi_constit_etajes_jet_phi", &jet_phi);
   ////tree->SetBranchAddress ("akt4hi_constit_etajes_jet_e", &jet_e);

   //tree->SetBranchAddress ("akt4hi_constit_jet_pt", &precalib_jet_pt);
   //tree->SetBranchAddress ("akt4hi_constit_jet_eta", &precalib_jet_eta);
   //tree->SetBranchAddress ("akt4hi_constit_jet_phi", &precalib_jet_phi);
   //tree->SetBranchAddress ("akt4hi_constit_jet_e", &precalib_jet_e);

   //tree->SetBranchAddress ("akt4emtopo_jet_n", &jet_n);

   //tree->SetBranchAddress ("akt4emtopo_em_jet_pt", &jet_pt);
   //tree->SetBranchAddress ("akt4emtopo_em_jet_eta", &jet_eta);
   //tree->SetBranchAddress ("akt4emtopo_em_jet_phi", &jet_phi);
   //tree->SetBranchAddress ("akt4emtopo_em_jet_e", &jet_e);

   //tree->SetBranchStatus ("akt4hi_jet_n", 0);

   //tree->SetBranchStatus ("akt4hi_em_xcalib_jet_pt", 0);
   //tree->SetBranchStatus ("akt4hi_em_xcalib_jet_eta", 0);
   //tree->SetBranchStatus ("akt4hi_em_xcalib_jet_phi", 0);
   //tree->SetBranchStatus ("akt4hi_em_xcalib_jet_e", 0);

   tree->SetBranchStatus ("akt4hi_em_etajes_jet_pt", 0);
   tree->SetBranchStatus ("akt4hi_em_etajes_jet_eta", 0);
   tree->SetBranchStatus ("akt4hi_em_etajes_jet_phi", 0);
   tree->SetBranchStatus ("akt4hi_em_etajes_jet_e", 0);

   tree->SetBranchStatus ("akt4hi_em_jet_pt", 0);
   tree->SetBranchStatus ("akt4hi_em_jet_eta", 0);
   tree->SetBranchStatus ("akt4hi_em_jet_phi", 0);
   tree->SetBranchStatus ("akt4hi_em_jet_e", 0);

   tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_pt", 0);
   tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_eta", 0);
   tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_phi", 0);
   tree->SetBranchStatus ("akt4hi_constit_xcalib_jet_e", 0);

   tree->SetBranchStatus ("akt4hi_constit_etajes_jet_pt", 0);
   tree->SetBranchStatus ("akt4hi_constit_etajes_jet_eta", 0);
   tree->SetBranchStatus ("akt4hi_constit_etajes_jet_phi", 0);
   tree->SetBranchStatus ("akt4hi_constit_etajes_jet_e", 0);

   tree->SetBranchStatus ("akt4hi_constit_jet_pt", 0);
   tree->SetBranchStatus ("akt4hi_constit_jet_eta", 0);
   tree->SetBranchStatus ("akt4hi_constit_jet_phi", 0);
   tree->SetBranchStatus ("akt4hi_constit_jet_e", 0);

   tree->SetBranchStatus ("akt4emtopo_jet_n", 0);

   tree->SetBranchStatus ("akt4emtopo_em_jet_pt", 0);
   tree->SetBranchStatus ("akt4emtopo_em_jet_eta", 0);
   tree->SetBranchStatus ("akt4emtopo_em_jet_phi", 0);
   tree->SetBranchStatus ("akt4emtopo_em_jet_e", 0);

   tree->SetBranchStatus ("akt4emtopo_calib_jet_pt", 0);
   tree->SetBranchStatus ("akt4emtopo_calib_jet_eta", 0);
   tree->SetBranchStatus ("akt4emtopo_calib_jet_phi", 0);
   tree->SetBranchStatus ("akt4emtopo_calib_jet_e", 0);
  }

  if (isMC) {
   tree->SetBranchAddress ("truth_jet_n", &truth_jet_n);
   tree->SetBranchAddress ("truth_jet_pt", &truth_jet_pt);
   tree->SetBranchAddress ("truth_jet_eta", &truth_jet_eta);
   tree->SetBranchAddress ("truth_jet_phi", &truth_jet_phi);
   tree->SetBranchAddress ("truth_jet_e", &truth_jet_e);

   if (getElectrons) {
    tree->SetBranchAddress ("truth_electron_n", &truth_electron_n);
    tree->SetBranchAddress ("truth_electron_pt", &truth_electron_pt);
    tree->SetBranchAddress ("truth_electron_eta", &truth_electron_eta);
    tree->SetBranchAddress ("truth_electron_phi", &truth_electron_phi);
    tree->SetBranchAddress ("truth_electron_charge", &truth_electron_charge);
   }
   else { 
    tree->SetBranchStatus ("truth_electron_n", 0);
    tree->SetBranchStatus ("truth_electron_pt", 0);
    tree->SetBranchStatus ("truth_electron_eta", 0);
    tree->SetBranchStatus ("truth_electron_phi", 0);
    tree->SetBranchStatus ("truth_electron_charge", 0);
   }

   if (getMuons) {
    tree->SetBranchAddress ("truth_muon_n", &truth_muon_n);
    tree->SetBranchAddress ("truth_muon_pt", &truth_muon_pt);
    tree->SetBranchAddress ("truth_muon_eta", &truth_muon_eta);
    tree->SetBranchAddress ("truth_muon_phi", &truth_muon_phi);
    tree->SetBranchAddress ("truth_muon_charge", &truth_muon_charge);
   }
   else {
    tree->SetBranchStatus ("truth_muon_n", 0);
    tree->SetBranchStatus ("truth_muon_pt", 0);
    tree->SetBranchStatus ("truth_muon_eta", 0);
    tree->SetBranchStatus ("truth_muon_phi", 0);
    tree->SetBranchStatus ("truth_muon_charge", 0);
   }

   if (getPhotons) {
    tree->SetBranchAddress ("truth_photon_n", &truth_photon_n);
    tree->SetBranchAddress ("truth_photon_pt", &truth_photon_pt);
    tree->SetBranchAddress ("truth_photon_eta", &truth_photon_eta);
    tree->SetBranchAddress ("truth_photon_phi", &truth_photon_phi);
   }
   else {
    tree->SetBranchStatus ("truth_photon_n", 0);
    tree->SetBranchStatus ("truth_photon_pt", 0);
    tree->SetBranchStatus ("truth_photon_eta", 0);
    tree->SetBranchStatus ("truth_photon_phi", 0);

   }
  }

  if (getElectrons) {
   tree->SetBranchAddress ("electron_n", &electron_n);
   tree->SetBranchAddress ("electron_pt", &electron_pt);
   tree->SetBranchAddress ("electron_eta", &electron_eta);
   tree->SetBranchAddress ("electron_phi", &electron_phi);
   tree->SetBranchAddress ("electron_charge", &electron_charge);
   tree->SetBranchAddress ("electron_tight", &electron_tight);
   tree->SetBranchAddress ("electron_loose", &electron_loose);
   tree->SetBranchAddress ("electron_d0sig", &electron_d0sig);
   tree->SetBranchAddress ("electron_delta_z0_sin_theta", &electron_delta_z0_sin_theta); 
  }
  else {
   tree->SetBranchStatus ("electron_n", 0);
   tree->SetBranchStatus ("electron_pt", 0);
   tree->SetBranchStatus ("electron_eta", 0);
   tree->SetBranchStatus ("electron_phi", 0);
   tree->SetBranchStatus ("electron_charge", 0);
   tree->SetBranchStatus ("electron_tight", 0);
   tree->SetBranchStatus ("electron_loose", 0);
   tree->SetBranchStatus ("electron_d0sig", 0);
   tree->SetBranchStatus ("electron_delta_z0_sin_theta", 0); 
  }

  if (getMuons) {
   tree->SetBranchAddress ("muon_n", &muon_n);
   tree->SetBranchAddress ("muon_pt", &muon_pt);
   tree->SetBranchAddress ("muon_eta", &muon_eta);
   tree->SetBranchAddress ("muon_phi", &muon_phi);
   tree->SetBranchAddress ("muon_charge", &muon_charge);
   tree->SetBranchAddress ("muon_quality", &muon_quality);
   tree->SetBranchAddress ("muon_tight", &muon_tight);
   tree->SetBranchAddress ("muon_loose", &muon_loose);
  }
  else {
   tree->SetBranchStatus ("muon_n", 0);
   tree->SetBranchStatus ("muon_pt", 0);
   tree->SetBranchStatus ("muon_eta", 0);
   tree->SetBranchStatus ("muon_phi", 0);
   tree->SetBranchStatus ("muon_charge", 0);
   tree->SetBranchStatus ("muon_quality", 0);
   tree->SetBranchStatus ("muon_tight", 0);
   tree->SetBranchStatus ("muon_loose", 0);
  }

  if (getPhotons) {
   tree->SetBranchAddress ("photon_n", &photon_n);
   tree->SetBranchAddress ("photon_pt", &photon_pt);
   tree->SetBranchAddress ("photon_eta", &photon_eta);
   tree->SetBranchAddress ("photon_phi", &photon_phi);
   tree->SetBranchAddress ("photon_tight", &photon_tight);
   tree->SetBranchAddress ("photon_isem", &photon_isem);
   tree->SetBranchAddress ("photon_convFlag", &photon_convFlag);
   tree->SetBranchAddress ("photon_Rconv", &photon_Rconv);
   tree->SetBranchAddress ("photon_topoetcone40", &photon_topoetcone40);
  }
  else {
   tree->SetBranchStatus ("photon_n", 0);
   tree->SetBranchStatus ("photon_pt", 0);
   tree->SetBranchStatus ("photon_eta", 0);
   tree->SetBranchStatus ("photon_phi", 0);
   tree->SetBranchStatus ("photon_tight", 0);
   tree->SetBranchStatus ("photon_isem", 0);
   tree->SetBranchStatus ("photon_convFlag", 0);
   tree->SetBranchStatus ("photon_Rconv", 0);
   tree->SetBranchStatus ("photon_topoetcone40", 0);
  }

  return;
}

void TreeVariables :: PrintAll (const long long entry) {
  tree->GetEntry (entry);
  cout << endl << "////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "// Getting event " << entry << "..." << endl;
  cout << "////////////////////////////////////////////////////////////////////////////////" << endl;

  if (!getSimpleJets) {

   if (getHIJets) {
    cout << endl << "Anti-Kt R=0.4 HI Insitu/xCalib + EtaJES jets (EM scale):" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4hi_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4hi_em_xcalib_jet_pt->at (j)
           << setw (10) << akt4hi_em_xcalib_jet_eta->at (j)
           << setw (10) << akt4hi_em_xcalib_jet_phi->at (j)
           << setw (10) << akt4hi_em_xcalib_jet_e->at (j) << endl;
    }

    cout << endl << "Anti-Kt R=0.4 HI EtaJES jets (EM scale):" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4hi_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4hi_em_etajes_jet_pt->at (j)
           << setw (10) << akt4hi_em_etajes_jet_eta->at (j)
           << setw (10) << akt4hi_em_etajes_jet_phi->at (j)
           << setw (10) << akt4hi_em_etajes_jet_e->at (j) << endl;
    }
 
    cout << endl << "Anti-Kt R=0.4 HI Reco jets (EM scale):" << endl; 
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4hi_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4hi_em_jet_pt->at (j)
           << setw (10) << akt4hi_em_jet_eta->at (j)
           << setw (10) << akt4hi_em_jet_phi->at (j)
           << setw (10) << akt4hi_em_jet_e->at (j) << endl;
    }

    cout << endl << "Anti-Kt R=0.4 HI Insitu/xCalib + EtaJES jets (Constit scale):" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4hi_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4hi_constit_xcalib_jet_pt->at (j)
           << setw (10) << akt4hi_constit_xcalib_jet_eta->at (j)
           << setw (10) << akt4hi_constit_xcalib_jet_phi->at (j)
           << setw (10) << akt4hi_constit_xcalib_jet_e->at (j) << endl;
    }

    cout << endl << "Anti-Kt R=0.4 HI EtaJES jets (Constit scale):" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4hi_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4hi_constit_etajes_jet_pt->at (j)
           << setw (10) << akt4hi_constit_etajes_jet_eta->at (j)
           << setw (10) << akt4hi_constit_etajes_jet_phi->at (j)
           << setw (10) << akt4hi_constit_etajes_jet_e->at (j) << endl;
    }
 
    cout << endl << "Anti-Kt R=0.4 HI Reco jets (Constit scale):" << endl; 
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4hi_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4hi_constit_jet_pt->at (j)
           << setw (10) << akt4hi_constit_jet_eta->at (j)
           << setw (10) << akt4hi_constit_jet_phi->at (j)
           << setw (10) << akt4hi_constit_jet_e->at (j) << endl;
    }
   } // end if getHIJets

   if (getEMTopoJets) {
    cout << endl << "Anti-Kt R=0.4 EMTopo Calibrated jets:" << endl;
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4emtopo_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4emtopo_calib_jet_pt->at (j)
           << setw (10) << akt4emtopo_calib_jet_eta->at (j)
           << setw (10) << akt4emtopo_calib_jet_phi->at (j)
           << setw (10) << akt4emtopo_calib_jet_e->at (j) << endl;
    }

    cout << endl << "Anti-Kt R=0.4 EMTopo Reco jets:" << endl; 
    cout << setw (2) << "#"
         << setw (10) << "Pt"
         << setw (10) << "Eta"
         << setw (10) << "Phi"
         << setw (10) << "E" << endl;
    for (int j = 0; j < akt4emtopo_jet_n; j++) {
      cout << setw (2) << j
           << setw (10) << akt4emtopo_em_jet_pt->at (j)
           << setw (10) << akt4emtopo_em_jet_eta->at (j)
           << setw (10) << akt4emtopo_em_jet_phi->at (j)
           << setw (10) << akt4emtopo_em_jet_e->at (j) << endl;
    }
   } // end if getEMTopoJets
  } // end if not getSimpleJets
  else { // if getSimpleJets
   cout << endl << "Calibrated jets:" << endl;
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

   cout << endl << "Reco jets:" << endl; 
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi"
        << setw (10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw (2) << j
          << setw (10) << precalib_jet_pt->at (j)
          << setw (10) << precalib_jet_eta->at (j)
          << setw (10) << precalib_jet_phi->at (j)
          << setw (10) << precalib_jet_e->at (j) << endl;
   }
  }

  if (isMC) {
   cout << endl << "Truth jets:" << endl;
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

   for (int j = 0; j < ntrk; j++) {
     cout << setw (10) << j
          << setw (10) << trk_quality_4->at (j)
          << setw (10) << trk_pt->at (j)
          << setw (10) << trk_eta->at (j)
          << setw (10) << trk_phi->at (j) << endl;
   }
  }

  if (getElectrons) {
   cout << endl << "Calibrated electrons:" << endl;
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi"
        << setw (10) << "Charge"
        << setw (10) << "LHLoose?"
        << setw (10) << "d0sig"
        << setw (15) << "dz0_sin_theta" << endl;
   for (int j = 0; j < electron_n; j++) {
     cout << setw (2) << j
          << setw (10) << electron_pt->at (j)
          << setw (10) << electron_eta->at (j)
          << setw (10) << electron_phi->at (j)
          << setw (10) << electron_charge->at (j)
          << setw (10) << electron_loose->at (j)
          << setw (10) << electron_d0sig->at (j)
          << setw (15) << electron_delta_z0_sin_theta->at (j) << endl;
   }

   cout << endl << "Truth electrons:" << endl;
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi"
        << setw (10) << "Charge" << endl;
   for (int j = 0; j < truth_electron_n; j++) {
     cout << setw (2) << j
          << setw (10) << truth_electron_pt->at (j)
          << setw (10) << truth_electron_eta->at (j)
          << setw (10) << truth_electron_phi->at (j)
          << setw (10) << truth_electron_charge->at (j) << endl;
   }
  } // end if getElectrons

  if (getMuons) {
   cout << endl << "Calibrated muons:" << endl;
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi"
        << setw (10) << "Charge"
        << setw (10) << "Loose?" << endl;

   for (int j = 0; j < muon_n; j++) {
     cout << setw (2) << j
          << setw (10) << muon_pt->at (j)
          << setw (10) << muon_eta->at (j)
          << setw (10) << muon_phi->at (j)
          << setw (10) << muon_charge->at (j)
          << setw (10) << muon_loose->at (j) << endl;
   }

   cout << endl << "Truth muons:" << endl;
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi"
        << setw (10) << "Charge" << endl;
   for (int j = 0; j < truth_muon_n; j++) {
     cout << setw (2) << j
          << setw (10) << truth_muon_pt->at (j)
          << setw (10) << truth_muon_eta->at (j)
          << setw (10) << truth_muon_phi->at (j)
          << setw (10) << truth_muon_charge->at (j) << endl;
   }
  } // end if getMuons

  if (getPhotons) {
   cout << endl << "Calibrated photons:" << endl;
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi"
        << setw (10) << "Tight?"
        << setw (10) << "IsoET" << endl;
   for (int j = 0; j < photon_n; j++) {
     cout << setw (2) << j
          << setw (10) << photon_pt->at (j)
          << setw (10) << photon_eta->at (j)
          << setw (10) << photon_phi->at (j)
          << setw (10) << photon_tight->at (j)
          << setw (10) << photon_topoetcone40->at (j) << endl;
   }

   cout << endl << "Truth photons:" << endl;
   cout << setw (2) << "#"
        << setw (10) << "Pt"
        << setw (10) << "Eta"
        << setw (10) << "Phi" << endl;
   for (int j = 0; j < truth_photon_n; j++) {
     cout << setw (2) << j
          << setw (10) << truth_photon_pt->at (j)
          << setw (10) << truth_photon_eta->at (j)
          << setw (10) << truth_photon_phi->at (j) << endl;
   }
  } // end if getPhotons

  return;
}

void TreeVariables :: SetGetMCInfo (const bool _getMCInfo, const double _crossSection_microbarns, const double _filterEfficiency, const int _numberEvents) {
  getMCInfo = _getMCInfo;
  crossSection_microbarns = _crossSection_microbarns;
  filterEfficiency = _filterEfficiency;
  numberEvents = _numberEvents;
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

void TreeVariables :: SetGetElectrons (const bool _getElectrons) {
  getElectrons = _getElectrons;
}

void TreeVariables :: SetGetPhotons (const bool _getPhotons) {
  getPhotons = _getPhotons;
}

void TreeVariables :: SetGetMuons (const bool _getMuons) {
  getMuons = _getMuons;
}

} // end namespace
