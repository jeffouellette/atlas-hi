#ifndef __TreeVariables_h__
#define __TreeVariables_h__

#include <TTree.h>
#include <vector>

/** 
 * This file defines the TreeVariables struct used by this analysis.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */

using namespace std;

namespace atlashi {

struct TreeVariables {
  private:
   TTree* tree;

   bool isMC = false;
   bool getMCInfo = true;
   bool getCollisionRateInfo = false;
   bool getVertices = false;
   bool getFCals = false;
   bool getTracks = false;
   bool getSimpleJets = false;
   bool getHIJets = false;
   bool getEMTopoJets = false;
   bool getElectrons = false;
   bool getPhotons = false;
   bool getMuons = false;
   bool getTruthJets = false;
   bool getTruthElectrons = false;
   bool getTruthPhotons = false;
   bool getTruthMuons = false;

  public:

   // public functions
   TreeVariables (TTree* t, const bool _isMC = false);
   ~TreeVariables ();    
   void SetBranchAddresses ();
   void PrintAll (const long long entry);

   // setter functions
   void SetGetMCInfo (const bool _getMCInfo = true, const double _crossSection_microbarns = 0, const double _filterEfficiency = 0, const int _numberEvents = 0);
   void SetGetCollisionRateInfo (const bool _getCollisionRateInfo = true);
   void SetGetVertices (const bool _getVertices = true);
   void SetGetFCals (const bool _getFCals = true);
   void SetGetTracks (const bool _getTracks = true);
   void SetGetSimpleJets (const bool _getSimpleJets = true);
   void SetGetHIJets (const bool _getHIJets = true);
   void SetGetEMTopoJets (const bool _getEMTopoJets = true);
   void SetGetElectrons (const bool _getElectrons = true);
   void SetGetPhotons (const bool _getPhotons = true);
   void SetGetMuons (const bool _getMuons = true);
   void SetGetTruthJets (const bool _getTruthJets = true);
   void SetGetTruthElectrons (const bool _getTruthElectrons = true);
   void SetGetTruthPhotons (const bool _getTruthPhotons = true);
   void SetGetTruthMuons (const bool _getTruthMuons = true);

   // public variables
   int eventNumber = 0;
   int runNumber = 0;
   unsigned int lumiBlock = 0;

   int numberEvents = 0;
   double crossSection_microbarns = 0;
   double filterEfficiency = 0;

   float actualInteractionsPerCrossing = 0;
   float averageInteractionsPerCrossing = 0;

   int nvert = 0;
   vector<float>* vert_x = NULL;
   vector<float>* vert_y = NULL;
   vector<float>* vert_z = NULL;
   vector<int>* vert_ntrk = NULL;
   vector<int>* vert_type = NULL;

   // fcal energy
   float fcalA_et;
   float fcalC_et;
 
   // tracking info (0th or primary vertex only)
   int ntrk = 0;
   vector<bool>* trk_quality_4;
   vector<float>* trk_d0;
   vector<float>* trk_z0;
   vector<float>* trk_theta;
   vector<float>* trk_charge;
   vector<float>* trk_pt;
   vector<float>* trk_eta;
   vector<float>* trk_phi;

   int clean_jet_n = 0;
   int total_jet_n = 0;
   
   int jet_n = 0; 

   vector<float>* jet_pt = NULL;
   vector<float>* jet_eta = NULL;
   vector<float>* jet_phi = NULL;
   vector<float>* jet_e = NULL;

   vector<float>* precalib_jet_pt = NULL;
   vector<float>* precalib_jet_eta = NULL;
   vector<float>* precalib_jet_phi = NULL;
   vector<float>* precalib_jet_e = NULL;

   int akt4hi_jet_n = 0;

   vector<float>* akt4hi_em_xcalib_jet_pt = NULL;
   vector<float>* akt4hi_em_xcalib_jet_eta = NULL;
   vector<float>* akt4hi_em_xcalib_jet_phi = NULL;
   vector<float>* akt4hi_em_xcalib_jet_e = NULL;

   vector<float>* akt4hi_em_etajes_jet_pt = NULL;
   vector<float>* akt4hi_em_etajes_jet_eta = NULL;
   vector<float>* akt4hi_em_etajes_jet_phi = NULL;
   vector<float>* akt4hi_em_etajes_jet_e = NULL;

   vector<float>* akt4hi_em_jet_pt = NULL;
   vector<float>* akt4hi_em_jet_eta = NULL;
   vector<float>* akt4hi_em_jet_phi = NULL;
   vector<float>* akt4hi_em_jet_e = NULL;

   vector<float>* akt4hi_constit_xcalib_jet_pt = NULL;
   vector<float>* akt4hi_constit_xcalib_jet_eta = NULL;
   vector<float>* akt4hi_constit_xcalib_jet_phi = NULL;
   vector<float>* akt4hi_constit_xcalib_jet_e = NULL;

   vector<float>* akt4hi_constit_etajes_jet_pt = NULL;
   vector<float>* akt4hi_constit_etajes_jet_eta = NULL;
   vector<float>* akt4hi_constit_etajes_jet_phi = NULL;
   vector<float>* akt4hi_constit_etajes_jet_e = NULL;

   vector<float>* akt4hi_constit_jet_pt = NULL;
   vector<float>* akt4hi_constit_jet_eta = NULL;
   vector<float>* akt4hi_constit_jet_phi = NULL;
   vector<float>* akt4hi_constit_jet_e = NULL;

   int akt4emtopo_jet_n = 0;

   vector<float>* akt4emtopo_em_jet_pt = NULL;
   vector<float>* akt4emtopo_em_jet_eta = NULL;
   vector<float>* akt4emtopo_em_jet_phi = NULL;
   vector<float>* akt4emtopo_em_jet_e = NULL;

   vector<float>* akt4emtopo_calib_jet_pt = NULL;
   vector<float>* akt4emtopo_calib_jet_eta = NULL;
   vector<float>* akt4emtopo_calib_jet_phi = NULL;
   vector<float>* akt4emtopo_calib_jet_e = NULL;

   int truth_jet_n = 0;
   vector<float>* truth_jet_pt = NULL;
   vector<float>* truth_jet_eta = NULL;
   vector<float>* truth_jet_phi = NULL;
   vector<float>* truth_jet_e = NULL;
   
   int muon_n = 0;
   vector<float>* muon_pt = NULL;
   vector<float>* muon_eta = NULL;
   vector<float>* muon_phi = NULL;
   vector<int>* muon_quality = NULL;
   vector<int>* muon_charge = NULL;
   vector<bool>* muon_tight = NULL;
   vector<bool>* muon_loose = NULL;

   int truth_muon_n = 0;
   vector<float>* truth_muon_pt = NULL;
   vector<float>* truth_muon_eta = NULL;
   vector<float>* truth_muon_phi = NULL;
   vector<int>* truth_muon_charge = NULL;
   
   int electron_n = 0;
   vector<float>* electron_pt = NULL;
   vector<float>* electron_eta = NULL;
   vector<float>* electron_phi = NULL;
   vector<int>* electron_charge = NULL;
   vector<bool>* electron_loose = NULL;
   vector<bool>* electron_tight = NULL;
   vector<float>* electron_d0sig = NULL;
   vector<float>* electron_delta_z0_sin_theta = NULL;

   int truth_electron_n = 0;
   vector<float>* truth_electron_pt = NULL;
   vector<float>* truth_electron_eta = NULL;
   vector<float>* truth_electron_phi = NULL;
   vector<int>* truth_electron_charge = NULL;
   
   int photon_n = 0;
   vector<float>* photon_pt = NULL;
   vector<float>* photon_eta = NULL;
   vector<float>* photon_phi = NULL;
   vector<bool>* photon_tight = NULL;
   vector<bool>* photon_loose = NULL;
   vector<unsigned int>* photon_isem = NULL;
   vector<int>* photon_convFlag = NULL;
   vector<float>* photon_Rconv = NULL;
   vector<float>* photon_topoetcone40 = NULL;

   int truth_photon_n = 0;
   vector<float>* truth_photon_pt = NULL;
   vector<float>* truth_photon_eta = NULL;
   vector<float>* truth_photon_phi = NULL;
};

} // end namespace

#endif
