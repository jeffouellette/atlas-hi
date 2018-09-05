#ifndef __TreeVariables_h__
#define __TreeVariables_h__

#include <set>
#include <TTree.h>

using namespace std;

namespace pPb8TeV2016DijetAnalysis {

struct TreeVariables {
  private:
   TTree* t;

   bool isMC = false;
   bool getCollisionRateInfo = false;
   bool getVertices = false;
   bool getFCals = false;
   bool getTracks = false;
   bool getSimpleJets = false;
   bool getHIJets = false;
   bool getEMTopoJets = false;

   void BranchHelper (const std::set<TString>& active, const char* name, int* ptr, const bool get=false);
   void BranchHelper (const std::set<TString>& active, const char* name, float* ptr, const bool get=false);
   void BranchHelper (const std::set<TString>& active, const char* name, double* ptr, const bool get=false);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<int>** ptr, const bool get=false);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<float>** ptr, const bool get=false);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<double>** ptr, const bool get=false);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<bool>** ptr, const bool get=false);

  public:
   TreeVariables (TTree* t, const bool _isMC);
   ~TreeVariables ();
   void SetBranchAddresses (const std::set<TString>& active);

   // setter functions
   void SetGetCollisionRateInfo (const bool _getCollisionRateInfo = true);
   void SetGetVertices (const bool _getVertices = true);
   void SetGetFCals (const bool _getFCals = true);
   void SetGetTracks (const bool _getTracks = true);
   void SetGetSimpleJets (const bool _getSimpleJets = true);
   void SetGetHIJets (const bool _getHIJets = true);
   void SetGetEMTopoJets (const bool _getEMTopoJets = true);

   int runNumber = 0;
   int eventNumber = 0;
   int lumiBlock = 0;
   double eventWeight = 0;

   float actualInteractionsPerCrossing = 0;
   float averageInteractionsPerCrossing = 0;

   int nvert = 0;
   vector<float>* vert_x = NULL;
   vector<float>* vert_y = NULL;
   vector<float>* vert_z = NULL;
   vector<int>* vert_ntrk = NULL;
   vector<int>* vert_type = NULL;

   int ntrk = 0;
   vector<bool>* trk_quality_4 = NULL;
   vector<float>* trk_d0 = NULL;
   vector<float>* trk_z0 = NULL;
   vector<float>* trk_theta = NULL;
   vector<float>* trk_charge = NULL;
   vector<float>* trk_pt = NULL;
   vector<float>* trk_eta = NULL;
   vector<float>* trk_phi = NULL;

   float fcalA_et = 0;
   float fcalC_et = 0;
  
   int total_jets = 0;
   int clean_jets = 0; 
   
   int njet = 0; 
   vector<float>* jet_pt = NULL;
   vector<float>* jet_eta = NULL;
   vector<float>* jet_phi = NULL;
   vector<float>* jet_e = NULL;

   int akt4hi_njet = 0;
   vector<float>* akt4hi_em_jet_pt = NULL;
   vector<float>* akt4hi_em_jet_eta = NULL;
   vector<float>* akt4hi_em_jet_phi = NULL;
   vector<float>* akt4hi_em_jet_e = NULL;

   vector<float>* akt4hi_etajes_jet_pt = NULL;
   vector<float>* akt4hi_etajes_jet_eta = NULL;
   vector<float>* akt4hi_etajes_jet_phi = NULL;
   vector<float>* akt4hi_etajes_jet_e = NULL;

   vector<float>* akt4hi_xcalib_jet_pt = NULL;
   vector<float>* akt4hi_xcalib_jet_eta = NULL;
   vector<float>* akt4hi_xcalib_jet_phi = NULL;
   vector<float>* akt4hi_xcalib_jet_e = NULL;

   int akt4truth_njet = 0;
   vector<float>* akt4truth_jet_pt = NULL;
   vector<float>* akt4truth_jet_eta = NULL;
   vector<float>* akt4truth_jet_phi = NULL;
   vector<float>* akt4truth_jet_e = NULL;

   int hlt_njet = 0;
   vector<float>* hlt_jet_pt = NULL;
   vector<float>* hlt_jet_eta = NULL;
   vector<float>* hlt_jet_phi = NULL;
   vector<float>* hlt_jet_e = NULL;
};

} // end namespace

#endif
