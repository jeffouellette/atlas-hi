#include <set>

struct TreeVariables {
  private:
   TTree* tree;
   void BranchHelper (const std::set<TString>& active, const char* name, int* ptr);
   void BranchHelper (const std::set<TString>& active, const char* name, float* ptr);
   void BranchHelper (const std::set<TString>& active, const char* name, double* ptr);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<int>** ptr);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<float>** ptr);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<double>** ptr);
   void BranchHelper (const std::set<TString>& active, const char* name, vector<bool>** ptr);

  public:
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
   vector<float>* init_jet_pt = NULL;
   vector<float>* init_jet_eta = NULL;
   vector<float>* init_jet_phi = NULL;
   vector<float>* init_jet_e = NULL;

   vector<float>* etajes_jet_pt = NULL;
   vector<float>* etajes_jet_eta = NULL;
   vector<float>* etajes_jet_phi = NULL;
   vector<float>* etajes_jet_e = NULL;

   vector<float>* xcalib_jet_pt = NULL;
   vector<float>* xcalib_jet_eta = NULL;
   vector<float>* xcalib_jet_phi = NULL;
   vector<float>* xcalib_jet_e = NULL;

   int truth_njet = 0;
   vector<float>* truth_jet_pt = NULL;
   vector<float>* truth_jet_eta = NULL;
   vector<float>* truth_jet_phi = NULL;
   vector<float>* truth_jet_e = NULL;

   int hlt_njet = 0;
   vector<float>* hlt_jet_pt = NULL;
   vector<float>* hlt_jet_eta = NULL;
   vector<float>* hlt_jet_phi = NULL;
   vector<float>* hlt_jet_e = NULL;
   
   TreeVariables (TTree* t);
   ~TreeVariables ();
   void SetBranchAddresses (const bool isMC, const std::set<TString>& active);
};


TreeVariables::TreeVariables(TTree* t) {
  tree = t;
}


TreeVariables::~TreeVariables() {
  if (vert_x) delete vert_x;
  if (vert_y) delete vert_y;
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

  if (init_jet_pt) delete init_jet_pt;
  if (init_jet_eta) delete init_jet_eta;
  if (init_jet_phi) delete init_jet_phi;
  if (init_jet_e) delete init_jet_e;
  if (etajes_jet_pt) delete etajes_jet_pt;
  if (etajes_jet_eta) delete etajes_jet_eta;
  if (etajes_jet_phi) delete etajes_jet_phi;
  if (etajes_jet_e) delete etajes_jet_e;
  if (xcalib_jet_pt) delete xcalib_jet_pt;
  if (xcalib_jet_eta) delete xcalib_jet_eta;
  if (xcalib_jet_phi) delete xcalib_jet_phi;
  if (xcalib_jet_e) delete xcalib_jet_e;

  if (truth_jet_pt) delete truth_jet_pt;
  if (truth_jet_eta) delete truth_jet_eta;
  if (truth_jet_phi) delete truth_jet_phi;
  if (truth_jet_e) delete truth_jet_e;

  if (hlt_jet_pt) delete hlt_jet_pt;
  if (hlt_jet_eta) delete hlt_jet_eta;
  if (hlt_jet_phi) delete hlt_jet_phi;
  if (hlt_jet_e) delete hlt_jet_e;
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, int* ptr) {
  if (active.count(TString(name)) == 1) tree->SetBranchAddress(name, ptr);
  else tree->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, float* ptr) {
  if (active.count(TString(name)) == 1) tree->SetBranchAddress(name, ptr);
  else tree->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, double* ptr) {
  if (active.count(TString(name)) == 1) tree->SetBranchAddress(name, ptr);
  else tree->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, vector<int>** ptr) {
  if (active.count(TString(name)) == 1) tree->SetBranchAddress(name, ptr);
  else tree->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, vector<float>** ptr) {
  if (active.count(TString(name)) == 1) tree->SetBranchAddress(name, ptr);
  else tree->SetBranchStatus(name, 0);
}


void TreeVariables::BranchHelper (const std::set<TString>& active, const char* name, vector<bool>** ptr) {
  if (active.count(TString(name)) == 1) tree->SetBranchAddress(name, ptr);
  else tree->SetBranchStatus(name, 0);
}


void TreeVariables::SetBranchAddresses (const bool isMC, const std::set<TString>& active) {
  BranchHelper (active, "runNumber", &runNumber);
  BranchHelper (active, "eventNumber", &eventNumber);
  BranchHelper (active, "lumiBlock", &lumiBlock);
  //if (isMC) BranchHelper (active, "EventWeight", &eventWeight);

  BranchHelper (active, "actualInteractionsPerCrossing", &actualInteractionsPerCrossing);
  BranchHelper (active, "averageInteractionsPerCrossing", &averageInteractionsPerCrossing);

  BranchHelper (active, "nvert", &nvert);
  //BranchHelper (active, "vert_x", &vert_x);
  //BranchHelper (active, "vert_y", &vert_y);
  BranchHelper (active, "vert_z", &vert_z);
  //BranchHelper (active, "vert_ntrk", &vert_ntrk);
  BranchHelper (active, "vert_type", &vert_type);

  //BranchHelper (active, "trk_quality_4", &trk_quality_4);
  //BranchHelper (active, "trk_d0", &trk_d0);
  //BranchHelper (active, "trk_z0", &trk_z0);
  //BranchHelper (active, "trk_theta", &trk_theta);
  //BranchHelper (active, "trk_charge", &trk_charge);
  //BranchHelper (active, "trk_pt", &trk_pt);
  //BranchHelper (active, "trk_eta", &trk_eta);
  //BranchHelper (active, "trk_phi", &trk_phi);

  BranchHelper (active, "fcalA_et", &fcalA_et);
  BranchHelper (active, "fcalC_et", &fcalC_et);

  BranchHelper (active, "total_jets", &total_jets);
  BranchHelper (active, "clean_jets", &clean_jets);

  BranchHelper (active, "njet", &njet);
  BranchHelper (active, "init_jet_pt", &init_jet_pt);
  BranchHelper (active, "init_jet_eta", &init_jet_eta);
  BranchHelper (active, "init_jet_phi", &init_jet_phi);
  BranchHelper (active, "init_jet_e", &init_jet_e);

  BranchHelper (active, "etajes_jet_pt", &etajes_jet_pt);
  BranchHelper (active, "etajes_jet_eta", &etajes_jet_eta);
  BranchHelper (active, "etajes_jet_phi", &etajes_jet_phi);
  BranchHelper (active, "etajes_jet_e", &etajes_jet_e);

  BranchHelper (active, "xcalib_jet_pt", &xcalib_jet_pt);
  BranchHelper (active, "xcalib_jet_eta", &xcalib_jet_eta);
  BranchHelper (active, "xcalib_jet_phi", &xcalib_jet_phi);
  BranchHelper (active, "xcalib_jet_e", &xcalib_jet_e);

  if (isMC) {
   BranchHelper (active, "truth_njet", &truth_njet);
   BranchHelper (active, "truth_jet_pt", &truth_jet_pt);
   BranchHelper (active, "truth_jet_eta", &truth_jet_eta);
   BranchHelper (active, "truth_jet_phi", &truth_jet_phi);
   BranchHelper (active, "truth_jet_e", &truth_jet_e);
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
