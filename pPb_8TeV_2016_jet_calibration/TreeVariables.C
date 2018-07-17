struct TreeVariables {
  private:
   TTree* tree;

  public:
   int eventNumber = 0;
   double eventWeight = 0;

   int nvert = 0;
   vector<float>* vert_z = NULL;
   vector<int>* vert_ntrk = NULL;
   vector<int>* vert_type = NULL;
   
   int jet_n = 0; 
   vector<float>* jet_pt = NULL;
   vector<float>* jet_eta = NULL;
   vector<float>* jet_phi = NULL;
   vector<float>* jet_e = NULL;
   vector<float>* init_jet_pt = NULL;
   
   int muon_n = 0;
   vector<float>* muon_pt = NULL;
   vector<float>* muon_eta = NULL;
   vector<float>* muon_phi = NULL;
   vector<int>* muon_quality = NULL;
   vector<int>* muon_charge = NULL;
   vector<bool>* muon_tight = NULL;
   vector<bool>* muon_loose = NULL;
   
   int electron_n = 0;
   vector<float>* electron_pt = NULL;
   vector<float>* electron_eta = NULL;
   vector<float>* electron_phi = NULL;
   vector<int>* electron_charge = NULL;
   vector<bool>* electron_loose = NULL;
   vector<bool>* electron_tight = NULL;
   
   int photon_n = 0;
   vector<float>* photon_pt = NULL;
   vector<float>* photon_eta = NULL;
   vector<float>* photon_phi = NULL;
   vector<bool>* photon_tight = NULL;
   vector<unsigned int>* photon_isem = NULL;
   vector<int>* photon_convFlag = NULL;
   vector<float>* photon_Rconv = NULL;
   vector<float>* photon_topoetcone40 = NULL;

   TreeVariables(TTree* t);
   ~TreeVariables();    
   void SetBranchAddresses(const bool isMC);
};


TreeVariables::TreeVariables(TTree* t) {
  tree = t;
}


TreeVariables::~TreeVariables() {
  if (vert_z) delete vert_z;
  if (vert_ntrk) delete vert_ntrk;
  if (vert_type) delete vert_type;

  if (jet_pt) delete jet_pt;
  if (jet_eta) delete jet_eta;
  if (jet_phi) delete jet_phi;
  if (jet_e) delete jet_e;
  if (init_jet_pt) delete init_jet_pt;

  if (muon_pt) delete muon_pt;
  if (muon_eta) delete muon_eta;
  if (muon_phi) delete muon_phi;
  if (muon_quality) delete muon_quality;
  if (muon_charge) delete muon_charge;
  if (muon_tight) delete muon_tight;
  if (muon_loose) delete muon_loose;

  if (electron_pt) delete electron_pt;
  if (electron_eta) delete electron_eta;
  if (electron_phi) delete electron_phi;
  if (electron_charge) delete electron_charge;
  if (electron_loose) delete electron_loose;
  if (electron_tight) delete electron_tight;

  if (photon_pt) delete photon_pt;
  if (photon_eta) delete photon_eta;
  if (photon_phi) delete photon_phi;
  if (photon_tight) delete photon_tight;
  if (photon_isem) delete photon_isem;
  if (photon_convFlag) delete photon_convFlag;
  if (photon_Rconv) delete photon_Rconv;
  if (photon_topoetcone40) delete photon_topoetcone40;
}


void TreeVariables::SetBranchAddresses(const bool isMC) {
  if (!isMC) tree->SetBranchAddress("eventNumber", &eventNumber);
  else tree->SetBranchAddress("EventWeight", &eventWeight);

  tree->SetBranchAddress("nvert", &nvert);
  tree->SetBranchAddress("vert_z", &vert_z);
  tree->SetBranchAddress("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress("vert_type", &vert_type);

  tree->SetBranchAddress("jet_n", &jet_n);
  tree->SetBranchAddress("xcalib_etajes_jet_pt", &jet_pt);
  tree->SetBranchAddress("xcalib_etajes_jet_eta", &jet_eta);
  tree->SetBranchAddress("xcalib_etajes_jet_phi", &jet_phi);
  tree->SetBranchAddress("xcalib_etajes_jet_e", &jet_e);
  //tree->SetBranchAddress("etajes_jet_pt", &jet_pt);
  //tree->SetBranchAddress("etajes_jet_eta", &jet_eta);
  //tree->SetBranchAddress("etajes_jet_phi", &jet_phi);
  //tree->SetBranchAddress("etajes_jet_e", &jet_e);
  //tree->SetBranchAddress("init_jet_pt", &jet_pt);
  //tree->SetBranchAddress("init_jet_eta", &jet_eta);
  //tree->SetBranchAddress("init_jet_phi", &jet_phi);
  //tree->SetBranchAddress("init_jet_e", &jet_e);

  tree->SetBranchAddress("electron_n", &electron_n);
  tree->SetBranchAddress("electron_pt", &electron_pt);
  tree->SetBranchAddress("electron_eta", &electron_eta);
  tree->SetBranchAddress("electron_phi", &electron_phi);
  tree->SetBranchAddress("electron_charge", &electron_charge);
  tree->SetBranchAddress("electron_tight", &electron_tight);
  tree->SetBranchAddress("electron_loose", &electron_loose);

  tree->SetBranchAddress("muon_n", &muon_n);
  tree->SetBranchAddress("muon_pt", &muon_pt);
  tree->SetBranchAddress("muon_eta", &muon_eta);
  tree->SetBranchAddress("muon_phi", &muon_phi);
  tree->SetBranchAddress("muon_charge", &muon_charge);
  tree->SetBranchAddress("muon_quality", &muon_quality);
  tree->SetBranchAddress("muon_tight", &muon_tight);
  tree->SetBranchAddress("muon_loose", &muon_loose);

  tree->SetBranchAddress("photon_n", &photon_n);
  tree->SetBranchAddress("photon_pt", &photon_pt);
  tree->SetBranchAddress("photon_eta", &photon_eta);
  tree->SetBranchAddress("photon_phi", &photon_phi);
  tree->SetBranchAddress("photon_tight", &photon_tight);
  tree->SetBranchAddress("photon_isem", &photon_isem);
  tree->SetBranchAddress("photon_convFlag", &photon_convFlag);
  tree->SetBranchAddress("photon_Rconv", &photon_Rconv);
  tree->SetBranchAddress("photon_topoetcone40", &photon_topoetcone40);

  return;
}
