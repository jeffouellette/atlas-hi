struct TreeVariables {
  private:
   TTree* tree;
   bool testTree;

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

   vector<float>* precalib_jet_pt = NULL;
   vector<float>* precalib_jet_eta = NULL;
   vector<float>* precalib_jet_phi = NULL;
   vector<float>* precalib_jet_e = NULL;

   vector<float>* em_xcalib_jet_pt = NULL;
   vector<float>* em_xcalib_jet_eta = NULL;
   vector<float>* em_xcalib_jet_phi = NULL;
   vector<float>* em_xcalib_jet_e = NULL;

   vector<float>* em_etajes_jet_pt = NULL;
   vector<float>* em_etajes_jet_eta = NULL;
   vector<float>* em_etajes_jet_phi = NULL;
   vector<float>* em_etajes_jet_e = NULL;

   vector<float>* emscale_jet_pt = NULL;
   vector<float>* emscale_jet_eta = NULL;
   vector<float>* emscale_jet_phi = NULL;
   vector<float>* emscale_jet_e = NULL;

   vector<float>* constit_xcalib_jet_pt = NULL;
   vector<float>* constit_xcalib_jet_eta = NULL;
   vector<float>* constit_xcalib_jet_phi = NULL;
   vector<float>* constit_xcalib_jet_e = NULL;

   vector<float>* constit_etajes_jet_pt = NULL;
   vector<float>* constit_etajes_jet_eta = NULL;
   vector<float>* constit_etajes_jet_phi = NULL;
   vector<float>* constit_etajes_jet_e = NULL;

   vector<float>* constit_jet_pt = NULL;
   vector<float>* constit_jet_eta = NULL;
   vector<float>* constit_jet_phi = NULL;
   vector<float>* constit_jet_e = NULL;

   int truth_jet_n = 0;
   vector<float>* truth_jet_pt = NULL;
   vector<float>* truth_jet_eta = NULL;
   vector<float>* truth_jet_phi = NULL;
   vector<float>* truth_jet_e = NULL;
   vector<unsigned int>* truth_jet_type = NULL;
   vector<unsigned int>* truth_jet_origin = NULL;
   
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
   vector<unsigned int>* truth_muon_type = NULL;
   vector<unsigned int>* truth_muon_origin = NULL;
   
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
   vector<unsigned int>* truth_electron_type = NULL;
   vector<unsigned int>* truth_electron_origin = NULL;
   
   int photon_n = 0;
   vector<float>* photon_pt = NULL;
   vector<float>* photon_eta = NULL;
   vector<float>* photon_phi = NULL;
   vector<bool>* photon_tight = NULL;
   vector<unsigned int>* photon_isem = NULL;
   vector<int>* photon_convFlag = NULL;
   vector<float>* photon_Rconv = NULL;
   vector<float>* photon_topoetcone40 = NULL;

   int truth_photon_n = 0;
   vector<float>* truth_photon_pt = NULL;
   vector<float>* truth_photon_eta = NULL;
   vector<float>* truth_photon_phi = NULL;
   vector<unsigned int>* truth_photon_type = NULL;
   vector<unsigned int>* truth_photon_origin = NULL;

   TreeVariables (TTree* t, const bool test = false);
   ~TreeVariables ();    
   void SetBranchAddresses (const bool isMC);
   void PrintAll (const long long entry);
};


TreeVariables::TreeVariables(TTree* t, const bool test = false) {
  tree = t;
  testTree = test;
}


TreeVariables::~TreeVariables() {
  if (vert_z) delete vert_z;
  if (vert_ntrk) delete vert_ntrk;
  if (vert_type) delete vert_type;

  if (jet_pt) delete jet_pt;
  if (jet_eta) delete jet_eta;
  if (jet_phi) delete jet_phi;
  if (jet_e) delete jet_e;

  if (emscale_jet_pt) delete emscale_jet_pt;
  if (emscale_jet_eta) delete emscale_jet_eta;
  if (emscale_jet_phi) delete emscale_jet_phi;
  if (emscale_jet_e) delete emscale_jet_e;

  if (constit_xcalib_jet_pt) delete constit_xcalib_jet_pt;
  if (constit_xcalib_jet_eta) delete constit_xcalib_jet_eta;
  if (constit_xcalib_jet_phi) delete constit_xcalib_jet_phi;
  if (constit_xcalib_jet_e) delete constit_xcalib_jet_e;

  if (constit_jet_pt) delete constit_jet_pt;
  if (constit_jet_eta) delete constit_jet_eta;
  if (constit_jet_phi) delete constit_jet_phi;
  if (constit_jet_e) delete constit_jet_e;

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
}


void TreeVariables::SetBranchAddresses(const bool isMC) {
  if (!isMC) tree->SetBranchAddress("eventNumber", &eventNumber);
  else tree->SetBranchAddress("EventWeight", &eventWeight);

  tree->SetBranchAddress("nvert", &nvert);
  tree->SetBranchAddress("vert_z", &vert_z);
  tree->SetBranchAddress("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress("vert_type", &vert_type);

  tree->SetBranchAddress("jet_n", &jet_n);

  if (testTree) {

   tree->SetBranchAddress("em_xcalib_jet_pt", &em_xcalib_jet_pt);
   tree->SetBranchAddress("em_xcalib_jet_eta", &em_xcalib_jet_eta);
   tree->SetBranchAddress("em_xcalib_jet_phi", &em_xcalib_jet_phi);
   tree->SetBranchAddress("em_xcalib_jet_e", &em_xcalib_jet_e);

   tree->SetBranchAddress("em_etajes_jet_pt", &em_etajes_jet_pt);
   tree->SetBranchAddress("em_etajes_jet_eta", &em_etajes_jet_eta);
   tree->SetBranchAddress("em_etajes_jet_phi", &em_etajes_jet_phi);
   tree->SetBranchAddress("em_etajes_jet_e", &em_etajes_jet_e);

   tree->SetBranchAddress("emscale_jet_pt", &emscale_jet_pt);
   tree->SetBranchAddress("emscale_jet_eta", &emscale_jet_eta);
   tree->SetBranchAddress("emscale_jet_phi", &emscale_jet_phi);
   tree->SetBranchAddress("emscale_jet_e", &emscale_jet_e);

   tree->SetBranchAddress("constit_xcalib_jet_pt", &constit_xcalib_jet_pt);
   tree->SetBranchAddress("constit_xcalib_jet_eta", &constit_xcalib_jet_eta);
   tree->SetBranchAddress("constit_xcalib_jet_phi", &constit_xcalib_jet_phi);
   tree->SetBranchAddress("constit_xcalib_jet_e", &constit_xcalib_jet_e);

   tree->SetBranchAddress("constit_etajes_jet_pt", &constit_etajes_jet_pt);
   tree->SetBranchAddress("constit_etajes_jet_eta", &constit_etajes_jet_eta);
   tree->SetBranchAddress("constit_etajes_jet_phi", &constit_etajes_jet_phi);
   tree->SetBranchAddress("constit_etajes_jet_e", &constit_etajes_jet_e);

   tree->SetBranchAddress("constit_jet_pt", &constit_jet_pt);
   tree->SetBranchAddress("constit_jet_eta", &constit_jet_eta);
   tree->SetBranchAddress("constit_jet_phi", &constit_jet_phi);
   tree->SetBranchAddress("constit_jet_e", &constit_jet_e);
  }
  else {
   tree->SetBranchAddress("em_xcalib_jet_pt", &jet_pt);
   tree->SetBranchAddress("em_xcalib_jet_eta", &jet_eta);
   tree->SetBranchAddress("em_xcalib_jet_phi", &jet_phi);
   tree->SetBranchAddress("em_xcalib_jet_e", &jet_e);

   //tree->SetBranchAddress("em_etajes_jet_pt", &jet_pt);
   //tree->SetBranchAddress("em_etajes_jet_eta", &jet_eta);
   //tree->SetBranchAddress("em_etajes_jet_phi", &jet_phi);
   //tree->SetBranchAddress("em_etajes_jet_e", &jet_e);

   tree->SetBranchAddress("emscale_jet_pt", &precalib_jet_pt);
   tree->SetBranchAddress("emscale_jet_eta", &precalib_jet_eta);
   tree->SetBranchAddress("emscale_jet_phi", &precalib_jet_phi);
   tree->SetBranchAddress("emscale_jet_e", &precalib_jet_e);

   //tree->SetBranchAddress("constit_xcalib_jet_pt", &jet_pt);
   //tree->SetBranchAddress("constit_xcalib_jet_eta", &jet_eta);
   //tree->SetBranchAddress("constit_xcalib_jet_phi", &jet_phi);
   //tree->SetBranchAddress("constit_xcalib_jet_e", &jet_e);

   //tree->SetBranchAddress("constit_etajes_jet_pt", &jet_pt);
   //tree->SetBranchAddress("constit_etajes_jet_eta", &jet_eta);
   //tree->SetBranchAddress("constit_etajes_jet_phi", &jet_phi);
   //tree->SetBranchAddress("constit_etajes_jet_e", &jet_e);

   //tree->SetBranchAddress("constit_jet_pt", &precalib_jet_pt);
   //tree->SetBranchAddress("constit_jet_eta", &precalib_jet_eta);
   //tree->SetBranchAddress("constit_jet_phi", &precalib_jet_phi);
   //tree->SetBranchAddress("constit_jet_e", &precalib_jet_e);
  }

  if (isMC) {
   tree->SetBranchAddress("truth_jet_n", &truth_jet_n);
   tree->SetBranchAddress("truth_jet_pt", &truth_jet_pt);
   tree->SetBranchAddress("truth_jet_eta", &truth_jet_eta);
   tree->SetBranchAddress("truth_jet_phi", &truth_jet_phi);
   tree->SetBranchAddress("truth_jet_e", &truth_jet_e);

   tree->SetBranchAddress("truth_electron_n", &truth_electron_n);
   tree->SetBranchAddress("truth_electron_pt", &truth_electron_pt);
   tree->SetBranchAddress("truth_electron_eta", &truth_electron_eta);
   tree->SetBranchAddress("truth_electron_phi", &truth_electron_phi);
   tree->SetBranchAddress("truth_electron_charge", &truth_electron_charge);

   tree->SetBranchAddress("truth_muon_n", &truth_muon_n);
   tree->SetBranchAddress("truth_muon_pt", &truth_muon_pt);
   tree->SetBranchAddress("truth_muon_eta", &truth_muon_eta);
   tree->SetBranchAddress("truth_muon_phi", &truth_muon_phi);
   tree->SetBranchAddress("truth_muon_charge", &truth_muon_charge);

   tree->SetBranchAddress("truth_photon_n", &truth_photon_n);
   tree->SetBranchAddress("truth_photon_pt", &truth_photon_pt);
   tree->SetBranchAddress("truth_photon_eta", &truth_photon_eta);
   tree->SetBranchAddress("truth_photon_phi", &truth_photon_phi);
  }

  tree->SetBranchAddress("electron_n", &electron_n);
  tree->SetBranchAddress("electron_pt", &electron_pt);
  tree->SetBranchAddress("electron_eta", &electron_eta);
  tree->SetBranchAddress("electron_phi", &electron_phi);
  tree->SetBranchAddress("electron_charge", &electron_charge);
  tree->SetBranchAddress("electron_tight", &electron_tight);
  tree->SetBranchAddress("electron_loose", &electron_loose);
  tree->SetBranchAddress("electron_d0sig", &electron_d0sig);
  tree->SetBranchAddress("electron_delta_z0_sin_theta", &electron_delta_z0_sin_theta); 

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

void TreeVariables::PrintAll (const long long entry) {
  tree->GetEntry(entry);
  cout << endl << "////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "// Getting event " << entry << "..." << endl;
  cout << "////////////////////////////////////////////////////////////////////////////////" << endl;

  if (testTree) {
   cout << endl << "Calibrated jets (EM scale):" << endl;
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << em_xcalib_jet_pt->at(j)
          << setw(10) << em_xcalib_jet_eta->at(j)
          << setw(10) << em_xcalib_jet_phi->at(j)
          << setw(10) << em_xcalib_jet_e->at(j) << endl;
   }

   cout << endl << "Calibrated jets (EtaJES scale):" << endl;
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << em_etajes_jet_pt->at(j)
          << setw(10) << em_etajes_jet_eta->at(j)
          << setw(10) << em_etajes_jet_phi->at(j)
          << setw(10) << em_etajes_jet_e->at(j) << endl;
   }
 
   cout << "Reco jets (EM scale):" << endl; 
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << emscale_jet_pt->at(j)
          << setw(10) << emscale_jet_eta->at(j)
          << setw(10) << emscale_jet_phi->at(j)
          << setw(10) << emscale_jet_e->at(j) << endl;
   }

   cout << "Calibrated jets (Constit scale):" << endl;
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << constit_xcalib_jet_pt->at(j)
          << setw(10) << constit_xcalib_jet_eta->at(j)
          << setw(10) << constit_xcalib_jet_phi->at(j)
          << setw(10) << constit_xcalib_jet_e->at(j) << endl;
   }

   cout << endl << "Calibrated jets (Constit-EtaJES scale):" << endl;
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << constit_etajes_jet_pt->at(j)
          << setw(10) << constit_etajes_jet_eta->at(j)
          << setw(10) << constit_etajes_jet_phi->at(j)
          << setw(10) << constit_etajes_jet_e->at(j) << endl;
   }
 
   cout << "Reco jets (Constit scale):" << endl; 
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << constit_jet_pt->at(j)
          << setw(10) << constit_jet_eta->at(j)
          << setw(10) << constit_jet_phi->at(j)
          << setw(10) << constit_jet_e->at(j) << endl;
   }
  }
  else {
   cout << endl << "Calibrated jets:" << endl;
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << jet_pt->at(j)
          << setw(10) << jet_eta->at(j)
          << setw(10) << jet_phi->at(j)
          << setw(10) << jet_e->at(j) << endl;
   }

   cout << "Reco jets:" << endl; 
   cout << setw(2) << "#"
        << setw(10) << "Pt"
        << setw(10) << "Eta"
        << setw(10) << "Phi"
        << setw(10) << "E" << endl;
   for (int j = 0; j < jet_n; j++) {
     cout << setw(2) << j
          << setw(10) << precalib_jet_pt->at(j)
          << setw(10) << precalib_jet_eta->at(j)
          << setw(10) << precalib_jet_phi->at(j)
          << setw(10) << precalib_jet_e->at(j) << endl;
   }
  }

  cout << "Truth jets:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi"
       << setw(10) << "E" << endl;
  for (int j = 0; j < truth_jet_n; j++) {
    cout << setw(2) << j
         << setw(10) << truth_jet_pt->at(j)
         << setw(10) << truth_jet_eta->at(j)
         << setw(10) << truth_jet_phi->at(j)
         << setw(10) << truth_jet_e->at(j) << endl;
  }

  cout << endl << "Calibrated electrons:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi" << endl;
  for (int j = 0; j < electron_n; j++) {
    cout << setw(2) << j
         << setw(10) << electron_pt->at(j)
         << setw(10) << electron_eta->at(j)
         << setw(10) << electron_phi->at(j) << endl;
  }

  cout << "Truth electrons:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi" << endl;
  for (int j = 0; j < truth_electron_n; j++) {
    cout << setw(2) << j
         << setw(10) << truth_electron_pt->at(j)
         << setw(10) << truth_electron_eta->at(j)
         << setw(10) << truth_electron_phi->at(j) << endl;
  }

  cout << endl << "Calibrated muons:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi" << endl;
  for (int j = 0; j < muon_n; j++) {
    cout << setw(2) << j
         << setw(10) << muon_pt->at(j)
         << setw(10) << muon_eta->at(j)
         << setw(10) << muon_phi->at(j) << endl;
  }

  cout << "Truth muons:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi" << endl;
  for (int j = 0; j < truth_muon_n; j++) {
    cout << setw(2) << j
         << setw(10) << truth_muon_pt->at(j)
         << setw(10) << truth_muon_eta->at(j)
         << setw(10) << truth_muon_phi->at(j) << endl;
  }

  cout << endl << "Calibrated photons:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi" << endl;
  for (int j = 0; j < photon_n; j++) {
    cout << setw(2) << j
         << setw(10) << photon_pt->at(j)
         << setw(10) << photon_eta->at(j)
         << setw(10) << photon_phi->at(j) << endl;
  }

  cout << "Truth photons:" << endl;
  cout << setw(2) << "#"
       << setw(10) << "Pt"
       << setw(10) << "Eta"
       << setw(10) << "Phi" << endl;
  for (int j = 0; j < truth_photon_n; j++) {
    cout << setw(2) << j
         << setw(10) << truth_photon_pt->at(j)
         << setw(10) << truth_photon_eta->at(j)
         << setw(10) << truth_photon_phi->at(j) << endl;
  }


  return;
}
