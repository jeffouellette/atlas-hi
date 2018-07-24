const float pi = TMath::Pi();

void jetEtaPhiCorr () {

  const string dataPath = "/Volumes/My Passport/Research/atlas-hi/pPb_8TeV_2016_dijetAnalysis/data/validation";
  const string fileName = "user.jeouelle.2.4.30hi.jetTreeMaker.valid.101.mc15_8TeV.420012.jetjet.JZ2R04.Both_myOutput.root";

  TFile* inFile = new TFile(Form("%s/%s", dataPath.c_str(), fileName.c_str()), "read");
  TTree* tree = (TTree*)inFile->Get("bush");

  vector<float>* jet_pt = NULL;
  vector<float>* jet_eta = NULL;
  vector<float>* jet_phi = NULL;
  vector<float>* jet_e = NULL;
  int njet = 0;

  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_phi", &jet_phi);
  tree->SetBranchAddress("jet_e", &jet_e);
  tree->SetBranchAddress("njet", &njet);

  TH2F* etaPhiHist;
  {
    etaPhiHist = new TH2F("etaPhiHist", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  }

  const int numentries = tree->GetEntries();
  double jpt, jeta, jphi, je;
  for (long long entry = 0; entry < numentries; entry++) {
    tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

    for (int j = 0; j < njet; j++) {
      jpt = (double)jet_pt->at(j);
      jeta = (double)jet_eta->at(j);
      jphi = (double)jet_phi->at(j);
      je = (double)jet_e->at(j);

      etaPhiHist->Fill(jeta, jphi);
    }
  }

  etaPhiHist->GetZaxis()->SetTitle("Counts");
  etaPhiHist->Draw("colz");

}
