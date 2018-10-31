#include <vector>
#include <algorithm>


void eventMismatchStudy() {
  TFile* xAODFile = new TFile("xAOD.root", "read");
  TFile* DAODFile = new TFile("DAOD.root", "read");

  TTree* xAODBush = (TTree*)xAODFile->Get("bush");
  TTree* DAODBush = (TTree*)DAODFile->Get("bush");

  std::vector<int> xAODEvents(0);
  std::vector<int> DAODEvents(0);

  const int xAODNumEntries = xAODBush->GetEntries();
  const int DAODNumEntries = DAODBush->GetEntries();

  int xAODEventNum = 0;
  int DAODEventNum = 0;

  std::vector<float>* xaod_jpt = NULL;
  std::vector<float>* xaod_jeta = NULL;
  int xaod_njet = 0;
  std::vector<float>* daod_jpt = NULL;
  std::vector<float>* daod_jeta = NULL;
  int daod_njet = 0;
  TH1F* xaod_jpt_spectrum = new TH1F("xaod_jpt_spectrum", "", 16, 20, 100);
  TH1F* daod_jpt_spectrum = new TH1F("daod_jpt_spectrum", "", 16, 20, 100);

  const double etabins[9] = {-4.9, -4.5, -3.2, -2.0, 0.0, 2.0, 3.2, 4.5, 4.9};
  const int numetabins = 8;
  const double pbins[13] = {20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
  const int numpbins = 12;

  TH2F* xaod_jeta_spectrum = new TH2F("xaod_jet_spectrum", "", numpbins, pbins, numetabins, etabins);
  TH2F* daod_jeta_spectrum = new TH2F("daod_jet_spectrum", "", numpbins, pbins, numetabins, etabins);

  xAODBush->SetBranchAddress("eventNumber", &xAODEventNum);
  xAODBush->SetBranchAddress("jet_pt", &xaod_jpt);
  xAODBush->SetBranchAddress("jet_eta", &xaod_jeta);
  xAODBush->SetBranchAddress("njet", &xaod_njet);
  DAODBush->SetBranchAddress("eventNumber", &DAODEventNum);
  DAODBush->SetBranchAddress("jet_pt", &daod_jpt);
  DAODBush->SetBranchAddress("jet_eta", &daod_jeta);
  DAODBush->SetBranchAddress("njet", &daod_njet);

  for (int entry = 0; entry < xAODNumEntries; entry++) {
   xAODBush->GetEntry(entry);
   xAODEvents.push_back(xAODEventNum);
   for (int jet = 0; jet < xaod_njet; jet++) {
    xaod_jpt_spectrum->Fill(xaod_jpt->at(jet));
    xaod_jeta_spectrum->Fill(xaod_jpt->at(jet), xaod_jeta->at(jet));
   }
  }

  for (int entry = 0; entry < DAODNumEntries; entry++) {
   DAODBush->GetEntry(entry);
   DAODEvents.push_back(DAODEventNum);
   for (int jet = 0; jet < daod_njet; jet++) {
    daod_jpt_spectrum->Fill(daod_jpt->at(jet));
    daod_jeta_spectrum->Fill(daod_jpt->at(jet), daod_jeta->at(jet));
   }
  }

  std::vector<int> mismatchedEvents(0);
  for (int entry = 0; entry < xAODNumEntries; entry++) {
   if (std::find(DAODEvents.begin(), DAODEvents.end(), xAODEvents.at(entry)) == DAODEvents.end()) {
    mismatchedEvents.push_back(entry);
   }
  }

  for (int eventNum : mismatchedEvents) {
   cout << eventNum << endl;
  }

}
