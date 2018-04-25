#include "Params.C"
#include "../../Initialization.C"

static const float muon_mass = 0.105658; // GeV

void ZGammaJetCrossCheck (const int dataSet, // Data set identifier. This should be a run number.
                          const double luminosity) { // Integrated luminosity for this run. Presumed constant over the run period.

  workPath = "/Users/jeffouellette/Research/atlas-hi/pPb_8TeV_2016_jet_calibration/";
  dataPath = "/Users/jeffouellette/Research/atlas-hi/pPb_8TeV_2016_jet_calibration/data/";
  // Setup trigger vectors
  initialize(dataSet, false);

  const double etacuts[7] = {-1.3, -0.8, -0.3, 0, 0.3, 0.8, 1.3};

  int eventNumber;
  //int ntrk;
  int nvert;
  vector<float>* vert_x;
  vector<float>* vert_y;
  vector<float>* vert_z;
  vector<int>* vert_ntrk;
  vector<int>* vert_type;
 
  int jet_n; 
  vector<float>* jet_pt;
  vector<float>* jet_eta;
  vector<float>* jet_phi;
  vector<float>* jet_e;
  
  int muon_n;
  vector<float>* muon_pt;
  vector<float>* muon_eta;
  vector<float>* muon_phi;
  vector<int>* muon_quality;
  vector<int>* muon_charge;
  vector<bool>* muon_tight;
  vector<bool>* muon_loose;

  int electron_n;
  vector<float>* electron_pt;
  vector<float>* electron_eta;
  vector<float>* electron_phi;
  vector<int>* electron_charge;
  vector<bool>* electron_loose;
  vector<bool>* electron_tight;

  int photon_n;
  vector<float>* photon_pt;
  vector<float>* photon_eta;
  vector<float>* photon_phi;
  vector<bool>* photon_tight;
  vector<unsigned int>* photon_isem;
  vector<int>* photon_convFlag;
  vector<float>* photon_Rconv;

  const int jetTrigLengthA = 27;
  const int jetTrigLengthB = 19;
  const int electronTrigLength = 5;
  const int muonTrigLength = 3;
  const int photonTrigLength = 6;

  const char* jetTriggerNamesA[jetTrigLengthA] = {
    "HLT_j15_ion_n320eta490_L1MBTS_1_1",
    "HLT_j25_ion_n320eta490_L1TE5",
    "HLT_j30_ion_0eta490_L1TE10",
    "HLT_j35_ion_n320eta490_L1TE10",
    "HLT_j40_ion_L1J5",
    "HLT_j50_ion_L1J10",
    "HLT_j60_ion_L1J20",
    "HLT_j75_ion_L1J20",
    "HLT_j90_ion_L1J20",
    "HLT_j100_ion_L1J20",
    "HLT_j45_ion_n200eta320",
    "HLT_j55_ion_n200eta320",
    "HLT_j65_ion_n200eta320",
    "HLT_j75_ion_n200eta320",
    "HLT_j45_ion_p200eta320",
    "HLT_j55_ion_p200eta320",
    "HLT_j65_ion_p200eta320",
    "HLT_j75_ion_p200eta320",
    "HLT_j15_ion_p320eta490_L1MBTS_1_1",
    "HLT_j25_ion_p320eta490_L1TE5",
    "HLT_j35_ion_p320eta490_L1TE10",
    "HLT_j45_ion_p320eta490",
    "HLT_j55_ion_p320eta490",
    "HLT_j65_ion_p320eta490",
    "HLT_j45_ion_n320eta490",
    "HLT_j55_ion_n320eta490",
    "HLT_j65_ion_n320eta490"
  };

  const char* jetTriggerNamesB[jetTrigLengthB] = {
    "HLT_j30_0eta490_L1TE10",
    "HLT_j30_L1J5",
    "HLT_j40_L1J5",
    "HLT_j50_L1J10",
    "HLT_j60",
    "HLT_j75_L1J20",
    "HLT_j90_L1J20",
    "HLT_j100_L1J20",
    "HLT_j45_p200eta320",
    "HLT_j55_p200eta320",
    "HLT_j65_p200eta320",
    "HLT_j75_p200eta320",
    "HLT_j15_p320eta490_L1MBTS_1_1",
    "HLT_j25_p320eta490_L1TE5",
    "HLT_j35_p320eta490_L1TE10",
    "HLT_j45_p320eta490",
    "HLT_j55_p320eta490",
    "HLT_j65_p320eta490",
    "HLT_j75_p320eta490"
  };

  const char* electronTriggerNames[electronTrigLength] = {
    "HLT_e10_lhloose_L1EM7",
    "HLT_e15_lhloose_L1EM7",
    "HLT_e15_lhloose_L1EM12",
    "HLT_e15_lhmedium_L1EM12",
    "HLT_e20_lhloose_L1EM15"
  };
  const float electronTriggerPtCuts[electronTrigLength] = {10, 15, 15, 15, 20};

  const char* muonTriggerNames[muonTrigLength] = {
    "HLT_mu15",
    "HLT_mu15_L1MU10",
    "HLT_mu15_L1MU6"
  };
  const float muonTriggerPtCuts[muonTrigLength] = {15, 15, 15};

  const char* photonTriggerNames[photonTrigLength] = {
    "HLT_g10_loose",
    "HLT_g15_loose",
    "HLT_g20_loose",
    "HLT_g25_loose",
    "HLT_g30_loose",
    "HLT_g35_loose"
  };
  const float photonTriggerPtCuts[photonTrigLength] = {10, 15, 20, 25, 30, 35};


  /**** Find the relevant TTree for this run ****/
  TTree* tree = NULL;
  {
    TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
    TList* sysfiles = dir.GetListOfFiles();
    if (sysfiles) {
      TSystemFile* sysfile;
      TString fname;
      TIter next(sysfiles);

      while ((sysfile = (TSystemFile*)next())) {
        fname = sysfile->GetName();
        if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
          if (debugStatements) cout << "Status: In DijetAnalysis.C (breakpoint B): Found " << fname.Data() << endl;
          if (fname.Contains(to_string(dataSet))) {
            tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
            break;
          }
        }
      }
    }
  }
  if (tree == NULL) {
      cout << "Error: In ZGammaJetCrossCheck.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
      return;
  }
  /**** End find TTree ****/

  tree->SetBranchAddress("eventNumber", &eventNumber);
  //tree->SetBranchAddress("ntrk", &ntrk);
  tree->SetBranchAddress("nvert", &nvert);
  tree->SetBranchAddress("vert_x", &vert_x);
  tree->SetBranchAddress("vert_y", &vert_y);
  tree->SetBranchAddress("vert_z", &vert_z);
  tree->SetBranchAddress("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress("vert_type", &vert_type);

  vector<Trigger*> jetTriggers = {};
  vector<Trigger*> electronTriggers = {};
  vector<Trigger*> muonTriggers = {};
  vector<Trigger*> photonTriggers = {};


  if (dataSet < 313629) { // ion triggers used up to run 313603
    for (int jetTriggerN = 0; jetTriggerN < jetTrigLengthA; jetTriggerN++) {
      for (Trigger* temp : triggerVec) {
        if (temp->name == jetTriggerNamesA[jetTriggerN]) {
          jetTriggers.push_back(temp);
          tree->SetBranchAddress(jetTriggerNamesA[jetTriggerN], &(temp->m_trig_bool));
          tree->SetBranchAddress(Form("%s_prescale", jetTriggerNamesA[jetTriggerN]), &(temp->m_trig_prescale));
          break;
        } else { continue; }
      }
    }
  } else {
    for (int jetTriggerN = 0; jetTriggerN < jetTrigLengthB; jetTriggerN++) {
      for (Trigger* temp : triggerVec) {
        if (temp->name == jetTriggerNamesB[jetTriggerN]) {
          jetTriggers.push_back(temp);
          tree->SetBranchAddress(jetTriggerNamesB[jetTriggerN], &(temp->m_trig_bool));
          tree->SetBranchAddress(Form("%s_prescale", jetTriggerNamesB[jetTriggerN]), &(temp->m_trig_prescale));
          break;
        } else { continue; }
      }
    }
  }

  for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* temp = new Trigger(electronTriggerNames[electronTriggerN], electronTriggerPtCuts[electronTriggerN], -2.47, 2.47);
    electronTriggers.push_back(temp);
  }

  for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* temp = new Trigger(muonTriggerNames[muonTriggerN], muonTriggerPtCuts[muonTriggerN], -2.47, 2.47);
    muonTriggers.push_back(temp);
  }

  for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* temp = new Trigger(photonTriggerNames[photonTriggerN], photonTriggerPtCuts[photonTriggerN], -2.47, 2.47);
    photonTriggers.push_back(temp);
  }

  tree->SetBranchAddress("jet_n", &jet_n);
  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_phi", &jet_phi);
  tree->SetBranchAddress("jet_e", &jet_e);

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

  const long long numEvents = tree->GetEntries();

  const int nhists = 6;
  TH2F* zjethists[nhists];
  for (int nhist = 0; nhist < nhists; nhist++) {
    zjethists[nhist] = new TH2F(Form("Zmumu%i", nhist), "", numpbins, pbins, 120, 0, 1.2);
    zjethists[nhist]->Sumw2();
  }

  // begin loop over events
  for (long long eventNum = 0; eventNum < numEvents; eventNum++) {
    tree->GetEntry(eventNum);

    // basic event selection.
    if ((nvert == 0) || (nvert > 0 && vert_type->at(0) != 1)) continue;

    if (muon_n > 2) continue; // reject events with more than 2 muons

    int leading_jet = -1;
    for (int jet = 0; jet < jet_n; jet++) {
      if (leading_jet == -1 || jet_pt->at(leading_jet) < jet_pt->at(jet)) leading_jet = jet;
    }
    float leading_jet_pt = jet_pt->at(leading_jet);
    float leading_jet_eta = jet_eta->at(leading_jet);
    float leading_jet_phi = jet_phi->at(leading_jet);
    float leading_jet_e = jet_e->at(leading_jet);

    if (leading_jet_eta < etacuts[0] || leading_jet_eta > etacuts[6]) continue;

    int etabin = 0;
    while (etacuts[etabin] < leading_jet_eta) etabin++;
    etabin--;
    
    // Z + jet -> mu- + mu+ + jet events
    if (muon_n == 2 && electron_n == 0 && photon_n == 0) {
      int leading_muon = -1;
      int subleading_muon = -1;
      for (int muon = 0; muon < muon_n; muon++) {
        if (leading_muon == -1 || muon_pt->at(leading_muon) < muon_pt->at(muon)) {
          subleading_muon = leading_muon;
          leading_muon = muon;
        } else if (subleading_muon == -1 || muon_pt->at(subleading_muon) < muon_pt->at(muon)) {
          subleading_muon = muon;
        }
      }
      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiM(muon_pt->at(leading_muon), muon_eta->at(leading_muon), muon_phi->at(leading_muon), muon_mass);
      mu2.SetPtEtaPhiM(muon_pt->at(subleading_muon), muon_eta->at(subleading_muon), muon_phi->at(subleading_muon), muon_mass);
      if (muon_charge[leading_muon] == muon_charge[subleading_muon]) continue; // require oppositely charged muons
      TLorentzVector Z = mu1 + mu2;
      float Z_pt = Z.Pt(); 
      float Z_eta = Z.Eta();
      float Z_phi = Z.Phi();
      float Z_m = Z.M();

      float ptratio = leading_jet_pt / Z_pt;
  
      zjethists[etabin]->Fill(leading_jet_pt, ptratio);
    } 
  }
  /** End event loop **/

  TFile* outfile = new TFile(Form("./rootFiles/run_%i.root", dataSet), "RECREATE");
  for (int nhist = 0; nhist < nhists; nhist++) {
    zjethists[nhist]->Write();
    delete zjethists[nhist];
  }
  TVectorD infoVec(2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write("infoVec");
  outfile->Close();
  delete outfile;
  return;
}
