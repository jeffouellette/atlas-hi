#include "../Params.C"
#include "../../Initialization.C"


double deltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {

  double deta = eta1 - eta2;
  double dphi = phi1 - phi2;

  if (dphi < -pi) dphi += 2*pi;
  if (dphi > pi) dphi -= 2*pi;

  return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );

}


double GetXCalibSystematicError(TFile* file, const double jpt, const double jeta) {

  short etabin = 0;
  while (xcalibEtabins[etabin] < TMath::Abs(jeta)) etabin++;
  etabin--;

  const TString hname = TString("fsys_rel_") + Form("%i", etabin);
  TH1D* fsys_rel = (TH1D*)file->Get(hname.Data());

  return (fsys_rel->GetBinContent(fsys_rel->FindBin(jpt)) - 1) * jpt;

}


void ZGammaJetCrossCheck (const int dataSet, // Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
                          const double luminosity = 0, // Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
                          const bool isMC = false, // is data/MC flag.
                          const bool isMCperiodAflag = false, // flag that is raised for MC (meaningless if isMC is false)
                          const string inFileName = "") // Input root file name where tree is stored; if == "" code will try to guess file name based on other info
{

  // Setup trigger vectors
  setupDirectories("", "pPb_8TeV_2016_jet_calibration/");
  if (!isMC) setupTriggers(dataSet);

  const bool isValidationSample = isMC && (TString(inFileName).Contains("valid") || TString(inFileName).Contains("Zee"));
  const string identifier = (isMC ? (string(isMCperiodAflag ? "pPb_":"Pbp_") + (dataSet > 0 ? (string(isValidationSample ? "Valid_":"Overlay_") + "GammaJet_Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet":(string("ZeeJet")+to_string(-dataSet))))) : to_string(dataSet));
  cout << "File Identifier: " << identifier << endl;

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

  const short jetTrigLengthA = 27;
  const short jetTrigLengthB = 19;
  const short electronTrigLength = 3;//5;
  const short muonTrigLength = 4;
  const short photonTrigLength = 7;

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
  //  "HLT_e10_lhloose", // these triggers don't really do anything since Z's only start appearing really at much higher electron pT (~40-50 GeV and up)
  //  "HLT_e15_lhloose",
    "HLT_e20_lhloose",
    "HLT_e22_lhloose",
    "HLT_e24_lhloose"
  };
  const float electronTriggerPtCuts[electronTrigLength] = {/*10, 15,*/ 20, 22, 24};

  const char* muonTriggerNames[muonTrigLength] = {
    "HLT_mu15",
    "HLT_mu18",
    "HLT_mu20",
    "HLT_mu20_L1MU15"
  };
  const float muonTriggerPtCuts[muonTrigLength] = {15, 18, 20, 20};

  const char* photonTriggerNames[photonTrigLength] = {
    "HLT_g10_loose",
    "HLT_g15_loose",
    "HLT_g20_loose",
    "HLT_g25_loose",
    "HLT_g30_loose",
    "HLT_g35_loose",
    "HLT_g60_loose"
  };
  const float photonTriggerPtCuts[photonTrigLength] = {10, 15, 20, 25, 30, 35, 60};

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
    string fileIdentifier;
    if (inFileName == "") {
      if (!isMC) fileIdentifier = to_string(dataSet);
      else fileIdentifier = string(dataSet > 0 ? ("Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet" : (string("ZeeJet")+to_string(-dataSet)))) + string(isMCperiodAflag ? ".pPb":".Pbp");
    } else fileIdentifier = inFileName;

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
          
          if (fname.Contains(fileIdentifier)) {
            file = new TFile(dataPath+fname, "READ");
            tree = (TTree*)file->Get("tree");
            break;
          }
        }
      }
    }
  }
  if (tree == NULL || file == NULL) {
    cout << "Error: In ZGammaJetCrossCheck.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
    return;
  }
  /**** End find TTree ****/

  if (!isMC) tree->SetBranchAddress("eventNumber", &eventNumber);
  tree->SetBranchAddress("nvert", &nvert);
  tree->SetBranchAddress("vert_z", &vert_z);
  tree->SetBranchAddress("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress("vert_type", &vert_type);

  vector<Trigger*> electronTriggers = {};
  vector<Trigger*> muonTriggers = {};
  vector<Trigger*> photonTriggers = {};

  if (!isMC) {
    if (dataSet < 313629) { // ion triggers used up to run 313603
      for (short jetTriggerN = 0; jetTriggerN < jetTrigLengthA; jetTriggerN++) {
        for (Trigger* temp : triggerVec) {
          if (temp->name == jetTriggerNamesA[jetTriggerN]) {
            tree->SetBranchStatus(jetTriggerNamesA[jetTriggerN], 0);
            tree->SetBranchStatus(Form("%s_prescale", jetTriggerNamesA[jetTriggerN]), 0);
            break;
          }
        }
      }
    } else {
      for (short jetTriggerN = 0; jetTriggerN < jetTrigLengthB; jetTriggerN++) {
        for (Trigger* temp : triggerVec) {
          if (temp->name == jetTriggerNamesB[jetTriggerN]) {
            tree->SetBranchStatus(jetTriggerNamesB[jetTriggerN], 0);
            tree->SetBranchStatus(Form("%s_prescale", jetTriggerNamesB[jetTriggerN]), 0);
            break;
          }
        }
      }
    }

    for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
      Trigger* temp = new Trigger(electronTriggerNames[electronTriggerN], electronTriggerPtCuts[electronTriggerN], -2.47, 2.47);
      electronTriggers.push_back(temp);
      tree->SetBranchAddress(electronTriggerNames[electronTriggerN], &(temp->m_trig_bool));
      tree->SetBranchAddress(Form("%s_prescale", electronTriggerNames[electronTriggerN]), &(temp->m_trig_prescale));
    }

    for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
      Trigger* temp = new Trigger(muonTriggerNames[muonTriggerN], muonTriggerPtCuts[muonTriggerN], -2.40, 2.40);
      muonTriggers.push_back(temp);
      tree->SetBranchAddress(muonTriggerNames[muonTriggerN], &(temp->m_trig_bool));
      tree->SetBranchAddress(Form("%s_prescale", muonTriggerNames[muonTriggerN]), &(temp->m_trig_prescale));
    }

    for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
      Trigger* temp = new Trigger(photonTriggerNames[photonTriggerN], photonTriggerPtCuts[photonTriggerN], -2.47, 2.47);
      photonTriggers.push_back(temp);
      tree->SetBranchAddress(photonTriggerNames[photonTriggerN], &(temp->m_trig_bool));
      tree->SetBranchAddress(Form("%s_prescale", photonTriggerNames[photonTriggerN]), &(temp->m_trig_prescale));
    }
  } // end branch triggers
  else {
    tree->SetBranchAddress("EventWeight", &eventWeight);
  }

  tree->SetBranchAddress("jet_n", &jet_n);
  tree->SetBranchAddress("xcalib_etajes_jet_pt", &jet_pt);
  tree->SetBranchAddress("xcalib_etajes_jet_eta", &jet_eta);
  tree->SetBranchAddress("xcalib_etajes_jet_phi", &jet_phi);
  tree->SetBranchAddress("xcalib_etajes_jet_e", &jet_e);

  tree->SetBranchAddress("init_jet_pt", &init_jet_pt);

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

  const long long numEntries = tree->GetEntries();

  TH2F* zeeJetHists[3][numetabins];
  TH2F* zmumuJetHists[3][numetabins];
  TH2F* gJetHists[3][numetabins];
  TH1F* zMassSpectra[2][numetabins+1];
  TH1F* zMassSpectra_AllSigns[2][numetabins+1];

  for (short etabin = 0; etabin < numetabins; etabin++) {
    for (short errType = 0; errType < 3; errType++) {
      string error = "sys_lo";
      if (errType == 1) error = "stat";
      else if (errType == 2) error = "sys_hi";

      zeeJetHists[errType][etabin] = new TH2F(Form("zeeJetPtRatio_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpbins, pbins, numxjrefbins, xjrefbins);
      zeeJetHists[errType][etabin]->Sumw2();

      zmumuJetHists[errType][etabin] = new TH2F(Form("zmumuJetPtRatio_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpbins, pbins, numxjrefbins, xjrefbins);
      zmumuJetHists[errType][etabin]->Sumw2();

      gJetHists[errType][etabin] = new TH2F(Form("gJetPtRatio_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpbins, pbins, numxjrefbins, xjrefbins);
      gJetHists[errType][etabin]->Sumw2();
    }

    for (short spcType = 0; spcType < 2; spcType++) {
      string species = "ee";
      if (spcType == 0) species = "mumu";

      zMassSpectra[spcType][etabin] = new TH1F(Form("z%sMassSpectrum_dataSet%s_%s_etabin%i", species.c_str(), identifier.c_str(), (isMC ? "mc":"data"), etabin), "", 50, 60, 110);
      zMassSpectra[spcType][etabin]->Sumw2();
      zMassSpectra_AllSigns[spcType][etabin] = new TH1F(Form("z%sMassSpectrum_AllSigns_dataSet%s_%s_etabin%i", species.c_str(), identifier.c_str(), (isMC ? "mc":"data"), etabin), "", 50, 60, 110);
      zMassSpectra_AllSigns[spcType][etabin]->Sumw2();

      if (etabin == 0) {
        zMassSpectra[spcType][numetabins] = new TH1F(Form("z%sMassSpectrum_dataSet%s_%s", species.c_str(), identifier.c_str(), (isMC ? "mc":"data")), "", 50, 60, 110);
        zMassSpectra[spcType][numetabins]->Sumw2();
        zMassSpectra_AllSigns[spcType][numetabins] = new TH1F(Form("z%sMassSpectrum_AllSigns_dataSet%s_%s", species.c_str(), identifier.c_str(), (isMC ? "mc":"data")), "", 50, 60, 110);
        zMassSpectra_AllSigns[spcType][numetabins]->Sumw2();
      }
    }
  }

  int Zee_n[numetabins] = {};
  int Zmumu_n[numetabins] = {};
  int g_n[numetabins] = {};

  TFile* xCalibSystematicsFile = new TFile((rootPath + "cc_sys_090816.root").c_str(), "READ");

  // begin loop over events
  for (long long entry = 0; entry < numEntries; entry++) {
    tree->GetEntry(entry);

    // basic event selection: require a primary vertex
    if ((nvert == 0) || (nvert > 0 && vert_type->at(0) != 1)) continue;

    if (photon_n < 1 && muon_n < 2 && electron_n < 2) continue; // reject events with less than 2 muons, 2 electrons AND 1 photon

    // Z + jet -> e- + e+ + jet events
    if (electron_n >= 2) {
      for (int e1 = 0; e1 < electron_n; e1++) {
        if (electron_pt->at(e1) < electron_pt_cut) continue;
        if ((1.37 < TMath::Abs(electron_eta->at(e1)) && TMath::Abs(electron_eta->at(e1)) < 1.52) || 2.47 < TMath::Abs(electron_eta->at(e1))) continue;
        for (int e2 = 0; e2 < e1; e2++) { 
          if (electron_pt->at(e2) < electron_pt_cut) continue;
          if ((1.37 < TMath::Abs(electron_eta->at(e2)) && TMath::Abs(electron_eta->at(e2)) < 1.52) || 2.47 < TMath::Abs(electron_eta->at(e2))) continue;

          TLorentzVector electron1, electron2;
          electron1.SetPtEtaPhiM(electron_pt->at(e1), electron_eta->at(e1), electron_phi->at(e1), electron_mass);
          electron2.SetPtEtaPhiM(electron_pt->at(e2), electron_eta->at(e2), electron_phi->at(e2), electron_mass);
          const int le = (electron_pt->at(e1) > electron_pt->at(e2) ? e1 : e2);
          const float leading_electron_pt = electron_pt->at(le);
          const float leading_electron_eta = electron_eta->at(le);

          float prescale = 1;
          if (!isMC) {
            Trigger* minPrescaleElectronTrigger = NULL;
            for (Trigger* trig : electronTriggers)
              if ( (minPrescaleElectronTrigger == NULL || trig->m_trig_prescale < minPrescaleElectronTrigger->m_trig_prescale)
                    && trig->m_trig_prescale > 0
                    && trig->min_pt <= leading_electron_pt
                    && trig->lower_eta <= leading_electron_eta 
                    && leading_electron_eta <= trig->upper_eta
                 )
                minPrescaleElectronTrigger = trig;
            if (minPrescaleElectronTrigger == NULL || !(minPrescaleElectronTrigger->m_trig_bool) || minPrescaleElectronTrigger->m_trig_prescale <= 0.) continue;
            prescale = minPrescaleElectronTrigger->m_trig_prescale;
          } else prescale = eventWeight;

          if (!(electron_loose->at(e1) && electron_loose->at(e2))) continue; // electron quality criteria
          const TLorentzVector Z = electron1 + electron2;
          const float Z_pt = Z.Pt(); 
          const float Z_eta = Z.Eta();
          const float Z_phi = Z.Phi();
          const float Z_m = Z.M();

          // find the 2 highest pt jets
          int leading_jet = -1;
          int subleading_jet = -1;
          TLorentzVector jet_tlv;
          for (int jet = 0; jet < jet_n; jet++) {
            jet_tlv.SetPtEtaPhiE(jet_pt->at(jet), jet_eta->at(jet), jet_phi->at(jet), jet_e->at(jet));
            if (jet_tlv.DeltaR(electron1) < 0.2 || jet_tlv.DeltaR(electron2) < 0.2) continue; // jet isolation criteria
            if (leading_jet == -1 || jet_pt->at(leading_jet) < jet_pt->at(jet)) {
              subleading_jet = leading_jet;
              leading_jet = jet;
            } else if (subleading_jet == -1 || jet_pt->at(subleading_jet) < jet_pt->at(jet)) {
              subleading_jet = jet;
            }
          }
          if (leading_jet == -1) continue; // true iff no jets are isolated enough from electrons (and there are jets)
          //if (init_jet_pt->at(leading_jet) < 20) continue; // exclude uncorrected jets with pt < 20 GeV

          const float leading_jet_pt = jet_pt->at(leading_jet);
          const float leading_jet_eta = jet_eta->at(leading_jet);
          const float leading_jet_phi = jet_phi->at(leading_jet);

          const float subleading_jet_pt = ((subleading_jet >= 0 && subleading_jet < jet_n) ? jet_pt->at(subleading_jet) : 0);
          const float subleading_jet_phi = ((subleading_jet >= 0 && subleading_jet < jet_n) ? jet_phi->at(subleading_jet) : 0);

          // Calculate opening angle in the transverse plane
          float dPhi = TMath::Abs(leading_jet_phi - Z_phi);
          while (dPhi > pi) dPhi = TMath::Abs(dPhi - 2*pi);
          if (dPhi < 7*pi/8) continue; // cut on Z+jet samples not back-to-back in the transverse plane

          // Determine if the opening angle between the subleading jet and the Z gives a non-vanishing pt ratio
          if (subleading_jet_pt > 12) {
            float subleading_dPhi = TMath::Abs(subleading_jet_phi - Z_phi);
            while (subleading_dPhi > pi) subleading_dPhi = TMath::Abs(subleading_dPhi - 2*pi);
            if (subleading_jet_pt / (Z_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2) continue; // suppress dijets
          }
  
          // Make sure the jet is within the relevant bounds
          if (leading_jet_eta < etabins[0] || leading_jet_eta > etabins[numetabins]) continue;
          // Put the jet in the right eta bin
          int etabin = 0;
          while (etabins[etabin] < leading_jet_eta) etabin++;
          etabin--;

          zMassSpectra_AllSigns[1][etabin]->Fill(Z_m, prescale);
          zMassSpectra_AllSigns[1][numetabins]->Fill(Z_m, prescale);

          if (electron_charge->at(e1) == electron_charge->at(e2)) continue;

          zMassSpectra[1][etabin]->Fill(Z_m, prescale);
          zMassSpectra[1][numetabins]->Fill(Z_m, prescale);
          Zee_n[etabin]++;

          if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m) continue; // cut on our sample Z boson mass

          const float leading_jet_pt_err = GetXCalibSystematicError(xCalibSystematicsFile, leading_jet_pt, leading_jet_eta);
          const float ptref = Z_pt * TMath::Cos(pi - dPhi);
          const float ptratio[3] = {leading_jet_pt-leading_jet_pt_err, leading_jet_pt, leading_jet_pt+leading_jet_pt_err};
    
          for (short errType = 0; errType < 3; errType++) zeeJetHists[errType][etabin]->Fill(Z_pt, ptratio[errType]/ptref, prescale);
        }
      } // end loop over electron pairs
    } // end Z->ee type events

    if (muon_n >= 2) {
      for (int m1 = 0; m1 < muon_n; m1++) {
        if (muon_pt->at(m1) < muon_pt_cut) continue;
        if (2.4 < TMath::Abs(muon_eta->at(m1))) continue;
        for (int m2 = 0; m2 < m1; m2++) {
          if (muon_pt->at(m2) < muon_pt_cut) continue; 
          if (2.4 < TMath::Abs(muon_eta->at(m2))) continue;

          TLorentzVector muon1, muon2;
          muon1.SetPtEtaPhiM(muon_pt->at(m1), muon_eta->at(m1), muon_phi->at(m1), muon_mass);
          muon2.SetPtEtaPhiM(muon_pt->at(m2), muon_eta->at(m2), muon_phi->at(m2), muon_mass);
          const int lm = (muon_pt->at(m1) > muon_pt->at(m2) ? m1 : m2);
          const float leading_muon_pt = muon_pt->at(lm);
          const float leading_muon_eta = muon_eta->at(lm);

          float prescale = 1;
          if (!isMC) {
            Trigger* minPrescaleMuonTrigger = NULL;
            for (Trigger* trig : muonTriggers)
              if ( (minPrescaleMuonTrigger == NULL || trig->m_trig_prescale < minPrescaleMuonTrigger->m_trig_prescale)
                     && trig->m_trig_prescale > 0
                     && trig->min_pt <= leading_muon_pt
                     && trig->lower_eta <= leading_muon_eta 
                     && leading_muon_eta <= trig->upper_eta
                 ) minPrescaleMuonTrigger = trig;
            if (minPrescaleMuonTrigger == NULL || !(minPrescaleMuonTrigger->m_trig_bool)) continue;
            prescale = minPrescaleMuonTrigger->m_trig_prescale;
          } else prescale = eventWeight;
          
          //if (!(muon_tight->at(m1) || muon_tight->at(m2))) continue;
          if (!(muon_loose->at(m1) && muon_loose->at(m2))) continue; // muon quality criteria

          const TLorentzVector Z = muon1 + muon2;
          const float Z_pt = Z.Pt(); 
          const float Z_eta = Z.Eta();
          const float Z_phi = Z.Phi();
          const float Z_m = Z.M();

          // find the 2 highest pt jets
          int leading_jet = -1;
          int subleading_jet = -1;
          TLorentzVector jet_tlv;
          for (int jet = 0; jet < jet_n; jet++) {
            jet_tlv.SetPtEtaPhiE(jet_pt->at(jet), jet_eta->at(jet), jet_phi->at(jet), jet_e->at(jet));
            if (jet_tlv.DeltaR(muon1) < 0.2 || jet_tlv.DeltaR(muon2) < 0.2) continue; // jet isolation criteria
            if (leading_jet == -1 || jet_pt->at(leading_jet) < jet_pt->at(jet)) {
              subleading_jet = leading_jet;
              leading_jet = jet;
            } else if (subleading_jet == -1 || jet_pt->at(subleading_jet) < jet_pt->at(jet)) {
              subleading_jet = jet;
            }
          }
          if (leading_jet == -1) continue; // true iff no jet is isolated enough from the muons
          //if (init_jet_pt->at(leading_jet) < 20) continue; // exclude uncorrected jets with pt < 20 GeV

          const float leading_jet_pt = jet_pt->at(leading_jet);
          const float leading_jet_eta = jet_eta->at(leading_jet);
          const float leading_jet_phi = jet_phi->at(leading_jet);

          const float subleading_jet_pt = ((subleading_jet >= 0 && subleading_jet < jet_n) ? jet_pt->at(subleading_jet) : 0);
          const float subleading_jet_phi = ((subleading_jet >= 0 && subleading_jet < jet_n) ? jet_phi->at(subleading_jet) : 0);

          // Calculate opening angle in the transverse plane
          float dPhi = TMath::Abs(leading_jet_phi - Z_phi);
          while (dPhi > pi) dPhi = TMath::Abs(dPhi - 2*pi);
          if (dPhi < 7*pi/8) continue; // cut on Z+jet samples not back-to-back in the transverse plane

          // Determine if the opening angle between the subleading jet and the Z gives a non-vanishing pt ratio
          if (subleading_jet_pt > 12) {
            float subleading_dPhi = TMath::Abs(subleading_jet_phi - Z_phi);
            while (subleading_dPhi > pi) subleading_dPhi = TMath::Abs(subleading_dPhi - 2*pi);
            if (subleading_jet_pt / (Z_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2) continue; // suppress dijets
          }

          // Make sure the jet is within the relevant bounds
          if (leading_jet_eta < etabins[0] || leading_jet_eta > etabins[numetabins]) continue;
          // Put the jet in the right eta bin
          short etabin = 0;
          while (etabins[etabin] < leading_jet_eta) etabin++;
          etabin--;

          zMassSpectra_AllSigns[0][etabin]->Fill(Z_m, prescale);
          zMassSpectra_AllSigns[0][numetabins]->Fill(Z_m, prescale);

          if (muon_charge->at(m1) == muon_charge->at(m2)) continue;

          zMassSpectra[0][etabin]->Fill(Z_m, prescale);
          zMassSpectra[0][numetabins]->Fill(Z_m, prescale);
          Zmumu_n[etabin]++;

          if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m) continue; // cut on our sample Z boson mass

          const float leading_jet_pt_err = GetXCalibSystematicError(xCalibSystematicsFile, leading_jet_pt, leading_jet_eta);
          const float ptref = Z_pt * TMath::Cos(pi - dPhi);
          const float ptratio[3] = {leading_jet_pt-leading_jet_pt_err, leading_jet_pt, leading_jet_pt+leading_jet_pt_err};

          for (short errType = 0; errType < 3; errType++) zmumuJetHists[errType][etabin]->Fill(Z_pt, ptratio[errType]/ptref, prescale);
        }
      } // end loop over muon pairs
    } // end Z->mumu type events

    if (photon_n >= 1) {
      for (int p = 0; p < photon_n; p++) {

        TLorentzVector photon;
        const float this_photon_pt = photon_pt->at(p);
        const float this_photon_eta = photon_eta->at(p);
        const float this_photon_phi = photon_phi->at(p);
        photon.SetPtEtaPhiM(this_photon_pt, this_photon_eta, this_photon_phi, 0);

        if ((1.37 < TMath::Abs(this_photon_eta) && TMath::Abs(this_photon_eta) < 1.52) || 2.37 < TMath::Abs(this_photon_eta)) continue;

        float prescale = 1;
        if (!isMC) {
          Trigger* minPrescalePhotonTrigger = NULL;
          for (Trigger* trig : photonTriggers)
            if ( (minPrescalePhotonTrigger == NULL || trig->m_trig_prescale < minPrescalePhotonTrigger->m_trig_prescale)
                   && trig->m_trig_prescale > 0
                   && trig->min_pt <= this_photon_pt
                   && trig->lower_eta <= this_photon_eta 
                   && this_photon_eta <= trig->upper_eta
               ) minPrescalePhotonTrigger = trig;
          if (minPrescalePhotonTrigger == NULL || !(minPrescalePhotonTrigger->m_trig_bool)) continue;
          prescale = minPrescalePhotonTrigger->m_trig_prescale;
        } else prescale = eventWeight;

        if (!photon_tight->at(p)) continue; // require tight cuts on photons

        if (photon_topoetcone40->at(p) > isolationEnergyCut) continue; // require maximum isolation energy on gammas

        int leading_jet = -1;
        int subleading_jet = -1;
        TLorentzVector jet_tlv;
        for (int jet = 0; jet < jet_n; jet++) {
          jet_tlv.SetPtEtaPhiE(jet_pt->at(jet), jet_eta->at(jet), jet_phi->at(jet), jet_e->at(jet));
          if (jet_tlv.DeltaR(photon) < 0.6) continue;

          if (leading_jet == -1 || jet_pt->at(leading_jet) < jet_pt->at(jet)) {
            subleading_jet = leading_jet;
            leading_jet = jet;
          } else if (subleading_jet == -1 || jet_pt->at(subleading_jet) < jet_pt->at(jet)) {
            subleading_jet = jet;
          }
        }

        if (leading_jet == -1) continue; // true iff there are no jets
        //if (init_jet_pt->at(leading_jet) < 20) continue; // exclude uncorrected jets with pt < 20 GeV

        const float leading_jet_pt = jet_pt->at(leading_jet);
        const float leading_jet_eta = jet_eta->at(leading_jet);
        const float leading_jet_phi = jet_phi->at(leading_jet);

        const float subleading_jet_pt = ((subleading_jet >= 0 && subleading_jet < jet_n) ? jet_pt->at(subleading_jet) : 0);
        const float subleading_jet_phi = ((subleading_jet >= 0 && subleading_jet < jet_n) ? jet_phi->at(subleading_jet) : 0);

        if (leading_jet_eta < etabins[0] || leading_jet_eta > etabins[numetabins]) continue;

        // Calculate opening angle in the transverse plane
        float dPhi = TMath::Abs(leading_jet_phi - this_photon_phi);
        while (dPhi > pi) dPhi = TMath::Abs(dPhi - 2*pi);
        if (dPhi < 7*pi/8) continue; // cut on gamma+jet samples not back-to-back in the transverse plane

        // Determine if the opening angle between the subleading jet and the Z gives a non-vanishing pt ratio
        if (subleading_jet_pt > 12) {
          float subleading_dPhi = TMath::Abs(subleading_jet_phi - this_photon_phi);
          while (subleading_dPhi > pi) subleading_dPhi = TMath::Abs(subleading_dPhi - 2*pi);
          if (subleading_jet_pt / (this_photon_pt * TMath::Cos(pi - subleading_dPhi)) > 0.3) continue; // suppress dijets
        }

        // Make sure the jet is within the relevant bounds
        if (leading_jet_eta < etabins[0] || leading_jet_eta > etabins[numetabins]) continue;
        // Put the jet in the right eta bin
        short etabin = 0;
        while (etabins[etabin] < leading_jet_eta) etabin++;
        etabin--;

        const float leading_jet_pt_err = GetXCalibSystematicError(xCalibSystematicsFile, leading_jet_pt, leading_jet_eta);
        const float ptref = this_photon_pt * TMath::Cos(pi - dPhi);
        const float ptratio[3] = {leading_jet_pt-leading_jet_pt_err, leading_jet_pt, leading_jet_pt+leading_jet_pt_err};

        for (short errType = 0; errType < 3; errType++) gJetHists[errType][etabin]->Fill(this_photon_pt, ptratio[errType]/ptref, prescale);
        g_n[etabin]++;
      }
    } 

  }

  // close root file with systematics
  xCalibSystematicsFile->Close();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  /** End event loop **/

  TFile* outfile = new TFile(Form("%sdataSet_%s.root", rootPath.c_str(), identifier.c_str()), "RECREATE");

  for (short etabin = 0; etabin < numetabins; etabin++) {
    for (short errType = 0; errType < 3; errType++) {
      zmumuJetHists[errType][etabin]->Write();
      if (zmumuJetHists[errType][etabin]) delete zmumuJetHists[errType][etabin];
      zeeJetHists[errType][etabin]->Write();
      if (zeeJetHists[errType][etabin]) delete zeeJetHists[errType][etabin];
      gJetHists[errType][etabin]->Write();
      if (gJetHists[errType][etabin]) delete gJetHists[errType][etabin];
    }
    for (short spcType = 0; spcType < 2; spcType++) {
      zMassSpectra[spcType][etabin]->Write();
      if (zMassSpectra[spcType][etabin]) delete zMassSpectra[spcType][etabin];
      zMassSpectra_AllSigns[spcType][etabin]->Write();
      if (zMassSpectra_AllSigns[spcType][etabin]) delete zMassSpectra_AllSigns[spcType][etabin];
      if (etabin == 0) {
        zMassSpectra[spcType][numetabins]->Write();
        if (zMassSpectra[spcType][numetabins]) delete zMassSpectra[spcType][numetabins];
        zMassSpectra_AllSigns[spcType][numetabins]->Write();
        if (zMassSpectra_AllSigns[spcType][numetabins]) delete zMassSpectra_AllSigns[spcType][numetabins];
      }
    }
  }

  TVectorD infoVec(2+3*numetabins);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  for (short etabin = 0; etabin < numetabins; etabin++) infoVec[2+etabin] = (double)Zee_n[etabin];
  for (short etabin = 0; etabin < numetabins; etabin++) infoVec[2+numetabins+etabin] = (double)Zmumu_n[etabin];
  for (short etabin = 0; etabin < numetabins; etabin++) infoVec[2+2*numetabins+etabin] = (double)g_n[etabin];
  infoVec.Write(Form("infoVec_%s", identifier.c_str()));
  outfile->Close();
  if (outfile) delete outfile;
  return;
}
