#include "ZGammaJetCrossCheck.h"

TFile* xCalibSystematicsFile = NULL;
TFile* dataOverMCFile = NULL;

vector<Trigger*> electronTriggers = {};
vector<Trigger*> muonTriggers = {};
vector<Trigger*> photonTriggers = {};

bool isPeriodA = false;

double deltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {

 double deta = eta1 - eta2;
 double dphi = phi1 - phi2;

 if (dphi < -pi) dphi += 2*pi;
 if (dphi > pi) dphi -= 2*pi;

 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );

}


double GetXCalibSystematicError(const double jpt, const double jeta) {
  TFile* file = xCalibSystematicsFile;
  if (!file || !file->IsOpen()) return 0;

  short etabin = 0;
  while (xcalibEtabins[etabin] < TMath::Abs(jeta)) etabin++;
  etabin--;

  const TString hname = TString("fsys_rel_") + Form("%i", etabin);
  TH1D* fsys_rel = (TH1D*)file->Get(hname.Data());

  return TMath::Abs(fsys_rel->GetBinContent(fsys_rel->FindBin(jpt)) - 1) * jpt;
}


double GetNewXCalibSystematicError(const double jeta, const double refpt) {
  TFile* file = dataOverMCFile;
  if (!file || !file->IsOpen()) return 0;

  short bin = 0;
  while (bin <= numetabins && etabins[bin] < jeta) bin++;
  bin--;

  const char* period = (isPeriodA ?  "periodA" : "periodB");
  const TString hname = TString(Form("gJetPtRatio_diff%i_stat_%s", bin, period));
  TH1D* hist = (TH1D*)file->Get(hname.Data());

  return hist->GetBinContent(hist->FindBin(refpt)) * refpt;
}


void ZGammaJetCrossCheck (const int dataSet,
                          const double luminosity = 0,
                          const bool isMC, 
                          const bool isMCperiodAflag,
                          const string inFileName)
{

  setupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  if (!isMC) isPeriodA = dataSet < 313500;
  else isPeriodA = isMCperiodAflag;

  const bool isValidationSample = isMC && (TString(inFileName).Contains("valid") || TString(inFileName).Contains("Zee"));
  const string identifier = (isMC ? (string(isMCperiodAflag ? "pPb_":"Pbp_") + (dataSet > 0 ? (string(isValidationSample ? "Valid_":"Overlay_") + "GammaJet_Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet":(string("ZeeJet")+to_string(-dataSet))))) : to_string(dataSet));
  cout << "File Identifier: " << identifier << endl;

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
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
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
  if (tree == NULL || file == NULL) {
   cout << "Error: In ZGammaJetCrossCheck.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables(tree);
  t->SetBranchAddresses(isMC);

  //if (!isMC) tree->SetBranchAddress("eventNumber", &eventNumber);
  //tree->SetBranchAddress("nvert", &nvert);
  //tree->SetBranchAddress("vert_z", &vert_z);
  //tree->SetBranchAddress("vert_ntrk", &vert_ntrk);
  //tree->SetBranchAddress("vert_type", &vert_type);

  if (!isMC) {
   for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* temp = new Trigger(electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
    temp->min_pt = electronTriggerMinPtCuts[electronTriggerN];
    temp->max_pt = electronTriggerMaxPtCuts[electronTriggerN];
    electronTriggers.push_back(temp);
    tree->SetBranchAddress(electronTriggerNames[electronTriggerN], &(temp->m_trig_bool));
    tree->SetBranchAddress(Form("%s_prescale", electronTriggerNames[electronTriggerN]), &(temp->m_trig_prescale));
   }

   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* temp = new Trigger(muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    temp->min_pt = muonTriggerMinPtCuts[muonTriggerN];
    temp->max_pt = muonTriggerMaxPtCuts[muonTriggerN];
    muonTriggers.push_back(temp);
    tree->SetBranchAddress(muonTriggerNames[muonTriggerN], &(temp->m_trig_bool));
    tree->SetBranchAddress(Form("%s_prescale", muonTriggerNames[muonTriggerN]), &(temp->m_trig_prescale));
   }

   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* temp = new Trigger(photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    temp->min_pt = photonTriggerMinPtCuts[photonTriggerN];
    temp->max_pt = photonTriggerMaxPtCuts[photonTriggerN];
    photonTriggers.push_back(temp);
    tree->SetBranchAddress(photonTriggerNames[photonTriggerN], &(temp->m_trig_bool));
    tree->SetBranchAddress(Form("%s_prescale", photonTriggerNames[photonTriggerN]), &(temp->m_trig_prescale));
   }
  } // end branch triggers
  //else {
  // tree->SetBranchAddress("EventWeight", &eventWeight);
  //}

  //tree->SetBranchAddress("jet_n", &jet_n);
  //tree->SetBranchAddress("xcalib_etajes_jet_pt", &jet_pt);
  //tree->SetBranchAddress("xcalib_etajes_jet_eta", &jet_eta);
  //tree->SetBranchAddress("xcalib_etajes_jet_phi", &jet_phi);
  //tree->SetBranchAddress("xcalib_etajes_jet_e", &jet_e);

  //tree->SetBranchAddress("electron_n", &electron_n);
  //tree->SetBranchAddress("electron_pt", &electron_pt);
  //tree->SetBranchAddress("electron_eta", &electron_eta);
  //tree->SetBranchAddress("electron_phi", &electron_phi);
  //tree->SetBranchAddress("electron_charge", &electron_charge);
  //tree->SetBranchAddress("electron_tight", &electron_tight);
  //tree->SetBranchAddress("electron_loose", &electron_loose);

  //tree->SetBranchAddress("muon_n", &muon_n);
  //tree->SetBranchAddress("muon_pt", &muon_pt);
  //tree->SetBranchAddress("muon_eta", &muon_eta);
  //tree->SetBranchAddress("muon_phi", &muon_phi);
  //tree->SetBranchAddress("muon_charge", &muon_charge);
  //tree->SetBranchAddress("muon_quality", &muon_quality);
  //tree->SetBranchAddress("muon_tight", &muon_tight);
  //tree->SetBranchAddress("muon_loose", &muon_loose);

  //tree->SetBranchAddress("photon_n", &photon_n);
  //tree->SetBranchAddress("photon_pt", &photon_pt);
  //tree->SetBranchAddress("photon_eta", &photon_eta);
  //tree->SetBranchAddress("photon_phi", &photon_phi);
  //tree->SetBranchAddress("photon_tight", &photon_tight);
  //tree->SetBranchAddress("photon_isem", &photon_isem);
  //tree->SetBranchAddress("photon_convFlag", &photon_convFlag);
  //tree->SetBranchAddress("photon_Rconv", &photon_Rconv);
  //tree->SetBranchAddress("photon_topoetcone40", &photon_topoetcone40);

  // initialize histograms
  TH2F* zeeJetHists[3][numetabins];
  //TH2F* zeeJetHistsSys[3][numetabins];
  TH2F* zmumuJetHists[3][numetabins];
  //TH2F* zmumuJetHistsSys[3][numetabins];
  TH2F* gJetHists[3][numetabins];
  TH2F* gJetHistsSys[3][numetabins];
  TH1F* zMassSpectra[2][numetabins+1];
  TH1F* zMassSpectra_AllSigns[2][numetabins+1];
  for (short etabin = 0; etabin < numetabins; etabin++) {
   for (short errType = 0; errType < 3; errType++) {
    string error = "sys_lo";
    if (errType == 1) error = "stat";
    else if (errType == 2) error = "sys_hi";

    zeeJetHists[errType][etabin] = new TH2F(Form("zeeJetPtRatio_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    zeeJetHists[errType][etabin]->Sumw2();
    //zeeJetHistsSys[errType][etabin] = new TH2F(Form("zeeJetPtRatioSys_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zeeJetHistsSys[errType][etabin]->Sumw2();

    zmumuJetHists[errType][etabin] = new TH2F(Form("zmumuJetPtRatio_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    zmumuJetHists[errType][etabin]->Sumw2();
    //zmumuJetHistsSys[errType][etabin] = new TH2F(Form("zmumuJetPtRatioSys_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zmumuJetHistsSys[errType][etabin]->Sumw2();

    gJetHists[errType][etabin] = new TH2F(Form("gJetPtRatio_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpgammabins, pgammabins, numxjrefbins, xjrefbins);
    gJetHists[errType][etabin]->Sumw2();
    if (errType == 1) {
     gJetHistsSys[errType][etabin] = new TH2F(Form("gJetPtRatioSys_dataSet%s_hist%i_%s_%s", identifier.c_str(), etabin, (isMC ? "mc":"data"), error.c_str()), "", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
     gJetHistsSys[errType][etabin]->Sumw2();
    }
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

  xCalibSystematicsFile = new TFile((rootPath + "cc_sys_090816.root").c_str(), "READ");
  dataOverMCFile = new TFile((rootPath + "cc_difference.root").c_str(), "READ");

  const long long numEntries = tree->GetEntries();

  // begin loop over events
  for (long long entry = 0; entry < numEntries; entry++) {
   tree->GetEntry(entry);

   // basic event selection: require a primary vertex
   if ((t->nvert == 0) || (t->nvert > 0 && t->vert_type->at(0) != 1)) continue;

   // reject events with less than 2 muons, 2 electrons AND 1 photon
   if (t->photon_n < 1 && t->muon_n < 2 && t->electron_n < 2) continue;

   // Z + jet -> e- + e+ + jet events
   for (int e1 = 0; e1 < t->electron_n; e1++) {
    if (t->electron_pt->at(e1) < electron_pt_cut) continue;
    if ((1.37 < TMath::Abs(t->electron_eta->at(e1)) &&
         TMath::Abs(t->electron_eta->at(e1)) < 1.52) ||
        2.47 < TMath::Abs(t->electron_eta->at(e1)))
     continue;
    for (int e2 = 0; e2 < e1; e2++) { 
     if (t->electron_pt->at(e2) < electron_pt_cut) continue;
     if ((1.37 < TMath::Abs(t->electron_eta->at(e2)) &&
          TMath::Abs(t->electron_eta->at(e2)) < 1.52) ||
         2.47 < TMath::Abs(t->electron_eta->at(e2))) continue;

     TLorentzVector electron1, electron2;
     electron1.SetPtEtaPhiM(t->electron_pt->at(e1), t->electron_eta->at(e1), t->electron_phi->at(e1), electron_mass);
     electron2.SetPtEtaPhiM(t->electron_pt->at(e2), t->electron_eta->at(e2), t->electron_phi->at(e2), electron_mass);
     const int le = (t->electron_pt->at(e1) > t->electron_pt->at(e2) ? e1 : e2);
     const float leading_electron_pt = t->electron_pt->at(le);
     const float leading_electron_eta = t->electron_eta->at(le);

     float prescale = 1;
     if (!isMC) {
      Trigger* electronTrigger = NULL;
      for (Trigger* trig : electronTriggers) {
       if (trig->m_trig_prescale > 0 &&
           trig->min_pt <= leading_electron_pt &&
           leading_electron_pt <= trig->max_pt &&
           trig->lower_eta <= leading_electron_eta &&
           leading_electron_eta <= trig->upper_eta)
        electronTrigger = trig;
      }
      if (electronTrigger == NULL ||
         !(electronTrigger->m_trig_bool) ||
         electronTrigger->m_trig_prescale <= 0.)
       continue;
      prescale = electronTrigger->m_trig_prescale;
     }
     else prescale = t->eventWeight;

     if (!(t->electron_loose->at(e1) && t->electron_loose->at(e2)))
      continue; // electron quality criteria
     const TLorentzVector Z = electron1 + electron2;
     const float Z_pt = Z.Pt(); 
     const float Z_eta = Z.Eta();
     const float Z_phi = Z.Phi();
     const float Z_m = Z.M();

     // find the 2 highest pt jets
     int leading_jet = -1;
     int subleading_jet = -1;
     TLorentzVector jet_tlv;
     for (int jet = 0; jet < t->jet_n; jet++) {
      jet_tlv.SetPtEtaPhiE(t->jet_pt->at(jet), t->jet_eta->at(jet), t->jet_phi->at(jet), t->jet_e->at(jet));
      if (jet_tlv.DeltaR(electron1) < 0.2 || jet_tlv.DeltaR(electron2) < 0.2)
       continue; // jet isolation criteria
      if (leading_jet == -1 ||
          t->jet_pt->at(leading_jet) < t->jet_pt->at(jet)) {
       subleading_jet = leading_jet;
       leading_jet = jet;
      }
      else if (subleading_jet == -1 ||
               t->jet_pt->at(subleading_jet) < t->jet_pt->at(jet)) {
       subleading_jet = jet;
      }
     }
     // true iff no jets are isolated enough from electrons (and there are
     // jets)
     if (leading_jet == -1) continue;

     // exclude uncorrected jets with pt < 20 GeV
     //if (init_jet_pt->at(leading_jet) < 20) continue;

     const float leading_jet_pt = t->jet_pt->at(leading_jet);
     const float leading_jet_eta = t->jet_eta->at(leading_jet);
     const float leading_jet_phi = t->jet_phi->at(leading_jet);

     const float subleading_jet_pt = ((subleading_jet >= 0 && subleading_jet < t->jet_n) ? t->jet_pt->at(subleading_jet) : 0);
     const float subleading_jet_phi = ((subleading_jet >= 0 && subleading_jet < t->jet_n) ? t->jet_phi->at(subleading_jet) : 0);

     // Calculate opening angle in the transverse plane
     float dPhi = TMath::Abs(leading_jet_phi - Z_phi);
     while (dPhi > pi) dPhi = TMath::Abs(dPhi - 2*pi);
     // cut on Z+jet samples not back-to-back in the transverse plane
     if (dPhi < 7*pi/8) continue;

     // Determine if the opening angle between the subleading jet and the Z
     // gives a non-vanishing pt ratio
     if (subleading_jet_pt > 12) {
      float subleading_dPhi = TMath::Abs(subleading_jet_phi - Z_phi);
      while (subleading_dPhi > pi) subleading_dPhi = TMath::Abs(subleading_dPhi - 2*pi);
      if (subleading_jet_pt / (Z_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2)
       continue; // suppress dijets
     }
  
     // Make sure the jet is within the relevant bounds
     if (leading_jet_eta < etabins[0] ||
         leading_jet_eta > etabins[numetabins])
      continue;
     // Put the jet in the right eta bin
     int etabin = 0;
     while (etabins[etabin] < leading_jet_eta) etabin++;
     etabin--;

     zMassSpectra_AllSigns[1][etabin]->Fill(Z_m, prescale);
     zMassSpectra_AllSigns[1][numetabins]->Fill(Z_m, prescale);

     if (t->electron_charge->at(e1) == t->electron_charge->at(e2)) continue;

     zMassSpectra[1][etabin]->Fill(Z_m, prescale);
     zMassSpectra[1][numetabins]->Fill(Z_m, prescale);
     Zee_n[etabin]++;

     if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m)
      continue; // cut on our sample Z boson mass

     const float leading_jet_pt_err = GetXCalibSystematicError(leading_jet_pt, leading_jet_eta);
     const float ptref = Z_pt * TMath::Cos(pi - dPhi);
     const float ptratio[3] = {leading_jet_pt-leading_jet_pt_err, leading_jet_pt, leading_jet_pt+leading_jet_pt_err};
   
     //for (short errType = 0; errType < 3; errType++) zeeJetHists[errType][etabin]->Fill(Z_pt, ptratio[errType]/ptref, prescale);
     for (short errType = 0; errType < 3; errType++) zeeJetHists[errType][etabin]->Fill(ptref, ptratio[errType]/ptref, prescale);
     //for (short errType = 0; errType < 3; errType++) zeeJetHistsSys[errType][etabin]->Fill(leading_jet_pt, ptratio[errType]/ptref, prescale);
    }
   } // end loop over electron pairs
    // end Z->ee type events

   for (int m1 = 0; m1 < t->muon_n; m1++) {
    if (t->muon_pt->at(m1) < muon_pt_cut) continue;
    if (2.4 < TMath::Abs(t->muon_eta->at(m1))) continue;
    for (int m2 = 0; m2 < m1; m2++) {
     if (t->muon_pt->at(m2) < muon_pt_cut) continue; 
     if (2.4 < TMath::Abs(t->muon_eta->at(m2))) continue;

     TLorentzVector muon1, muon2;
     muon1.SetPtEtaPhiM(t->muon_pt->at(m1), t->muon_eta->at(m1), t->muon_phi->at(m1), muon_mass);
     muon2.SetPtEtaPhiM(t->muon_pt->at(m2), t->muon_eta->at(m2), t->muon_phi->at(m2), muon_mass);
     const int lm = (t->muon_pt->at(m1) > t->muon_pt->at(m2) ? m1 : m2);
     const float leading_muon_pt = t->muon_pt->at(lm);
     const float leading_muon_eta = t->muon_eta->at(lm);

     float prescale = 1;
     if (!isMC) {
      Trigger* muonTrigger = NULL;
      for (Trigger* trig : muonTriggers) {
       if (trig->m_trig_prescale > 0 &&
           trig->min_pt <= leading_muon_pt &&
           leading_muon_pt <= trig->max_pt &&
           trig->lower_eta <= leading_muon_eta &&
           leading_muon_eta <= trig->upper_eta)
        muonTrigger = trig;
      }
      if (muonTrigger == NULL ||
         !(muonTrigger->m_trig_bool))
       continue;
      prescale = muonTrigger->m_trig_prescale;
     }
     else prescale = t->eventWeight;

     //if (!(t->muon_tight->at(m1) || t->muon_tight->at(m2))) continue;
     if (!(t->muon_loose->at(m1) && t->muon_loose->at(m2)))
      continue; // muon quality criteria

     const TLorentzVector Z = muon1 + muon2;
     const float Z_pt = Z.Pt(); 
     const float Z_eta = Z.Eta();
     const float Z_phi = Z.Phi();
     const float Z_m = Z.M();

     // find the 2 highest pt jets
     int leading_jet = -1;
     int subleading_jet = -1;
     TLorentzVector jet_tlv;
     for (int jet = 0; jet < t->jet_n; jet++) {
      jet_tlv.SetPtEtaPhiE(t->jet_pt->at(jet), t->jet_eta->at(jet), t->jet_phi->at(jet), t->jet_e->at(jet));
      if (jet_tlv.DeltaR(muon1) < 0.2 || jet_tlv.DeltaR(muon2) < 0.2)
       continue; // jet isolation criteria
      if (leading_jet == -1 ||
          t->jet_pt->at(leading_jet) < t->jet_pt->at(jet)) {
       subleading_jet = leading_jet;
       leading_jet = jet;
      }
      else if (subleading_jet == -1 ||
               t->jet_pt->at(subleading_jet) < t->jet_pt->at(jet)) {
       subleading_jet = jet;
      }
     }
     // true iff no jet is isolated enough from the muons
     if (leading_jet == -1) continue;

     // exclude uncorrected jets with pt < 20 GeV
     //if (init_jet_pt->at(leading_jet) < 20) continue;

     const float leading_jet_pt = t->jet_pt->at(leading_jet);
     const float leading_jet_eta = t->jet_eta->at(leading_jet);
     const float leading_jet_phi = t->jet_phi->at(leading_jet);

     const float subleading_jet_pt = ((subleading_jet >= 0 && subleading_jet < t->jet_n) ? t->jet_pt->at(subleading_jet) : 0);
     const float subleading_jet_phi = ((subleading_jet >= 0 && subleading_jet < t->jet_n) ? t->jet_phi->at(subleading_jet) : 0);

     // Calculate opening angle in the transverse plane
     float dPhi = TMath::Abs(leading_jet_phi - Z_phi);
     while (dPhi > pi) dPhi = TMath::Abs(dPhi - 2*pi);
     // cut on Z+jet samples not back-to-back in the transverse plane
     if (dPhi < 7*pi/8) continue;

     // Determine if the opening angle between the subleading jet and the Z
     // gives a non-vanishing pt ratio
     if (subleading_jet_pt > 12) {
      float subleading_dPhi = TMath::Abs(subleading_jet_phi - Z_phi);
      while (subleading_dPhi > pi) subleading_dPhi = TMath::Abs(subleading_dPhi - 2*pi);
      if (subleading_jet_pt / (Z_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2)
       continue; // suppress dijets
     }

     // Make sure the jet is within the relevant bounds
     if (leading_jet_eta < etabins[0] ||
         leading_jet_eta > etabins[numetabins])
      continue;
     // Put the jet in the right eta bin
     short etabin = 0;
     while (etabins[etabin] < leading_jet_eta) etabin++;
     etabin--;

     zMassSpectra_AllSigns[0][etabin]->Fill(Z_m, prescale);
     zMassSpectra_AllSigns[0][numetabins]->Fill(Z_m, prescale);

     if (t->muon_charge->at(m1) == t->muon_charge->at(m2)) continue;

     zMassSpectra[0][etabin]->Fill(Z_m, prescale);
     zMassSpectra[0][numetabins]->Fill(Z_m, prescale);
     Zmumu_n[etabin]++;

     if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m)
      continue; // cut on our sample Z boson mass

     const float leading_jet_pt_err = GetXCalibSystematicError(leading_jet_pt, leading_jet_eta);
     const float ptref = Z_pt * TMath::Cos(pi - dPhi);
     const float ptratio[3] = {leading_jet_pt-leading_jet_pt_err, leading_jet_pt, leading_jet_pt+leading_jet_pt_err};

     //for (short errType = 0; errType < 3; errType++)
     // zmumuJetHists[errType][etabin]->Fill(Z_pt, ptratio[errType]/ptref, prescale);
     for (short errType = 0; errType < 3; errType++)
      zmumuJetHists[errType][etabin]->Fill(ptref, ptratio[errType]/ptref, prescale);
     //for (short errType = 0; errType < 3; errType++)
     // zmumuJetHistsSys[errType][etabin]->Fill(leading_jet_pt, ptratio[errType]/ptref, prescale);
    }
   } // end loop over muon pairs
    // end Z->mumu type events

   for (int p = 0; p < t->photon_n; p++) {

    TLorentzVector photon;
    const float this_photon_pt = t->photon_pt->at(p);
    const float this_photon_eta = t->photon_eta->at(p);
    const float this_photon_phi = t->photon_phi->at(p);
    photon.SetPtEtaPhiM(this_photon_pt, this_photon_eta, this_photon_phi, 0);

    if ((1.37 < TMath::Abs(this_photon_eta) &&
         TMath::Abs(this_photon_eta) < 1.52) ||
        2.37 < TMath::Abs(this_photon_eta))
     continue;
    //if (1.37 < TMath::Abs(this_photon_eta))
    // continue;

    float prescale = 1;
    if (!isMC) {
     Trigger* photonTrigger = NULL;
     for (Trigger* trig : photonTriggers) {
      if (trig->m_trig_prescale > 0 &&
          trig->min_pt <= this_photon_pt &&
          this_photon_pt <= trig->max_pt)
       photonTrigger = trig;
     }
     if (photonTrigger == NULL ||
         !(photonTrigger->m_trig_bool))
      continue;
     prescale = photonTrigger->m_trig_prescale;
    }
    else prescale = t->eventWeight;

    if (!t->photon_tight->at(p))
     continue; // require tight cuts on photons

    if (t->photon_topoetcone40->at(p) > isolationEnergyIntercept + isolationEnergySlope*this_photon_pt)
     continue; // require maximum isolation energy on gammas

    int leading_jet = -1;
    int subleading_jet = -1;
    TLorentzVector jet_tlv;
    for (int jet = 0; jet < t->jet_n; jet++) {
      jet_tlv.SetPtEtaPhiE(t->jet_pt->at(jet), t->jet_eta->at(jet), t->jet_phi->at(jet), t->jet_e->at(jet));
      if (jet_tlv.DeltaR(photon) < 0.6) continue;

      if (leading_jet == -1 ||
          t->jet_pt->at(leading_jet) < t->jet_pt->at(jet)) {
       subleading_jet = leading_jet;
       leading_jet = jet;
      }
      else if (subleading_jet == -1 ||
               t->jet_pt->at(subleading_jet) < t->jet_pt->at(jet)) {
       subleading_jet = jet;
      }
    }

    // true iff there are no jets
    if (leading_jet == -1) continue;
    // exclude uncorrected jets with pt < 20 GeV
    //if (init_jet_pt->at(leading_jet) < 20) continue;

    const float leading_jet_pt = t->jet_pt->at(leading_jet);
    const float leading_jet_eta = t->jet_eta->at(leading_jet);
    const float leading_jet_phi = t->jet_phi->at(leading_jet);

    const float subleading_jet_pt = ((subleading_jet >= 0 && subleading_jet < t->jet_n) ? t->jet_pt->at(subleading_jet) : 0);
    const float subleading_jet_phi = ((subleading_jet >= 0 && subleading_jet < t->jet_n) ? t->jet_phi->at(subleading_jet) : 0);

    if (leading_jet_eta < etabins[0] ||
        leading_jet_eta > etabins[numetabins])
     continue;

    // Calculate opening angle in the transverse plane
    float dPhi = TMath::Abs(leading_jet_phi - this_photon_phi);
    while (dPhi > pi) dPhi = TMath::Abs(dPhi - 2*pi);
    // cut on gamma+jet samples not back-to-back in the transverse plane
    if (dPhi < 7*pi/8) continue;

    // Determine if the opening angle between the subleading jet and the Z
    // gives a non-vanishing pt ratio
    if (subleading_jet_pt > 12) {
     float subleading_dPhi = TMath::Abs(subleading_jet_phi - this_photon_phi);
     while (subleading_dPhi > pi)
      subleading_dPhi = TMath::Abs(subleading_dPhi - 2*pi);
     if (subleading_jet_pt / (this_photon_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2)
      continue; // suppress dijets
    }

    // Make sure the jet is within the relevant bounds
    if (leading_jet_eta < etabins[0] ||
        leading_jet_eta > etabins[numetabins])
     continue;
    // Put the jet in the right eta bin
    short etabin = 0;
    while (etabins[etabin] < leading_jet_eta) etabin++;
    etabin--;

    const float leading_jet_pt_err = GetXCalibSystematicError(leading_jet_pt, leading_jet_eta);
    const float ptref = this_photon_pt * TMath::Cos(pi - dPhi);
    const float ptratio[3] = {leading_jet_pt-leading_jet_pt_err, leading_jet_pt, leading_jet_pt+leading_jet_pt_err};

    //for (short errType = 0; errType < 3; errType++)
    // gJetHists[errType][etabin]->Fill(this_photon_pt, ptratio[errType]/ptref, prescale);
    for (short errType = 0; errType < 3; errType++)
     gJetHists[errType][etabin]->Fill(ptref, ptratio[errType]/ptref, prescale);

    // if data, calculate additional systematics for this jet
    if (!isMC)
     gJetHistsSys[1][etabin]->Fill(leading_jet_pt, GetNewXCalibSystematicError(leading_jet_eta, ptref)/leading_jet_pt, prescale);
    g_n[etabin]++;
   }
    
  } // end loop over events

  // close root files with systematics
  xCalibSystematicsFile->Close();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  dataOverMCFile->Close();
  if (dataOverMCFile) delete dataOverMCFile;
  
  /** End event loop **/

  TFile* outFile = new TFile(Form("%sdataSet_%s.root", rootPath.c_str(), identifier.c_str()), "RECREATE");


  // Write histograms to output and clean memory
  for (short etabin = 0; etabin < numetabins; etabin++) {
   for (short errType = 0; errType < 3; errType++) {
    zmumuJetHists[errType][etabin]->Write();
    if (zmumuJetHists[errType][etabin])
     delete zmumuJetHists[errType][etabin];
    //zmumuJetHistsSys[errType][etabin]->Write();
    //if (zmumuJetHistsSys[errType][etabin])
    // delete zmumuJetHistsSys[errType][etabin];
    zeeJetHists[errType][etabin]->Write();
    if (zeeJetHists[errType][etabin])
     delete zeeJetHists[errType][etabin];
    //zeeJetHistsSys[errType][etabin]->Write();
    //if (zeeJetHistsSys[errType][etabin])
    // delete zeeJetHistsSys[errType][etabin];
    gJetHists[errType][etabin]->Write();
    if (gJetHists[errType][etabin])
     delete gJetHists[errType][etabin];
    if (errType == 1) {
     gJetHistsSys[errType][etabin]->Write();
     if (gJetHistsSys[errType][etabin])
      delete gJetHistsSys[errType][etabin];
    }
   }
   for (short spcType = 0; spcType < 2; spcType++) {
    zMassSpectra[spcType][etabin]->Write();
    if (zMassSpectra[spcType][etabin])
     delete zMassSpectra[spcType][etabin];
    zMassSpectra_AllSigns[spcType][etabin]->Write();
    if (zMassSpectra_AllSigns[spcType][etabin])
     delete zMassSpectra_AllSigns[spcType][etabin];
    if (etabin == 0) {
     zMassSpectra[spcType][numetabins]->Write();
     if (zMassSpectra[spcType][numetabins])
      delete zMassSpectra[spcType][numetabins];
     zMassSpectra_AllSigns[spcType][numetabins]->Write();
     if (zMassSpectra_AllSigns[spcType][numetabins])
      delete zMassSpectra_AllSigns[spcType][numetabins];
    }
   }
  }

  TVectorD infoVec(2+3*numetabins);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  for (short etabin = 0; etabin < numetabins; etabin++)
   infoVec[2+etabin] = (double)Zee_n[etabin];
  for (short etabin = 0; etabin < numetabins; etabin++)
   infoVec[2+numetabins+etabin] = (double)Zmumu_n[etabin];
  for (short etabin = 0; etabin < numetabins; etabin++)
   infoVec[2+2*numetabins+etabin] = (double)g_n[etabin];
  infoVec.Write(Form("infoVec_%s", identifier.c_str()));
  outFile->Close();
  if (outFile) delete outFile;
  return;
}
