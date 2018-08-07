#include "ZGammaJetCrossCheck.h"

TFile* xCalibSystematicsFile = NULL;
TFile* dataOverMCFile = NULL;

vector<Trigger*> electronTriggers = {};
vector<Trigger*> muonTriggers = {};
vector<Trigger*> photonTriggers = {};

bool isPeriodA = false;


double GetXCalibSystematicError(const double jpt, const double jeta) {
  TFile* file = xCalibSystematicsFile;
  if (!file || !file->IsOpen())
   return 0;

  if (TMath::Abs(jeta) < xcalibEtabins[0] ||
      xcalibEtabins[sizeof(xcalibEtabins)/sizeof(xcalibEtabins[0]) -1] < TMath::Abs(jeta)) {
   return 0;
  }

  short etabin = 0;
  while (xcalibEtabins[etabin] < TMath::Abs(jeta)) etabin++;
  etabin--;

  const TString hname = TString("fsys_rel_") + Form("%i", etabin);
  TH1D* fsys_rel = (TH1D*)file->Get(hname.Data());

  return TMath::Abs(fsys_rel->GetBinContent(fsys_rel->FindBin(jpt)) - 1) * jpt;
}


double GetNewXCalibSystematicError(const double jeta, const double refpt) {
  TFile* file = dataOverMCFile;
  if (!file || !file->IsOpen())
   return 0;

  if (jeta < etabins[0] ||
      etabins[numetabins] < jeta) {
   return 0;
  }

  short bin = 0;
  while (bin <= numetabins && etabins[bin] < jeta) bin++;
  bin--;

  const char* period = (isPeriodA ?  "periodA" : "periodB");
  const TString hname = TString(Form("gJetPtRatio_diff%i_stat_%s", bin, period));
  TH1D* hist = (TH1D*)file->Get(hname.Data());

  return hist->GetBinContent(hist->FindBin(refpt)) * refpt;
}


TString GetIdentifier (const int dataSet, const bool isMC, const bool isValidationSample, const bool periodA) {
  if (!isMC) return to_string(dataSet);
  TString id = "";
  if (isPeriodA) id = "pPb_";
  else id = "Pbp_";
  if (dataSet > 0) { // true for GammaJet samples
   if (isValidationSample) id = id + "Valid_";
   else id = id + "Overlay_";
   id = id + "GammaJet_Slice" + to_string(dataSet);
  }
  else {
   if (dataSet == 0) { // true for Zmumu samples
    id = id + "ZmumuJet";
   }
   else if (dataSet == -6) { // true for Zee overlay samples
    id = id + "ZeeJet_Overlay";
   }
   else { // true for Zee signal-only samples
    id = id + "ZeeJet_Slice" + to_string(-dataSet);
   }
  }
  return id;
}


void ZGammaJetCrossCheck (const int dataSet,
                          const double luminosity,
                          const double weight,
                          const bool isMC, 
                          const bool isMCperiodAflag,
                          const TString inFileName)
{

  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  if (!isMC) isPeriodA = dataSet < 313500;
  else isPeriodA = isMCperiodAflag;

  const bool isValidationSample = isMC && (TString(inFileName).Contains("valid") || TString(inFileName).Contains("Zee"));
  const TString identifier = GetIdentifier(dataSet, isMC, isValidationSample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   if (inFileName == "") {
    if (!isMC) fileIdentifier = to_string(dataSet);
    else fileIdentifier = TString(dataSet > 0 ? ("Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet" : ("ZeeJet"+to_string(-dataSet)))) + (isMCperiodAflag ? ".pPb":".Pbp");
   } else fileIdentifier = inFileName;

   TSystemDirectory dir(dataPath.Data(), dataPath.Data());
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

  if (!isMC) {
   for (int electronTriggerN = 0; electronTriggerN < electronTrigLength; electronTriggerN++) {
    Trigger* temp = new Trigger(electronTriggerNames[electronTriggerN], electronTriggerMinPtCuts[electronTriggerN], -2.47, 2.47);
    temp->min_pt = electronTriggerMinPtCuts[electronTriggerN];
    temp->max_pt = electronTriggerMaxPtCuts[electronTriggerN];
    electronTriggers.push_back(temp);
    tree->SetBranchAddress(electronTriggerNames[electronTriggerN], &(temp->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", electronTriggerNames[electronTriggerN]), &(temp->trigPrescale));
   }

   for (int muonTriggerN = 0; muonTriggerN < muonTrigLength; muonTriggerN++) {
    Trigger* temp = new Trigger(muonTriggerNames[muonTriggerN], muonTriggerMinPtCuts[muonTriggerN], -2.40, 2.40);
    temp->min_pt = muonTriggerMinPtCuts[muonTriggerN];
    temp->max_pt = muonTriggerMaxPtCuts[muonTriggerN];
    muonTriggers.push_back(temp);
    tree->SetBranchAddress(muonTriggerNames[muonTriggerN], &(temp->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", muonTriggerNames[muonTriggerN]), &(temp->trigPrescale));
   }

   for (int photonTriggerN = 0; photonTriggerN < photonTrigLength; photonTriggerN++) {
    Trigger* temp = new Trigger(photonTriggerNames[photonTriggerN], photonTriggerMinPtCuts[photonTriggerN], -2.47, 2.47);
    temp->min_pt = photonTriggerMinPtCuts[photonTriggerN];
    temp->max_pt = photonTriggerMaxPtCuts[photonTriggerN];
    photonTriggers.push_back(temp);
    tree->SetBranchAddress(photonTriggerNames[photonTriggerN], &(temp->trigBool));
    tree->SetBranchAddress(Form("%s_prescale", photonTriggerNames[photonTriggerN]), &(temp->trigPrescale));
   }
  } // end branch triggers

  // initialize histograms
  TH2D* zeeJetHists[3][numetabins+1];
  //TH2D* zeeJetHistsSys[3][numetabins];
  TH2D* zmumuJetHists[3][numetabins+1];
  //TH2D* zmumuJetHistsSys[3][numetabins];
  TH2D* gJetHists[3][numetabins+1];
  TH2D* gJetHistsSys[3][numetabins+1];
  TH1D* zMassSpectra[2][numetabins+1];
  TH2D* lowResponseEtaPhi = new TH2D(Form("lowResponseEtaPhi_dataSet%s", identifier.Data()), ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  TH2D* highResponseEtaPhi = new TH2D(Form("highResponseEtaPhi_dataSet%s", identifier.Data()), ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  TH2D* lowResponseEtaPt = new TH2D(Form("lowResponseEtaPt_dataSet%s", identifier.Data()), ";#eta;#it{p}_{T}^{reco};", 25, -2.5, 2.5, 66, 20, 350);
  TH2D* highResponseEtaPt = new TH2D(Form("highResponseEtaPt_dataSet%s", identifier.Data()), ";#eta;#phi;", 25, -2.5, 2.5, 66, 20, 350);
  TH2D* responsePt = new TH2D(Form("responsePt_dataSet%s", identifier.Data()), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};#it{p}_{T}^{reco} #left[GeV#right];", 200, 0, 4.0, 66, 20, 350);
  TH1D* jetSpectrum = new TH1D(Form("jetPt_dataSet%s", identifier.Data()), ";#it{p}_{T}^{jet} #left[GeV#right];", 33, 20, 350);
  TH1D* jetEnergyResponseCalib[numpbins+1][numetabins+1];
  TH1D* jetEnergyResponseReco[numpbins+1][numetabins+1];
  TH1D* photonEnergyResponse[numpbins+1][numetabins+1];
  TH1D* electronEnergyScale = new TH1D(Form("electronEnergyScale_dataSet%s", identifier.Data()), "", 200, 0, 2);
  electronEnergyScale->Sumw2();


  for (short etabin = 0; etabin <= numetabins; etabin++) {
   for (short errType = 0; errType < 3; errType++) {
    TString error = "sys_lo";
    if (errType == 1) error = "stat";
    else if (errType == 2) error = "sys_hi";

    zeeJetHists[errType][etabin] = new TH2D(Form("zeeJetPtRatio_dataSet%s_etabin%i_%s_%s", identifier.Data(), etabin, (isMC ? "mc":"data"), error.Data()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    zeeJetHists[errType][etabin]->Sumw2();
    //zeeJetHistsSys[errType][etabin] = new TH2D(Form("zeeJetPtRatioSys_dataSet%s_etabin%i_%s_%s", identifier.Data(), etabin, (isMC ? "mc":"data"), error.Data()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zeeJetHistsSys[errType][etabin]->Sumw2();

    zmumuJetHists[errType][etabin] = new TH2D(Form("zmumuJetPtRatio_dataSet%s_etabin%i_%s_%s", identifier.Data(), etabin, (isMC ? "mc":"data"), error.Data()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    zmumuJetHists[errType][etabin]->Sumw2();
    //zmumuJetHistsSys[errType][etabin] = new TH2D(Form("zmumuJetPtRatioSys_dataSet%s_etabin%i_%s_%s", identifier.Data(), etabin, (isMC ? "mc":"data"), error.Data()), "", numpzbins, pzbins, numxjrefbins, xjrefbins);
    //zmumuJetHistsSys[errType][etabin]->Sumw2();

    gJetHists[errType][etabin] = new TH2D(Form("gJetPtRatio_dataSet%s_etabin%i_%s_%s", identifier.Data(), etabin, (isMC ? "mc":"data"), error.Data()), "", numpbins, pbins, numxjrefbins, xjrefbins);
    gJetHists[errType][etabin]->Sumw2();
    if (errType == 1) {
     gJetHistsSys[errType][etabin] = new TH2D(Form("gJetPtRatioSys_dataSet%s_etabin%i_%s_%s", identifier.Data(), etabin, (isMC ? "mc":"data"), error.Data()), "", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
     gJetHistsSys[errType][etabin]->Sumw2();
    }
   }

   for (short spcType = 0; spcType < 2; spcType++) {
    const TString species = (spcType == 0 ? "mumu":"ee");

    zMassSpectra[spcType][etabin] = new TH1D(Form("z%sMassSpectrum_dataSet%s_%s_etabin%i", species.Data(), identifier.Data(), (isMC ? "mc":"data"), etabin), "", 50, 60, 110);
    zMassSpectra[spcType][etabin]->Sumw2();
   }
  }

  for (short pbin = 0; pbin <= numpbins; pbin++) {
   for (short etabin = 0; etabin <= numetabins; etabin++) {
    jetEnergyResponseCalib[pbin][etabin] = new TH1D(Form("jetEnergyResponseCalib_dataSet%s_pbin%i_etabin%i", identifier.Data(), pbin, etabin), "", 50, 0, 2);
    jetEnergyResponseCalib[pbin][etabin]->Sumw2();
    jetEnergyResponseReco[pbin][etabin] = new TH1D(Form("jetEnergyResponseReco_dataSet%s_pbin%i_etabin%i", identifier.Data(), pbin, etabin), "", 50, 0, 2);
    jetEnergyResponseReco[pbin][etabin]->Sumw2();
    photonEnergyResponse[pbin][etabin] = new TH1D(Form("photonEnergyResponse_dataSet%s_pbin%i_etabin%i", identifier.Data(), pbin, etabin), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 50, 0, 2);
    photonEnergyResponse[pbin][etabin]->Sumw2();
   }
  }

  int* nJet[numpbins+1] = {};
  int* nGamma[numpbins+1] = {};
  for (short pbin = 0; pbin <= numpbins; pbin++) {
   nJet[pbin] = new int[numetabins+1];
   nGamma[pbin] = new int[numetabins+1];
   for (short etabin = 0; etabin <= numetabins; etabin++) {
    nJet[pbin][etabin] = 0;
    nGamma[pbin][etabin] = 0;
   }
  }
  int nZeeMass[numetabins+1] = {};
  int nZeeJet[numetabins+1] = {};
  int nZmumuMass[numetabins+1] = {};
  int nZmumuJet[numetabins+1] = {};
  int nGammaJet[numetabins+1] = {};

  xCalibSystematicsFile = new TFile(rootPath + "cc_sys_090816.root", "READ");
  dataOverMCFile = new TFile(rootPath + "cc_difference.root", "READ");

  const long long numEntries = tree->GetEntries();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (long long entry = 0; entry < numEntries; entry++) {
   tree->GetEntry(entry);


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if ((t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at(0) != 1)) continue;


   /////////////////////////////////////////////////////////////////////////////
   // if MC, fill jet energy response
   /////////////////////////////////////////////////////////////////////////////
   if (isMC) {
    for (int j = 0; j < t->jet_n; j++) {
     if (t->jet_pt->at(j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (InDisabledHEC (t->jet_eta->at(j), t->jet_phi->at(j)))
      continue; // only truth match jets outside disabled HEC
     if (!InHadCal (t->jet_eta->at(j)))
      continue; // reject jets reconstructed outside reasonable HCal bounds.

     double minDeltaR = 1000;
     // loop over truth+reco electrons, muons, and photons
     for (int e = 0; e < t->truth_electron_n; e++) { // truth electrons
      //if (t->truth_electron_pt->at(e) < electron_pt_cut) continue;
      const double dR = DeltaR (t->jet_eta->at(j), t->truth_electron_eta->at(e), t->jet_phi->at(j), t->truth_electron_phi->at(e));
      if (dR < minDeltaR) minDeltaR = dR;
     }
     //for (int e = 0; e < t->electron_n; e++) { // reco electrons
     // if (t->electron_pt->at(e) < electron_pt_cut) continue;
     // const double dR = DeltaR (t->jet_eta->at(j), t->electron_eta->at(e), t->jet_phi->at(j), t->electron_phi->at(e));
     // if (dR < minDeltaR) minDeltaR = dR;
     //}
     //for (int m = 0; m < t->truth_muon_n; m++) { // truth muons
     // const double dR = DeltaR (t->jet_eta->at(j), t->truth_muon_eta->at(m), t->jet_phi->at(j), t->truth_muon_phi->at(m));
     // if (dR < minDeltaR) minDeltaR = dR;
     //}
     //for (int m = 0; m < t->muon_n; m++) { // reco muons
     // const double dR = DeltaR (t->jet_eta->at(j), t->muon_eta->at(m), t->jet_phi->at(j), t->muon_phi->at(m));
     // if (dR < minDeltaR) minDeltaR = dR;
     //}
     for (int p = 0; p < t->truth_photon_n; p++) { // truth photons
      //if (t->truth_photon_pt->at(p) < photon_pt_cut) continue;
      const double dR = DeltaR (t->jet_eta->at(j), t->truth_photon_eta->at(p), t->jet_phi->at(j), t->truth_photon_phi->at(p));
      if (dR < minDeltaR) minDeltaR = dR;
     }
     //for (int p = 0; p < t->photon_n; p++) { // photons
     // if (t->photon_pt->at(p) < photon_pt_cut) continue;
     // const double dR = DeltaR (t->jet_eta->at(j), t->photon_eta->at(p), t->jet_phi->at(j), t->photon_phi->at(p));
     // if (dR < minDeltaR) minDeltaR = dR;
     //}
     if (minDeltaR < 0.4)
      continue; // reject jets close to some other lepton or photon

     jetSpectrum->Fill (t->jet_pt->at(j), t->eventWeight / weight / weight);

     minDeltaR = 1000;
     int truth_jet = -1;
     for (int tj = 0; tj < t->truth_jet_n; tj++) {
      double dR = DeltaR (t->jet_eta->at(j), t->truth_jet_eta->at(tj), t->jet_phi->at(j), t->truth_jet_phi->at(tj));
      if (dR < minDeltaR) {
       minDeltaR = dR;
       truth_jet = tj;
      }
     }
     if (0 <= truth_jet && truth_jet < t->truth_jet_n && minDeltaR < 0.4) {

      // Put the jet in the right eta bin
      short etabin = 0;
      // Make sure jet is in eta bounds
      if (etabins[0] < t->jet_eta->at(j) ||
          t->jet_eta->at(j) < etabins[numetabins]) {
       while (etabins[etabin] < t->jet_eta->at(j)) etabin++;
      }
      etabin--;
      short pbin = 0;
      if (pbins[0] < t->jet_pt->at(j) ||
          t->jet_pt->at(j) < pbins[numpbins]) {
       while (pbins[pbin] < t->jet_pt->at(j)) pbin++;
      }
      pbin--;

      const double jer = t->jet_pt->at(j) / t->truth_jet_pt->at(truth_jet);
      const double pjer = t->precalib_jet_pt->at(j) / t->truth_jet_pt->at(truth_jet);
      if (0 <= pbin && pbin < numpbins &&
          0 <= etabin && etabin < numetabins) {
       nJet[pbin][etabin]++;
       jetEnergyResponseCalib[pbin][etabin]->Fill (jer, t->eventWeight / weight / weight);
       jetEnergyResponseReco[pbin][etabin]->Fill (pjer, t->eventWeight / weight / weight);
      }
      if (0 <= pbin && pbin < numpbins) {
       nJet[pbin][numetabins]++;
       jetEnergyResponseCalib[pbin][numetabins]->Fill (jer, t->eventWeight / weight / weight);
       jetEnergyResponseReco[pbin][numetabins]->Fill (pjer, t->eventWeight / weight / weight);
      }
      if (0 <= etabin && etabin < numetabins) {
       nJet[numpbins][etabin]++;
       jetEnergyResponseCalib[numpbins][etabin]->Fill (jer, t->eventWeight / weight / weight);
       jetEnergyResponseReco[numpbins][etabin]->Fill (pjer, t->eventWeight / weight / weight);
      }
      nJet[numpbins][numetabins]++;
      jetEnergyResponseCalib[numpbins][numetabins]->Fill (jer, t->eventWeight / weight / weight);
      jetEnergyResponseReco[numpbins][numetabins]->Fill (pjer, t->eventWeight / weight / weight);

      // Extra plots
      if (t->jet_pt->at(j) / t->truth_jet_pt->at(truth_jet) < 1.18) {
       lowResponseEtaPhi->Fill (t->jet_eta->at(j), t->jet_phi->at(j), t->eventWeight / weight / weight);
       lowResponseEtaPt->Fill (t->jet_eta->at(j), t->jet_pt->at(j), t->eventWeight / weight / weight);
      }
      else {
       highResponseEtaPhi->Fill (t->jet_eta->at(j), t->jet_phi->at(j), t->eventWeight / weight / weight);
       highResponseEtaPt->Fill (t->jet_eta->at(j), t->jet_pt->at(j), t->eventWeight / weight / weight);
      }
      responsePt->Fill (jer, t->jet_pt->at(j), t->eventWeight / weight / weight);

     }
    }
   }// end jet energy response


   /////////////////////////////////////////////////////////////////////////////
   // now reject events with less than 2 muons, 2 electrons AND 1 photon
   /////////////////////////////////////////////////////////////////////////////
   if (t->photon_n < 1 && t->muon_n < 2 && t->electron_n < 2) continue;


   /////////////////////////////////////////////////////////////////////////////
   // Z (ee) + jet events
   /////////////////////////////////////////////////////////////////////////////
   for (int e1 = 0; e1 < t->electron_n; e1++) { // loop over primary electron

    // electron cuts
    if (t->electron_pt->at(e1) < electron_pt_cut)
     continue; // basic electron pT cuts
    if (!InEMCal (t->electron_eta->at(e1)))
     continue; // reject electrons reconstructed outside EMCal
    if (!t->electron_loose->at(e1))
     continue; // reject non-loose electrons
    if (t->electron_d0sig->at(e1) > 5)
     continue; // d0 (transverse impact parameter) significance cut
    if (t->electron_delta_z0_sin_theta->at(e1) > 0.5)
     continue; // z0 (longitudinal impact parameter) vertex compatibility cut

    // fill the electron energy scale plot first, finding the truth-reco pair
    if (isMC && t->truth_electron_n > 0) {
     int closestEl = (e1 < t->truth_electron_n ? e1:t->truth_electron_n-1); // best initial guess
     
     double minDeltaR = DeltaR (t->electron_eta->at(e1), t->truth_electron_eta->at(closestEl), t->electron_phi->at(e1), t->truth_electron_phi->at(closestEl));
     for (int te = 0; te < t->truth_electron_n; te++) {
      if (t->truth_electron_pt->at(te) < electron_pt_cut) continue;
      if (t->electron_charge->at(e1) != t->truth_electron_charge->at(te)) continue;
      const double dR = DeltaR (t->electron_eta->at(e1), t->truth_electron_eta->at(te), t->electron_phi->at(e1), t->truth_electron_phi->at(te));
      if (dR < minDeltaR) {
       closestEl = te;
       minDeltaR = dR;
      }
     }
     if (minDeltaR < 0.4)
      electronEnergyScale->Fill (t->electron_pt->at(e1) / t->truth_electron_pt->at(closestEl));//, t->eventWeight);
    }

    for (int e2 = 0; e2 < e1; e2++) { // loop over secondary electron
     // electron cuts
     if (t->electron_pt->at(e2) < electron_pt_cut)
      continue; // basic electron pT cut
     if (!InEMCal (t->electron_eta->at(e2)))
      continue; // reject electrons reconstructed outside the EMCal
     if (!t->electron_loose->at(e2))
      continue; // reject non-loose electrons
     if (t->electron_d0sig->at(e2) > 5)
      continue; // d0 (transverse impact parameter) significance cut
     if (t->electron_delta_z0_sin_theta->at(e2) > 0.5)
      continue; // z0 (longitudinal impact parameter) vertex compatibility cut

     // relevant electron kinematic data
     TLorentzVector electron1, electron2;
     electron1.SetPtEtaPhiM (t->electron_pt->at(e1), t->electron_eta->at(e1), t->electron_phi->at(e1), 0);//electron_mass);
     electron2.SetPtEtaPhiM (t->electron_pt->at(e2), t->electron_eta->at(e2), t->electron_phi->at(e2), 0);//electron_mass);
     const int le = (t->electron_pt->at(e1) > t->electron_pt->at(e2) ? e1 : e2);
     const double leading_electron_pt = t->electron_pt->at(le);
     const double leading_electron_eta = t->electron_eta->at(le);

     // triggering and event weighting
     float prescale = 1;
     if (!isMC) {
      Trigger* electronTrigger = NULL;
      for (Trigger* trig : electronTriggers) {
       if (trig->trigPrescale > 0 &&
           (electronTrigger == NULL ||
            (trig->trigPrescale < electronTrigger->trigPrescale &&
             trig->min_pt <= leading_electron_pt &&
             leading_electron_pt <= trig->max_pt)))
        electronTrigger = trig;
      }
      if (electronTrigger == NULL ||
          !electronTrigger->trigBool ||
          electronTrigger->trigPrescale <= 0.)
       continue;
      //bool nonElectronTrigger = false;
      //for (Trigger* trig : photonTriggers) {
      // if (trig->trigBool) nonElectronTrigger = true;
      //}
      //for (Trigger* trig : muonTriggers) {
      // if (trig->trigBool) nonElectronTrigger = true;
      //}
      //if (nonElectronTrigger)
      // continue; // reject events where a non-electron trigger also fired.
      prescale = electronTrigger->trigPrescale;
     }
     else prescale = t->eventWeight;
     prescale = prescale / weight / weight;

     // Reco ee invariant 4-momentum (best guess for Z boson)
     const TLorentzVector Z = electron1 + electron2;
     const double Z_pt = Z.Pt(); 
     const double Z_eta = Z.Eta();
     const double Z_phi = Z.Phi();
     const double Z_m = Z.M();

     // Z boson, dielectron cuts
     if (Z_pt < Z_pt_cut)
      continue; // pt cut on Z bosons
     if (t->electron_charge->at(e1) == t->electron_charge->at(e2))
      continue; // opposite charge requirement

     // Fill dielectron mass spectra
     zMassSpectra[1][numetabins]->Fill (Z_m, prescale);

     // Increment total Z count
     nZeeMass[numetabins]++;

     // Put the Z in the right eta bin
     short etabin = 0;
     if (etabins[0] < Z_eta ||
         Z_eta < etabins[numetabins]) {
      while (etabins[etabin] < Z_eta) etabin++;
     }
     etabin--;

     if (etabin != -1) {
      // Fill additional dielectron mass spectrum
      zMassSpectra[1][etabin]->Fill (Z_m, prescale);

      // Fill appropriate Z counter
      nZeeMass[etabin]++;
     }

     // additional Z boson cuts
     if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m)
      continue; // cut on our sample Z boson mass

     // Jet finding
     int lj = -1; // allowed to be in disabled HEC while finding
     int sj = -1; // not allowed to be in disabled HEC while finding
     for (int j = 0; j < t->jet_n; j++) {
      if (t->jet_pt->at(j) < jet_pt_cut)
       continue; // basic jet pT cut
      if (!InHadCal (t->jet_eta->at(j), 0.4))
       continue; // require jets inside hadronic calorimeter
      // Require jets to be isolated from both electrons
      if (DeltaR (t->electron_eta->at(e1), t->jet_eta->at(j), t->electron_phi->at(e1), t->jet_phi->at(j)) < 0.2 ||
          DeltaR (t->electron_eta->at(e2), t->jet_eta->at(j), t->electron_phi->at(e2), t->jet_phi->at(j)) < 0.2)
       continue; // require jets to be isolated from both electrons
      if (DeltaPhi (t->jet_phi->at(j), Z_phi) < 3*pi/4)
       continue; // require jet to be back-to-back with Z in transverse plane

      // compare to leading jet
      else if (lj == -1 ||
               t->jet_pt->at(lj) < t->jet_pt->at(j)) {
       if (lj != -1 && !InDisabledHEC (t->jet_eta->at(lj), t->jet_phi->at(lj)))
        sj = lj;
       lj = j;
      }

      // compare to subleading jet
      else if ((sj == -1 ||
                t->jet_pt->at(sj) < t->jet_pt->at(j)) &&
               !InDisabledHEC (t->jet_eta->at(j), t->jet_phi->at(j))) {
       sj = j;
      }
     } // end jet finding loop
     if (lj == -1) // true iff no candidate jet is found
      continue; // reject on no candidate jet

     // exclude uncorrected jets with pt < 20 GeV
     //if (precalib_jet_pt->at(lj) < 20) continue;

     // relevant jet kinematic data
     const double ljet_pt = t->jet_pt->at(lj);
     const double ljet_eta = t->jet_eta->at(lj);
     const double ljet_phi = t->jet_phi->at(lj);
     const double sjet_pt = ((0 <= sj && sj < t->jet_n) ? t->jet_pt->at(sj) : 0);
     const double sjet_phi = ((0 <= sj && sj < t->jet_n) ? t->jet_phi->at(sj) : 0);

     // jet cuts
     if (InDisabledHEC (ljet_eta, ljet_phi))
      continue; // Reject event on additional HEC cuts
     //if (ljet_eta < etabins[0] ||
     //    ljet_eta > etabins[numetabins])
     // continue; // Make sure jet is in eta bounds
     if (sjet_pt > 12) {
      const double subleading_dPhi = DeltaPhi (sjet_phi, Z_phi);
      //if (sjet_pt / Z_pt > 0.3)
      if (sjet_pt / (Z_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2)
      //if (sjet_pt * TMath::Cos(pi - subleading_dPhi) / Z_pt > 0.2)
       continue; // suppress dijets by requiring leading jet to dominate ptref
     }

     // Calculate opening angle in the transverse plane
     const double dPhi = DeltaPhi (ljet_phi, Z_phi);

     // Calculate systematics on jet pT
     const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError (ljet_pt, ljet_eta));

     // Calculate ptref and xjrefs
     const double ptref = Z_pt * TMath::Cos(pi - dPhi);
     const double ptratio[3] = {(ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

     // Fill xjref histograms 
     //for (short errType = 0; errType < 3; errType++) zeeJetHists[errType][etabin]->Fill (Z_pt, ptratio[errType]/ptref, prescale);
     for (short errType = 0; errType < 3; errType++) {
      if (etabin != -1) zeeJetHists[errType][etabin]->Fill (Z_pt, ptratio[errType], prescale);
      zeeJetHists[errType][numetabins]->Fill (Z_pt, ptratio[errType], prescale);
     }

     // Increment Z+jet counters
     if (etabin != -1) nZeeJet[etabin]++;
     nZeeJet[numetabins]++;

     // Fill additional systematics histograms
     //for (short errType = 0; errType < 3; errType++) zeeJetHistsSys[errType][etabin]->Fill (ljet_pt, ptratio[errType]/ptref, prescale);
    }
   } // end loop over electron pairs
   // end Z->ee type events


   /////////////////////////////////////////////////////////////////////////////
   // Z->mumu + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int m1 = 0; m1 < t->muon_n; m1++) { // loop over primary muon
    // primary muon cuts
    if (t->muon_pt->at(m1) < muon_pt_cut)
     continue; // basic muon pT cuts
    if (!t->muon_loose->at(m1))
     continue; // require loose muons
    if (2.4 < abs(t->muon_eta->at(m1)))
     continue; // reject muons reconstructed outside muon spectrometer

    for (int m2 = 0; m2 < m1; m2++) { // loop over secondary muon
     // secondary muon cuts
     if (t->muon_pt->at(m2) < muon_pt_cut)
      continue; // basic muon pT cuts
     if (!t->muon_loose->at(m2))
      continue; // require loose muons
     if (2.4 < abs(t->muon_eta->at(m2))) 
      continue; // reject muons reconstructed outside muon spectrometer

     // relevant muon kinematic data
     TLorentzVector muon1, muon2;
     muon1.SetPtEtaPhiM (t->muon_pt->at(m1), t->muon_eta->at(m1), t->muon_phi->at(m1), muon_mass);
     muon2.SetPtEtaPhiM (t->muon_pt->at(m2), t->muon_eta->at(m2), t->muon_phi->at(m2), muon_mass);
     const int lm = (t->muon_pt->at(m1) > t->muon_pt->at(m2) ? m1 : m2);
     const double leading_muon_pt = t->muon_pt->at(lm);
     const double leading_muon_eta = t->muon_eta->at(lm);

     // triggering and event weighting
     double prescale = 1;
     if (!isMC) {
      Trigger* muonTrigger = NULL;
      for (Trigger* trig : muonTriggers) {
       if (trig->trigPrescale > 0 &&
           (muonTrigger == NULL ||
            (trig->trigPrescale < muonTrigger->trigPrescale &&
             trig->min_pt <= leading_muon_pt &&
             leading_muon_pt <= trig->max_pt)))
        muonTrigger = trig;
      }
      if (muonTrigger == NULL ||
          !muonTrigger->trigBool ||
          muonTrigger->trigPrescale <= 0.)
       continue;
      bool nonMuonTrigger = false;
      for (Trigger* trig : electronTriggers) {
       if (trig->trigBool) nonMuonTrigger = true;
      }
      for (Trigger* trig : photonTriggers) {
       if (trig->trigBool) nonMuonTrigger = true;
      }
      if (nonMuonTrigger)
       continue; // reject events where a non-electron trigger also fired.
      prescale = muonTrigger->trigPrescale;
     }
     else prescale = t->eventWeight;
     prescale = prescale / weight / weight;

     // Reco mumu invariant 4-momentum (best guess for Z boson)
     const TLorentzVector Z = muon1 + muon2;
     const double Z_pt = Z.Pt(); 
     const double Z_eta = Z.Eta();
     const double Z_phi = Z.Phi();
     const double Z_m = Z.M();

     // Z boson, dimuon cuts
     if (Z_pt < Z_pt_cut)
      continue; // pt cut on Z bosons
     if (t->muon_charge->at(m1) == t->muon_charge->at(m2))
      continue; // require oppositely charged muons

     // Fill dimuon mass spectrum
     zMassSpectra[0][numetabins]->Fill (Z_m, prescale);

     // Increment total Z counter
     nZmumuMass[numetabins]++;

     // Put the Z boson in the right eta bin
     short etabin = 0;
     if (etabins[0] < Z_eta ||
         Z_eta < etabins[numetabins]) {
      while (etabins[etabin] < Z_eta) etabin++;
     }
     etabin--;

     if (etabin != -1) {
      // Fill additional dimuon mass spectrum specific to opposite jet
      zMassSpectra[0][etabin]->Fill (Z_m, prescale);

      // Increment appropriate Z counter
      nZmumuMass[etabin]++;
     }

     // additional Z boson cuts
     if (Z_m < Z_mass - Z_mass_lower_cut || Z_mass + Z_mass_upper_cut < Z_m)
      continue; // cut on our sample Z boson mass

     // jet finding
     int lj = -1;
     int sj = -1;
     for (int j = 0; j < t->jet_n; j++) {
      if (t->jet_pt->at(j) < jet_pt_cut)
       continue; // basic jet pT cut
      if (!InHadCal (t->jet_eta->at(j), 0.4))
       continue; // require jets inside hadronic calorimeter
      if (DeltaR (t->muon_eta->at(m1), t->jet_eta->at(j), t->muon_phi->at(m1), t->jet_phi->at(j)) < 0.2 ||
          DeltaR (t->muon_eta->at(m2), t->jet_eta->at(j), t->muon_phi->at(m2), t->jet_phi->at(j)) < 0.2)
       continue; // require jets to be isolated from both muons
      if (DeltaPhi (t->jet_phi->at(j), Z_phi) < 3*pi/4)
       continue; // require jet to be back-to-back with Z in transverse plane

      // compare to leading jet
      else if (lj == -1 ||
               t->jet_pt->at(lj) < t->jet_pt->at(j)) {
       if (lj == -1 || !InDisabledHEC (t->jet_eta->at(lj), t->jet_phi->at(lj)))
        sj = lj;
       lj = j;
      }

      // compare to subleading jet
      else if ((sj == -1 ||
                t->jet_pt->at(sj) < t->jet_pt->at(j)) &&
               !InDisabledHEC (t->jet_eta->at(j), t->jet_phi->at(j))) {
       sj = j;
      }
     } // end jet finding loop
     if (lj == -1) // true iff no candidate jet is found
      continue; // reject on no candidate jet

     // relevant jet kinematic data
     const double ljet_pt = t->jet_pt->at(lj);
     const double ljet_eta = t->jet_eta->at(lj);
     const double ljet_phi = t->jet_phi->at(lj);
     const double sjet_pt = ((0 <= sj && sj < t->jet_n) ? t->jet_pt->at(sj) : 0);
     const double sjet_phi = ((0 <= sj && sj < t->jet_n) ? t->jet_phi->at(sj) : 0);

     // jet cuts
     if (InDisabledHEC (ljet_eta, ljet_phi))
      continue; // Reject event on additional HEC cuts
     //if (ljet_eta < etabins[0] ||
     //    ljet_eta > etabins[numetabins])
     // continue; // Make sure the jet is within the relevant bounds
     if (sjet_pt > 12) {
      const double subleading_dPhi = DeltaPhi (sjet_phi, Z_phi);
      //if (sjet_pt / Z_pt > 0.3)
      if (sjet_pt / (Z_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2)
      //if (sjet_pt * TMath::Cos(pi - subleading_dPhi) / Z_pt > 0.2)
       continue; // suppress dijets by requiring leading jet to dominate ptref
     }

     // Calculate opening angle in the transverse plane
     const double dPhi = DeltaPhi (ljet_phi, Z_phi);

     // Calculate systematics on jet pT
     const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError (ljet_pt, ljet_eta));

     // Calculate ptref and xjrefs
     const double ptref = Z_pt * TMath::Cos(pi - dPhi);
     const double ptratio[3] = {(ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

     // Fill dimuon xjref histograms
     //for (short errType = 0; errType < 3; errType++)
     // zmumuJetHists[errType][etabin]->Fill (Z_pt, ptratio[errType]/ptref, prescale);
     for (short errType = 0; errType < 3; errType++) {
      if (etabin != -1) zmumuJetHists[errType][etabin]->Fill (Z_pt, ptratio[errType], prescale);
      zmumuJetHists[errType][numetabins]->Fill (Z_pt, ptratio[errType], prescale);
     }

     // Increment Z+jet counters
     if (etabin != -1) nZmumuJet[etabin]++;
     nZmumuJet[numetabins]++;

     // Fill additional systematics histograms
     //for (short errType = 0; errType < 3; errType++)
     // zmumuJetHistsSys[errType][etabin]->Fill (ljet_pt, ptratio[errType]/ptref, prescale);
    }
   } // end loop over muon pairs
   // end Z->mumu type events


   /////////////////////////////////////////////////////////////////////////////
   // photon + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    TLorentzVector photon;
    const double this_photon_pt = t->photon_pt->at(p);
    const double this_photon_eta = t->photon_eta->at(p);
    const double this_photon_phi = t->photon_phi->at(p);
    photon.SetPtEtaPhiM (this_photon_pt, this_photon_eta, this_photon_phi, 0);

    // photon cuts
    if (this_photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at(p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at(p) > isolationEnergyIntercept + isolationEnergySlope*this_photon_pt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (this_photon_eta) || InDisabledHEC (this_photon_eta, this_photon_phi))
     continue; // require photon to be in EMCal

    // do photon truth matching to estimate photon energy scale
    if (isMC) {
     int truth_photon = (p < t->truth_photon_n ? p : t->truth_photon_n-1); // best initial guess
     double minDeltaR = 1000;
     for (int tp = 0; tp < t->truth_photon_n; tp++) {
      const double dR = DeltaR (this_photon_eta, t->truth_photon_eta->at(tp), this_photon_phi, t->truth_photon_phi->at(tp));
      if (dR < minDeltaR) {
       truth_photon = tp;
       minDeltaR = tp;
      }
     }
     if (minDeltaR >= 0.4)
      continue; // reco photons not matched to a truth photon within dR=0.4 are skipped

     // Put the jet in the right eta, pt bin
     short etabin = 0;
     if (etabins[0] < this_photon_eta ||
         this_photon_eta < etabins[numetabins]) {
      while (etabins[etabin] < this_photon_eta) etabin++;
     }
     etabin--;
     short pbin = 0;
     if (pbins[0] < this_photon_pt ||
         this_photon_pt < pbins[numpbins]) {
      while (pbins[pbin] < this_photon_pt) pbin++;
     }
     pbin--;

     const double per = this_photon_pt / t->truth_photon_pt->at(truth_photon);
     if (0 <= pbin && pbin < numpbins &&
         0 <= etabin && etabin < numetabins) {
      photonEnergyResponse[pbin][etabin]->Fill (per, t->eventWeight / weight / weight);
      nGamma[pbin][etabin]++;
     }
     if (0 <= pbin && pbin < numpbins) {
      photonEnergyResponse[pbin][numetabins]->Fill (per, t->eventWeight / weight / weight);
      nGamma[pbin][numetabins]++;
     }
     if (0 <= etabin && etabin < numetabins) {
      photonEnergyResponse[numpbins][etabin]->Fill (per, t->eventWeight / weight / weight);
      nGamma[numpbins][etabin]++;
     }
     photonEnergyResponse[numpbins][numetabins]->Fill (per, t->eventWeight / weight / weight);
     nGamma[numpbins][numetabins]++;
    }

    // triggering and event weighting
    float prescale = 1;
    if (!isMC) {
     Trigger* photonTrigger = NULL;
     for (Trigger* trig : photonTriggers) {
      if (trig->trigPrescale > 0 &&
          (photonTrigger == NULL ||
           (trig->trigPrescale < photonTrigger->trigPrescale &&
            trig->min_pt <= this_photon_pt &&
            this_photon_pt <= trig->max_pt)))
       photonTrigger = trig;
      //if (trig->trigPrescale > 0 &&
      //    trig->min_pt <= this_photon_pt &&
      //    this_photon_pt <= trig->max_pt)
      // photonTrigger = trig;
     }
     if (photonTrigger == NULL ||
         !photonTrigger->trigBool)
      continue;
     bool nonPhotonTrigger = false;
     for (Trigger* trig : electronTriggers) {
      if (trig->trigBool) nonPhotonTrigger = true;
     }
     for (Trigger* trig : muonTriggers) {
      if (trig->trigBool) nonPhotonTrigger = true;
     }
     if (nonPhotonTrigger)
      continue; // reject events where a non-electron trigger also fired.
     prescale = photonTrigger->trigPrescale;
    }
    else prescale = t->eventWeight;
    prescale = prescale / weight / weight;

    // Put the jet in the right eta bin
    short etabin = 0;
    if (etabins[0] < this_photon_eta ||
        this_photon_eta < etabins[numetabins]) {
     while (etabins[etabin] < this_photon_eta) etabin++;
    }
    etabin--;

    // jet finding
    int lj = -1; // allowed to be in HEC when finding
    int sj = -1; // not allowed to be in HEC when finding
    for (int j = 0; j < t->jet_n; j++) {
     // cuts on leading jet
     if (t->jet_pt->at(j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (!InHadCal (t->jet_eta->at(j), 0.4))
      continue; // require jets inside hadronic calorimeter
     double minDeltaR = 1000;
     for (int p2 = 0; p2 < t->photon_n; p2++) { // find minimum dR to a reco photon
      if (!t->photon_tight->at(p2))
       continue; // only compare with tight photons
      const double dR = DeltaR (t->photon_eta->at(p2), t->jet_eta->at(j), t->photon_phi->at(p2), t->jet_phi->at(j));
      if (dR < minDeltaR)
       minDeltaR = dR;
     }
     if (minDeltaR < 0.4)
      continue; // require jet is not reconstructed as a photon
     if (DeltaPhi (this_photon_phi, t->jet_phi->at(j)) < 3*pi/4)
      continue; // cut on gamma+jet samples not back-to-back in the transverse plane

     // compare to the leading jet
     else if (lj == -1 ||
              t->jet_pt->at(lj) < t->jet_pt->at(j)) {
      if (lj != -1 &&
          !InDisabledHEC (t->jet_eta->at(lj), t->jet_phi->at(lj)))
       sj = lj;
      lj = j;
     }

     // compare to the subleading jet
     else if ((sj == -1 ||
               t->jet_pt->at(sj) < t->jet_pt->at(j)) &&
              !InDisabledHEC (t->jet_eta->at(j), t->jet_phi->at(j))) {
      sj = j;
     }
    } // end jet finding loop
    if (lj == -1) // true iff there are no jets opposite photon
     continue; // reject on no jets

    // store relevant jet kinematics
    const double ljet_pt = t->jet_pt->at(lj);
    const double ljet_eta = t->jet_eta->at(lj);
    const double ljet_phi = t->jet_phi->at(lj);
    const double sjet_pt = ((0 <= sj && sj < t->jet_n) ? t->jet_pt->at(sj) : 0);
    const double sjet_phi = ((0 <= sj && sj < t->jet_n) ? t->jet_phi->at(sj) : 0);

    // cuts on jets
    if (InDisabledHEC (ljet_eta, ljet_phi))
     continue; // require jet to be outside of disabled HEC
    //if (ljet_eta < etabins[0] ||
    //    ljet_eta > etabins[numetabins])
    // continue; // Make sure jet is in eta bounds
    bool hasOtherJet = false;
    for (int j = 0; j < t->jet_n; j++) {
     if (j == lj)
      continue; // don't look at the leading jet, its our candidate :)
     if (DeltaR (t->jet_eta->at(j), this_photon_eta, t->jet_phi->at(j), this_photon_phi) < 0.4)
      continue; // reject jets that are just this photon
     if (t->jet_pt->at(j) < 12 || InDisabledHEC (t->jet_eta->at(j), t->jet_phi->at(j)))
      continue; // basic jet pT cut, also reject on the disabled HEC
     const double s_dphi = DeltaPhi (t->jet_phi->at(j), this_photon_phi);
     const double s_xjref = t->jet_pt->at(j) / (this_photon_pt);// * TMath::Cos(pi - s_dphi));
     if (0.1 < s_xjref) {
      hasOtherJet = true;
      break;
     }
    }
    if (hasOtherJet)
     continue; // cut on other jets that look back-to-back with gamma
    //if (sjet_pt > 12) {
    // const double subleading_dPhi = DeltaPhi (sjet_phi, this_photon_phi);
    // //if (sjet_pt / this_photon_pt > 0.1)
    // if (sjet_pt / (this_photon_pt * TMath::Cos(pi - subleading_dPhi)) > 0.2)
    // //if (sjet_pt * TMath::Cos(pi - subleading_dPhi) / this_photon_pt > 0.2)
    // //if (sjet_pt / this_photon_pt > 0.02)
    //  continue; // suppress dijets by requiring leading jet to dominate ptref
    //}

    // Calculate opening angle in the transverse plane
    const double dPhi = DeltaPhi (ljet_phi, this_photon_phi);

    // Calculate systematics
    const double ljet_pt_err = (isMC ? 0:GetXCalibSystematicError(ljet_pt, ljet_eta));

    // Calculate ptref and xjrefs
    const double ptref = this_photon_pt * TMath::Cos(pi - dPhi);
    const double ptratio[3] = {(ljet_pt-ljet_pt_err)/ptref, ljet_pt/ptref, (ljet_pt+ljet_pt_err)/ptref};

    // Fill xjref histograms
    for (short errType = 0; errType < 3; errType++) {
     if (etabin != -1) gJetHists[errType][etabin]->Fill (this_photon_pt, ptratio[errType], prescale);
     gJetHists[errType][numetabins]->Fill (this_photon_pt, ptratio[errType], prescale);
    }

    // if data, calculate additional systematics for this jet
    if (!isMC) {
     const double newJetPtSys = GetNewXCalibSystematicError (ljet_eta, this_photon_pt) / ljet_pt;
     if (etabin != -1) gJetHistsSys[1][etabin]->Fill (ljet_pt, newJetPtSys, prescale);
     gJetHistsSys[1][numetabins]->Fill (ljet_pt, newJetPtSys, prescale);
    }

    // Increment gamma+jet counters
    if (etabin != -1)
     nGammaJet[etabin]++;
    nGammaJet[numetabins]++;

    //if (90 < ptref && ptref < 110 && 0.6 < ptratio[1] && ptratio[1] < 0.8) {
    // cout << "Found a weird event! Dataset: " << identifier.Data() << ", entry: " << entry << endl;
    //}
   }
    
  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  dataOverMCFile->Close();
  if (dataOverMCFile) delete dataOverMCFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form("%sdataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");


  // Write histograms to output and clean memory
  for (short etabin = 0; etabin <= numetabins; etabin++) {
   for (short errType = 0; errType < 3; errType++) {
    zmumuJetHists[errType][etabin]->Write();
    if (zmumuJetHists[errType][etabin]) delete zmumuJetHists[errType][etabin];
    //zmumuJetHistsSys[errType][etabin]->Write();
    //if (zmumuJetHistsSys[errType][etabin]) delete zmumuJetHistsSys[errType][etabin];

    zeeJetHists[errType][etabin]->Write();
    if (zeeJetHists[errType][etabin]) delete zeeJetHists[errType][etabin];
    //zeeJetHistsSys[errType][etabin]->Write();
    //if (zeeJetHistsSys[errType][etabin]) delete zeeJetHistsSys[errType][etabin];

    gJetHists[errType][etabin]->Write();
    if (gJetHists[errType][etabin]) delete gJetHists[errType][etabin];
   }

   gJetHistsSys[1][etabin]->Write();
   if (gJetHistsSys[1][etabin]) delete gJetHistsSys[1][etabin];

   for (short spcType = 0; spcType < 2; spcType++) {
    zMassSpectra[spcType][etabin]->Write();
    if (zMassSpectra[spcType][etabin]) delete zMassSpectra[spcType][etabin];
   }
  }
  electronEnergyScale->Write();
  if (electronEnergyScale) delete electronEnergyScale;
  lowResponseEtaPhi->Write();
  if (lowResponseEtaPhi) delete lowResponseEtaPhi;
  highResponseEtaPhi->Write();
  if (highResponseEtaPhi) delete highResponseEtaPhi;
  lowResponseEtaPt->Write();
  if (lowResponseEtaPt) delete lowResponseEtaPt;
  highResponseEtaPt->Write();
  if (highResponseEtaPt) delete highResponseEtaPt;
  responsePt->Write();
  if (responsePt) delete responsePt;
  jetSpectrum->Write();
  if (jetSpectrum) delete jetSpectrum;

  for (short pbin = 0; pbin <= numpbins; pbin++) {
   for (short etabin = 0; etabin <= numetabins; etabin++) {
    jetEnergyResponseCalib[pbin][etabin]->Write();
    if (jetEnergyResponseCalib[pbin][etabin]) delete jetEnergyResponseCalib[pbin][etabin];
    jetEnergyResponseReco[pbin][etabin]->Write();
    if (jetEnergyResponseReco[pbin][etabin]) delete jetEnergyResponseReco[pbin][etabin];
    photonEnergyResponse[pbin][etabin]->Write();
    if (photonEnergyResponse[pbin][etabin]) delete photonEnergyResponse[pbin][etabin];
   }
  }

  TVectorD infoVec(2);
  infoVec[0] = luminosity;
  infoVec[1] = dataSet;
  infoVec.Write(Form("infoVec_%s", identifier.Data()));

  TVectorD nJetVec((numpbins+1)*(numetabins+1));
  TVectorD nGammaVec((numpbins+1)*(numetabins+1));
  TVectorD nZeeMassVec(numetabins+1);
  TVectorD nZeeJetVec(numetabins+1);
  TVectorD nZmumuMassVec(numetabins+1);
  TVectorD nZmumuJetVec(numetabins+1);
  TVectorD nGammaJetVec(numetabins+1);
  for (short etabin = 0; etabin <= numetabins; etabin++) {
   for (short pbin = 0; pbin <= numpbins; pbin++) {
    nJetVec[pbin + (numpbins+1)*etabin] = (double)nJet[pbin][etabin];
    nGammaVec[pbin + (numpbins+1)*etabin] = (double)nGamma[pbin][etabin];
   }
   nZeeMassVec[etabin] = (double)nZeeMass[etabin];
   nZeeJetVec[etabin] = (double)nZeeJet[etabin];
   nZmumuMassVec[etabin] = (double)nZmumuMass[etabin];
   nZmumuJetVec[etabin] = (double)nZmumuJet[etabin];
   nGammaJetVec[etabin] = (double)nGammaJet[etabin];
  }
  nJetVec.Write(Form("nJetVec_%s", identifier.Data()));
  nGammaVec.Write(Form("nGammaVec_%s", identifier.Data()));
  nZeeMassVec.Write(Form("nZeeMassVec_%s", identifier.Data()));
  nZeeJetVec.Write(Form("nZeeJetVec_%s", identifier.Data()));
  nZmumuMassVec.Write(Form("nZmumuMassVec_%s", identifier.Data()));
  nZmumuJetVec.Write(Form("nZmumuJetVec_%s", identifier.Data()));
  nGammaJetVec.Write(Form("nGammaJetVec_%s", identifier.Data()));

  outFile->Close();
  if (outFile) delete outFile;
  return;
}
