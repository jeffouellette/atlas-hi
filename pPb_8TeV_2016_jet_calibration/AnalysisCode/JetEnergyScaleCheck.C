#include "JetEnergyScaleCheck.h"

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


TString GetIdentifier (const int dataSet, const bool isValidationSample, const bool periodA) {
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


void JetEnergyScaleCheck (const int dataSet,
                          const bool isPeriodAflag,
                          const TString inFileName)
{

  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  isPeriodA = isPeriodAflag;

  const bool isValidationSample = TString(inFileName).Contains("valid");
  const TString identifier = GetIdentifier(dataSet, isMC, isValidationSample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   if (inFileName == "") {
    fileIdentifier = TString(dataSet > 0 ? ("Slice" + to_string(dataSet)) : (dataSet==0 ? "ZmumuJet" : ("ZeeJet"+to_string(-dataSet)))) + (isMCperiodAflag ? ".pPb":".Pbp");
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

  TreeVariables* t = new TreeVariables(tree, true);
  t->SetBranchAddresses(true);

  // initialize histograms
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

  xCalibSystematicsFile = new TFile(rootPath + "cc_sys_090816.root", "READ");

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

     jetSpectrum->Fill (t->jet_pt->at(j), evtWeight);

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
       jetEnergyResponseCalib[pbin][etabin]->Fill (jer, evtWeight);
       jetEnergyResponseReco[pbin][etabin]->Fill (pjer, evtWeight);
      }
      if (0 <= pbin && pbin < numpbins) {
       nJet[pbin][numetabins]++;
       jetEnergyResponseCalib[pbin][numetabins]->Fill (jer, evtWeight);
       jetEnergyResponseReco[pbin][numetabins]->Fill (pjer, evtWeight);
      }
      if (0 <= etabin && etabin < numetabins) {
       nJet[numpbins][etabin]++;
       jetEnergyResponseCalib[numpbins][etabin]->Fill (jer, evtWeight);
       jetEnergyResponseReco[numpbins][etabin]->Fill (pjer, evtWeight);
      }
      nJet[numpbins][numetabins]++;
      jetEnergyResponseCalib[numpbins][numetabins]->Fill (jer, evtWeight);
      jetEnergyResponseReco[numpbins][numetabins]->Fill (pjer, evtWeight);

      // Extra plots
      if (t->jet_pt->at(j) / t->truth_jet_pt->at(truth_jet) < 1.18) {
       lowResponseEtaPhi->Fill (t->jet_eta->at(j), t->jet_phi->at(j), evtWeight);
       lowResponseEtaPt->Fill (t->jet_eta->at(j), t->jet_pt->at(j), evtWeight);
      }
      else {
       highResponseEtaPhi->Fill (t->jet_eta->at(j), t->jet_phi->at(j), evtWeight);
       highResponseEtaPt->Fill (t->jet_eta->at(j), t->jet_pt->at(j), evtWeight);
      }
      responsePt->Fill (jer, t->jet_pt->at(j), evtWeight);

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


   /////////////////////////////////////////////////////////////////////////////
   // photon + jet type events
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    const double photon_pt = t->photon_pt->at(p);
    const double photon_eta = t->photon_eta->at(p);
    const double photon_phi = t->photon_phi->at(p);

    // photon cuts
    if (photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at(p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at(p) > isolationEnergyIntercept + isolationEnergySlope*photon_pt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (photon_eta) || InDisabledHEC (photon_eta, photon_phi))
     continue; // require photon to be in EMCal

    // do photon truth matching to estimate photon energy scale
    if (isMC) {
     int truth_photon = (p < t->truth_photon_n ? p : t->truth_photon_n-1); // best initial guess
     double minDeltaR = 1000;
     for (int tp = 0; tp < t->truth_photon_n; tp++) {
      const double dR = DeltaR (photon_eta, t->truth_photon_eta->at(tp), photon_phi, t->truth_photon_phi->at(tp));
      if (dR < minDeltaR) {
       truth_photon = tp;
       minDeltaR = tp;
      }
     }
     if (minDeltaR >= 0.4)
      continue; // reco photons not matched to a truth photon within dR=0.4 are skipped

     // Put the photon in the right eta, pt bin
     short etabin = 0;
     if (etabins[0] < photon_eta &&
         photon_eta < etabins[numetabins]) {
      while (etabins[etabin] < photon_eta) etabin++;
     }
     etabin--;
     short pbin = 0;
     if (pbins[0] < photon_pt &&
         photon_pt < pbins[numpbins]) {
      while (pbins[pbin] < photon_pt) pbin++;
     }
     pbin--;

     const double per = photon_pt / t->truth_photon_pt->at(truth_photon);
     if (0 <= pbin && pbin < numpbins &&
         0 <= etabin && etabin < numetabins) {
      photonEnergyResponse[pbin][etabin]->Fill (per, evtWeight);
      nGamma[pbin][etabin]++;
     }
     if (0 <= pbin && pbin < numpbins) {
      photonEnergyResponse[pbin][numetabins]->Fill (per, evtWeight);
      nGamma[pbin][numetabins]++;
     }
     if (0 <= etabin && etabin < numetabins) {
      photonEnergyResponse[numpbins][etabin]->Fill (per, evtWeight);
      nGamma[numpbins][etabin]++;
     }
     photonEnergyResponse[numpbins][numetabins]->Fill (per, evtWeight);
     nGamma[numpbins][numetabins]++;
    }
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

  const char* outFileName = Form("%s/JetEnergyScaleCheck/dataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");


  // Write histograms to output and clean memory
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
  for (short etabin = 0; etabin <= numetabins; etabin++) {
   for (short pbin = 0; pbin <= numpbins; pbin++) {
    nJetVec[pbin + (numpbins+1)*etabin] = (double)nJet[pbin][etabin];
    nGammaVec[pbin + (numpbins+1)*etabin] = (double)nGamma[pbin][etabin];
   }
  }
  nJetVec.Write(Form("nJetVec_%s", identifier.Data()));
  nGammaVec.Write(Form("nGammaVec_%s", identifier.Data()));

  outFile->Close();
  if (outFile) delete outFile;
  return;
}
