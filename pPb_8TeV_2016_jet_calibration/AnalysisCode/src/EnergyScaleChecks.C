#include "EnergyScaleChecks.h"

#include <ArrayTemplates.h>

#include <TFile.h>
#include <TSystemDirectory.h>
#include <TH2D.h>
#include <TVectorT.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

TFile* xCalibSystematicsFile = NULL;


double GetXCalibSystematicError(const double jpt, const double jeta) {
  TFile* file = xCalibSystematicsFile;
  if (!file || !file->IsOpen())
   return 0;

  if (TMath::Abs(jeta) < xcalibEtabins[0] ||
      xcalibEtabins[sizeof(xcalibEtabins)/sizeof(xcalibEtabins[0]) -1] < TMath::Abs(jeta)) {
   return 0;
  }

  short iEta = 0;
  while (xcalibEtabins[iEta] < TMath::Abs(jeta)) iEta++;
  iEta--;

  const TString hname = TString("fsys_rel_") + Form ("%i", iEta);
  TH1D* fsys_rel = (TH1D*)file->Get(hname.Data());

  return TMath::Abs(fsys_rel->GetBinContent(fsys_rel->FindBin(jpt)) - 1) * jpt;
}


TString GetIdentifier (const int dataSet, const bool isPeriodA, const TString inFileName) {

  TString id = TString (isPeriodA ? "pPb_" : "Pbp_") + TString (inFileName.Contains ("signalonly") ? "Signal" : "Overlay");

  if (inFileName.Contains ("valid")) id = id + "_Valid";

  if (inFileName.Contains ("jetjet")) { // dijet
   if (dataSet <= 0) return "";
   id = id + "_Dijet_Slice" + to_string (dataSet);
  }
  else if (inFileName.Contains ("42310") && inFileName.Contains ("Slice")) { // gamma+jet
   if (dataSet <= 0) return "";
   id = id + "_GammaJet_Slice" + to_string (dataSet);
  }
  else if (inFileName.Contains ("ZeeJet")) { // Zee+jet
   if (dataSet < 0) return "";
   id = id + "_ZeeJet" + (dataSet == 0 ? "" : "_Slice" + to_string (dataSet));
  }
  else if (inFileName.Contains ("ZmumuJet")) { // Zmumu+jet
   if (dataSet != 0) return "";
   id = id + "_ZmumuJet";
  }

  return id;
}


void EnergyScaleChecks (const int dataSet,
                        const bool isPeriodA,
                        const TString inFileName)
{

  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  const TString identifier = GetIdentifier(dataSet, isPeriodA, inFileName);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier = inFileName;
   if (fileIdentifier == "") {
    cout << "File name not provided! Exiting." << endl;
    return;
   }

   const TString dataPathTemp = dataPath;// + "/rtrk_data/";
   TSystemDirectory dir(dataPathTemp.Data(), dataPathTemp.Data());
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
     if (debugStatements) cout << "Status: In EnergyScaleChecks.C (breakpoint B): Found " << fname.Data() << endl;
     
     if (fname.Contains(fileIdentifier)) {
      file = new TFile(dataPathTemp+fname, "READ");
      tree = (TTree*)file->Get("tree");
      break;
     }
    }
   }
  }
  if (tree == NULL || file == NULL) {
   cout << "Error: In EnergyScaleChecks.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables(tree, true);
  t->SetGetVertices ();
  t->SetGetHIJets ();
  t->SetGetEMTopoJets ();
  t->SetGetElectrons ();
  t->SetGetPhotons ();
  t->SetBranchAddresses();

  TH1D* electronEnergyScale = new TH1D (Form ("electronEnergyScale_dataSet%s", identifier.Data()), "", 200, 0, 2);
  electronEnergyScale->Sumw2();

  TH1D**** jetEnergyResponseCalib = Get3DArray <TH1D*> (2, numpbins+1, numetabins+1);
  TH1D**** jetEnergyResponseReco = Get3DArray <TH1D*> (2, numpbins+1, numetabins+1);
  TH1D*** photonEnergyResponse = Get2DArray <TH1D*> (numpbins+1, numetabins+1);

  for (short iP = 0; iP <= numpbins; iP++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     jetEnergyResponseCalib[iAlgo][iP][iEta] = new TH1D (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), identifier.Data(), iP, iEta), "", 100, 0, 2);
     jetEnergyResponseCalib[iAlgo][iP][iEta]->Sumw2();
     jetEnergyResponseReco[iAlgo][iP][iEta] = new TH1D (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), identifier.Data(), iP, iEta), "", 100, 0, 2);
     jetEnergyResponseReco[iAlgo][iP][iEta]->Sumw2();
    }
    photonEnergyResponse[iP][iEta] = new TH1D (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", identifier.Data(), iP, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 100, 0, 2);
    photonEnergyResponse[iP][iEta]->Sumw2();
   }
  }

  int*** nJet = Get3DArray <int> (2, numpbins+1, numetabins+1);
  int** nGamma = Get2DArray <int> (numpbins+1, numetabins+1);

  int* jet_n;
  vector<float>* jet_pt;
  vector<float>* jet_eta;
  vector<float>* jet_phi;
  vector<float>* jet_e;
  vector<float>* precalib_jet_pt;
  vector<float>* precalib_jet_e;

  xCalibSystematicsFile = new TFile(rootPath + "cc_sys_090816.root", "READ");

  const long long numEntries = tree->GetEntries();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (long long entry = 0; entry < numEntries; entry++) {
   tree->GetEntry(entry);

   const double evtWeight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;
   //const double evtWeight = 1.2910E+03 / 0.0056462 / 3997692; // for dijet signal only sample
   //const double evtWeight = 1; // since AMI is down currently


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if ((t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at(0) != 1)) continue;

   // basically just do the jet part twice, once for each algorithm
   for (short iAlgo = 0; iAlgo < 2; iAlgo++) {

    if (iAlgo == 0) { // HI algorithm at the EM scale
     jet_n = &(t->akt4hi_jet_n);
     jet_pt = t->akt4hi_em_xcalib_jet_pt;
     jet_eta = t->akt4hi_em_xcalib_jet_eta;
     jet_phi = t->akt4hi_em_xcalib_jet_phi;
     jet_e = t->akt4hi_em_xcalib_jet_e;
     precalib_jet_pt = t->akt4hi_em_jet_pt;
     precalib_jet_e = t->akt4hi_em_jet_e;
    }
    else { // EMTopo algorithm at the EM scale
     jet_n = &(t->akt4emtopo_jet_n);
     jet_pt = t->akt4emtopo_calib_jet_pt;
     jet_eta = t->akt4emtopo_calib_jet_eta;
     jet_phi = t->akt4emtopo_calib_jet_phi;
     jet_e = t->akt4emtopo_calib_jet_e;
     precalib_jet_pt = t->akt4emtopo_em_jet_pt;
     precalib_jet_e = t->akt4emtopo_em_jet_e;
    }


    /////////////////////////////////////////////////////////////////////////////
    // fill jet energy response
    /////////////////////////////////////////////////////////////////////////////
    for (int j = 0; j < *jet_n; j++) {
     if (jet_pt->at(j) < jet_pt_cut)
      continue; // basic jet pT cut
     if (InDisabledHEC (jet_eta->at(j), jet_phi->at(j)))
      continue; // only truth match jets outside disabled HEC
     if (!InHadCal (jet_eta->at(j)))
      continue; // reject jets reconstructed outside reasonable HCal bounds.

     double minDeltaR = 1000;
     // loop over truth electrons and photons
     for (int e = 0; e < t->truth_electron_n; e++) { // truth electrons
      const double dR = DeltaR (jet_eta->at(j), t->truth_electron_eta->at(e), jet_phi->at(j), t->truth_electron_phi->at(e));
      if (dR < minDeltaR) minDeltaR = dR;
     }
     for (int p = 0; p < t->truth_photon_n; p++) { // truth photons
      const double dR = DeltaR (jet_eta->at(j), t->truth_photon_eta->at(p), jet_phi->at(j), t->truth_photon_phi->at(p));
      if (dR < minDeltaR) minDeltaR = dR;
     }
     if (minDeltaR < 0.6)
      continue; // reject jets close to some other lepton or photon

     minDeltaR = 1000;
     int truth_jet = -1;
     for (int tj = 0; tj < t->truth_jet_n; tj++) {
      double dR = DeltaR (jet_eta->at(j), t->truth_jet_eta->at(tj), jet_phi->at(j), t->truth_jet_phi->at(tj));
      if (dR < minDeltaR) {
       minDeltaR = dR;
       truth_jet = tj;
      }
     }
     if (0 <= truth_jet && truth_jet < t->truth_jet_n && minDeltaR < 0.4) {

      // Put the jet in the right eta bin
      short iEta = 0;
      // Make sure jet is in eta bounds
      if (etabins[0] < jet_eta->at(j) ||
          jet_eta->at(j) < etabins[numetabins]) {
       while (etabins[iEta] < jet_eta->at(j)) iEta++;
      }
      iEta--;
      short iP = 0;
      double jer, pjer;
      if (!calcPtClosure) {
       if (pbins[0] < t->truth_jet_e->at(truth_jet) ||
           t->truth_jet_e->at(truth_jet) < pbins[numpbins]) {
        while (pbins[iP] < t->truth_jet_e->at(truth_jet)) iP++;
       }
       jer = jet_e->at(j) / t->truth_jet_e->at(truth_jet);
       pjer = precalib_jet_e->at(j) / t->truth_jet_e->at(truth_jet);
      }
      else if (pbins[0] < t->truth_jet_pt->at(truth_jet) ||
               t->truth_jet_pt->at(truth_jet) < pbins[numpbins]) {
       while (pbins[iP] < t->truth_jet_pt->at(truth_jet)) iP++;
       jer = jet_pt->at(j) / t->truth_jet_pt->at(truth_jet);
       pjer = precalib_jet_pt->at(j) / t->truth_jet_pt->at(truth_jet);
      }
      else continue;
      iP--;

      if (0 <= iP && iP < numpbins &&
          0 <= iEta && iEta < numetabins) {
       nJet[iAlgo][iP][iEta]++;
       jetEnergyResponseCalib[iAlgo][iP][iEta]->Fill (jer, evtWeight);
       jetEnergyResponseReco[iAlgo][iP][iEta]->Fill (pjer, evtWeight);
      }
      if (0 <= iP && iP < numpbins) {
       nJet[iAlgo][iP][numetabins]++;
       jetEnergyResponseCalib[iAlgo][iP][numetabins]->Fill (jer, evtWeight);
       jetEnergyResponseReco[iAlgo][iP][numetabins]->Fill (pjer, evtWeight);
      }
      if (0 <= iEta && iEta < numetabins) {
       nJet[iAlgo][numpbins][iEta]++;
       jetEnergyResponseCalib[iAlgo][numpbins][iEta]->Fill (jer, evtWeight);
       jetEnergyResponseReco[iAlgo][numpbins][iEta]->Fill (pjer, evtWeight);
      }
      nJet[iAlgo][numpbins][numetabins]++;
      jetEnergyResponseCalib[iAlgo][numpbins][numetabins]->Fill (jer, evtWeight);
      jetEnergyResponseReco[iAlgo][numpbins][numetabins]->Fill (pjer, evtWeight);
     }
    }
   } // end jet algorithm loop


   /////////////////////////////////////////////////////////////////////////////
   // electron containing events
   /////////////////////////////////////////////////////////////////////////////
   for (int e = 0; e < t->electron_n; e++) { // loop over primary electron

    // electron cuts
    if (t->electron_pt->at(e) < electron_pt_cut)
     continue; // basic electron pT cuts
    if (!InEMCal (t->electron_eta->at(e)))
     continue; // reject electrons reconstructed outside EMCal
    if (!t->electron_loose->at(e))
     continue; // reject non-loose electrons
    if (t->electron_d0sig->at(e) > 5)
     continue; // d0 (transverse impact parameter) significance cut
    if (t->electron_delta_z0_sin_theta->at(e) > 0.5)
     continue; // z0 (longitudinal impact parameter) vertex compatibility cut

    // fill the electron energy scale plot first, finding the truth-reco pair
    if (t->truth_electron_n > 0) {
     int closestEl = (e < t->truth_electron_n ? e:t->truth_electron_n-1); // best initial guess
     
     double minDeltaR = DeltaR (t->electron_eta->at(e), t->truth_electron_eta->at(closestEl), t->electron_phi->at(e), t->truth_electron_phi->at(closestEl));
     for (int te = 0; te < t->truth_electron_n; te++) {
      if (t->truth_electron_pt->at(te) < electron_pt_cut) continue;
      if (t->electron_charge->at(e) != t->truth_electron_charge->at(te)) continue;
      const double dR = DeltaR (t->electron_eta->at(e), t->truth_electron_eta->at(te), t->electron_phi->at(e), t->truth_electron_phi->at(te));
      if (dR < minDeltaR) {
       closestEl = te;
       minDeltaR = dR;
      }
     }
     if (minDeltaR < 0.4)
      electronEnergyScale->Fill (t->electron_pt->at(e) / t->truth_electron_pt->at(closestEl));//, t->eventWeight);
    }
   }


   /////////////////////////////////////////////////////////////////////////////
   // photon containing events
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
    short iEta = 0;
    if (etabins[0] < photon_eta &&
        photon_eta < etabins[numetabins]) {
     while (etabins[iEta] < photon_eta) iEta++;
    }
    iEta--;
    short iP = 0;
    if (pbins[0] < photon_pt &&
        photon_pt < pbins[numpbins]) {
     while (pbins[iP] < photon_pt) iP++;
    }
    iP--;

    const double per = photon_pt / t->truth_photon_pt->at(truth_photon);
    if (0 <= iP && iP < numpbins &&
        0 <= iEta && iEta < numetabins) {
     photonEnergyResponse[iP][iEta]->Fill (per, evtWeight);
     nGamma[iP][iEta]++;
    }
    if (0 <= iP && iP < numpbins) {
     photonEnergyResponse[iP][numetabins]->Fill (per, evtWeight);
     nGamma[iP][numetabins]++;
    }
    if (0 <= iEta && iEta < numetabins) {
     photonEnergyResponse[numpbins][iEta]->Fill (per, evtWeight);
     nGamma[numpbins][iEta]++;
    }
    photonEnergyResponse[numpbins][numetabins]->Fill (per, evtWeight);
    nGamma[numpbins][numetabins]++;
   }

  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form ("%s/EnergyScaleChecks/dataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");


  // Write histograms to output and clean memory
  electronEnergyScale->Write();
  if (electronEnergyScale) delete electronEnergyScale;

  for (short iP = 0; iP <= numpbins; iP++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     jetEnergyResponseCalib[iAlgo][iP][iEta]->Write();
     jetEnergyResponseReco[iAlgo][iP][iEta]->Write();
    }
    photonEnergyResponse[iP][iEta]->Write();
   }
  }

  Delete3DArray (jetEnergyResponseCalib, 2, numpbins+1, numetabins+1);
  Delete3DArray (jetEnergyResponseReco, 2, numpbins+1, numetabins+1);
  Delete2DArray (photonEnergyResponse, numpbins+1, numetabins+1);

  TVectorD nJetVec(2*(numpbins+1)*(numetabins+1));
  TVectorD nGammaVec((numpbins+1)*(numetabins+1));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (short iP = 0; iP <= numpbins; iP++) {
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     nJetVec[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)] = (double)nJet[iAlgo][iP][iEta];
    }
    nGammaVec[iP + (numpbins+1)*iEta] = (double)nGamma[iP][iEta];
   }
  }
  nJetVec.Write(Form ("nJetVec_%s", identifier.Data()));
  nGammaVec.Write(Form ("nGammaVec_%s", identifier.Data()));

  Delete3DArray (nJet, 2, numpbins+1, numetabins+1);
  Delete2DArray (nGamma, numpbins+1, numetabins+1);

  outFile->Close();
  if (outFile) delete outFile;
  return;
}

} // end namespace
