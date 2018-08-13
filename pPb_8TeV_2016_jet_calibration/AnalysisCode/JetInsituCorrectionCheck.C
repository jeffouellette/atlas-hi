#include "JetInsituCorrectionCheck.h"

TFile* xCalibSystematicsFile = NULL;

bool isPeriodA = false;


double GetXCalibSystematicError(const double jpt, const double jeta) {
  TFile* file = xCalibSystematicsFile;
  if (!file || !file->IsOpen()) {
   cout << "Cannot find xCalib systematics file!" << endl;
   return 0;
  }

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


void JetInsituCorrectionCheck (const int dataSet,
                               const double luminosity = 0, 
                               const bool periodA = false)
{
  isPeriodA = periodA;
  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  const TString identifier = to_string(dataSet);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   const TString fileIdentifier = to_string(dataSet);

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
     if (debugStatements) cout << "Status: In JetInsituCorrectionCheck.C (breakpoint B): Found " << fname.Data() << endl;
     
     if (fname.Contains(fileIdentifier)) {
      file = new TFile(dataPath+fname, "READ");
      tree = (TTree*)file->Get("tree");
      break;
     }
    }
   }
  }
  if (tree == NULL || file == NULL) {
   cout << "Error: In JetInsituCorrectionCheck.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables(tree, true);
  t->SetBranchAddresses(false);

  // initialize histograms
  TH2D* jetPreInsituSpectrum = new TH2D (Form ("jetPreInsituSpectrum_dataSet%s", identifier.Data()), ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPreInsituSpectrum->Sumw2();
  TH2D* jetPostInsituSpectrum = new TH2D (Form ("jetPostInsituSpectrum_dataSet%s", identifier.Data()), ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPostInsituSpectrum->Sumw2();
  TH2D* jetPostInsituSpectrumSysHi = new TH2D (Form ("jetPostInsituSpectrumSysHi_dataSet%s", identifier.Data()), ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPostInsituSpectrumSysHi->Sumw2();
  TH2D* jetPostInsituSpectrumSysLo = new TH2D (Form ("jetPostInsituSpectrumSysLo_dataSet%s", identifier.Data()), ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPostInsituSpectrumSysLo->Sumw2();

  TH2D* jetInsituResponse[numpbins];
  TH2D* jetInsituResponseSysHi[numpbins];
  TH2D* jetInsituResponseSysLo[numpbins];

  for (int pbin = 0; pbin < numpbins; pbin++) {

   jetInsituResponse[pbin] = new TH2D (Form ("jetInsituResponse_dataSet%s_pbin%i", identifier.Data(), pbin), ";#it{p}_{T}^{Calib} / #it{p}_{T}^{EtaJES};#eta;", 200, 0.9, 1.1, numetabins, etabins);
   jetInsituResponse[pbin]->Sumw2();
   jetInsituResponseSysHi[pbin] = new TH2D (Form ("jetInsituResponseSysHi_dataSet%s_pbin%i", identifier.Data(), pbin), ";#it{p}_{T}^{Calib} / #it{p}_{T}^{EtaJES};#eta;", 200, 0.9, 1.1, numetabins, etabins);
   jetInsituResponseSysHi[pbin]->Sumw2();
   jetInsituResponseSysLo[pbin] = new TH2D (Form ("jetInsituResponseSysLo_dataSet%s_pbin%i", identifier.Data(), pbin), ";#it{p}_{T}^{Calib} / #it{p}_{T}^{EtaJES};#eta;", 200, 0.9, 1.1, numetabins, etabins);
   jetInsituResponseSysLo[pbin]->Sumw2();
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
   // fill jet energy response from insitu corrections
   /////////////////////////////////////////////////////////////////////////////
   for (int j = 0; j < t->jet_n; j++) {

    if (t->em_xcalib_jet_pt->at(j) < jet_pt_cut)
     continue; // basic jet pT cut
    if (InDisabledHEC (t->em_xcalib_jet_eta->at(j), t->em_xcalib_jet_phi->at(j)))
     continue; // only use jets outside disabled HEC
    if (!InHadCal (t->em_xcalib_jet_eta->at(j)))
     continue; // reject jets reconstructed outside reasonable HCal bounds.

    double minDeltaR = 1000;
    for (int p = 0; p < t->photon_n; p++) {
     // photon cuts
     //if (t->photon_pt->at(p) < photon_pt_cut)
     // continue; // basic pT cut on photons
     //if (!t->photon_tight->at(p))
     // continue; // require tight cuts on photons
     //if (t->photon_topoetcone40->at(p) > isolationEnergyIntercept + isolationEnergySlope*t->photon_pt->at(p))
     // continue; // require maximum isolation energy on gammas
     //if (!InEMCal (t->photon_eta->at(p)) || InDisabledHEC (t->photon_eta->at(p), t->photon_phi->at(p)))
     // continue; // require photon to be in EMCal

     const double dR = DeltaR (t->em_xcalib_jet_eta->at(j), t->photon_eta->at(p), t->em_xcalib_jet_phi->at(j), t->photon_phi->at(p));
     if (dR < minDeltaR) minDeltaR = dR;
    }
    for (int e = 0; e < t->electron_n; e++) {
     // electron cuts
     //if (t->electron_pt->at(e) < electron_pt_cut)
     // continue; // basic electron pT cuts
     //if (!InEMCal (t->electron_eta->at(e)))
     // continue; // reject electrons reconstructed outside EMCal
     //if (!t->electron_loose->at(e))
     // continue; // reject non-loose electrons
     //if (t->electron_d0sig->at(e) > 5)
     // continue; // d0 (transverse impact parameter) significance cut
     //if (t->electron_delta_z0_sin_theta->at(e) > 0.5)
     // continue; // z0 (longitudinal impact parameter) vertex compatibility cut

     const double dR = DeltaR (t->em_xcalib_jet_eta->at(j), t->electron_eta->at(e), t->em_xcalib_jet_phi->at(j), t->electron_phi->at(e));
     if (dR < minDeltaR) minDeltaR = dR;
    }

    if (minDeltaR < 0.4)
     continue; // reject jets reconstructed next to a photon or electron

    const double insituSys = GetXCalibSystematicError (t->em_xcalib_jet_pt->at(j), t->em_xcalib_jet_eta->at(j));

    const double etaPlot = t->em_xcalib_jet_eta->at(j);//(isPeriodA ? -t->em_xcalib_jet_eta->at(j) : t->em_xcalib_jet_eta->at(j));

    jetPreInsituSpectrum->Fill (t->em_etajes_jet_pt->at(j), etaPlot);
    jetPostInsituSpectrum->Fill (t->em_xcalib_jet_pt->at(j), etaPlot);
    jetPostInsituSpectrumSysHi->Fill (t->em_xcalib_jet_pt->at(j) + insituSys, etaPlot);
    jetPostInsituSpectrumSysLo->Fill (t->em_xcalib_jet_pt->at(j) - insituSys, etaPlot);

    short pbin = 0;
    while (pbin < numpbins && pbins[pbin] < t->em_etajes_jet_pt->at(j)) pbin++;
    pbin--;

    if (0 <= pbin && pbin < numpbins) {
     const double ratio = t->em_xcalib_jet_pt->at(j) / t->em_etajes_jet_pt->at(j);
     const double ratiohi = (t->em_xcalib_jet_pt->at(j) + insituSys) / t->em_etajes_jet_pt->at(j);
     const double ratiolo = (t->em_xcalib_jet_pt->at(j) - insituSys) / t->em_etajes_jet_pt->at(j);
     jetInsituResponse[pbin]->Fill (ratio, etaPlot);
     jetInsituResponseSysHi[pbin]->Fill (ratiohi, etaPlot);
     jetInsituResponseSysLo[pbin]->Fill (ratiolo, etaPlot);
     if (-1.3 < etaPlot && etaPlot < -0.5 &&
         0.94 < ratio && ratio < 0.97) {
      cout << "Found a weird event! Run: " << dataSet << ", Entry: " << entry << endl;
     }
    }

   } // end jet energy response

  } // end loop over events


  // close root files with systematics
  xCalibSystematicsFile->Close();
  if (xCalibSystematicsFile) delete xCalibSystematicsFile;
  
  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  const char* outFileName = Form("%s/JetInsituCorrectionCheck/dataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");

  // Write histograms to output and clean memory
  jetPreInsituSpectrum->Write();
  if (jetPreInsituSpectrum) delete jetPreInsituSpectrum;
  jetPostInsituSpectrum->Write();
  if (jetPostInsituSpectrum) delete jetPostInsituSpectrum;
  jetPostInsituSpectrumSysHi->Write();
  if (jetPostInsituSpectrumSysHi) delete jetPostInsituSpectrumSysHi;
  jetPostInsituSpectrumSysLo->Write();
  if (jetPostInsituSpectrumSysLo) delete jetPostInsituSpectrumSysLo;

  for (short pbin = 0; pbin < numpbins; pbin++) {
   jetInsituResponse[pbin]->Write();
   if (jetInsituResponse[pbin]) delete jetInsituResponse[pbin];
   jetInsituResponseSysHi[pbin]->Write();
   if (jetInsituResponseSysHi[pbin]) delete jetInsituResponseSysHi[pbin];
   jetInsituResponseSysLo[pbin]->Write();
   if (jetInsituResponseSysLo[pbin]) delete jetInsituResponseSysLo[pbin];
  }

  outFile->Close();
  if (outFile) delete outFile;
  return;
}
