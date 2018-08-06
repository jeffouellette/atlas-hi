#include "electronContaminationStudy.h"

static const double pebins[18] = {20, 25, 35, 45, 55, 65, 75, 85, 105, 125, 150, 175, 200, 250, 300, 350, 400, 550};
static const short numpebins = sizeof(pebins)/sizeof(pebins[0]) - 1;

static const double eetabins[6] = { -2.37, -1.56, -1.37, 1.37, 1.56, 2.37};
static const short numeetabins = sizeof(eetabins)/sizeof(eetabins[0]) - 1;

TString GetIdentifier (const bool periodA) {
  TString id = "";
  if (periodA) id = "pPb_";
  else id = "Pbp_";
  id = id + "ZeeJet_Overlay";
  return id;
}


void electronContaminationStudy (const int dataSet,
                                 const bool isPeriodA,
                                 const TString inFileName)
{

  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  const TString identifier = GetIdentifier(isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   fileIdentifier = inFileName;

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
  t->SetBranchAddresses(true);

  // initialize histograms
  TH2D* electronContamination = new TH2D(Form("electronContamination_dataSet%s", identifier.Data()), ";#it{p}_{T}^{e} #left[GeV#right];#eta;", numpebins, pebins,  numeetabins, eetabins);
  electronContamination->Sumw2();
  TH2D* electronSpectrum = new TH2D (Form("electronSpectrum_dataSet%s", identifier.Data()), ";#it{p}_{T}^{e} #left[GeV#right];#eta;", numpebins, pebins, numeetabins, eetabins);
  electronSpectrum->Sumw2();


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
   // events with photons
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    TLorentzVector photon;
    const double photon_pt = t->photon_pt->at(p);
    const double photon_eta = t->photon_eta->at(p);
    const double photon_phi = t->photon_phi->at(p);
    photon.SetPtEtaPhiM (photon_pt, photon_eta, photon_phi, 0);

    // photon cuts
    if (photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at(p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at(p) > isolationEnergyIntercept + isolationEnergySlope*photon_pt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (photon_eta) || InDisabledHEC (photon_eta, photon_phi))
     continue; // require photon to be in EMCal

    double minDeltaR = 1000;
    int truth_electron = -1;
    for (int te = 0; te < t->truth_electron_n; te++) {
     const double dR = DeltaR (photon_eta, t->truth_electron_eta->at(te), photon_phi, t->truth_electron_phi->at(te));
     if (dR < minDeltaR) {
      truth_electron = te;
      minDeltaR = dR;
     }
    }
    if (truth_electron == -1 || minDeltaR > 0.4)
     continue; // unable to truth match to an electron

    electronContamination->Fill (t->truth_electron_pt->at(truth_electron), t->truth_electron_eta->at(truth_electron));

   }

   /////////////////////////////////////////////////////////////////////////////
   // events with electrons
   /////////////////////////////////////////////////////////////////////////////
   for (int e = 0; e < t->electron_n; e++) {
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

    electronSpectrum->Fill (t->electron_pt->at(e), t->electron_eta->at(e), t->eventWeight);
   }
    
  } // end loop over events

  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////


  const char* outFileName = Form("%s/electronContamination/dataSet_%s.root", rootPath.Data(), identifier.Data());
  TFile* outFile = new TFile(outFileName, "RECREATE");

  electronContamination->Write();
  if (electronContamination) delete electronContamination;
  electronSpectrum->Write();
  if (electronSpectrum) delete electronSpectrum;

  // Write histograms to output and clean memory

  outFile->Close();
  if (outFile) delete outFile;
  return;
}
