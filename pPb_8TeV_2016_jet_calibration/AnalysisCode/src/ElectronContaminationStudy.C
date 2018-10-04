#include "ElectronContaminationStudy.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Utils.h"

#include <TSystemDirectory.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <iostream>

namespace pPb8TeV2016JetCalibration {

static const double pebins[18] = {20, 25, 35, 45, 55, 65, 75, 85, 105, 125, 150, 175, 200, 250, 300, 350, 400, 550};
static const short numpebins = sizeof (pebins)/sizeof (pebins[0]) - 1;

static const double eetabins[6] = { -2.37, -1.56, -1.37, 1.37, 1.56, 2.37};
static const short numeetabins = sizeof (eetabins)/sizeof (eetabins[0]) - 1;

void ElectronContaminationStudy (const bool isPeriodA,
                                 const TString inFileName)
{

  SetupDirectories ("", "pPb_8TeV_2016_jet_calibration/");

  const TString identifier = GetIdentifier (0, inFileName, true, false, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = NULL;
  TTree* tree = NULL;
  {
   TString fileIdentifier;
   fileIdentifier = inFileName;

   TSystemDirectory dir (dataPath.Data (), dataPath.Data ());
   TList* sysfiles = dir.GetListOfFiles ();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile* sysfile;
   TString fname;
   TIter next (sysfiles);

   while ( (sysfile = (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In ElectronContaminationStudy.C (breakpoint B): Found " << fname.Data () << endl;
     
     if (fname.Contains (fileIdentifier)) {
      file = new TFile (dataPath+fname, "READ");
      tree = (TTree*)file->Get ("tree");
      break;
     }
    }
   }
  }
  if (tree == NULL || file == NULL) {
   cout << "Error: In ElectronContaminationStudy.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, true);
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetPhotons ();
  t->SetBranchAddresses ();

  // initialize histograms
  TH2D* fakePhotonSpectrum = new TH2D (Form ("fakePhotonSpectrum_dataSet%s", identifier.Data ()), "", numpebins, pebins,  numeetabins, eetabins);
  fakePhotonSpectrum->Sumw2 ();

  TH2D* fakePhotonCounts = new TH2D (Form ("fakePhotonCounts_dataSet%s", identifier.Data ()), "", numpebins, pebins, numeetabins, eetabins);
  fakePhotonCounts->Sumw2 ();

  TH2D* allElectronCounts = new TH2D (Form ("allElectronCounts_dataSet%s", identifier.Data ()), "", numpebins, pebins, numeetabins, eetabins);
  allElectronCounts->Sumw2 ();

  TH2D* allGammaCounts = new TH2D (Form ("allGammaCounts_dataSet%s", identifier.Data ()), "", numpebins, pebins, numeetabins, eetabins);
  allGammaCounts->Sumw2 ();

  TH2D* electronSpectrum = new TH2D (Form ("electronSpectrum_dataSet%s", identifier.Data ()), "", numpebins, pebins, numeetabins, eetabins);
  electronSpectrum->Sumw2 ();

  TH2D* truthElectronRecoPhotonCounts = new TH2D (Form ("truthElectronRecoPhotonCounts_dataSet%s", identifier.Data ()), "", numpebins, pebins, numeetabins, eetabins);
  truthElectronRecoPhotonCounts->Sumw2 ();

  TH2D* truthElectronRecoElectronCounts = new TH2D (Form ("truthElectronRecoElectronCounts_dataSet%s", identifier.Data ()), "", numpebins, pebins, numeetabins, eetabins);
  truthElectronRecoElectronCounts->Sumw2 ();


  const long long numEntries = tree->GetEntries ();

  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (long long entry = 0; entry < numEntries; entry++) {
   tree->GetEntry (entry);


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;

   const double weight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;


   /////////////////////////////////////////////////////////////////////////////
   // events with photons
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    TLorentzVector photon;
    const double photon_pt = t->photon_pt->at (p);
    const double photon_eta = t->photon_eta->at (p);
    const double photon_phi = t->photon_phi->at (p);
    photon.SetPtEtaPhiM (photon_pt, photon_eta, photon_phi, 0);

    // photon cuts
    if (photon_pt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at (p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*photon_pt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (photon_eta) || InDisabledHEC (photon_eta, photon_phi))
     continue; // require photon to be in EMCal

    const double eta = (!isPeriodA ? photon_eta : -photon_eta);
    allGammaCounts->Fill (photon_pt, eta);

    double minDeltaR = 1000;
    int truth_electron = -1;
    for (int te = 0; te < t->truth_electron_n; te++) {
     const double dR = DeltaR (photon_eta, t->truth_electron_eta->at (te), photon_phi, t->truth_electron_phi->at (te));
     if (dR < minDeltaR) {
      truth_electron = te;
      minDeltaR = dR;
     }
    }
    if (truth_electron == -1 || minDeltaR > 0.4)
     continue; // unable to truth match to an electron

    fakePhotonSpectrum->Fill (photon_pt, eta, weight);
    fakePhotonCounts->Fill (photon_pt, eta);

   }

   /////////////////////////////////////////////////////////////////////////////
   // events with truth electrons
   /////////////////////////////////////////////////////////////////////////////
   for (int e = 0; e < t->truth_electron_n; e++) {

    const double eta = (!isPeriodA ? t->truth_electron_eta->at (e) : -t->truth_electron_eta->at (e));

    allElectronCounts->Fill (t->truth_electron_pt->at (e), eta);
    electronSpectrum->Fill (t->truth_electron_pt->at (e), eta, weight);

    // try to reco match to a photon
    double minDeltaR = 1000;
    int best_p = -1;
    for (int p = 0; p < t->photon_n; p++) { // loop over reco photons
     // photon cuts
     if (t->photon_pt->at (p) < photon_pt_cut)
      continue; // basic pT cut on photons
     if (!t->photon_tight->at (p))
      continue; // require tight cuts on photons
     if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*t->photon_pt->at (p))
      continue; // require maximum isolation energy on gammas
     if (!InEMCal (t->photon_eta->at (p)) || InDisabledHEC (t->photon_eta->at (p), t->photon_phi->at (p)))
      continue; // require photon to be in EMCal

     const double dR = DeltaR (t->truth_electron_eta->at (e), t->photon_eta->at (p), t->truth_electron_phi->at (e), t->photon_phi->at (p));
     if (dR < minDeltaR) {
      best_p = p;
      minDeltaR = dR;
     } 
    }
    const double pDeltaR = minDeltaR;

    // try to reco match to an electron
    minDeltaR = 1000;
    int best_e = -1;
    for (int re = 0; re < t->electron_n; re++) { // loop over reco electrons
     // electron cuts
     if (t->electron_pt->at (re) < electron_pt_cut)
      continue; // basic electron pT cuts
     if (!InEMCal (t->electron_eta->at (re)))
      continue; // reject electrons reconstructed outside EMCal
     if (!t->electron_loose->at (re))
      continue; // reject non-loose electrons
     if (t->electron_d0sig->at (re) > 5)
      continue; // d0 (transverse impact parameter) significance cut
     if (t->electron_delta_z0_sin_theta->at (re) > 0.5)
      continue; // z0 (longitudinal impact parameter) vertex compatibility cut

     const double dR = DeltaR (t->truth_electron_eta->at (e), t->electron_eta->at (re), t->truth_electron_phi->at (e), t->electron_phi->at (re));
     if (dR < minDeltaR) {
      best_e = re;
      minDeltaR = dR;
     }
    }
    const double eDeltaR = minDeltaR;

    // if reco matched to a photon
    if (pDeltaR < eDeltaR && pDeltaR < 0.4) {
     truthElectronRecoPhotonCounts->Fill (t->truth_electron_pt->at (e), eta);
    }
    else if (eDeltaR < pDeltaR && eDeltaR < 0.4) {
     truthElectronRecoElectronCounts->Fill (t->truth_electron_pt->at (e), eta);
    }

   }

    
  } // end loop over events

  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////


  const char* outFileName = Form ("%s/ElectronContamination/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  fakePhotonSpectrum->Write ();
  if (fakePhotonSpectrum) delete fakePhotonSpectrum;
  fakePhotonCounts->Write ();
  if (fakePhotonCounts) delete fakePhotonCounts;
  allElectronCounts->Write ();
  if (allElectronCounts) delete allElectronCounts;
  allGammaCounts->Write ();
  if (allGammaCounts) delete allGammaCounts;
  electronSpectrum->Write ();
  if (electronSpectrum) delete electronSpectrum;
  truthElectronRecoPhotonCounts->Write ();
  if (truthElectronRecoPhotonCounts) delete truthElectronRecoPhotonCounts;
  truthElectronRecoElectronCounts->Write ();
  if (truthElectronRecoElectronCounts) delete truthElectronRecoElectronCounts;

  // Write histograms to output and clean memory

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
