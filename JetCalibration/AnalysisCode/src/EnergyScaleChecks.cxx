#include "EnergyScaleChecks.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>
#include <TreeVariables.h>

#include <TH2D.h>
#include <TVectorT.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {

void EnergyScaleChecks (const char* directory,
                        const int dataSet,
                        const bool isPeriodA,
                        const char* inFileName,
                        const double crossSection_microbarns,
                        const double filterEfficiency,
                        const int numberEvents)
{

  SetupDirectories ("", "JetCalibration/");

  const bool isSignalOnlySample = TString (inFileName).Contains ("signalonly");
  const TString identifier = GetIdentifier (dataSet, inFileName, true, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, true, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("tree");
  if (tree == NULL || file == NULL) {
    cout << "Error: In EnergyScaleChecks.C: TTree not obtained for given data set. Quitting." << endl;
    return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, true);
  if (crossSection_microbarns != 0.)
    t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetVertices ();
  t->SetGetHIJets ();
  t->SetGetEMTopoJets ();
  t->SetGetTruthJets ();
  t->SetGetElectrons ();
  t->SetGetPhotons ();
  t->SetGetTruthElectrons ();
  t->SetGetTruthPhotons ();
  t->SetBranchAddresses ();

  int jet_n;
  vector<float>* jet_pt = NULL, *jet_eta = NULL, *jet_phi = NULL, *jet_e = NULL, *precalib_jet_pt = NULL, *precalib_jet_e = NULL;

  const long long numEntries = tree->GetEntries ();
  
  const char* outFileName = Form ("%s/EnergyScaleChecks/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  TFile* outFile = new TFile (outFileName, "RECREATE");

  float evtWeight;
  float jpt = 0, jeta = 0, jphi = 0, je = 0, jes = 0, pjes = 0, jpts = 0, pjpts = 0;
  float ppt = 0, peta = 0, pphi = 0, pes = 0;

  TTree* outJetTree = new TTree ("jeffsjets", "jeffsjets");
  TTree* outPhotonTree = new TTree ("jeffsphotons", "jeffsphotons");

  outJetTree->SetDirectory (outFile);
  outPhotonTree->SetDirectory (outFile);

  outJetTree->Branch ("evt_weight", &evtWeight, "evt_weight/F");
  outJetTree->Branch ("jet_pt", &jpt, "jet_pt/F");
  outJetTree->Branch ("jet_eta", &jeta, "jet_eta/F");
  outJetTree->Branch ("jet_phi", &jphi, "jet_phi/F");
  outJetTree->Branch ("jet_e", &je, "jet_e/F");
  outJetTree->Branch ("jet_energy_scale", &jes, "jet_energy_scale/F");
  outJetTree->Branch ("precalib_jet_energy_scale", &pjes, "precalib_jet_energy_scale/F");
  outJetTree->Branch ("jet_pt_scale", &jpts, "jet_pt_scale/F");
  outJetTree->Branch ("precalib_jet_pt_scale", &pjpts, "precalib_jet_pt_scale/F");

  outPhotonTree->Branch ("evt_weight", &evtWeight, "evt_weight/F");
  outPhotonTree->Branch ("photon_pt", &ppt, "photon_pt/F");
  outPhotonTree->Branch ("photon_eta", &peta, "photon_eta/F");
  outPhotonTree->Branch ("photon_phi", &pphi, "photon_phi/F");
  outPhotonTree->Branch ("photon_energy_scale", &pes, "photon_energy_scale/F");


  //////////////////////////////////////////////////////////////////////////////
  // begin loop over events
  //////////////////////////////////////////////////////////////////////////////
  for (long long entry = 0; entry < numEntries; entry++) {
    tree->GetEntry (entry);

    evtWeight = (double)t->crossSection_microbarns * (double)t->filterEfficiency / (double)numEntries;


    /////////////////////////////////////////////////////////////////////////////
    // basic event selection: e.g., require a primary vertex
    /////////////////////////////////////////////////////////////////////////////
    if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1))
      continue;

    // basically just do the jet part twice, once for each algorithm
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {

      if (iAlgo == 0) { // HI algorithm at the EM scale
        jet_n = * (& (t->akt4hi_jet_n));
        jet_pt = t->akt4hi_em_xcalib_jet_pt;
        jet_eta = t->akt4hi_em_xcalib_jet_eta;
        jet_phi = t->akt4hi_em_xcalib_jet_phi;
        jet_e = t->akt4hi_em_xcalib_jet_e;
        precalib_jet_pt = t->akt4hi_em_jet_pt;
        precalib_jet_e = t->akt4hi_em_jet_e;
      }
      else { // EMTopo algorithm at the EM scale
        jet_n = * (& (t->akt4emtopo_jet_n));
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
      for (int j = 0; j < jet_n; j++) {
        jpt = jet_pt->at (j);
        jeta = jet_eta->at (j);
        jphi = jet_phi->at (j);
        je = jet_e->at (j);

        if (jpt < jet_pt_cut)
          continue; // basic jet pT cut
        if (InDisabledHEC (jeta, jphi, 0.4))
          continue; // only truth match jets outside disabled HEC
        if (!InHadCal (jeta))
          continue; // reject jets reconstructed outside reasonable HCal bounds.

        double minDeltaR = 1000;
        // loop over truth electrons and photons
        for (int e = 0; e < t->truth_electron_n; e++) { // truth electrons
          //if (t->truth_electron_pt->at (e) < electron_pt_cut) continue; // pt cut on truth electrons

          const double dR = DeltaR (jeta, t->truth_electron_eta->at (e), jphi, t->truth_electron_phi->at (e));
          if (dR < minDeltaR)
            minDeltaR = dR;
          //if (dR < 0.6)
          // truth_electron_matched_n[iAlgo]++;
        }
        for (int p = 0; p < t->truth_photon_n; p++) { // truth photons
          //if (t->truth_photon_pt->at (p) < photon_pt_cut) continue; // pt cut on truth photons

          const double dR = DeltaR (jeta, t->truth_photon_eta->at (p), jphi, t->truth_photon_phi->at (p));
          if (dR < minDeltaR)
           minDeltaR = dR;
          //if (dR < 0.6)
          // truth_photon_matched_n[iAlgo]++;
        }

        if (minDeltaR < 0.6)
          continue; // reject jets close to some other lepton or photon
      

        minDeltaR = 1000;
        int truth_jet = -1;
        for (int tj = 0; tj < t->truth_jet_n; tj++) {
          double dR = DeltaR (jeta, t->truth_jet_eta->at (tj), jphi, t->truth_jet_phi->at (tj));
          if (dR < minDeltaR) {
            minDeltaR = dR;
            truth_jet = tj;
          }
        }

        if (1 <= dataSet && dataSet <= numdpbins && (t->truth_jet_pt->at (truth_jet) < dpbins[dataSet-1] || dpbins[dataSet] < t->truth_jet_pt->at (truth_jet)))
          continue;

        if (!(0 <= truth_jet && truth_jet < t->truth_jet_n && minDeltaR < 0.2))
          continue;

        jes = je / t->truth_jet_e->at (truth_jet);
        pjes = precalib_jet_e->at (j) / t->truth_jet_e->at (truth_jet);
        jpts = jpt / t->truth_jet_pt->at (truth_jet);
        pjpts = precalib_jet_pt->at (j) / t->truth_jet_pt->at (truth_jet);

        outJetTree->Fill ();
       
      }
    } // end jet algorithm loop


    /////////////////////////////////////////////////////////////////////////////
    // photon containing events
    /////////////////////////////////////////////////////////////////////////////
    for (int p = 0; p < t->photon_n; p++) { // loop over all photons

     /////////////////////////////////////////////////////////////////////////////
     // relevant photon kinematic data
     /////////////////////////////////////////////////////////////////////////////
     ppt = t->photon_pt->at (p);
     peta = t->photon_eta->at (p);
     pphi = t->photon_phi->at (p);

     /////////////////////////////////////////////////////////////////////////////
     // photon cuts
     /////////////////////////////////////////////////////////////////////////////
     if (ppt < photon_pt_cut)
       continue; // basic pT cut on photons
     if (!t->photon_tight->at (p))
       continue; // require tight cuts on photons
     if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*ppt)
       continue; // require maximum isolation energy on gammas
     if (!InEMCal (peta) || InDisabledHEC (peta, pphi, 0.1))
       continue; // require photon to be in EMCal

     /////////////////////////////////////////////////////////////////////////////
     // do photon truth matching to estimate photon energy scale
     /////////////////////////////////////////////////////////////////////////////
     int truth_photon = (p < t->truth_photon_n ? p : t->truth_photon_n-1); // best initial guess
     double minDeltaR = 1000;
     for (int tp = 0; tp < t->truth_photon_n; tp++) {
       const double dR = DeltaR (peta, t->truth_photon_eta->at (tp), pphi, t->truth_photon_phi->at (tp));
       if (dR < minDeltaR) {
         truth_photon = tp;
         minDeltaR = tp;
       }
     }
     if (minDeltaR >= 0.2)
       continue; // reco photons not matched to a truth photon within dR=0.2 are skipped

     /////////////////////////////////////////////////////////////////////////////
     // reject truth matched photons outside the DP range
     /////////////////////////////////////////////////////////////////////////////
     if (1 <= dataSet && dataSet <= numdpbins && (t->truth_photon_pt->at (truth_photon) < dpbins[dataSet-1] || dpbins[dataSet] < t->truth_photon_pt->at (truth_photon)))
       continue;

     pes = ppt / t->truth_photon_pt->at (truth_photon);
     outPhotonTree->Fill ();
    }

  } // end loop over events


  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////

  outFile->Write ();

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
