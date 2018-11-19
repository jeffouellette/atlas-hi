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

  //TH1D* electronEnergyScale = new TH1D ("electronEnergyScale", "", 200, 0, 2);
  //electronEnergyScale->Sumw2 ();

  //TH1D***** jetEnergyResponseCalib = Get4DArray <TH1D*> (2, numpbins+1, numetabins+1, numphibins);
  //TH1D***** jetEnergyResponseReco = Get4DArray <TH1D*> (2, numpbins+1, numetabins+1, numphibins);
  //TH1D*** photonEnergyResponse = Get2DArray <TH1D*> (numpbins+1, numetabins+1);

  //for (short iP = 0; iP <= numpbins; iP++) {
  // for (short iEta = 0; iEta <= numetabins; iEta++) {
  //  for (short iPhi = 0; iPhi < numphibins; iPhi++) {
  //   for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
  //    const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
  //    jetEnergyResponseCalib[iAlgo][iP][iEta][iPhi] = new TH1D (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi), "", 100, 0, 2);
  //    jetEnergyResponseCalib[iAlgo][iP][iEta][iPhi]->Sumw2 ();
  //    jetEnergyResponseReco[iAlgo][iP][iEta][iPhi] = new TH1D (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi), "", 100, 0, 2);
  //    jetEnergyResponseReco[iAlgo][iP][iEta][iPhi]->Sumw2 ();
  //   }
  //  }
  //  photonEnergyResponse[iP][iEta] = new TH1D (Form ("photonEnergyResponse_iP%i_iEta%i", iP, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 100, 0, 2);
  //  photonEnergyResponse[iP][iEta]->Sumw2 ();
  // }
  //}

  //int**** nJet = Get4DArray <int> (2, numpbins+1, numetabins+1, numphibins);
  //int** nGamma = Get2DArray <int> (numpbins+1, numetabins+1);

  int jet_n;
  vector<float>* jet_pt = NULL, *jet_eta = NULL, *jet_phi = NULL, *jet_e = NULL, *precalib_jet_pt = NULL, *precalib_jet_e = NULL;

  const long long numEntries = tree->GetEntries ();

  //int* reco_jet_n = Get1DArray <int> (2);
  //int* truth_electron_matched_n = Get1DArray <int> (2);
  //int* truth_photon_matched_n = Get1DArray <int> (2);

  
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

   evtWeight = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;


   /////////////////////////////////////////////////////////////////////////////
   // basic event selection: e.g., require a primary vertex
   /////////////////////////////////////////////////////////////////////////////
   if ( (t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at (0) != 1)) continue;

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
     if (InDisabledHEC (jeta, jphi, 0.2))
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
     //reco_jet_n[iAlgo]++;

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

     //// Put the jet in the right eta bin
     //short iEta = 0;
     //// Make sure jet is in eta bounds
     //if (etabins[0] < t->truth_jet_eta->at (truth_jet) ||
     //    t->truth_jet_eta->at (truth_jet) < etabins[numetabins]) {
     // while (etabins[iEta] < t->truth_jet_eta->at (truth_jet)) iEta++;
     //}
     //iEta--;

     //// Put the jet in the right phi bin
     //short iPhi = 0;
     //// Make sure jet is in phi bounds
     //if (phibins[0] <= t->truth_jet_phi->at (truth_jet) ||
     //    t->truth_jet_phi->at (truth_jet) < phibins[numphibins]) {
     // while (phibins[iPhi] < t->truth_jet_phi->at (truth_jet)) iPhi++;
     //}
     //iPhi--;
     //if (!(0 <= iPhi && iPhi < numphibins))
     // continue;

     //short iP = 0;
     //if (!calcPtClosure) {
     // if (pbins[0] < t->truth_jet_e->at (truth_jet) ||
     //     t->truth_jet_e->at (truth_jet) < pbins[numpbins]) {
     //  while (pbins[iP] < t->truth_jet_e->at (truth_jet)) iP++;
     // }
     //}
     //else if (pbins[0] < t->truth_jet_pt->at (truth_jet) ||
     //         t->truth_jet_pt->at (truth_jet) < pbins[numpbins]) {
     // while (pbins[iP] < t->truth_jet_pt->at (truth_jet)) iP++;
     //}
     //iP--;

     jes = je / t->truth_jet_e->at (truth_jet);
     pjes = precalib_jet_e->at (j) / t->truth_jet_e->at (truth_jet);
     jpts = jpt / t->truth_jet_pt->at (truth_jet);
     pjpts = precalib_jet_pt->at (j) / t->truth_jet_pt->at (truth_jet);

     //if (0 <= iP && iP < numpbins &&
     //    0 <= iEta && iEta < numetabins) {
     // nJet[iAlgo][iP][iEta][iPhi]++;
     // jetEnergyResponseCalib[iAlgo][iP][iEta][iPhi]->Fill (jer, evtWeight);
     // jetEnergyResponseReco[iAlgo][iP][iEta][iPhi]->Fill (pjer, evtWeight);
     //}
     //if (0 <= iP && iP < numpbins) {
     // nJet[iAlgo][iP][numetabins][iPhi]++;
     // jetEnergyResponseCalib[iAlgo][iP][numetabins][iPhi]->Fill (jer, evtWeight);
     // jetEnergyResponseReco[iAlgo][iP][numetabins][iPhi]->Fill (pjer, evtWeight);
     //}
     //if (0 <= iEta && iEta < numetabins) {
     // nJet[iAlgo][numpbins][iEta][iPhi]++;
     // jetEnergyResponseCalib[iAlgo][numpbins][iEta][iPhi]->Fill (jer, evtWeight);
     // jetEnergyResponseReco[iAlgo][numpbins][iEta][iPhi]->Fill (pjer, evtWeight);
     //}
     //nJet[iAlgo][numpbins][numetabins][iPhi]++;
     //jetEnergyResponseCalib[iAlgo][numpbins][numetabins][iPhi]->Fill (jer, evtWeight);
     //jetEnergyResponseReco[iAlgo][numpbins][numetabins][iPhi]->Fill (pjer, evtWeight);
     outJetTree->Fill ();
     
    }
   } // end jet algorithm loop


   ///////////////////////////////////////////////////////////////////////////////
   //// electron containing events
   ///////////////////////////////////////////////////////////////////////////////
   //for (int e = 0; e < t->electron_n; e++) { // loop over primary electron

   // // electron cuts
   // if (t->electron_pt->at (e) < electron_pt_cut)
   //  continue; // basic electron pT cuts
   // if (!InEMCal (t->electron_eta->at (e)))
   //  continue; // reject electrons reconstructed outside EMCal
   // if (!t->electron_loose->at (e))
   //  continue; // reject non-loose electrons
   // if (t->electron_d0sig->at (e) > 5)
   //  continue; // d0 (transverse impact parameter) significance cut
   // if (t->electron_delta_z0_sin_theta->at (e) > 0.5)
   //  continue; // z0 (longitudinal impact parameter) vertex compatibility cut

   // // fill the electron energy scale plot first, finding the truth-reco pair
   // if (t->truth_electron_n > 0) {
   //  int closestEl = (e < t->truth_electron_n ? e:t->truth_electron_n-1); // best initial guess
   //  
   //  double minDeltaR = DeltaR (t->electron_eta->at (e), t->truth_electron_eta->at (closestEl), t->electron_phi->at (e), t->truth_electron_phi->at (closestEl));
   //  for (int te = 0; te < t->truth_electron_n; te++) {
   //   if (t->truth_electron_pt->at (te) < electron_pt_cut) continue;
   //   if (t->electron_charge->at (e) != t->truth_electron_charge->at (te)) continue;
   //   const double dR = DeltaR (t->electron_eta->at (e), t->truth_electron_eta->at (te), t->electron_phi->at (e), t->truth_electron_phi->at (te));
   //   if (dR < minDeltaR) {
   //    closestEl = te;
   //    minDeltaR = dR;
   //   }
   //  }
   //  if (minDeltaR < 0.4)
   //   electronEnergyScale->Fill (t->electron_pt->at (e) / t->truth_electron_pt->at (closestEl));//, t->eventWeight);
   // }
   //}


   /////////////////////////////////////////////////////////////////////////////
   // photon containing events
   /////////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < t->photon_n; p++) { // loop over all photons
    // relevant photon kinematic data
    ppt = t->photon_pt->at (p);
    peta = t->photon_eta->at (p);
    pphi = t->photon_phi->at (p);

    // photon cuts
    if (ppt < photon_pt_cut)
     continue; // basic pT cut on photons
    if (!t->photon_tight->at (p))
     continue; // require tight cuts on photons
    if (t->photon_topoetcone40->at (p) > isolationEnergyIntercept + isolationEnergySlope*ppt)
     continue; // require maximum isolation energy on gammas
    if (!InEMCal (peta) || InDisabledHEC (peta, pphi))
     continue; // require photon to be in EMCal

    // do photon truth matching to estimate photon energy scale
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

    // reject truth matched photons outside the DP range
    if (1 <= dataSet && dataSet <= numdpbins && (t->truth_photon_pt->at (truth_photon) < dpbins[dataSet-1] || dpbins[dataSet] < t->truth_photon_pt->at (truth_photon)))
     continue;

    //// Put the photon in the right eta, pt bin
    //short iEta = 0;
    //if (etabins[0] < t->truth_photon_eta->at (truth_photon) &&
    //    t->truth_photon_eta->at (truth_photon) < etabins[numetabins]) {
    // while (etabins[iEta] < t->truth_photon_eta->at (truth_photon)) iEta++;
    //}
    //iEta--;
    //short iP = 0;
    //if (pbins[0] < t->truth_photon_pt->at (truth_photon) &&
    //    t->truth_photon_pt->at (truth_photon) < pbins[numpbins]) {
    // while (pbins[iP] < t->truth_photon_pt->at (truth_photon)) iP++;
    //}
    //iP--;

    pes = ppt / t->truth_photon_pt->at (truth_photon);
    //if (0 <= iP && iP < numpbins &&
    //    0 <= iEta && iEta < numetabins) {
    // photonEnergyResponse[iP][iEta]->Fill (per, evtWeight);
    // nGamma[iP][iEta]++;
    //}
    //if (0 <= iP && iP < numpbins) {
    // photonEnergyResponse[iP][numetabins]->Fill (per, evtWeight);
    // nGamma[iP][numetabins]++;
    //}
    //if (0 <= iEta && iEta < numetabins) {
    // photonEnergyResponse[numpbins][iEta]->Fill (per, evtWeight);
    // nGamma[numpbins][iEta]++;
    //}
    //photonEnergyResponse[numpbins][numetabins]->Fill (per, evtWeight);
    //nGamma[numpbins][numetabins]++;
    outPhotonTree->Fill ();
   }

  } // end loop over events


  //////////////////////////////////////////////////////////////////////////////
  // End event loop
  //////////////////////////////////////////////////////////////////////////////


  //for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
  // if (iAlgo == 0) cout << "Anti-kt 4 HI Jets:     ";
  // else cout << "Anti-kt 4 EMTopo Jets: ";
  // cout << "truth-matched to electrons > " << electron_pt_cut << " GeV / total jets = " << truth_electron_matched_n[iAlgo] << " / " << reco_jet_n[iAlgo] << " = " << 100. * truth_electron_matched_n[iAlgo] / reco_jet_n[iAlgo] << "\%" << endl;
  // if (iAlgo == 0) cout << "Anti-kt 4 HI Jets:     ";
  // else cout << "Anti-kt 4 EMTopo Jets: ";
  // cout << "truth-matched to photons > " << photon_pt_cut << " GeV / total jets =   " << truth_photon_matched_n[iAlgo] << " / " << reco_jet_n[iAlgo] << " = " << 100. * truth_photon_matched_n[iAlgo] / reco_jet_n[iAlgo] << "\%" << endl;
  //}


  //const char* outFileName = Form ("%s/EnergyScaleChecks/dataSet_%s.root", rootPath.Data (), identifier.Data ());
  //TFile* outFile = new TFile (outFileName, "RECREATE");


  //// Write histograms to output and clean memory
  //electronEnergyScale->Write ();
  //if (electronEnergyScale) delete electronEnergyScale;

  //for (short iP = 0; iP <= numpbins; iP++) {
  // for (short iEta = 0; iEta <= numetabins; iEta++) {
  //  for (short iPhi = 0; iPhi < numphibins; iPhi++) {
  //   for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
  //    jetEnergyResponseCalib[iAlgo][iP][iEta][iPhi]->Write ();
  //    jetEnergyResponseReco[iAlgo][iP][iEta][iPhi]->Write ();
  //   }
  //  }
  //  photonEnergyResponse[iP][iEta]->Write ();
  // }
  //}

  //Delete4DArray (jetEnergyResponseCalib, 2, numpbins+1, numetabins+1, numphibins);
  //Delete4DArray (jetEnergyResponseReco, 2, numpbins+1, numetabins+1, numphibins);
  //Delete2DArray (photonEnergyResponse, numpbins+1, numetabins+1);

  //TVectorD nJetVec (2* (numpbins+1)* (numetabins+1));
  //TVectorD nGammaVec ( (numpbins+1)* (numetabins+1));
  //for (short iEta = 0; iEta <= numetabins; iEta++) {
  // for (short iP = 0; iP <= numpbins; iP++) {
  //  for (short iPhi = 0; iPhi < numphibins; iPhi++) {
  //   for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
  //    nJetVec[iPhi + iP*(numphibins+1) + iEta*(numpbins+1)*(numphibins+1) + iAlgo*(numpbins+1)*(numetabins+1)*(numphibins+1)] = (double)nJet[iAlgo][iP][iEta][iPhi];
  //   }
  //  }
  //  nGammaVec[iP + (numpbins+1)*iEta] = (double)nGamma[iP][iEta];
  // }
  //}
  //nJetVec.Write (Form ("nJetVec_%s", identifier.Data ()));
  //nGammaVec.Write (Form ("nGammaVec_%s", identifier.Data ()));

  //Delete4DArray (nJet, 2, numpbins+1, numetabins+1, numphibins);
  //Delete2DArray (nGamma, numpbins+1, numetabins+1);

  outFile->Write ();

  outFile->Close ();
  if (outFile) delete outFile;
  return;
}

} // end namespace
