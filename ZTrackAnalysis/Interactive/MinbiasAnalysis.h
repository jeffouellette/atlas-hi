#ifndef __MinbiasAnalysis_h__
#define __MinbiasAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace atlashi;

class MinbiasAnalysis : public FullAnalysis {

  private:
  TFile* zMixFile = nullptr;

  TTree* LoadEventMixingTree (const char* _inFile, const char* _treeName);

  public:
  MinbiasAnalysis (const char* _name = "bkg") : FullAnalysis () {
    name = _name;
    //eventWeightsExt = _name;
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    backgroundSubtracted = true;
    iaaCalculated = true;
    icpCalculated = true;
  }

  void Execute (const char* inFileName = "outFile.root", const char* outFileName = "savedHists.root") override;

  void CombineHists () override;
  void ScaleHists () override;

  //void GenerateWeights ();
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Loads the file for mixing Z's into minimum bias events
////////////////////////////////////////////////////////////////////////////////////////////////
TTree* MinbiasAnalysis :: LoadEventMixingTree (const char* _inFile, const char* _treeName) {
  if (zMixFile && zMixFile->IsOpen ())
    zMixFile->Close ();

  SetupDirectories ("", "ZTrackAnalysis/");

  TString inFile = TString (_inFile);
  if (!inFile.Contains (".root"))
    inFile = inFile + ".root";
  inFile.Replace (0, 7, "Data");

  zMixFile = new TFile (Form ("%s/%s", rootPath.Data (), inFile.Data ()), "read");

  if (!zMixFile || !zMixFile->IsOpen ()) {
    cout << "Cannot find zMixFile! Will return null tree!" << endl;
    return nullptr;
  }

  return (TTree*) zMixFile->Get (_treeName);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*) inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*) inFile->Get ("ppZTrackTree");

  CreateHists ();

  int run_number = 0, ntrk = 0;
  unsigned int event_number = 0;
  float event_weight = 1;//, fcal_weight = 1, q2_weight = 1, psi2_weight = 1, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_fcal_et = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  float z_pt = 0, z_phi = 0;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("run_number",   &run_number);
    PbPbTree->SetBranchAddress ("event_number", &event_number);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);
    PbPbTree->SetBranchAddress ("ntrk",         &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",       trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",      trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",      trk_phi);
    {
      const long memory = 8000000000;
      long memoryLoaded = PbPbTree->LoadBaskets (memory);
      cout << "Loaded tree, total space = " << memoryLoaded / 32000 << " GB" << endl;
    }
                        //2000000000 = 2GB

    int iMixEvt = 0;
    const int nMixEvts = PbPbTree->GetEntries ();
    PbPbTree->GetEntry (iMixEvt);

    std::vector <int> eventOrder = {};
    for (int i = 0; i < nMixEvts; i++) eventOrder.push_back (i);
    std::random_shuffle (eventOrder.begin (), eventOrder.end ());

    TTree* zTree = LoadEventMixingTree (inFileName, "PbPbZTrackTree");
    if (!zTree)
      cout << "Got a null mixing tree!" << endl;
    const int nZEvts = zTree->GetEntries ();
    zTree->SetBranchAddress ("z_pt",          &z_pt);
    zTree->SetBranchAddress ("z_phi",         &z_phi);
    zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
    zTree->SetBranchAddress ("event_weight",  &event_weight);

    if (nZEvts == 0)
      cout << "Warning! No Z's to mix with in this run!" << endl;
    cout << "For PbPb tree, maximum mixing fraction = " << nMixEvts / nZEvts << endl;
    if (mixingFraction * nZEvts > nMixEvts)
      cout << "Warning! Mixing fraction too high, will use all minimum bias events, mixing fraction = " << nMixEvts / nZEvts << endl;

    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (iZEvt % nZEvts);

      const short iPtZ = GetPtZBin (z_pt);

      const short iFCalEt = GetSuperFineCentBin (z_fcal_et);
      if (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1)
        continue;
      
      bool goodMixEvent = false;
      int _iMixEvt = iMixEvt;
      do {
        PbPbTree->GetEntry (eventOrder[iMixEvt]);
        iMixEvt = (iMixEvt+1) % nMixEvts;
        goodMixEvent = (iFCalEt == GetSuperFineCentBin (fcal_et));
      } while (!goodMixEvent && iMixEvt != _iMixEvt);
      if (_iMixEvt == iMixEvt)
        cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;

      const short iSpc = 0;
      const short iPhi = 0;

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      const short iFinerCent = GetFinerCentBin (fcal_et);
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      //fcal_weight = h_PbPbFCal_weights[iPtZ]->GetBinContent (h_PbPbFCal_weights[iPtZ]->FindBin (fcal_et));
      //q2_weight = h_PbPbQ2_weights[iFinerCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iFinerCent][iPtZ]->FindBin (q2));
      //psi2_weight = h_PbPbPsi2_weights[iFinerCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iFinerCent][iPtZ]->FindBin (psi2));

      //event_weight = fcal_weight * q2_weight * psi2_weight * vz_weight;

      if (event_weight == 0)
        continue;

      h_fcal_et->Fill (fcal_et);
      //h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);
      //h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

      h_q2[iFinerCent]->Fill (q2);
      h_q2_reweighted[iFinerCent]->Fill (q2, event_weight);
      h_psi2[iFinerCent]->Fill (psi2);
      h_psi2_reweighted[iFinerCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < ptTrkBins[iPtZ][0] || trkpt >= ptTrkBins[iPtZ][nPtTrkBins[iPtZ]])
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const float trkPur = doTrackPurVar ? GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true) : 1.;
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = event_weight * trkPur / trkEff;

        TH1D* h = h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trkpt, trkWeight);

        h = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trkpt, trkWeight);

        short iPtTrk = 0;
        while (iPtTrk < nPtTrkBins[iPtZ]) {
          if (trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            break;
          iPtTrk++;
        }
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;
        h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
        h->Fill (dphi, trkWeight);

        const float xHZ = trkpt / z_pt;
        if (xHZ < xHZBins[iPtZ][0] || xHZ >= xHZBins[iPtZ][nXHZBins[iPtZ]])
          continue;

        h = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
        h->Fill (xHZ, trkWeight);
      }

    } // end loop over Pb+Pb tree
    cout << "Done minbias Pb+Pb loop." << endl;
  }


  if (ppTree) {
    ppTree->SetBranchAddress ("run_number",   &run_number);
    ppTree->SetBranchAddress ("event_number", &event_number);
    ppTree->SetBranchAddress ("vz",           &vz);
    ppTree->SetBranchAddress ("ntrk",         &ntrk);
    ppTree->SetBranchAddress ("trk_pt",       trk_pt);
    ppTree->SetBranchAddress ("trk_eta",      trk_eta);
    ppTree->SetBranchAddress ("trk_phi",      trk_phi);
    ppTree->LoadBaskets (2000000000);

    int iMixEvt = 0;
    const int nMixEvts = ppTree->GetEntries ();
    ppTree->GetEntry (iMixEvt);

    std::vector <int> eventOrder = {};
    for (int i = 0; i < nMixEvts; i++) eventOrder.push_back (i);
    std::random_shuffle (eventOrder.begin (), eventOrder.end ());

    TTree* zTree = LoadEventMixingTree (Form ("%i", run_number), "ppZTrackTree");
    if (!zTree)
      cout << "Got a null mixing tree!" << endl;
    const int nZEvts = zTree->GetEntries ();
    zTree->SetBranchAddress ("z_pt",   &z_pt);
    zTree->SetBranchAddress ("z_phi",  &z_phi);
    zTree->SetBranchAddress ("event_weight", &event_weight);

    if (nZEvts == 0)
      cout << "Warning! No Z's to mix with in this run!" << endl;
    cout << "For pp tree, maximum mixing fraction = " << nMixEvts / nZEvts << endl;
    if (mixingFraction * nZEvts > nMixEvts)
      cout << "Warning! Mixing fraction too high, will use all minimum bias events, mixing fraction = " << nMixEvts / nZEvts << endl;

    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (iZEvt % nZEvts);

      const short iPtZ = GetPtZBin (z_pt);

      bool goodMixEvent = false;
      int _iMixEvt = iMixEvt;
      do {
        ppTree->GetEntry (eventOrder[iMixEvt]);
        iMixEvt = (iMixEvt+1) % nMixEvts;
        goodMixEvent = true;
      } while (!goodMixEvent && iMixEvt != _iMixEvt);
      if (_iMixEvt == iMixEvt)
        cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;

    //for (int iMixEvt = 0; iMixEvt < nMixEvts; iMixEvt++) {
    //  if (nMixEvts > 0 && iMixEvt % (nMixEvts / 100) == 0)
    //    cout << iMixEvt / (nMixEvts / 100) << "\% done...\r" << flush;
    //  ppTree->GetEntry (iMixEvt);

    //  zTree->GetEntry (iZEvt % nZEvts);

      const short iSpc = 0;
      const short iPhi = 0;
      const short iCent = 0; // iCent = 0 for pp

      //nch_weight = h_ppNch_weights->GetBinContent (h_ppNch_weights->FindBin (ntrk));

      //event_weight = nch_weight;

      if (event_weight == 0)
        continue;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < ptTrkBins[iPtZ][0] || trkpt >= ptTrkBins[iPtZ][nPtTrkBins[iPtZ]])
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const float trkPur = doTrackPurVar ? GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false) : 1.;
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = event_weight * trkPur / trkEff;

        TH1D* h = h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trkpt, trkWeight);

        h = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trkpt, trkWeight);

        short iPtTrk = 0;
        while (iPtTrk < nPtTrkBins[iPtZ]) {
          if (trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            break;
          iPtTrk++;
        }
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;
        h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
        h->Fill (dphi, trkWeight);

        const float xHZ = trkpt / z_pt;
        if (xHZ < xHZBins[iPtZ][0] || xHZ >= xHZBins[iPtZ][nXHZBins[iPtZ]])
          continue;

        h = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
        h->Fill (xHZ, trkWeight);
      }
    } // end loop over pp tree
    cout << "Done minbias pp loop." << endl;
  }

  //CombineHists ();
  //ScaleHists ();

  SaveHists (outFileName);

  inFile->Close ();
  SaferDelete (inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (iSpc == 0 && iPhi == 0)
            continue;
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_raw_pt[0][iPtZ][0][iCent]);
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_pt[0][iPtZ][0][iCent]);
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_xzh[0][iPtZ][0][iCent]);
        } // end loop over phi
        h_z_trk_zpt[iSpc][iPtZ][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][0][iCent]);
        h_z_trk_zxzh[iSpc][iPtZ][iCent]->Add (h_z_trk_xzh[iSpc][iPtZ][0][iCent]);

        if (iSpc == 0)
          continue;

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Add (h_z_trk_phi[0][iPtZ][iPtTrk][iCent]);
        }
      } // end loop over species
    } // end loop over pT^Z

    for (short iSpc = 1; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_counts[iSpc][iPtZ][iCent]->Add (h_z_counts[0][iPtZ][iCent]);
      } // end loop over pT^Z
    } // end loop over species
  } // end loop over centralities
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
        const float counts = countsHist->GetBinContent (1);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          const double countsdPhi = counts * (pi);
          if (countsdPhi > 0) {
            h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1./ counts);
            h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
            h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
          }
        } // end loop over phi

        if (counts > 0) {
          h_z_trk_zpt[iSpc][iPtZ][iCent]->Scale (1./ (counts * (pi)), "width");
          h_z_trk_zxzh[iSpc][iPtZ][iCent]->Scale (1./ (counts * (pi)), "width");
        }

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          TH1D* h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
          h->Rebin (2);
          if (iPtTrk > 3)
            h->Rebin (2);
          if (iCent != 0)
            h->Rebin (2);
          if (counts > 0)
            h->Scale (1. / counts);
        }
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  histsScaled = true;
  return;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Main macro. Loops over minbias trees and fills histograms appropriately. (OLD VERSION)
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: GenerateWeights () {
//  const int gw_nFCalBins = 100;
//  const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
//  const int gw_nQ2Bins = 20;
//  const double* gw_q2Bins = linspace (0, 0.3, gw_nQ2Bins);
//
//  const int gw_nNchBins = 160;
//  const double* gw_nchBins = linspace (-0.5, 160.5, gw_nNchBins);
//  
//  TH1D* _h_PbPbFCalDist = nullptr;
//  TH1D* _h_PbPbFCal_weights = nullptr;
//  TH1D* _h_PbPbQ2Dist[numFinerCentBins];
//  TH1D* _h_PbPbQ2_weights[numFinerCentBins];
//
//  TH1D* _h_ppNchDist = nullptr;
//  TH1D* _h_ppNch_weights = nullptr;
//  
//
//  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
//  TFile* inFile = new TFile (Form ("%s/Nominal/eventWeightsTree.root", rootPath.Data ()), "read");
//
//  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
//  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");
//
//  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "recreate");
//
//  _h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
//  _h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
//  _h_PbPbFCalDist->Sumw2 ();
//  _h_PbPbFCal_weights->Sumw2 ();
//  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
//    _h_PbPbQ2Dist[iCent] = new TH1D (Form ("h_PbPbQ2Dist_iCent%i_%s", iCent, name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
//    _h_PbPbQ2_weights[iCent] = new TH1D (Form ("h_PbPbQ2_weights_iCent%i_%s", iCent, name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
//    _h_PbPbQ2Dist[iCent]->Sumw2 ();
//    _h_PbPbQ2_weights[iCent]->Sumw2 ();
//  }
//  _h_ppNchDist = new TH1D (Form ("h_ppNchDist_%s", name.c_str ()), "", gw_nNchBins, gw_nchBins);
//  _h_ppNchDist->Sumw2 ();
//  _h_ppNch_weights = new TH1D (Form ("h_ppNch_weights_%s", name.c_str ()), "", gw_nNchBins, gw_nchBins);
//  _h_ppNch_weights->Sumw2 ();
//
//  //bool passes_toroid = true;
//  float fcal_et = 0, q2 = 0, event_weight = 1;
//
//  ////////////////////////////////////////////////////////////////////////////////////////////////
//  // Loop over PbPb tree
//  ////////////////////////////////////////////////////////////////////////////////////////////////
//  if (PbPbTree) {
//    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
//    PbPbTree->SetBranchAddress ("q2",         &q2);
//    PbPbTree->SetBranchStatus ("vz",          0);
//    PbPbTree->SetBranchStatus ("psi2",        0);
//
//    const int nEvts = PbPbTree->GetEntries ();
//    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
//      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
//        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
//      PbPbTree->GetEntry (iEvt);
//      _h_PbPbFCalDist->Fill (fcal_et);
//    }
//    cout << "Done 1st Pb+Pb loop." << endl;
//
//    if (_h_PbPbFCalDist->Integral () > 0)
//      _h_PbPbFCalDist->Scale (1./_h_PbPbFCalDist->Integral ());
//
//    SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
//    TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
//    TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get (Form ("h_PbPbFCalDist_iPtZ%i_data", nPtZBins));
//
//    TH1D* referenceQ2Dist[numFinerCentBins];
//    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
//      referenceQ2Dist[iCent] = (TH1D*)ztrackFile->Get (Form ("h_PbPbQ2Dist_iCent%i_iPtZ%i_data", iCent, nPtZBins));
//    }
//
//    for (int ix = 1; ix <= _h_PbPbFCal_weights->GetNbinsX (); ix++) {
//      const double fcal_weight = (_h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / _h_PbPbFCalDist->GetBinContent (ix) : 0);
//      _h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
//    }
//
//    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
//      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
//        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
//      PbPbTree->GetEntry (iEvt);
//
//      short iCent = 0;
//      while (iCent < numFinerCentBins) {
//        if (fcal_et < finerCentBins[iCent])
//          break;
//        else
//          iCent++;
//      }
//      if (iCent < 1 || iCent > numFinerCentBins-1)
//        continue;
//
//      event_weight = _h_PbPbFCal_weights->GetBinContent (_h_PbPbFCal_weights->FindBin (fcal_et));
//
//      _h_PbPbQ2Dist[iCent]->Fill (q2, event_weight);
//    }
//    cout << "Done 2nd Pb+Pb loop." << endl;
//
//    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
//      if (_h_PbPbQ2Dist[iCent]->Integral () > 0)
//        _h_PbPbQ2Dist[iCent]->Scale (1./_h_PbPbQ2Dist[iCent]->Integral ());
//    }
//
//    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
//      for (int ix = 1; ix <= _h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
//        const double q2_weight = (_h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? referenceQ2Dist[iCent]->GetBinContent (ix) / _h_PbPbQ2Dist[iCent]->GetBinContent (ix) : 0);
//        const double q2_weight_err = (_h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? q2_weight * sqrt (pow (_h_PbPbQ2Dist[iCent]->GetBinError (ix) / _h_PbPbQ2Dist[iCent]->GetBinContent (ix), 2) + pow (referenceQ2Dist[iCent]->GetBinError (ix) / referenceQ2Dist[iCent]->GetBinContent (ix), 2)) : 0);
//        _h_PbPbQ2_weights[iCent]->SetBinContent (ix, q2_weight);
//        _h_PbPbQ2_weights[iCent]->SetBinError (ix, q2_weight_err);
//      }
//    }
//
//    ztrackFile->Close ();
//    SaferDelete (ztrackFile);
//
//    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
//  }
//
//  if (ppTree) {
//    int ntrk = 0;
//    //ppTree->SetBranchAddress ("ntrk", &ntrk);
//
//    const int nEvts = ppTree->GetEntries ();
//    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
//      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
//        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
//      ppTree->GetEntry (iEvt);
//
//      _h_ppNchDist->Fill (ntrk);//, event_weight);
//    }
//    cout << "Done pp loop." << endl;
//
//    if (_h_ppNchDist->Integral () > 0)
//      _h_ppNchDist->Scale (1./_h_ppNchDist->Integral ());
//
//    if (name == "data") {
//      for (int ix = 1; ix <= _h_ppNch_weights->GetNbinsX (); ix++) {
//        _h_ppNch_weights->SetBinContent (ix, 1);
//      }
//    } else {
//      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
//      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
//      TH1D* referenceNchDist = (TH1D*)ztrackFile->Get ("h_ppNchDist_data");
//
//      for (int ix = 1; ix <= _h_ppNch_weights->GetNbinsX (); ix++) {
//        const double nch_weight = (_h_ppNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / _h_ppNchDist->GetBinContent (ix) : 0);
//        _h_ppNch_weights->SetBinContent (ix, nch_weight);
//      }
//
//      ztrackFile->Close ();
//
//      if (name == "data")
//        SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
//      else if (name == "mc")
//        SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
//      else if (name == "minbias")
//        SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
//      else if (name == "truth")
//        SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
//    }
//  }
//
//
//  inFile->Close ();
//  SaferDelete (inFile);
//
//  eventWeightsFile->cd ();
//
//  SafeWrite (_h_PbPbFCalDist);
//  SafeWrite (_h_PbPbFCal_weights);
//  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
//    SafeWrite (_h_PbPbQ2Dist[iCent]);
//    SafeWrite (_h_PbPbQ2_weights[iCent]);
//  }
//  SafeWrite (_h_ppNchDist);
//  SafeWrite (_h_ppNch_weights);
//
//  eventWeightsFile->Close ();
//}


#endif
