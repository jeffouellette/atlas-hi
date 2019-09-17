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

  //TH1D*****   h_z_trk_pt_corr      = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);   // iSpc, iPtZ, iPhi, iCent
  //TH1D*****   h_z_trk_xzh_corr     = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);   // iSpc, iPtZ, iPhi, iCent
  //TH1D****    h_z_trk_zpt_corr     = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);               // iSpc, iPtZ, iCent
  //TH1D****    h_z_trk_zxzh_corr    = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);               // iSpc, iPtZ, iCent

  MinbiasAnalysis (const char* _name = "bkg") : FullAnalysis () {
    name = _name;
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    backgroundSubtracted = true;
    iaaCalculated = true;
    icpCalculated = true;
  }

  //void CreateHists () override;
  //void CopyAnalysis (PhysicsAnalysis* a, const bool copyBkgs = false) override;
  //void CombineHists ();
  //void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true);
  //void SaveHists (const char* histFileName = "savedHists.root");
  //void ScaleHists ();
  void Execute (const char* inFileName = "outFile.root", const char* outFileName = "savedHists.root") override;
};




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Create new histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: CreateHists () {
//
//  FullAnalysis :: CreateHists ();
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_pt_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_xzh_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
//        }
//      }
//    }
//  }
//
//  histsLoaded = true;
//  return;
//}
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
//// Fill combined species histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: CombineHists () {
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 2; iSpc++) {
//      for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          if (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_pt_corr[2][iPtZ][iPhi][iCent]->Add (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]);
//          if (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_xzh_corr[2][iPtZ][iPhi][iCent]->Add (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]);
//
//          if (iPhi != 0) {
//            if (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zpt_corr[iSpc][iPtZ][iCent]->Add (h_z_trk_zpt_corr[iSpc][iPtZ][iPhi][iCent]);
//            if (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zpt_corr[2][iPtZ][iCent]->Add (h_z_trk_zpt_corr[iSpc][iPtZ][iPhi][iCent]);
//            if (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zxzh_corr[iSpc][iPtZ][iCent]->Add (h_z_trk_zxzh_corr[iSpc][iPtZ][iPhi][iCent]);
//            if (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zxzh_corr[2][iPtZ][iCent]->Add (h_z_trk_zxzh_corr[iSpc][iPtZ][iPhi][iCent]);
//          }
//        } // end loop over phi
//      } // end loop over pT^Z
//    } // end loop over species
//  } // end loop over centralities
//  return;
//}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Scale histograms for plotting, calculating signals, etc.
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: ScaleHists () {
//  if (histsScaled || !histsLoaded)
//    return;
//
//  FullAnalysis :: ScaleHists ();
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
//        const float counts = countsHist->GetBinContent (1);
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          const double countsdPhi = counts * (phiHighBins[iPhi]-phiLowBins[iPhi]);
//
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
//        }
//      }
//    }
//  }
//}
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load pre-filled histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
//  if (histsLoaded)
//    return;
//
//  FullAnalysis :: LoadHists (histFileName, false);
//
//  if (!histFile) {
//    SetupDirectories ("", "ZTrackAnalysis/");
//    histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");
//  }
//  TDirectory* _gDirectory = gDirectory;
//  if (!histFile->IsOpen ()) {
//    cout << "Error in MinbiasAnalysis :: LoadHists: histFile not open after calling parent function, exiting." << endl;
//    return;
//  }
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        }
//      }
//    }
//  }
//
//  _gDirectory->cd ();
//
//  histsLoaded = true;
//
//  if (_finishHists) {
//    MinbiasAnalysis :: CombineHists ();
//    MinbiasAnalysis :: ScaleHists ();
//  }
//
//  return;
//}
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
//// Save histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: SaveHists (const char* histFileName) {
//  FullAnalysis :: SaveHists (histFileName);
//
//  if (!histFile) {
//    SetupDirectories ("", "ZTrackAnalysis/");
//    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "update");
//    histFile->cd ();
//  }
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          SafeWrite (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]);
//          SafeWrite (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]);
//        }
//      }
//    }
//  }
//
//  histFile->Close ();
//  histFile = nullptr;
//  histsLoaded = false;
//  return;
//}




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

  int run_number = 0, ntrk = 0, z_ntrk = 0;
  bool isEE = false;
  unsigned int event_number = 0;
  float event_weight = 1;//, fcal_weight = 1, q2_weight = 1, psi2_weight = 1, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_fcal_et = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  //float z_trk_pt[10000], z_trk_eta[10000], z_trk_phi[10000];
  //float*** yield_pt = Get3DArray <float> (2, maxNXHZBins, numPhiBins);
  //float*** yield_xhz =  Get3DArray <float> (2, maxNXHZBins, numPhiBins);
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
    std::vector <bool> eventsMixed = {};
    for (int i = 0; i < nMixEvts; i++) {
      eventOrder.push_back (i);
      eventsMixed.push_back (false);
    }
    std::random_shuffle (eventOrder.begin (), eventOrder.end ());

    TTree* zTree = LoadEventMixingTree (inFileName, "PbPbZTrackTree");
    if (!zTree)
      cout << "Got a null mixing tree!" << endl;
    const int nZEvts = zTree->GetEntries ();
    zTree->SetBranchAddress ("isEE",          &isEE);
    zTree->SetBranchAddress ("z_pt",          &z_pt);
    zTree->SetBranchAddress ("z_phi",         &z_phi);
    zTree->SetBranchAddress ("ntrk",          &z_ntrk);
    //zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
    //zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
    //zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
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

      {
        const short iFCalEt = GetSuperFineCentBin (z_fcal_et);
        if (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1)
          continue;
        
        bool goodMixEvent = false;
        const int _iMixEvt = iMixEvt;
        do {
          iMixEvt = (iMixEvt+1) % nMixEvts;
          PbPbTree->GetEntry (eventOrder[iMixEvt]);
          goodMixEvent = (iFCalEt == GetSuperFineCentBin (fcal_et) && !eventsMixed[iMixEvt]);
        } while (!goodMixEvent && iMixEvt != _iMixEvt);
        if (_iMixEvt == iMixEvt) {
          cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          continue;
        }
        eventsMixed[iMixEvt] = true;
      }

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
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

      //for (short iData : {0, 1}) {
      //  for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
      //    for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++)
      //      yield_pt[iData][iPtTrk][idPhi] = 0;
      //    for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++)
      //      yield_xhz[iData][iXHZ][idPhi] = 0;
      //  }
      //}

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const float trkPur = doTrackPurVar ? 1. : GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = event_weight * trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_z_trk_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
          }
        }

        //for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
        //  if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
        //    for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
        //      if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
        //        yield_pt[0][iPtTrk][idPhi] += trkWeight;
        //    }
        //    for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++) {
        //      if (xHZBins[iPtZ][iXHZ] <= trkpt / z_pt && trkpt / z_pt < ptTrkBins[iPtZ][iXHZ+1])
        //        yield_xhz[0][iXHZ][idPhi] += trkWeight;
        //    }
        //  }
        //}
      }

      //for (int iTrk = 0; iTrk < z_ntrk; iTrk++) {
      //  const float trkpt = z_trk_pt[iTrk];

      //  if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
      //    continue;

      //  const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, z_trk_eta[iTrk], true);
      //  const float trkPur = doTrackPurVar ? 1. : GetTrackingPurity (fcal_et, trkpt, z_trk_eta[iTrk], true);
      //  if (trkEff == 0 || trkPur == 0)
      //    continue;
      //  const float trkWeight = event_weight * trkPur / trkEff;

      //  float dphi = DeltaPhi (z_phi, z_trk_phi[iTrk], false);
      //  for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
      //    if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
      //      for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
      //        if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
      //          yield_pt[1][iPtTrk][idPhi] += trkWeight;
      //      }
      //      for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++) {
      //        if (xHZBins[iPtZ][iXHZ] <= trkpt / z_pt && trkpt / z_pt < ptTrkBins[iPtZ][iXHZ+1])
      //          yield_xhz[1][iXHZ][idPhi] += trkWeight;
      //      }
      //    }
      //  }
      //}

      //for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
      //  for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++)
      //    h_z_trk_raw_pt_corr[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtTrk, h_z_trk_raw_pt_corr[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtTrk) + yield_pt[0][iPtTrk][idPhi]*yield_pt[1][iPtTrk][idPhi]);
      //  for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++)
      //    h_z_trk_xzh_corr[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXHZ, h_z_trk_xzh_corr[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXHZ) + yield_xhz[0][iXHZ][idPhi]*yield_xhz[1][iXHZ][idPhi]);
      //}

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
    std::vector <bool> eventsMixed = {};
    for (int i = 0; i < nMixEvts; i++) {
      eventOrder.push_back (i);
      eventsMixed.push_back (false);
    }
    std::random_shuffle (eventOrder.begin (), eventOrder.end ());

    TTree* zTree = LoadEventMixingTree (inFileName, "ppZTrackTree");
    if (!zTree)
      cout << "Got a null mixing tree!" << endl;
    const int nZEvts = zTree->GetEntries ();
    zTree->SetBranchAddress ("isEE",          &isEE);
    zTree->SetBranchAddress ("z_pt",          &z_pt);
    zTree->SetBranchAddress ("z_phi",         &z_phi);
    zTree->SetBranchAddress ("ntrk",          &z_ntrk);
    //zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
    //zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
    //zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    zTree->SetBranchAddress ("event_weight",  &event_weight);

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

      {
        bool goodMixEvent = false;
        const int _iMixEvt = iMixEvt;
        do {
          iMixEvt = (iMixEvt+1) % nMixEvts;
          ppTree->GetEntry (eventOrder[iMixEvt]);
          goodMixEvent = (!eventsMixed[iMixEvt]);
        } while (!goodMixEvent && iMixEvt != _iMixEvt);
        if (_iMixEvt == iMixEvt) {
          cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          continue;
        }
        eventsMixed[iMixEvt] = true;
      }

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
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

      //for (short iData : {0, 1}) {
      //  for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
      //    for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++)
      //      yield_pt[iData][iPtTrk][idPhi] = 0;
      //    for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++)
      //      yield_xhz[iData][iXHZ][idPhi] = 0;
      //  }
      //}

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const float trkPur = doTrackPurVar ? 1. : GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = event_weight * trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_z_trk_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
          }
        }

        //for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
        //  if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
        //    for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
        //      if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
        //        yield_pt[0][iPtTrk][idPhi] += trkWeight;
        //    }
        //    for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++) {
        //      if (xHZBins[iPtZ][iXHZ] <= trkpt / z_pt && trkpt / z_pt < ptTrkBins[iPtZ][iXHZ+1])
        //        yield_xhz[0][iXHZ][idPhi] += trkWeight;
        //    }
        //  }
        //}
      }

      //for (int iTrk = 0; iTrk < z_ntrk; iTrk++) {
      //  const float trkpt = z_trk_pt[iTrk];

      //  if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
      //    continue;

      //  const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, z_trk_eta[iTrk], true);
      //  const float trkPur = doTrackPurVar ? 1. : GetTrackingPurity (fcal_et, trkpt, z_trk_eta[iTrk], true);
      //  if (trkEff == 0 || trkPur == 0)
      //    continue;
      //  const float trkWeight = event_weight * trkPur / trkEff;

      //  float dphi = DeltaPhi (z_phi, z_trk_phi[iTrk], false);
      //  for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
      //    if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
      //      for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
      //        if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
      //          yield_pt[1][iPtTrk][idPhi] += trkWeight;
      //      }
      //      for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++) {
      //        if (xHZBins[iPtZ][iXHZ] <= trkpt / z_pt && trkpt / z_pt < ptTrkBins[iPtZ][iXHZ+1])
      //          yield_xhz[1][iXHZ][idPhi] += trkWeight;
      //      }
      //    }
      //  }
      //}

      //for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
      //  for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++)
      //    h_z_trk_raw_pt_corr[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtTrk, h_z_trk_pt_corr[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtTrk) + yield_pt[0][iPtTrk][idPhi]*yield_pt[1][iPtTrk][idPhi]);
      //  for (short iXHZ = 0; iXHZ < nXHZBins[iPtZ]; iXHZ++)
      //    h_z_trk_xzh_corr[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXHZ, h_z_xzh_corr[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXHZ) + yield_xhz[0][iXHZ][idPhi]*yield_xhz[1][iXHZ][idPhi]);
      //}

    } // end loop over pp tree
    cout << "Done minbias pp loop." << endl;

    //Delete3DArray (yield_pt, 2, maxNPtTrkBins, numPhiBins);
    //Delete3DArray (yield_xhz, 2, maxNXHZBins, numPhiBins);
  }

  SaveHists (outFileName);

  inFile->Close ();
  SaferDelete (inFile);
}


#endif
