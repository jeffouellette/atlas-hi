#ifndef __MinbiasAnalysis_h__
#define __MinbiasAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class MinbiasAnalysis : public FullAnalysis {

  public:
  MinbiasAnalysis (const char* _name = "minbias", const char* subDir = "Nominal") : FullAnalysis () {
    name = _name;
    directory = Form ("MinbiasAnalysis/%s/", subDir);
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    backgroundSubtracted = true;
    LoadTrackingEfficiencies ();
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute () override;

  void CombineHists () override;
  void ScaleHists () override;

  void GenerateWeights ();
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: Execute () {
  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
  h_PbPbFCal_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbFCal_weights_minbias");
  for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
    h_PbPbQ2_weights[iCent] = (TH1D*) eventWeightsFile->Get (Form ("h_PbPbQ2_weights_iCent%i_minbias", iCent));
  }

  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*) inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*) inFile->Get ("ppZTrackTree");

  CreateHists ();

  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, event_weight = 1, fcal_weight = 1, q2_weight = 1, vz_weight = 1;
  vector<float>* trk_yield = nullptr;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchAddress ("vz",         &vz);
    PbPbTree->SetBranchAddress ("trk_yield",  &trk_yield);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iSpc = 0;
      const short iPtZ = nPtZBins-1;
      const short iXZTrk = 0;
      short iCent = 0;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      short iFinerCent = 0;
      while (iFinerCent < numFinerCentBins) {
        if (fcal_et < finerCentBins[iFinerCent])
          break;
        else
          iFinerCent++;
      }
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      fcal_weight = h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));
      q2_weight = h_PbPbQ2_weights[iFinerCent]->GetBinContent (h_PbPbQ2_weights[iFinerCent]->FindBin (q2));

      event_weight = fcal_weight * q2_weight * vz_weight;
      //event_weight = h_PbPb_event_reweights->GetBinContent (h_PbPb_event_reweights->FindBin (fcal_et, q2, vz));

      if (event_weight == 0)
        continue;

      h_fcal_et->Fill (fcal_et);
      //h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);
      //h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

      h_q2[iFinerCent]->Fill (q2);
      h_q2_reweighted[iFinerCent]->Fill (q2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        TH1D* h = h_z_trk_pt[iSpc][iPtZ][idPhi][iCent];
        //for (int iPhiTrk = 0; iPhiTrk < 3; iPhiTrk++) {
          for (int iPtTrk = 0; iPtTrk < 7; iPtTrk++) {
            const float trkEff = 1;
            //const float trkEff = GetTrackingEfficiency (fcal_et, h->GetBinCenter (iPtTrk), 0, true); // TODO temporary -- will need to update grid code when tracking efficiencies are available!
            //if (trkEff == 0)
            //  continue;
            h->SetBinContent (iPtTrk+1, h->GetBinContent (iPtTrk+1) + trk_yield->at (iPtTrk+7*idPhi) * event_weight / trkEff);
            h->SetBinError (iPtTrk+1, sqrt (pow (h->GetBinError (iPtTrk+1), 2) + pow (sqrt (trk_yield->at (iPtTrk+7*idPhi)) * event_weight / trkEff, 2))); 
          } // end loop over PtTrk bins
        //} // end loop over dPhiTrk (from grid)
      } // end loop over dPhi bins

    } // end loop over Pb+Pb tree
    cout << "Done minbias Pb+Pb loop." << endl;
  }


  if (ppTree) {
    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchAddress ("trk_yield",  &trk_yield);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      const short iSpc = 0;
      const short iPtZ = nPtZBins-1;
      const short iCent = 0; // iCent = 0 for pp

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        TH1D* h = h_z_trk_pt[iSpc][iPtZ][idPhi][iCent];
        //for (int iPhiTrk = 0; iPhiTrk < 3; iPhiTrk++) {
          for (int iPtTrk = 0; iPtTrk < 7; iPtTrk++) {
            const float trkEff = 1;
            //const float trkEff = GetTrackingEfficiency (fcal_et, h->GetBinCenter (iPtTrk), 0, false); // TODO temporary -- will need to update grid code when tracking efficiencies are available!
            //if (trkEff == 0)
            //  continue;
            h->SetBinContent (iPtTrk+1, h->GetBinContent (iPtTrk+1) + trk_yield->at (iPtTrk+7*idPhi) * event_weight / trkEff);
            h->SetBinError (iPtTrk+1, sqrt (pow (h->GetBinError (iPtTrk+1), 2) + pow (sqrt (trk_yield->at (iPtTrk+7*idPhi)) * event_weight / trkEff, 2))); 
          } // end loop over PtTrk bins
        //} // end loop over dPhiTrk (from grid)
      } // end loop over dPhi bins
    } // end loop over pp tree
    cout << "Done minbias pp loop." << endl;
  }

  CombineHists ();
  ScaleHists ();

  SaveHists ();

  inFile->Close ();
  SaferDelete (inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        if (iSpc == 0 && iPtZ == nPtZBins-1)
          continue;

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_pt[0][nPtZBins-1][iPhi][iCent]);
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_xzh[0][nPtZBins-1][iPhi][iCent]);
        } // end loop over phi
        
        h_z_counts[iSpc][iPtZ][iCent]->Add (h_z_counts[0][nPtZBins-1][iCent]);
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
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
          const double yieldNormFactor = countsHist->GetBinContent (1) * (pi/3.);
          //RescaleWithError (h, yieldNormFactor, yieldNormFactorError);
          if (yieldNormFactor > 0) {
            h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1. / yieldNormFactor);
            h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Scale (1. / yieldNormFactor);
          }
        } // end loop over phi
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (OLD VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: GenerateWeights () {
  const int gw_nFCalBins = 100;
  const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
  const int gw_nQ2Bins = 50;
  const double* gw_q2Bins = linspace (0, 0.3, gw_nQ2Bins);
  
  TH1D* _h_PbPbFCalDist = nullptr;
  TH1D* _h_PbPbFCal_weights = nullptr;
  TH1D* _h_PbPbQ2Dist[numFinerCentBins];
  TH1D* _h_PbPbQ2_weights[numFinerCentBins];
  

  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  TFile* inFile = new TFile (Form ("%s/Nominal/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "recreate");

  _h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
  _h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
  _h_PbPbFCalDist->Sumw2 ();
  _h_PbPbFCal_weights->Sumw2 ();
  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    _h_PbPbQ2Dist[iCent] = new TH1D (Form ("h_PbPbQ2Dist_iCent%i_%s", iCent, name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
    _h_PbPbQ2_weights[iCent] = new TH1D (Form ("h_PbPbQ2_weights_iCent%i_%s", iCent, name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
    _h_PbPbQ2Dist[iCent]->Sumw2 ();
    _h_PbPbQ2_weights[iCent]->Sumw2 ();
  }

  bool passes_toroid = true;
  float fcal_et = 0, q2 = 0, event_weight = 1;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchStatus ("vz",          0);
    PbPbTree->SetBranchStatus ("trk_yield",   0);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);
      _h_PbPbFCalDist->Fill (fcal_et);
    }
    cout << "Done 1st Pb+Pb loop." << endl;

    if (_h_PbPbFCalDist->Integral () > 0)
      _h_PbPbFCalDist->Scale (1./_h_PbPbFCalDist->Integral ());

    SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
    TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
    TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("h_PbPbFCalDist_data");

    TH1D* referenceQ2Dist[numFinerCentBins];
    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      referenceQ2Dist[iCent] = (TH1D*)ztrackFile->Get (Form ("h_PbPbQ2Dist_iCent%i_data", iCent));
    }

    for (int ix = 1; ix <= _h_PbPbFCal_weights->GetNbinsX (); ix++) {
      const double fcal_weight = (_h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / _h_PbPbFCalDist->GetBinContent (ix) : 0);
      _h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
    }

    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      short iCent = 0;
      while (iCent < numFinerCentBins) {
        if (fcal_et < finerCentBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numFinerCentBins-1)
        continue;

      event_weight = _h_PbPbFCal_weights->GetBinContent (_h_PbPbFCal_weights->FindBin (fcal_et));

      _h_PbPbQ2Dist[iCent]->Fill (q2, event_weight);
    }
    cout << "Done 2nd Pb+Pb loop." << endl;

    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      if (_h_PbPbQ2Dist[iCent]->Integral () > 0)
        _h_PbPbQ2Dist[iCent]->Scale (1./_h_PbPbQ2Dist[iCent]->Integral ());
    }

    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      for (int ix = 1; ix <= _h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
        const double q2_weight = (_h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? referenceQ2Dist[iCent]->GetBinContent (ix) / _h_PbPbQ2Dist[iCent]->GetBinContent (ix) : 0);
        _h_PbPbQ2_weights[iCent]->SetBinContent (ix, q2_weight);
      }
    }

    ztrackFile->Close ();
    SaferDelete (ztrackFile);

    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  }


  inFile->Close ();
  SaferDelete (inFile);

  eventWeightsFile->cd ();

  SafeWrite (_h_PbPbFCalDist);
  SafeWrite (_h_PbPbFCal_weights);
  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    SafeWrite (_h_PbPbQ2Dist[iCent]);
    SafeWrite (_h_PbPbQ2_weights[iCent]);
  }

  eventWeightsFile->Close ();
}


#endif
