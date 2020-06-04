#ifndef __HijingMixedAnalysis_h__
#define __HijingMixedAnalysis_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_set>

using namespace std;
using namespace atlashi;

class HijingMixedAnalysis : public PhysicsAnalysis {

  public:

  bool takeNonTruthTracks = false;
  bool doCentMixing = true;
  bool doQ2Mixing = false;
  bool doPsi2Mixing = true;
  bool doPsi3Mixing = false;
  bool doVZMixing = false;

  int nQ2MixBins = 1;
  double* q2MixBins = nullptr;
  int nPsi2MixBins = 16;
  double* psi2MixBins = nullptr;
  int nPsi3MixBins = 1;
  double* psi3MixBins = nullptr;
  const int nIPBins = 156;//2640;//508;//298;//478;//298;//156;
  double* ipBins = nullptr;

  short GetQ2MixBin (const float q2) {
    if (!q2MixBins)
      return -1;
    short i = 0;
    while (i < nQ2MixBins) {
      if (q2 < q2MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetPsi2MixBin (const float psi2) {
    if (!psi2MixBins)
      return -1;
    short i = 0;
    while (i < nPsi2MixBins) {
      if (psi2 < psi2MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetPsi3MixBin (const float psi3) {
    if (!psi3MixBins)
      return -1;
    short i = 0;
    while (i < nPsi3MixBins) {
      if (psi3 < psi3MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetIPBin (const float ip) {
    if (!useImpactParameter)
      return -1;
    short i = 0;
    while (i < nIPBins-1) {
      if (ip < ipBins[i+1])
        break;
      i++;
    }
    return i;
  }

  HijingMixedAnalysis (const char* _name = "bkg") : PhysicsAnalysis () {
    name = _name;
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    isBkg = true;
    hasBkg = false;
    backgroundSubtracted = true;
    histsUnfolded = true;
    iaaCalculated = true;
    //icpCalculated = true;
    useHijingEffs = true;
    useImpactParameter = true;
    doPPMBMixing = false;
  }

  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  void Execute (const char* inFileName, const char* outFileName) override;


  private:
  // impact parameter centralities
  // 0%     1%      10%     20%     30%     40%     50%     60%     70%     80%
  // 0      1.564   4.953   7.009   8.582   9.911   11.083  12.136  13.111  14.032
  void GenerateIPMixBins () {
    int i = 0;
    ipBins = new double[nIPBins+1];

    for (i = 0; i < 10; i++)
      ipBins[i] = i*(1.564)/10.; // 0-1% central
    for (i = 10; i < 36; i++)
      ipBins[i] = 1.564 + (i-10)*(4.953-1.564)/36.; // 1-10% central
    for (i = 36; i < 126; i++)
      ipBins[i] = 4.953 + (i-36)*(11.083-4.953)/80.; // 10-50% central
    for (i = 126; i < 156; i++)
      ipBins[i] = 11.083 + (i-126)*(14.032-11.083)/30.; // 50-80% central
    ipBins[156] = 14.032;

    //for (i = 0; i < 160; i++)
    //  ipBins[i] = i*(1.564)/160.; // 0-1% central
    //for (i = 160; i < 880; i++)
    //  ipBins[i] = 1.564 + (i-160)*(4.953-1.564)/720.; // 1-10% central
    //for (i = 880; i < 2160; i++)
    //  ipBins[i] = 4.953 + (i-880)*(11.083-4.953)/1280.; // 10-50% central
    //for (i = 2160; i < 2640; i++)
    //  ipBins[i] = 11.083 + (i-2160)*(14.032-11.083)/480.; // 50-80% central
    //ipBins[2640] = 14.032;

    //for (i = 0; i < 40; i++)
    //  ipBins[i] = i*(1.564)/40.;
    //for (i = 40; i < 188; i++)
    //  ipBins[i] = 1.564 + (i-40)*(4.953-1.564)/144.;
    //for (i = 188; i < 268; i++)
    //  ipBins[i] = 4.953 + (i-188)*(11.083-4.953)/80.;
    //for (i = 268; i < 298; i++)
    //  ipBins[i] = 11.083 + (i-268)*(14.032-11.083)/30.;
    //ipBins[298] = 14.032;

    //for (i = 0; i < 80; i++)
    //  ipBins[i] = i*(1.564)/80.;
    //for (i = 80; i < 368; i++)
    //  ipBins[i] = 1.564 + (i-80)*(4.953-1.564)/288.;
    //for (i = 368; i < 448; i++)
    //  ipBins[i] = 4.953 + (i-368)*(11.083-4.953)/80.;
    //for (i = 448; i < 508; i++)
    //  ipBins[i] = 11.083 + (i-448)*(14.032-11.083)/60.;
    //ipBins[508] = 14.032;
  }
};




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load histograms into memory, then combine channels.
//////////////////////////////////////////////////////////////////////////////////////////////////
void HijingMixedAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);

  PhysicsAnalysis :: CombineHists ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void HijingMixedAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  LoadEventWeights ();
  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));

  SetupDirectories ("", "ZTrackAnalysis/");

  cout << "Arguments provided: " << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "outFileName = " << outFileName << endl;


  CreateHists ();


  // Get TTree with mixed events
  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;
  TTree* inTree = (TTree*) inFile->Get ("PbPbMixedTree");
  if (!inTree) {
    cout << "Could not find TTree, exiting" << endl;
    return;
  }
  const int nMixEvts = inTree->GetEntries ();
  cout << "N events in mixed event tree = " << nMixEvts << endl;


  // variables for branches
  int event_number = 0, run_number = 0;
  bool isEE = false;//, passes_toroid = false;
  float event_weight = 1., fcal_weight = 1., q2_weight = 1., psi2_weight = 1.;
  float fcal_et = 0, vz = 0, zdcEnergy = 0, ntrk_perp = 0;
  float ip = 0, eventPlane = 0;
  float phi_transmin = 0, phi_transmax = 0;
  float q2 = 0, q3 = 0, q4 = 0;
  float psi2 = 0, psi3 = 0, psi4 = 0;
  float z_pt = 0, z_phi = 0, z_y = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  int ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  int***    trks_counts   = Get3DArray <int> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (maxNPtchBins, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (maxNPtchBins, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (maxNPtchBins, 40);

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is found
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (inTree) {
    inTree->SetBranchAddress ("z_run_number",       &run_number);
    inTree->SetBranchAddress ("z_event_number",     &event_number);
    inTree->SetBranchAddress ("z_event_weight",     &event_weight);
    inTree->SetBranchAddress ("z_fcal_et",          &fcal_et);
    inTree->SetBranchAddress ("z_zdcEnergy",        &zdcEnergy);
    inTree->SetBranchAddress ("z_ip",               &ip);
    inTree->SetBranchAddress ("z_eventPlane",       &eventPlane);
    inTree->SetBranchAddress ("z_q2",               &q2);
    inTree->SetBranchAddress ("z_psi2",             &psi2);
    inTree->SetBranchAddress ("z_q3",               &q3);
    inTree->SetBranchAddress ("z_psi3",             &psi3);
    inTree->SetBranchAddress ("z_q4",               &q4);
    inTree->SetBranchAddress ("z_psi4",             &psi4);
    inTree->SetBranchAddress ("z_vz",               &vz);
    inTree->SetBranchAddress ("ntrk",               &ntrk);
    inTree->SetBranchAddress ("trk_pt",             &trk_pt);
    inTree->SetBranchAddress ("trk_eta",            &trk_eta);
    inTree->SetBranchAddress ("trk_phi",            &trk_phi);
    inTree->SetBranchAddress ("isEE",               &isEE);
    inTree->SetBranchAddress ("z_pt",               &z_pt);
    inTree->SetBranchAddress ("z_phi",              &z_phi);
    inTree->SetBranchAddress ("z_y",                &z_y);
    inTree->SetBranchAddress ("z_m",                &z_m);
    inTree->SetBranchAddress ("l1_pt",              &l1_pt);
    inTree->SetBranchAddress ("l1_eta",             &l1_eta);
    inTree->SetBranchAddress ("l1_phi",             &l1_phi);
    inTree->SetBranchAddress ("l1_trk_pt",          &l1_trk_pt);
    inTree->SetBranchAddress ("l1_trk_eta",         &l1_trk_eta);
    inTree->SetBranchAddress ("l1_trk_phi",         &l1_trk_phi);
    inTree->SetBranchAddress ("l1_charge",          &l1_charge);
    inTree->SetBranchAddress ("l2_pt",              &l2_pt);
    inTree->SetBranchAddress ("l2_eta",             &l2_eta);
    inTree->SetBranchAddress ("l2_phi",             &l2_phi);
    inTree->SetBranchAddress ("l2_trk_pt",          &l2_trk_pt);
    inTree->SetBranchAddress ("l2_trk_eta",         &l2_trk_eta);
    inTree->SetBranchAddress ("l2_trk_phi",         &l2_trk_phi);
    inTree->SetBranchAddress ("l2_charge",          &l2_charge);

    for (int iMixEvt = 0; iMixEvt < nMixEvts; iMixEvt++) {
      if (nMixEvts > 100 && iMixEvt % (nMixEvts / 100) == 0)
        cout << iMixEvt / (nMixEvts / 100) << "\% done...\r" << flush;

      inTree->GetEntry (iMixEvt);

      if (fabs (vz) > 150) continue;

      if (event_weight == 0) continue;
      if (isnan (event_weight))
        cout << "Warning: event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      if (iPtZ < 2) continue; // skip unneeded events

      const short iSpc = (isEE ? 0 : 1); // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = (useImpactParameter ? GetIPCentBin (ip) : GetCentBin (fcal_et));
      //const short iCent = GetCentBin (z_fcal_et);
      if (iCent < 1 || iCent > numCentBins-1) continue;
      //const short iFineCent = GetFineCentBin (z_fcal_et);
      //if (iFineCent < 1 || iFineCent > numFineCentBins-1) continue;
      const short iFineCent = 1;

      // do a reweighting procedure
      // only uses weights if the appropriate flags have been set. In most cases these default to false (exception is centrality weighting which can be default true).
      {
        if (useCentWgts)  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
        if (useQ2Wgts)    q2_weight   = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
        if (usePsi2Wgts) {
          float dphi = DeltaPhi (z_phi, psi2, false);
          if (dphi > pi/2)  dphi = pi - dphi;
          psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));
        }

        event_weight *= fcal_weight * q2_weight * psi2_weight;
      }

      // fill event category histograms
      {
        h_fcal_et->Fill (fcal_et);
        h_fcal_et_reweighted->Fill (fcal_et, event_weight);
        h_q2[iFineCent]->Fill (q2);
        h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
        h_psi2[iFineCent]->Fill (psi2);
        h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
        h_PbPb_vz->Fill (vz);
        h_PbPb_vz_reweighted->Fill (vz, event_weight);
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        const double trkEff = GetTrackingEfficiency (useImpactParameter ? ip : fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (useImpactParameter ? ip : fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkPur == 0 || trkEff == 0) continue;
        const float trkWeight = trkPur / trkEff;
        //const float trkWeight = 1;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {   
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {   
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks
      
      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (iPtch, iCent); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1
      
      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over Pb+Pb tree

    cout << endl << "Done minbias Pb+Pb loop." << endl;
  }


  Delete3DArray (&trks_counts, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights1, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights2, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete2DArray (&trks_counts_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights1_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights2_inPhi, maxNPtchBins, 40);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (&inFile);
}


#endif
