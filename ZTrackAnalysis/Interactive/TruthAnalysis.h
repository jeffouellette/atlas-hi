#ifndef __TruthAnalysis_h__
#define __TruthAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <TRandom3.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class TruthAnalysis : public FullAnalysis {

  public:

  bool doSmearing = false;
  TFile* trkMomResFile = nullptr;
  bool tmrsLoaded   = false;

  TH2D** h2_tmr_factors = Get1DArray <TH2D*> (numCentBins);
  TRandom3* tmr_rndm = nullptr;

  TruthAnalysis (const char* _name = "truth") : FullAnalysis () {
    name = _name;
    plotFill = false;
    useAltMarker = false;
    hasBkg = false;
    histsUnfolded = true;
    isMC = true;
  }

  void Execute (const char* inFileName, const char* outFileName) override;
  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;

  void LoadTrackMomentumResolutions ();
  double SmearTrackMomentum (const float fcal_et, const double truth_pt, const double truth_eta, const bool isPbPb);
};




//////////////////////////////////////////////////////////////////////////////////////////////////
// Load track momentum resolution factors into memory for truth-level smearing.
//////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: LoadTrackMomentumResolutions () {
  if (tmrsLoaded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory; 
  TString _tmrDir = "Nominal";
  if (useHITight)
    _tmrDir = "Variations/TrackHITightWPVariation";
  else if (doTrackEffVar)
    _tmrDir = "Variations/TrackEffPionsVariation";
  cout << Form ("Reading tracking efficiencies from %s/TrackingMomentumResolution/%s/trackingMomentumResolutionFactors.root", rootPath.Data (), _tmrDir.Data ()) << endl;
  trkMomResFile = new TFile (Form ("%s/TrackingMomentumResolution/%s/trackingMomentumResolutionFactors.root", rootPath.Data (), _tmrDir.Data ()), "read");

  for (int iCent = 0; iCent < numTrkCorrCentBins; iCent++) {
    h2_tmr_factors[iCent] = (TH2D*) trkMomResFile->Get (Form ("h2_avg_tmr_iCent%i", iCent));
  }

  tmrsLoaded = true;

  tmr_rndm = new TRandom3 ();
  tmr_rndm->SetSeed (8141995);

  _gDirectory->cd ();
  return;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
// Smears track momentum by the resolution factors to identify any unfolding effect.
//////////////////////////////////////////////////////////////////////////////////////////////////
double TruthAnalysis :: SmearTrackMomentum (const float fcal_et, const double truth_pt, const double truth_eta, const bool isPbPb) {
  if (!tmrsLoaded)
    LoadTrackMomentumResolutions ();

  short iCent = 0;
  if (isPbPb && !useImpactParameter) {
    while (iCent < numTrkCorrCentBins) {
      if (fcal_et < trkCorrCentBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent == numTrkCorrCentBins) // force ultra-central events to behave like 0-10% Pb+Pb
      iCent--;
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }
  else if (isPbPb && useImpactParameter) {
    iCent = GetIPCentBin (fcal_et); // fcal_et variable should actually be impact parameter if useImpactParameter is true
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }

  TH2D* h2 = h2_tmr_factors[iCent];

  const int xbin = h2->GetXaxis ()->FindFixBin (truth_eta);
  const int ybin = h2->GetYaxis ()->FindFixBin (truth_pt);
  if (xbin < 1 || h2->GetXaxis ()->GetNbins () < xbin)
    return 0;
  if (ybin < 1 || h2->GetYaxis ()->GetNbins () < ybin)
    return 0;

  const double tmr = h2->GetBinContent (xbin, ybin);

  return tmr_rndm->Gaus (1, tmr/100.) * truth_pt;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
// Load histograms into memory, then combine channels.
//////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);
  PhysicsAnalysis :: CombineHists ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1;//, fcal_weight = 1, q2_weight = 1, psi2_weight = 1, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, ip = 0, eventPlane = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_eta = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;//, njet = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  //
  int***    trks_counts   = Get3DArray <int> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (maxNPtchBins, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (maxNPtchBins, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (maxNPtchBins, 40);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight",     &event_weight);
    PbPbTree->SetBranchAddress ("isEE",             &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",          &fcal_et);
    PbPbTree->SetBranchAddress ("impactParameter",  &ip);
    PbPbTree->SetBranchAddress ("eventPlane",       &eventPlane);
    PbPbTree->SetBranchAddress ("q2",               &q2);
    PbPbTree->SetBranchAddress ("psi2",             &psi2);
    PbPbTree->SetBranchAddress ("vz",               &vz);
    PbPbTree->SetBranchAddress ("z_pt",             &z_pt);
    PbPbTree->SetBranchAddress ("z_y",              &z_y);
    PbPbTree->SetBranchAddress ("z_phi",            &z_phi);
    PbPbTree->SetBranchAddress ("z_m",              &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",            &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",           &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",           &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",        &l1_charge);
    PbPbTree->SetBranchAddress ("l2_pt",            &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",           &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",           &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",        &l2_charge);
    PbPbTree->SetBranchAddress ("ntrk",             &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",           &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",          &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",          &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = (useImpactParameter ? GetIPCentBin (ip) : GetCentBin (fcal_et));
      //const short iCent = GetCentBin (z_fcal_et);
      if (iCent < 1 || iCent > numCentBins-1) continue;
      //const short iFineCent = GetFineCentBin (z_fcal_et);
      //if (iFineCent < 1 || iFineCent > numFineCentBins-1) continue;
      const short iFineCent = 1;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin

      if (event_weight == 0) continue;

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2) dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      if (iPtZ < 2) continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFineCent]->Fill (q2);
      h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
      h_psi2[iFineCent]->Fill (psi2);
      h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = (doSmearing ? SmearTrackMomentum (fcal_et, trk_pt[iTrk], trk_eta[iTrk], true) : trk_pt[iTrk]);
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        const float trkWeight = 1;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

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
        if (7*pi/8 <= dPhi) {
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
    cout << "Done truth-level Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight",  &event_weight);
    ppTree->SetBranchAddress ("isEE",          &isEE);
    ppTree->SetBranchAddress ("vz",            &vz);
    ppTree->SetBranchAddress ("z_pt",          &z_pt);
    ppTree->SetBranchAddress ("z_y",           &z_y);
    ppTree->SetBranchAddress ("z_phi",         &z_phi);
    ppTree->SetBranchAddress ("z_m",           &z_m);
    ppTree->SetBranchAddress ("l1_pt",         &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",        &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",        &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",     &l1_charge);
    ppTree->SetBranchAddress ("l2_pt",         &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",        &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",        &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",     &l2_charge);
    ppTree->SetBranchAddress ("ntrk",          &ntrk);
    ppTree->SetBranchAddress ("trk_pt",        &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",       &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",       &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin

      TLorentzVector ztlv;
      ztlv.SetPxPyPzE (z_pt * cos (z_phi), z_pt * sin (z_phi), sqrt (z_pt*z_pt + z_m*z_m) * sinh (z_y), sqrt (z_pt*z_pt + z_m*z_m) * cosh (z_y));
      z_eta = ztlv.Eta ();

      //nch_weight = h_ppNch_weights->GetBinContent (h_ppNch_weights->FindBin (ntrk));

      //event_weight = event_weight * vz_weight * nch_weight;
      if (event_weight == 0)
        continue;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = (doSmearing ? SmearTrackMomentum (fcal_et, trk_pt[iTrk], trk_eta[iTrk], false) : trk_pt[iTrk]);
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        const float trkWeight = 1;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

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

    } // end loop over pp tree
    cout << "Done truth-level pp loop." << endl;
  }

  Delete3DArray (&trks_counts, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights1, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights2, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete2DArray (&trks_counts_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights1_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights2_inPhi, maxNPtchBins, 40);

  SaveHists (outFileName);

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }
}


#endif
