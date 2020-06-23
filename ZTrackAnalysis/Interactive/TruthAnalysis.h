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

  void ComparePPYields (TruthAnalysis* ta);
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

  const int xbin = h2->GetXaxis ()->FindFixBin (truth_pt);
  const int ybin = h2->GetYaxis ()->FindFixBin (truth_eta);
  if (xbin < 1 || h2->GetXaxis ()->GetNbins () < xbin)
    return 0;
  if (ybin < 1 || h2->GetYaxis ()->GetNbins () < ybin)
    return 0;

  const double tmr = h2->GetBinContent (xbin, ybin);

  if (tmr == 0)
    return truth_pt;

  return tmr_rndm->Gaus (1, tmr/100.) * truth_pt;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
// Load histograms into memory, then combine channels.
//////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);
  PhysicsAnalysis :: CombineHists ();
  PhysicsAnalysis :: TruncatePhysicsPlots ();
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




//////////////////////////////////////////////////////////////////////////////////////////////////
// Compare pp yields with a variation of the truth analysis.
// Designed for studying potential impact of hadron yield unfolding.
//////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: ComparePPYields (TruthAnalysis* ta) {
  TLatex* tl = new TLatex ();

  TCanvas* c_truth_pp_compare = new TCanvas ("c_truth_pp_compare", "", 800, 880);

  const double lMargin = 0.15;
  const double rMargin = 0.04;
  const double utMargin = 0.04;
  const double ubMargin = 0;//0.01;
  const double dtMargin = 0;//0.02;
  const double dbMargin = 0.35;

  TPad* uPad = new TPad ("c_truth_pp_compare_uPad", "", 0, 0.3, 1, 1);
  TPad* dPad = new TPad ("c_truth_pp_compare_dPad", "", 0, 0, 1, 0.3);

  uPad->SetLeftMargin (lMargin);
  uPad->SetRightMargin (rMargin);
  dPad->SetLeftMargin (lMargin);
  dPad->SetRightMargin (rMargin);
  uPad->SetTopMargin (utMargin);
  uPad->SetBottomMargin (ubMargin);
  dPad->SetTopMargin (dtMargin);
  dPad->SetBottomMargin (dbMargin);

  uPad->Draw ();
  dPad->Draw ();

  uPad->cd ();
  uPad->SetLogx ();
  uPad->SetLogy ();

  {
    TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
    xax->SetTitleSize (0);
    xax->SetTickLength (0.03 * (1.-utMargin-ubMargin) / uPad->GetHNDC ());
    xax->SetTickLength (0.03);
    xax->SetLabelSize (0);

    xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
    yax->SetTitleFont (43);
    yax->SetTitleSize (36);
    yax->SetTickLength (0.02 * (1.-utMargin-ubMargin) / uPad->GetHNDC ());
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    double ymin = 2e-3;
    double ymax = 8e3;
    yax->SetRangeUser (ymin, ymax);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);
  }

  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
    const int iCent = 0;

    TGAE* g = make_graph (h_trk_pt_ptz[2][iPtZ][iCent]);

    Style_t markerStyle = kOpenCircle;
    float markerSize = 1.6;

    OffsetYAxis (g, pow (10, iPtZ-2), true);
    RecenterGraph (g);
    ResetXErrors (g);
    deltaize (g, 0.95, true);
    ResetXErrors (g);

    g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (markerSize);
    g->SetLineWidth (3);
    g->SetMarkerColor (finalColors[iPtZ-1]);
    g->SetLineColor (finalColors[iPtZ-1]);

    ((TGAE*) g->Clone ())->Draw ("P");

    //markerStyle = kDot;

    //g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (0);
    //g->SetMarkerColor (kBlack);

    ((TGAE*) g->Clone ())->Draw ("P");

    SaferDelete (&g);
  } // end loop over iPtZ


  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
    const int iCent = 0;

    TGAE* g = make_graph (ta->h_trk_pt_ptz[2][iPtZ][iCent]);

    Style_t markerStyle = kOpenSquare;
    float markerSize = 1.6;

    OffsetYAxis (g, pow (10, iPtZ-2), true);
    RecenterGraph (g);
    ResetXErrors (g);
    deltaize (g, 1.05, true);
    ResetXErrors (g);

    g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (markerSize);
    g->SetLineWidth (3);
    g->SetMarkerColor (finalColors[iPtZ-1]);
    g->SetLineColor (finalColors[iPtZ-1]);

    ((TGAE*) g->Clone ())->Draw ("P");

    //markerStyle = kDot;

    //g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (0);
    //g->SetMarkerColor (kBlack);

    ((TGAE*) g->Clone ())->Draw ("P");

    SaferDelete (&g);
  } // end loop over iPtZ

  tl->SetTextFont (43);
  tl->SetTextColor (kBlack);
  tl->SetTextAlign (11);

  tl->SetTextSize (28);
  tl->DrawLatexNDC (0.20, 0.14, "#bf{#it{ATLAS}} Internal");
  tl->SetTextSize (28);
  tl->DrawLatexNDC (0.20, 0.09, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  tl->DrawLatexNDC (0.20, 0.04, "Powheg + Pythia 8.186");

  tl->SetTextSize (22);
  tl->DrawLatexNDC (0.730, 0.885, "Truth #it{p}_{T}^{Z} [GeV]");
  tl->DrawLatexNDC (0.730, 0.835, "15-30 (#times 1)");
  tl->DrawLatexNDC (0.730, 0.785, "30-60 (#times 10)");
  tl->DrawLatexNDC (0.730, 0.735, "60+ (#times 10^{2})");
  MakeDataBox   (0.56, 0.840, finalFillColors[1], 0.00, kOpenCircle, 1.6);
  MakeDataBox   (0.56, 0.790, finalFillColors[2], 0.00, kOpenCircle, 1.6);
  MakeDataBox   (0.56, 0.740, finalFillColors[3], 0.00, kOpenCircle, 1.6);
  MakeDataBox   (0.70, 0.840, finalFillColors[1], 0.00, kOpenSquare, 1.6);
  MakeDataBox   (0.70, 0.790, finalFillColors[2], 0.00, kOpenSquare, 1.6);
  MakeDataBox   (0.70, 0.740, finalFillColors[3], 0.00, kOpenSquare, 1.6);
  //myMarkerAndBoxAndLineText (0.72, 0.840, 1.4, 1001, finalFillColors[1], 0.00, finalColors[1], kOpenSquare, 1.6, "", 0.036);
  //myMarkerAndBoxAndLineText (0.63, 0.840, 1.4, 1001, finalFillColors[1], 0.00, finalColors[1], kOpenCircle, 1.6, "", 0.036);
  //myMarkerAndBoxAndLineText (0.72, 0.790, 1.4, 1001, finalFillColors[2], 0.00, finalColors[2], kOpenSquare, 1.6, "", 0.036);
  //myMarkerAndBoxAndLineText (0.63, 0.790, 1.4, 1001, finalFillColors[2], 0.00, finalColors[2], kOpenCircle, 1.6, "", 0.036);
  //myMarkerAndBoxAndLineText (0.72, 0.740, 1.4, 1001, finalFillColors[3], 0.00, finalColors[3], kOpenSquare, 1.6, "", 0.036);
  //myMarkerAndBoxAndLineText (0.63, 0.740, 1.4, 1001, finalFillColors[3], 0.00, finalColors[3], kOpenCircle, 1.6, "", 0.036);

  tl->SetTextSize (22);
  tl->DrawLatexNDC (0.450, 0.885, "Nominal");
  tl->SetTextSize (22);
  tl->DrawLatexNDC (0.585, 0.885, "Smeared");
  dPad->cd ();
  dPad->SetLogx ();

  {
    TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

    TAxis* xax = h->GetXaxis ();
    TAxis* yax = h->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetTitleFont (43);
    xax->SetTitleSize (36);
    xax->SetTitleOffset (3.6);
    xax->SetTickLength (0.03 * (1.-dtMargin-dbMargin) / dPad->GetHNDC ());
    xax->SetLabelSize (0);
    xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

    yax->SetTitle ("Ratio");
    yax->SetTitleFont (43);
    yax->SetTitleSize (36);
    yax->CenterTitle (true);
    yax->SetTickLength (0.02 * (1.-dtMargin-dbMargin) / dPad->GetHNDC ());
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    double ymin = 0.993;
    double ymax = 1.007;
    yax->SetRangeUser (ymin, ymax);
    yax->SetNdivisions (504);

    h->SetLineWidth (0);

    h->DrawCopy ("");
    SaferDelete (&h);

    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextAlign (21);
    const double yoff = ymin - (0.12 * (ymax - ymin) / (1.-dtMargin-dbMargin));
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");

    TLine* line = new TLine (pTchBins[nPtZBins-1][0], 1., pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]], 1.);
    line->SetLineColor (kBlack);
    line->SetLineStyle (2);
    line->Draw ("same");
  }

  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
    const int iCent = 0;
    TGAE* this_graph = make_graph (h_trk_pt_ptz[2][iPtZ][iCent]);
    TGAE* that_graph = make_graph (ta->h_trk_pt_ptz[2][iPtZ][iCent]);

    RecenterGraph (this_graph);
    RecenterGraph (that_graph);

    TGAE* g_ratio = new TGAE ();

    double x, y_num, y_den;
    double y_this_hi, y_this_lo, y_that_hi, y_that_lo;
    for (int iX = 0; iX < this_graph->GetN (); iX++) {
      that_graph->GetPoint (iX, x, y_num);
      this_graph->GetPoint (iX, x, y_den);

      y_this_hi = this_graph->GetErrorYhigh (iX);
      y_this_lo = this_graph->GetErrorYlow (iX);
      y_that_hi = that_graph->GetErrorYhigh (iX);
      y_that_lo = that_graph->GetErrorYlow (iX);

      const double ratio = y_num / y_den;

      g_ratio->SetPoint (g_ratio->GetN (), x, ratio);
      g_ratio->SetPointEYhigh (g_ratio->GetN () - 1, fabs (ratio) * sqrt (fabs (pow (y_that_hi/y_num, 2) - pow (y_this_hi/y_den, 2))));
      g_ratio->SetPointEYlow (g_ratio->GetN () - 1, fabs (ratio) * sqrt (fabs (pow (y_that_lo/y_num, 2) - pow (y_this_lo/y_den, 2))));
    }

    Style_t markerStyle = kOpenCircle;
    float markerSize = 1.6;

    g_ratio->SetMarkerStyle (markerStyle);
    g_ratio->SetMarkerSize (markerSize);
    g_ratio->SetLineWidth (3);
    g_ratio->SetMarkerColor (finalColors[iPtZ-1]);
    g_ratio->SetLineColor (finalColors[iPtZ-1]);

    ((TGAE*) g_ratio->Clone ())->Draw ("P");

    markerStyle = kDot;

    g_ratio->SetMarkerStyle (markerStyle);
    g_ratio->SetMarkerSize (markerSize);
    g_ratio->SetMarkerColor (finalColors[iPtZ-1]);

    ((TGAE*) g_ratio->Clone ())->Draw ("P");

    SaferDelete (&g_ratio);

    SaferDelete (&this_graph);
    SaferDelete (&that_graph);
  }

  c_truth_pp_compare->SaveAs ("../Plots/TrackMomentumResolutionStudy/truth_comparison_pTch.pdf");

}


#endif
