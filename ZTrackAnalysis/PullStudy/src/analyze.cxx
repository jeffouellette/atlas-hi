#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <string>

#include <math.h>

//#include "eps09_cxx/eps09.h"

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <GlobalParams.h>
#include <Utilities.h>

using namespace atlashi;

string FormatCounts (int counts) {
  if (counts < 1000) return "";
  else if (1000 <= counts && counts < 10000) {
    string countsStr = FormatMeasurement (counts, 0, 1);
    countsStr = countsStr.substr(0, 1) + "k";
    return countsStr;
  }
  else if (10000 <= counts && counts < 100000) {
    string countsStr = FormatMeasurement (counts, 0, 2);
    countsStr = countsStr.substr(0, 2) + "k";
    return countsStr;
  }
  else if (100000 <= counts && counts < 1000000) {
    string countsStr = FormatMeasurement (counts, 0, 3);
    countsStr = countsStr.substr(0, 3) + "k";
    return countsStr;
  }
  else return "";
}


int GetNumInGroup1 (short iPtZ, short iCent) {
  if (iPtZ == 2) {
    if (iCent == 1)       return 260;
    else if (iCent == 2)  return 673;
    else if (iCent == 3)  return 782;
  }
  else if (iPtZ == 3) {
    if (iCent == 1)       return 148;
    else if (iCent == 2)  return 364;
    else if (iCent == 3)  return 443;
  }
  else if (iPtZ == 4) {
    if (iCent == 1)       return 52;
    else if (iCent == 2)  return 147;
    else if (iCent == 3)  return 145;
  }
  return 0;
}


int GetNumInGroup2 (short iPtZ, short iCent) {
  if (iPtZ == 2) {
    if (iCent == 1)       return 397;
    else if (iCent == 2)  return 828;
    else if (iCent == 3)  return 805;
  }
  else if (iPtZ == 3) {
    if (iCent == 1)       return 206;
    else if (iCent == 2)  return 484;
    else if (iCent == 3)  return 406;
  }
  else if (iPtZ == 4) {
    if (iCent == 1)       return 89;
    else if (iCent == 2)  return 171;
    else if (iCent == 3)  return 173;
  }
  return 0;
}


int main () {

  SetAtlasStyle();

  const int nSeeds = 80;
  TFile* inFile;
  TTree* inTree;

  //int code = 0;
  //int id1 = 0;
  //int id2 = 0;
  //float x1pdf = 0;
  //float x2pdf = 0;
  //float Q = 0;
  //bool isValence1 = false;
  //bool isValence2 = false;

  float z_pt = 0;
  float z_eta = 0;
  float z_phi = 0;
  float z_m = 0;
  int part_n = 0;
  float part_pt[10000];
  float part_eta[10000];
  float part_phi[10000];


  float trackPtYields[2][6];
  float trackPtYieldWeightsSq[2][6];
  float trackXYields[2][6];
  float trackXYieldWeightsSq[2][6];

  TH1D* h_trk_pth_pull[3][3][6];
  TH1D* h_trk_xhz_pull[3][3][6];

  for (int iPtZ : {2, 3, 4}) {
    for (int iCent : {1, 2, 3}) {

      for (int iX = 0; iX < 6; iX++) {
        h_trk_pth_pull[iPtZ-2][iCent-1][iX] = new TH1D (Form ("h_trk_pth_pull_iPtZ%i_iCent%i_iPtTrk%i", iPtZ, iCent, iX), ";Pull;Counts", 12, -4, 4);
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX] = new TH1D (Form ("h_trk_xhz_pull_iPtZ%i_iCent%i_iXHZ%i", iPtZ, iCent, iX), ";Pull;Counts", 12, -4, 4);
      }

      const int nGroup1 = GetNumInGroup1 (iPtZ, iCent);
      const int nGroup2 = GetNumInGroup2 (iPtZ, iCent);

      for (int iSeed = 0; iSeed < nSeeds; iSeed++) {

        inFile = new TFile (Form ("output/seed%i_iPtZ%i_iCent%i.root", iSeed, iPtZ, iCent), "read");
        inTree = (TTree*) inFile->Get ("tree");

        //inTree->SetBranchAddress ("code",       &code);
        //inTree->SetBranchAddress ("id1",        &id1);
        //inTree->SetBranchAddress ("id2",        &id2);
        //inTree->SetBranchAddress ("x1pdf",      &x1pdf);
        //inTree->SetBranchAddress ("x2pdf",      &x2pdf);
        //inTree->SetBranchAddress ("Q",          &Q);
        //inTree->SetBranchAddress ("isValence1", &isValence1);
        //inTree->SetBranchAddress ("isValence2", &isValence2);
        inTree->SetBranchAddress ("z_pt",       &z_pt);
        inTree->SetBranchAddress ("z_eta",      &z_eta);
        inTree->SetBranchAddress ("z_phi",      &z_phi);
        inTree->SetBranchAddress ("z_m",        &z_m);
        inTree->SetBranchAddress ("part_n",     &part_n);
        inTree->SetBranchAddress ("part_pt",    &part_pt);
        inTree->SetBranchAddress ("part_eta",   &part_eta);
        inTree->SetBranchAddress ("part_phi",   &part_phi);

        const int nEvents = inTree->GetEntries ();
        assert (nEvents == nGroup1+nGroup2);

        for (int iEvent = 0; iEvent < nEvents; iEvent++) {
          inTree->GetEntry (iEvent);

          const int iGroup = (iEvent < nGroup1 ? 0 : 1); // 0 corresponds to group 1, 1 to group 2

          float tracks_1to2GeV = 0;
          float tracks_2to4GeV = 0;
          float tracks_4to8GeV = 0;
          float tracks_8to15GeV = 0;
          float tracks_15to30GeV = 0;
          float tracks_30to60GeV = 0;
          float tracks_1_60to1_30 = 0;
          float tracks_1_30to1_15 = 0;
          float tracks_1_15to1_8 = 0;
          float tracks_1_8to1_4 = 0;
          float tracks_1_4to1_2 = 0;
          float tracks_1_2to1 = 0;
      
          for (int iPart = 0; iPart < part_n; iPart++) {
            if (fabs (part_eta[iPart]) > 2.5)
              continue;
            if (DeltaPhi (part_phi[iPart], z_phi) < 3*pi/4)
              continue;

            const float trkpt = part_pt[iPart];
            const float xhz = trkpt / z_pt;

            if (30 <= trkpt && trkpt < 60) tracks_30to60GeV += 1.;
            else if (15 <= trkpt) tracks_15to30GeV += 1.;
            else if (8 <= trkpt) tracks_8to15GeV += 1.;
            else if (4 <= trkpt) tracks_4to8GeV += 1.;
            else if (2 <= trkpt) tracks_2to4GeV += 1.;
            else if (1 <= trkpt) tracks_1to2GeV += 1.;

            if (1./60. <= xhz && xhz <= 1./30.) tracks_1_60to1_30 += 1.;
            else if (xhz <= 1./15.) tracks_1_30to1_15 += 1.;
            else if (xhz <= 1./8.) tracks_1_15to1_8 += 1.;
            else if (xhz <= 1./4.) tracks_1_8to1_4 += 1.;
            else if (xhz <= 1./2.) tracks_1_4to1_2 += 1.;
            else if (xhz <= 1.) tracks_1_2to1 += 1.;
          } // end loop over iPart

          trackPtYields[iGroup][0] += tracks_1to2GeV;
          trackPtYields[iGroup][1] += tracks_2to4GeV;
          trackPtYields[iGroup][2] += tracks_4to8GeV;
          trackPtYields[iGroup][3] += tracks_8to15GeV;
          trackPtYields[iGroup][4] += tracks_15to30GeV;
          trackPtYields[iGroup][5] += tracks_30to60GeV;
          trackPtYieldWeightsSq[iGroup][0] += pow (tracks_1to2GeV, 2);
          trackPtYieldWeightsSq[iGroup][1] += pow (tracks_2to4GeV, 2);
          trackPtYieldWeightsSq[iGroup][2] += pow (tracks_4to8GeV, 2);
          trackPtYieldWeightsSq[iGroup][3] += pow (tracks_8to15GeV, 2);
          trackPtYieldWeightsSq[iGroup][4] += pow (tracks_15to30GeV, 2);
          trackPtYieldWeightsSq[iGroup][5] += pow (tracks_30to60GeV, 2);

          trackXYields[iGroup][0] += tracks_1_60to1_30;
          trackXYields[iGroup][1] += tracks_1_30to1_15;
          trackXYields[iGroup][2] += tracks_1_15to1_8;
          trackXYields[iGroup][3] += tracks_1_8to1_4;
          trackXYields[iGroup][4] += tracks_1_4to1_2;
          trackXYields[iGroup][5] += tracks_1_2to1;
          trackXYieldWeightsSq[iGroup][0] += pow (tracks_1_60to1_30, 2);
          trackXYieldWeightsSq[iGroup][1] += pow (tracks_1_30to1_15, 2);
          trackXYieldWeightsSq[iGroup][2] += pow (tracks_1_15to1_8, 2);
          trackXYieldWeightsSq[iGroup][3] += pow (tracks_1_8to1_4, 2);
          trackXYieldWeightsSq[iGroup][4] += pow (tracks_1_4to1_2, 2);
          trackXYieldWeightsSq[iGroup][5] += pow (tracks_1_2to1, 2);
        } // end loop over iEvent

        inFile->Close ();
        SaferDelete (inFile);

        // calculate pull on per-Z yields
        float perZYield1, perZYield2, sigmaPerZYield1, sigmaPerZYield2;
        for (int iX = 0; iX < 6; iX++) {
          perZYield1 = trackPtYields[0][iX] / nGroup1;
          sigmaPerZYield1 = sqrt (trackPtYieldWeightsSq[0][iX]) / nGroup1;
          perZYield2 = trackPtYields[1][iX] / nGroup2;
          sigmaPerZYield2 = sqrt (trackPtYieldWeightsSq[1][iX]) / nGroup2;

          const float pullPt = (perZYield1 - perZYield2) / sqrt (sigmaPerZYield1*sigmaPerZYield1 + sigmaPerZYield2*sigmaPerZYield2);
          h_trk_pth_pull[iPtZ-2][iCent-1][iX]->Fill (pullPt);

          perZYield1 = trackXYields[0][iX] / nGroup1;
          sigmaPerZYield1 = sqrt (trackXYieldWeightsSq[0][iX]) / nGroup1;
          perZYield2 = trackXYields[1][iX] / nGroup2;
          sigmaPerZYield2 = sqrt (trackXYieldWeightsSq[1][iX]) / nGroup2;

          const float pullX = (perZYield1 - perZYield2) / sqrt (sigmaPerZYield1*sigmaPerZYield1 + sigmaPerZYield2*sigmaPerZYield2);
          h_trk_xhz_pull[iPtZ-2][iCent-1][iX]->Fill (pullX);
        }

        for (int iX = 0; iX < 6; iX++) {
          trackPtYields[0][iX]  = 0.;
          trackPtYields[1][iX]  = 0.;
          trackPtYieldWeightsSq[0][iX]  = 0.;
          trackPtYieldWeightsSq[1][iX]  = 0.;
          trackXYields[0][iX]   = 0.;
          trackXYields[1][iX]   = 0.;
          trackXYieldWeightsSq[0][iX]   = 0.;
          trackXYieldWeightsSq[1][iX]   = 0.;
        }

      } // end loop over seeds

    } // end loop over iCent
  } // end loop over iPtZ


  TCanvas* c_pull_pth_iPtZ2 = new TCanvas ("c_pull_pth_iPtZ2", "", 800, 600);
  c_pull_pth_iPtZ2->Divide (4, 3);
  //TCanvas* c_pull_xhz_iPtZ2 = new TCanvas ("c_pull_xhz_iPtZ2", "", 800, 600);
  //c_pull_xhz_iPtZ2->Divide (4, 3);
  TCanvas* c_pull_pth_iPtZ3 = new TCanvas ("c_pull_pth_iPtZ3", "", 1000, 600);
  c_pull_pth_iPtZ3->Divide (5, 3);
  //TCanvas* c_pull_xhz_iPtZ3 = new TCanvas ("c_pull_xhz_iPtZ3", "", 1000, 600);
  //c_pull_xhz_iPtZ3->Divide (5, 3);
  TCanvas* c_pull_pth_iPtZ4 = new TCanvas ("c_pull_pth_iPtZ4", "", 1200, 600);
  c_pull_pth_iPtZ4->Divide (6, 3);
  //TCanvas* c_pull_xhz_iPtZ4 = new TCanvas ("c_pull_xhz_iPtZ4", "", 1200, 600);
  //c_pull_xhz_iPtZ4->Divide (6, 3);

  for (int iX = 0; iX < 4; iX++) {
    for (int iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ2->cd (4*(iCent-1)+iX+1);
      h_trk_pth_pull[0][iCent-1][iX]->Draw ("hist");

      myText (0.5, 0.8, kBlack, Form ("#sigma = %.2f", h_trk_pth_pull[0][iCent-1][iX]->GetStdDev ()), 0.1);

      //c_pull_xhz_iPtZ2->cd (4*(iCent-1)+iX+1);
      //h_trk_xhz_pull[0][iCent-1][iX+2]->Draw ("hist");
    } // end loop over iCent
  } // end loop over x-values
  for (int iX = 0; iX < 5; iX++) {
    for (int iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ3->cd (5*(iCent-1)+iX+1);
      h_trk_pth_pull[1][iCent-1][iX]->Draw ("hist");

      myText (0.5, 0.8, kBlack, Form ("#sigma = %.2f", h_trk_pth_pull[1][iCent-1][iX]->GetStdDev ()), 0.1);

//      c_pull_xhz_iPtZ3->cd (5*(iCent-1)+iX+1);
//      h_trk_xhz_pull[1][iCent-1][iX+1]->Draw ("hist");
    } // end loop over iCent
  } // end loop over x-values
  for (int iX = 0; iX < 6; iX++) {
    for (int iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ4->cd (6*(iCent-1)+iX+1);
      h_trk_pth_pull[2][iCent-1][iX]->Draw ("hist");

      myText (0.5, 0.8, kBlack, Form ("#sigma = %.2f", h_trk_pth_pull[2][iCent-1][iX]->GetStdDev ()), 0.1);

//      c_pull_xhz_iPtZ4->cd (6*(iCent-1)+iX+1);
//      h_trk_xhz_pull[2][iCent-1][iX]->Draw ("hist");
    } // end loop over iCent
  } // end loop over x-values



  return 0;
}

void analyze () {
  main ();
}
