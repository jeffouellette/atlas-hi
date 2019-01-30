#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

void ZTrackHistMaker () {
  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* tree = (TTree*)inFile->Get ("ZTrackTree");

  float fcal_et = 0, z_pt = 0, z_eta = 0, z_phi = 0, z_m = 0, event_weight = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;
  tree->SetBranchAddress ("fcal_et", &fcal_et);
  tree->SetBranchAddress ("l1_pt", &l1_pt);
  tree->SetBranchAddress ("l1_eta", &l1_eta);
  tree->SetBranchAddress ("l1_phi", &l1_phi);
  tree->SetBranchAddress ("l2_pt", &l2_pt);
  tree->SetBranchAddress ("l2_eta", &l2_eta);
  tree->SetBranchAddress ("l2_phi", &l2_phi);
  tree->SetBranchAddress ("z_pt", &z_pt);
  tree->SetBranchAddress ("z_eta", &z_eta);
  tree->SetBranchAddress ("z_phi", &z_phi);
  tree->SetBranchAddress ("z_m", &z_m);
  tree->SetBranchAddress ("trk_pt", &trk_pt);
  tree->SetBranchAddress ("trk_eta", &trk_eta);
  tree->SetBranchAddress ("trk_phi", &trk_phi);
  tree->SetBranchAddress ("event_weight", &event_weight);

  const int nXZTrkBins = 50;
  const double* xZTrkBins = logspace (2e-3, 2e1, nXZTrkBins);

  const int nPtTrkBins = 50;
  const double* ptTrkBins = logspace (1, 100, nPtTrkBins);

  const int numPhiBins = 5;
  const double phiLowBins[5] = {0, 3*pi/16, 7*pi/16, 11*pi/16, 15*pi/16};
  const double phiHighBins[5] = {pi/16, 5*pi/16, 9*pi/16, 13*pi/16, pi};

  TH2D* ZTrackXPhi = new TH2D ("ZTrackXPhi", "", 80, -pi/2, 3*pi/2, nXZTrkBins, xZTrkBins);
  TH2D* ZTrackPtPhi = new TH2D ("ZTrackPtPhi", "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
  ZTrackXPhi->Sumw2 ();

  TH1D* ZPtSpec = new TH1D ("ZPtSpec", "", 50, 0, 250);
  ZPtSpec->Sumw2 ();
  TH1D* ZMSpec = new TH1D ("ZMSpec", "", 50, 66, 116);
  ZMSpec->Sumw2 ();

  TH1D** ZTracksPt = Get1DArray<TH1D*> (5);
  TH1D** ZTracksXZ = Get1DArray<TH1D*> (5);
  for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
    ZTracksPt[iPhi] = new TH1D (Form ("ZTracksPt_iPhi%i", iPhi), "", nPtTrkBins, ptTrkBins);
    ZTracksPt[iPhi]->Sumw2 ();
    ZTracksXZ[iPhi] = new TH1D (Form ("ZTracksXZ_iPhi%i", iPhi), "", nXZTrkBins, xZTrkBins);
    ZTracksXZ[iPhi]->Sumw2 ();
  }

  const int nEvts = tree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    tree->GetEntry (iEvt);

    ZPtSpec->Fill (z_pt, event_weight);
    ZMSpec->Fill (z_m, event_weight);

    if (z_pt < 5)
      continue;

    for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
      const float trkpt = 1e-3 * trk_pt->at (iTrk);

      if (trkpt < 1)
        continue;

      const float deltaR1 = DeltaR (l1_eta, trk_eta->at (iTrk), l1_phi, trk_phi->at (iTrk));
      const float deltaPt1 = fabs (l1_pt - trkpt);

      const float deltaR2 = DeltaR (l2_eta, trk_eta->at (iTrk), l2_phi, trk_phi->at (iTrk));
      const float deltaPt2 = fabs (l2_pt - trkpt);

      //if ((deltaR1 < 0.2 && deltaPt1 < 15) || (deltaR2 < 0.2 && deltaPt2 < 15)) {
      if (deltaR1 < 0.2 || deltaR2 < 0.2) {
      //  cout << "Trk pt: " << trk_pt->at (iTrk) << ", l1_pt: " << l1_pt << ", l2_pt: " << l2_pt << endl;
        continue;
      }

      const float xZTrk = trkpt / z_pt;

      float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);

      short idPhi = 0;
      while (idPhi < numPhiBins) {
        if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
          break;
        else
          idPhi++;
      }

      if (idPhi >= 0 && idPhi < numPhiBins) {
        ZTracksPt[idPhi]->Fill (trkpt, event_weight);
        ZTracksXZ[idPhi]->Fill (xZTrk, event_weight);
      }

      dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
      if (dphi < -pi/2)
        dphi = dphi + 2*pi;

      ZTrackXPhi->Fill (dphi, xZTrk, event_weight);
      ZTrackPtPhi->Fill (dphi, trkpt, event_weight);
    }
  }

  TCanvas* c1 = new TCanvas ("c1", "", 800, 600);
  c1->cd ();

  c1->SetLogy ();

  ZTrackXPhi->GetXaxis ()->SetTitle ("#phi_{Z} - #phi_{Trk}");
  ZTrackXPhi->GetYaxis ()->SetTitle ("#it{p}_{T}^{trk} / #it{p}_{T}^{Z}");
  //ZTrackXPhi->GetYaxis ()->SetTitle ("#it{p}_{T}^{trk} cos |#Delta#phi|#left[GeV#right]");

  ZTrackXPhi->RebinY (2);
  ZTrackXPhi->Draw ("lego2");

  c1->SaveAs (Form ("%s/ZTrackCorr.pdf", plotPath.Data ()));


  TCanvas* c2 = new TCanvas ("c2", "", 800, 600);
  c2->cd ();
  c2->SetLogy ();

  const float ptBinNums[7] = {1, 5, 10, 15, 22, 36, 50};
  const int numPtBins = sizeof (ptBinNums) / sizeof (ptBinNums[0]) - 1;

  const Color_t colors[9] = {kBlack, kRed, kBlue, kMagenta, 8, kCyan+1, kOrange, kViolet, kGray};

  TH1D** ZTrackPhiPtBins = Get1DArray <TH1D*> (numPtBins);
  double min = 1e30, max = 0;
  for (int iPt = 0; iPt < numPtBins; iPt++) {
    ZTrackPhiPtBins[iPt] = ZTrackPtPhi->ProjectionX (Form ("ZTrackPhi_iPt%i", iPt), ptBinNums[iPt], ptBinNums[iPt+1]-1);
    if (ZTrackPhiPtBins[iPt]->GetMinimum () < min) min = ZTrackPhiPtBins[iPt]->GetMinimum ();
    if (ZTrackPhiPtBins[iPt]->GetMaximum () > max) max = ZTrackPhiPtBins[iPt]->GetMaximum ();
  }

  for (int iPt = 0; iPt < numPtBins; iPt++) {
    ZTrackPhiPtBins[iPt]->GetYaxis ()->SetRangeUser (0.4*min, 2.2*max);
    ZTrackPhiPtBins[iPt]->GetXaxis ()->SetTitle ("#Delta#phi");
    ZTrackPhiPtBins[iPt]->GetYaxis ()->SetTitle ("Pairs");
    ZTrackPhiPtBins[iPt]->SetLineColor (colors[iPt]);
    ZTrackPhiPtBins[iPt]->SetMarkerColor (colors[iPt]);
    if (iPt == 0)
      ZTrackPhiPtBins[iPt]->Draw ("e1");
    else
      ZTrackPhiPtBins[iPt]->Draw ("same e1");

    const float pt_lo = ZTrackPtPhi->GetYaxis ()->GetBinLowEdge (ptBinNums[iPt]);
    const float pt_hi = ZTrackPtPhi->GetYaxis ()->GetBinLowEdge (ptBinNums[iPt+1]);
    myText (0.7, 0.88-0.07*iPt, colors[iPt], Form ("%.1f < #it{p}_{T}^{trk} < %.1f GeV", pt_lo, pt_hi), 0.04);
  }

  c2->SaveAs (Form ("%s/dPhi_pTtrk.pdf", plotPath.Data ()));


  TCanvas* c3 = new TCanvas ("c3", "", 800, 600);
  c3->cd ();
  c3->SetLogy ();

  ZPtSpec->GetXaxis ()->SetTitle ("#it{p}_{T}^{Z} #left[GeV#right]");

  TH1D* ZPtSpecAboveX = new TH1D ("ZPtSpecAboveX", "", 50, 0, 250);
  for (int ix = 1; ix <= ZPtSpecAboveX->GetNbinsX (); ix++) {
    double content = 0;
    double varsq = 0;
    for (int ixprime = ix; ixprime <= ZPtSpec->GetNbinsX (); ixprime++) {
      content += ZPtSpec->GetBinContent (ixprime);
      varsq += pow (ZPtSpec->GetBinError (ixprime), 2);
    }
    ZPtSpecAboveX->SetBinContent (ix, content);
    ZPtSpecAboveX->SetBinError (ix, sqrt (varsq));
  }

  ZPtSpecAboveX->SetLineColor (kBlue);
  ZPtSpecAboveX->SetMarkerColor (kBlue);

  ZPtSpecAboveX->Draw ("e1");
  ZPtSpec->Draw ("same e1");

  myText (0.5, 0.88, kBlack, "n(#it{p}_{T}) #equiv #Z's in #it{p}_{T} bin", 0.04);
  myText (0.5, 0.78, kBlue, "N(#it{p}_{T}) #equiv #Z's at or above #it{p}_{T} bin = #int_{#it{p}_{T}}^{#infty}n(x)dx", 0.04);

  c3->SaveAs (Form ("%s/z_pt_spectrum.pdf", plotPath.Data ()));


  TCanvas* c4 = new TCanvas ("c4", "", 800, 600);
  c4->cd ();

  ZMSpec->GetXaxis ()->SetTitle ("m_{Z} #left[GeV#right]");
  ZMSpec->Draw ();

  c4->SaveAs (Form ("%s/z_mass_spectrum.pdf", plotPath.Data ()));


  TCanvas* c5 = new TCanvas ("c5", "", 800, 600);
  c5->cd ();

  const double padRatio = 1.2; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetTopMargin (0.04);
  topPad->SetBottomMargin (0);
  topPad->SetLeftMargin (0.15);
  topPad->SetRightMargin (0.04);
  bottomPad->SetTopMargin (0);
  bottomPad->SetBottomMargin (0.20);
  bottomPad->SetLeftMargin (0.15);
  bottomPad->SetRightMargin (0.04);
  topPad->Draw ();
  bottomPad->Draw ();

  topPad->cd ();
  gPad->SetLogx ();
  gPad->SetLogy ();

  bottomPad->cd ();
  gPad->SetLogx ();

  for (int iPhi = 0; iPhi < 5; iPhi++) {
    topPad->cd ();

    ZTracksPt[iPhi]->SetMarkerColor (colors[iPhi]);
    ZTracksPt[iPhi]->SetLineColor (colors[iPhi]);
    ZTracksPt[iPhi]->SetMarkerSize (0.75);
    //ZTracksPt[iPhi]->GetXaxis ()->SetTitle ("#x_{Z}^{trk}");
    ZTracksPt[iPhi]->GetXaxis ()->SetTitle ("#it{p}_{T}^{trk} #left[GeV#right]");
    ZTracksPt[iPhi]->GetYaxis ()->SetTitle ("Pairs / d#Delta#phi");
    ZTracksPt[iPhi]->GetXaxis ()->SetTitleSize (0.04 / uPadY);
    ZTracksPt[iPhi]->GetYaxis ()->SetTitleSize (0.04 / uPadY);
    ZTracksPt[iPhi]->GetXaxis ()->SetLabelSize (0.04 / uPadY);
    ZTracksPt[iPhi]->GetYaxis ()->SetLabelSize (0.04 / uPadY);
    ZTracksPt[iPhi]->GetYaxis ()->SetTitleOffset (1.5 * uPadY);
    ZTracksPt[iPhi]->Scale (1. / (phiHighBins[iPhi] - phiLowBins[iPhi]));

    if (iPhi == 0)
      ZTracksPt[iPhi]->Draw ("e1");
    else
      ZTracksPt[iPhi]->Draw ("same e1");

    const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/16 ? "#pi/16" : Form ("%.0f#pi/16", phiLowBins[iPhi]*16/pi)) : "0";
    const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/16 ? "#pi/16" : Form ("%.0f#pi/16", phiHighBins[iPhi]*16/pi)) : "#pi";

    myText (0.7, 0.88-0.07*iPhi, colors[iPhi], Form ("%s < #Delta#phi < %s", lo, hi), 0.04 / uPadY);

    bottomPad->cd ();

    TH1D* nearSideRatio_hist = new TH1D (Form ("nearSideRatioPt_iPhi%i", iPhi), "", nPtTrkBins, ptTrkBins);
    nearSideRatio_hist->Sumw2 ();
    nearSideRatio_hist->Add (ZTracksPt[iPhi]);
    nearSideRatio_hist->Divide (ZTracksPt[0]);

    TGraphAsymmErrors* nearSideRatio = make_graph (nearSideRatio_hist);
    deltaize (nearSideRatio, 1+(2.5-iPhi)*0.01, true);

    nearSideRatio->SetLineColor (colors[iPhi]);
    nearSideRatio->SetMarkerColor (colors[iPhi]);
    nearSideRatio->SetMarkerSize (0.75);
    nearSideRatio->GetXaxis ()->SetTitle ("#it{p}_{T}^{trk} #left[GeV#right]");
    nearSideRatio->GetYaxis ()->SetTitle ("Ratio to #color[1]{Black}");
    nearSideRatio->GetXaxis ()->SetTitleSize (0.04 / dPadY);
    nearSideRatio->GetYaxis ()->SetTitleSize (0.04 / dPadY);
    nearSideRatio->GetXaxis ()->SetLabelSize (0.04 / dPadY);
    nearSideRatio->GetYaxis ()->SetLabelSize (0.04 / dPadY);
    nearSideRatio->GetXaxis ()->SetTitleOffset (1);
    nearSideRatio->GetYaxis ()->SetTitleOffset (1.5 * dPadY);

    if (iPhi == 0)
      nearSideRatio->Draw ("e1");
    else
      nearSideRatio->Draw ("same e1");
    
  }

  topPad->cd ();
  myText (0.22, 0.28, kBlack, "dR_{1,2} > 0.2", 0.04 / uPadY);// or d#it{p}_{T}^{1,2} > 5 GeV", 0.04 / uPadY);

  c5->SaveAs (Form ("%s/pTtrk_dPhi.pdf", plotPath.Data ()));


  TCanvas* c6 = new TCanvas ("c6", "", 800, 600);
  c6->cd ();

  topPad->Draw ();
  bottomPad->Draw ();

  for (int iPhi = 0; iPhi < 5; iPhi++) {
    topPad->cd ();

    ZTracksXZ[iPhi]->SetMarkerColor (colors[iPhi]);
    ZTracksXZ[iPhi]->SetLineColor (colors[iPhi]);
    ZTracksXZ[iPhi]->SetMarkerSize (0.75);
    //ZTracksXZ[iPhi]->GetXaxis ()->SetTitle ("#x_{Z}^{trk}");
    ZTracksXZ[iPhi]->GetXaxis ()->SetTitle ("#it{p}_{T}^{trk} #left[GeV#right]");
    ZTracksXZ[iPhi]->GetYaxis ()->SetTitle ("Pairs / d#Delta#phi");
    ZTracksXZ[iPhi]->GetXaxis ()->SetTitleSize (0.04 / uPadY);
    ZTracksXZ[iPhi]->GetYaxis ()->SetTitleSize (0.04 / uPadY);
    ZTracksXZ[iPhi]->GetXaxis ()->SetLabelSize (0.04 / uPadY);
    ZTracksXZ[iPhi]->GetYaxis ()->SetLabelSize (0.04 / uPadY);
    ZTracksXZ[iPhi]->GetYaxis ()->SetTitleOffset (1.5 * uPadY);
    ZTracksXZ[iPhi]->Scale (1. / (phiHighBins[iPhi] - phiLowBins[iPhi]));

    if (iPhi == 0)
      ZTracksXZ[iPhi]->Draw ("e1");
    else
      ZTracksXZ[iPhi]->Draw ("same e1");

    const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/16 ? "#pi/16" : Form ("%.0f#pi/16", phiLowBins[iPhi]*16/pi)) : "0";
    const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/16 ? "#pi/16" : Form ("%.0f#pi/16", phiHighBins[iPhi]*16/pi)) : "#pi";

    myText (0.7, 0.88-0.07*iPhi, colors[iPhi], Form ("%s < #Delta#phi < %s", lo, hi), 0.04 / uPadY);

    bottomPad->cd ();

    TH1D* nearSideRatio_hist = new TH1D (Form ("nearSideRatioXZ_iPhi%i", iPhi), "", nXZTrkBins, xZTrkBins);
    nearSideRatio_hist->Sumw2 ();
    nearSideRatio_hist->Add (ZTracksXZ[iPhi]);
    nearSideRatio_hist->Divide (ZTracksXZ[0]);

    TGraphAsymmErrors* nearSideRatio = make_graph (nearSideRatio_hist);
    deltaize (nearSideRatio, 1+(2.5-iPhi)*0.01, true);

    nearSideRatio->SetLineColor (colors[iPhi]);
    nearSideRatio->SetMarkerColor (colors[iPhi]);
    nearSideRatio->SetMarkerSize (0.75);
    nearSideRatio->GetXaxis ()->SetTitle ("x_{Z}^{trk}");
    nearSideRatio->GetYaxis ()->SetTitle ("Ratio to #color[1]{Black}");
    nearSideRatio->GetXaxis ()->SetTitleSize (0.04 / dPadY);
    nearSideRatio->GetYaxis ()->SetTitleSize (0.04 / dPadY);
    nearSideRatio->GetXaxis ()->SetLabelSize (0.04 / dPadY);
    nearSideRatio->GetYaxis ()->SetLabelSize (0.04 / dPadY);
    nearSideRatio->GetXaxis ()->SetTitleOffset (1);
    nearSideRatio->GetYaxis ()->SetTitleOffset (1.5 * dPadY);

    if (iPhi == 0)
      nearSideRatio->Draw ("e1");
    else
      nearSideRatio->Draw ("same e1");
    
  }
  c6->SaveAs (Form ("%s/xZTrk_dPhi.pdf", plotPath.Data ()));


}
