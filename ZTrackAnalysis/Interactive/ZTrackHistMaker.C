#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

const Color_t colors[10] =    {kBlack, kRed, kBlue, kGreen+2, kMagenta, 8, kCyan+1, kOrange, kViolet, kGray};
const Color_t fillColors[10] = {kGray, kRed-10, kBlue-10, kGreen-10, kMagenta-10, kGreen-2, kCyan-6, kOrange, kViolet-2, kGray};
const float fillAlpha = 1;

const int nXZTrkBins = 75;
const double* xZTrkBins = logspace (2e-3, 2e1, nXZTrkBins);

const int nPtTrkBins = 6;
const double* ptTrkBins = logspace (3, 80, nPtTrkBins);

double* doubleLogSpace (const double lo, const double hi, const int num) {
  double* arr = Get1DArray <double> (2*num+1);
  double* logarr = logspace (lo, hi, num-1);
  for (int i = 0; i < num; i++) {
    arr[num+i+1] = logarr[i];
    arr[i] = -logarr[num-1-i];
  }
  arr[num] = 0;
  delete logarr;
  return arr;
}

const int numZMissingPtBins = 16;
const double* zMissingPtBins = doubleLogSpace (1, 150, 8);
//const double* zMissingPtBins = linspace (-150, 150, numZMissingPtBins);
const double phiTrkBins[3] = {0, pi/8, pi/2};
const int numPhiTrkBins = sizeof (phiTrkBins) / sizeof (phiTrkBins[0]) - 1;

// Centrality cuts in GeV:  80%  -  30%  -   0%
//const double centBins[3] = {63.719, 1368.75, 5000};
//const int centCuts[3] =    {80,     30,      0};

// Centrality cuts in GeV:  80%  -  40%  -  15%  -   0%
//const double centBins[4] = {63.719, 875.41, 2476.58, 5000};
//const int centCuts[4] =    {80,     40,     15,      0};

// Centrality cuts in GeV:  80%  -  30%  -  10%  -   0%
const double centBins[4] = {63.719, 1368.75, 2989.31, 5000};
const int centCuts[4] =    {80,     30,     10,      0};

// Centrality cuts in GeV:    80%   -   70%   -   60%  -   50%  -   40%  -   30%  -   20%  -   10%  -   0%
//const double centBins[9]   = {63.719,   144.14,   289.595, 525.092, 875.41,  1368.75, 2046.51, 2989.31, 5000};
//const double centCuts[9]   = {80,       70,       60,      50,      40,      30,      20,      10,      0};
//const double centThicks[8] = {0.223102, 0.564885, 1.2811,  2.63435, 4.94611, 8.63844, 14.3345, 23.3523};
const int numCentBins = sizeof (centBins) / sizeof (centBins[0]); // no minus 1 to include pp bin

const double phiLowBins[4] = {0, pi/2, 7*pi/8, 0};
const double phiHighBins[4] = {pi/2, 7*pi/8, pi, pi};
const int numPhiBins = sizeof (phiLowBins) / sizeof (phiLowBins[0]);

const double zPtBins[3] = {0, 5, 10000};
const int nZPtBins = sizeof (zPtBins) / sizeof (zPtBins[0]) - 1;

// Analysis checks
TH1D** FCalSpec           = Get1DArray <TH1D*> (2);                             // iData
TH2D** FCalQ2Corr         = Get1DArray <TH2D*> (2);                             // iData
TH1D*** ZPhiYields        = Get2DArray <TH1D*> (2, numCentBins);                // iData, iCent
TH1D**** ZPtSpecs         = Get3DArray <TH1D*> (2, numCentBins, 3);             // iData, iCent, iSpc
TH1D**** ZMYields         = Get3DArray <TH1D*> (2, numCentBins, 3);             // iData, iCent, iSpc
TH1D*** ElectronSpec      = Get2DArray <TH1D*> (2, numCentBins);                // iData, iCent
TH1D*** MuonSpec          = Get2DArray <TH1D*> (2, numCentBins);                // iData, iCent
TH1D**** TrackSpec        = Get3DArray <TH1D*> (2, numCentBins, 3);             // iData, iCent, iSpc
//TH2D**** DRDists          = Get3DArray <TH2D*> (2, numCentBins, 3);           // iData, iCent, iSpc
TH2D**** ZTrackPtPhi      = Get3DArray <TH2D*> (2, numCentBins, 3);             // iData, iCent, iSpc (0=ee, 1=mumu, 2=combined)
TH1D***** ZTrackPhiPtBins = Get4DArray <TH1D*> (2, nPtTrkBins, numCentBins, 3); // iData, iPtTrk, iCent, iSpc

// Physics plots
TH2D****** ZMissingPt         = Get5DArray <TH2D*> (3, nZPtBins, 2, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iData, iPhi, iCent
TH1D****** ZMissingPtAvgs     = Get5DArray <TH1D*> (3, nZPtBins, 2, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iData, iPhi, iCent
TH1D****** ZMissingPtInts     = Get5DArray <TH1D*> (3, nZPtBins, 2, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iData, iPhi, iCent
TH1D****** ZTracksPt          = Get5DArray <TH1D*> (3, nZPtBins, 2, numPhiBins, numCentBins);    // iSpc, iPtZ, iData, iPhi, iCent
TH1D******* ZTracksSubYields  = Get6DArray <TH1D*> (3, nZPtBins, 2, 2, numPhiBins, numCentBins); // iSpc, iPtZ, iData, iBkg, iPhi, iCent
TH1D******* ZTracksIAARatios  = Get6DArray <TH1D*> (3, nZPtBins, 2, 2, numPhiBins, numCentBins); // iSpc, iPtZ, iData, iBkg, iPhi, iCent
TH1D******* ZTracksICPRatios  = Get6DArray <TH1D*> (3, nZPtBins, 2, 2, numPhiBins, numCentBins); // iSpc, iPtZ, iData, iBkg, iPhi, iCent

TFile* histFile = nullptr;


TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2) {
  TPad* p;
  if (gDirectory->Get ("box_pad"))
    p = (TPad*)gDirectory->Get ("box_pad");
  else {
    p = new TPad ("box_pad", "box_pad", 0., 0., 1., 1.);
    p->SetFillStyle (0);
  } 
  p->Draw ();
  p->cd ();
  TBox* b = new TBox (x1, y1, x2, y2);
  return b;
}


void ResetXErrors (TGraphAsymmErrors* tg) {
  for (int ix = 0; ix < tg->GetN (); ix++) {
    tg->SetPointEXlow (ix, 0);
    tg->SetPointEXhigh (ix, 0);
  }
  return;
}


void RescaleWithError (TH1* h, double sf, double err) {
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    double y = h->GetBinContent (ix), yerr = h->GetBinError (ix);
    h->SetBinContent (ix, y / sf);
    h->SetBinError (ix, sqrt (pow (yerr / sf, 2) + pow (y * err / (sf * sf), 2)));
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void LoadHists () {
  if (histFile && histFile->IsOpen ())
    return;

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");
  histFile = new TFile ("savedHists.root", "read");
  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        ZTrackPtPhi[iData][iCent][iSpc] = (TH2D*)histFile->Get (Form ("ZTrackPtPhi_%s_%s_iCent%i", data, spc, iCent));
        ZPtSpecs[iData][iCent][iSpc] = (TH1D*)histFile->Get (Form ("ZPtSpec_%s_%s_iCent%i", data, spc, iCent));
        ZMYields[iData][iCent][iSpc] = (TH1D*)histFile->Get (Form ("ZMSpec_%s_%s_iCent%i", data, spc, iCent));
        TrackSpec[iData][iCent][iSpc] = (TH1D*)histFile->Get (Form ("TrackSpec_%s_%s_iCent%i", data, spc, iCent));
        //DRDists[iData][iCent][iSpc] = (TH2D*)histFile->Get (Form ("DRDist_%s_%s_iCent%i", data, spc, iCent));
      }
      ZPhiYields[iData][iCent] = (TH1D*)histFile->Get (Form ("ZPhiYield_%s_iCent%i", data, iCent));
      ElectronSpec[iData][iCent] = (TH1D*)histFile->Get (Form ("ElectronSpec_%s_iCent%i", data, iCent));
      MuonSpec[iData][iCent] = (TH1D*)histFile->Get (Form ("MuonSpec_%s_iCent%i", data, iCent));
      
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent] = (TH2D*)histFile->Get (Form ("ZMissingPt_%s_%s_iPtZ%i_iPhi%i_iCent%i", data, spc, iPtZ, iPhi, iCent));
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent] = (TH1D*)histFile->Get (Form ("ZTracksPt_%s_%s_iPtZ%i_iPhi%i_iCent%i", data, spc, iPtZ, iPhi, iCent));
          }
        }
      }
    }
    FCalSpec[iData] = (TH1D*)histFile->Get (Form ("FCalSpec_%s", data));
    FCalQ2Corr[iData] = (TH2D*)histFile->Get (Form ("FCalQ2Corr_%s", data));
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void SaveHists () {
  TFile* histFile = new TFile ("savedHists.root", "recreate");
  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        ZTrackPtPhi[iData][iCent][iSpc]->Write ();
        ZPtSpecs[iData][iCent][iSpc]->Write ();
        ZMYields[iData][iCent][iSpc]->Write ();
        TrackSpec[iData][iCent][iSpc]->Write ();
        //DRDists[iData][iCent][iSpc]->Write ();
      }
      ZPhiYields[iData][iCent]->Write ();
      ElectronSpec[iData][iCent]->Write ();
      MuonSpec[iData][iCent]->Write ();
      
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Write ();
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]->Write ();
          }
        }
      }
    }
    FCalSpec[iData]->Write ();
    FCalQ2Corr[iData]->Write ();
  }
  histFile->Close ();
  histFile = nullptr;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - xZTrk 3d distribution
////////////////////////////////////////////////////////////////////////////////////////////////
void Plot3DDist () {
  TCanvas* c = new TCanvas ("3dCanvas", "", 800, 600);
  c->cd ();

  c->SetLogy ();

  ZTrackPtPhi[0][numCentBins-1][2]->GetXaxis ()->SetTitle ("#phi_{Z} - #phi_{Trk}");
  ZTrackPtPhi[0][numCentBins-1][2]->GetYaxis ()->SetTitle ("#it{p}_{T}^{ ch} / #it{p}_{T}^{ Z}");

  //ZTrackPtPhi[0][numCentBins-1][2]->RebinY (2);
  ZTrackPtPhi[0][numCentBins-1][2]->Draw ("lego2");

  c->SaveAs (Form ("%s/ZTrackCorr.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotFCalDists () {
  TCanvas* c = new TCanvas ("FCalCanvas", "", 800, 600);
  c->cd ();

  c->SetLogy ();

  FCalSpec[0]->Scale (1./FCalSpec[0]->Integral ());
  FCalSpec[1]->Scale (1./FCalSpec[1]->Integral ());

  FCalSpec[0]->SetLineColor (kBlack);
  FCalSpec[1]->SetLineColor (kBlue);

  FCalSpec[0]->GetYaxis ()->SetRangeUser (1.5e-4, 1e-1);

  FCalSpec[0]->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
  FCalSpec[0]->GetYaxis ()->SetTitle ("A.U.");
  FCalSpec[1]->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
  FCalSpec[1]->GetYaxis ()->SetTitle ("A.U.");

  FCalSpec[0]->Draw ("hist");
  FCalSpec[1]->Draw ("hist same");

  myText (0.65, 0.88, kBlack, "2018 Pb+Pb", 0.04);
  myText (0.65, 0.81, kBlue, "Pythia8 + Hijing", 0.04);

  c->SaveAs (Form ("%s/FCalDist.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - pTTrk 2d projections
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotdPhiPtTrk () {
  TCanvas* c = new TCanvas ("dPhiPtTrkCanvas", "", 1600, 1200);
  c->cd ();
  c->Divide (2, 2);
  //c->SetLogy ();

  TF1***** Fits = Get4DArray <TF1*> (2, nPtTrkBins, numCentBins, 3); // iData, iPtTrk, iCent, iSpc
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      c->cd (iCent+1);

      gPad->SetTopMargin (0.01);
      gPad->SetBottomMargin (0.12);
      gPad->SetRightMargin (0.01);
      gPad->SetLeftMargin (0.12);

      for (short iData = 0; iData < 2; iData++) {
        const char* data = iData == 0 ? "data" : "mc";
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          ZTrackPhiPtBins[iData][iPtTrk][iCent][iSpc] = ZTrackPtPhi[iData][iCent][iSpc]->ProjectionX (Form ("ZTrackPhi_%s_iPtTrk%i_iCent%i", data, iPtTrk, iCent), iPtTrk+1, iPtTrk+1);
          TH1D* thisHist = ZTrackPhiPtBins[iData][iPtTrk][iCent][iSpc];
          thisHist->Rebin (2);
          if (iPtTrk > 3)
            thisHist->Rebin (2);
          if (iCent != 0)
            thisHist->Rebin (2);
          if (ZPtSpecs[iData][iCent][iSpc]->Integral () > 0)
            thisHist->Scale (1./ZPtSpecs[iData][iCent][iSpc]->Integral (), "width");

          TF1* constFit = new TF1 (Form ("fit_%s_iPtTrk%i_iCent%i", data, iPtTrk, iCent), "[0]+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x)", -pi/2, 3*pi/2);
          thisHist->Fit (constFit, "RN0");
          //const float min = constFit->GetMinimum (-pi/2, 3*pi/2);
          const double min = thisHist->Integral (thisHist->GetXaxis ()->FindBin (-pi/3), thisHist->GetXaxis ()->FindBin (pi/3)) / (thisHist->GetXaxis      ()->FindBin (pi/3) - thisHist->GetXaxis ()->FindBin (-pi/3) + 1);

          for (int ix = 1; ix <= thisHist->GetNbinsX (); ix++)
            thisHist->SetBinContent (ix, thisHist->GetBinContent (ix) - min); 

          constFit->SetParameter (0, constFit->GetParameter (0) - min);
          Fits[iData][iPtTrk][iCent][iSpc] = constFit;
        } // end loop over pT^trk bins
      } // end loop over data/MC
      
      double min = 1e30, max = 0;
      for (short iData = 0; iData < 2; iData++) {
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
          TH1D* thisHist = ZTrackPhiPtBins[iData][iPtTrk][iCent][iSpc];
          if (thisHist->GetMinimum () < min) min = thisHist->GetMinimum ();
          if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
        } // end loop over pT^trk bins
      } // end loop over data/MC

      for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
        TH1D* thisHist = (TH1D*)ZTrackPhiPtBins[1][iPtTrk][iCent][iSpc]->Clone ();

        thisHist->GetYaxis ()->SetRangeUser (1.3*min, 1.3*max);

        thisHist->SetLineColor (kBlack);
        thisHist->SetLineWidth (0);
        thisHist->SetFillColorAlpha (fillColors[iPtTrk], fillAlpha);
        thisHist->SetMarkerSize (0);

        thisHist->GetXaxis ()->SetTitle ("#Delta#phi");
        thisHist->GetYaxis ()->SetTitle ("1/N_{Z} dN_{ch}/d#Delta#phi (\"ZYAM\")");
        thisHist->GetXaxis ()->SetTitleOffset (0.6);
        thisHist->GetYaxis ()->SetTitleOffset (0.8);
        thisHist->GetXaxis ()->SetTitleSize (0.08);
        thisHist->GetYaxis ()->SetTitleSize (0.06);
        thisHist->GetXaxis ()->SetLabelSize (0.06);
        thisHist->GetYaxis ()->SetLabelSize (0.06);

        thisHist->DrawCopy (iPtTrk == 0 ? "b" : "same b");
        thisHist->SetLineWidth (1);
        thisHist->Draw ("hist same");
      }

      for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
        TGraphAsymmErrors* thisGraph = make_graph (ZTrackPhiPtBins[0][iPtTrk][iCent][iSpc]);

        thisGraph->GetYaxis ()->SetRangeUser (0, 1.3*max);

        thisGraph->SetLineColor (colors[iPtTrk]);
        thisGraph->SetMarkerColor (colors[iPtTrk]);

        thisGraph->Draw ("P");

        Fits[0][iPtTrk][iCent][iSpc]->SetLineColor (colors[iPtTrk]);
        Fits[0][iPtTrk][iCent][iSpc]->Draw ("same");

        if (iCent == 0) {
          myText (0.2, 0.9, kBlack, "#it{pp}", 0.06);
          const float pt_lo = ptTrkBins[iPtTrk];
          const float pt_hi = ptTrkBins[iPtTrk+1];
          myText (0.2, 0.80-0.075*iPtTrk, colors[iPtTrk], Form ("%.1f < #it{p}_{T}^{ ch} < %.1f GeV", pt_lo, pt_hi), 0.06);
        }
        else {
          myText (0.2, 0.9, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
        }
      } // end loop over pT^trk bins
    } // end loop over centrality

    c->SaveAs (Form ("%s/dPhi_pTtrk_%s.pdf", plotPath.Data (), spc));

  }
  Delete4DArray (Fits, 2, nPtTrkBins, numCentBins, 3);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot lepton Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotLeptonPtSpectra () {
  TCanvas* c = new TCanvas ("LeptonPtCanvas", "", 1000, 600);
  c->cd ();
  c->SetLogy ();
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TH1D* thisHist = (TH1D*)ElectronSpec[0][iCent]->Clone ();
    thisHist->Rebin (5);
    thisHist->Scale (1./ZPtSpecs[0][iCent][0]->Integral (), "width");

    thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ e} [GeV]");
    thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrowee} dN_{e}/d#it{p}_{T} [GeV^{-1}]");

    thisHist->SetLineColor (colors[iCent]);
    thisHist->SetMarkerColor (colors[iCent]);
    thisHist->SetMarkerStyle (kFullCircle);
    thisHist->SetMarkerSize (0.75);

    thisHist->Draw (iCent == 0 ? "e1" : "e1 same");
  }
  myText (0.76, 0.88, colors[0], "#it{pp}", 0.04);
  for (short iCent = 1; iCent < numCentBins; iCent++) {
    myText (0.76, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
  }
  c->SaveAs (Form ("%s/ElectronPtSpectra.pdf", plotPath.Data ()));

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TH1D* thisHist = (TH1D*)MuonSpec[0][iCent]->Clone ();
    thisHist->Scale (1./ZPtSpecs[0][iCent][1]->Integral (), "width");
    thisHist->Rebin (5);

    thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ #mu} [GeV]");
    thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrow#mu#mu} dN_{#mu}/d#it{p}_{T} [GeV^{-1}]");

    thisHist->SetLineColor (colors[iCent]);
    thisHist->SetMarkerColor (colors[iCent]);
    thisHist->SetMarkerStyle (kFullCircle);
    thisHist->SetMarkerSize (0.75);

    thisHist->Draw (iCent == 0 ? "e1" : "e1 same");
  }
  gPad->RedrawAxis ();

  myText (0.76, 0.88, colors[0], "#it{pp}", 0.04);
  for (short iCent = 1; iCent < numCentBins; iCent++) {
    myText (0.76, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
  }
  c->SaveAs (Form ("%s/MuonPtSpectra.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot track Pt spectra for each lepton species
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotLeptonTrackPtSpectra () {
  TCanvas* c = new TCanvas ("LeptonTrackPtCanvas", "", 800, 600);
  c->cd ();
  c->SetLogy ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TH1D* thisHist = (TH1D*)TrackSpec[0][iCent][iSpc]->Clone (Form("test_%i_%i", iSpc, iCent));
      thisHist->Rebin (5);
      thisHist->Scale (1./ZPtSpecs[0][iCent][iSpc]->Integral (ZPtSpecs[0][iCent][iSpc]->GetXaxis ()->FindBin (5), ZPtSpecs[0][iCent][iSpc]->GetNbinsX ()), "width");
      cout << iSpc << ", " << iCent << ", " << ZPtSpecs[0][iCent][iSpc]->Integral (ZPtSpecs[0][iCent][iSpc]->GetXaxis ()->FindBin (5), ZPtSpecs[0][iCent][iSpc]->GetNbinsX ()) << endl;

      thisHist->GetYaxis ()->SetRangeUser (6e-6, 450);

      if (iSpc == 0) {
        thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
        thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrowee} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }
      else {
        thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
        thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrow#mu#mu} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }

      thisHist->SetLineColor (colors[iCent]);
      thisHist->SetMarkerColor (colors[iCent]);
      thisHist->SetMarkerStyle (kFullCircle);
      thisHist->SetMarkerSize (0.75);

      thisHist->Draw (iCent == 0 ? "e1" : "e1 same");
    }
    gPad->RedrawAxis ();

    myText (0.66, 0.82, colors[0], "#it{pp}", 0.04);
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      myText (0.66, 0.82-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
    }

    myText (0.66, 0.88, kBlack, iSpc == 0 ? "Z#rightarrowee Events" : "Z#rightarrow#mu#mu Events", 0.04);
    c->SaveAs (Form ("%s/%sTrackPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : "Muon"));
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotZPtSpectra () {
  TCanvas* c = new TCanvas ("ZPtCanvas", "", 800, 600);
  c->cd ();
  c->SetLogy ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    //TH1D* thisHist = (TH1D*)ZPtSpecs[0][iCent][0]->Clone ();
    TH1D* thisHist = (TH1D*)ZPtSpecs[0][iCent][1]->Clone ();
    //thisHist->Add (ZPtSpecs[0][iCent][1]);
    //thisHist->Rebin (5);
    thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");

    TH1D* integral = new TH1D (Form ("ZPtSpec_integral_iCent%i", iCent), "", 300, 0, 300);
    //integral->Rebin (5);

    if (iCent == 0)
      cout << "Centrality bin: pp" << endl;
    else
      cout << "Centrality bin: " << Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]) << endl;

    for (int ix = 1; ix <= integral->GetNbinsX (); ix++) {
      double content = 0, varsq = 0;
      content = thisHist->IntegralAndError (ix, thisHist->GetNbinsX (), varsq);
      integral->SetBinContent (ix, content);
      integral->SetBinError (ix, sqrt (varsq));
    }

    cout << "\t#Z's >= 0 GeV: " << integral->GetBinContent (1) << endl;
    cout << "\t#Z's <  3 GeV: " << integral->GetBinContent (1) - integral->GetBinContent (4) << endl;
    cout << "\t#Z's >= 5 GeV: " << integral->GetBinContent (6) << endl;
    cout << "\t#Z's >= 10 GeV: " << integral->GetBinContent (11) << endl;

    thisHist->SetLineColor (colors[iCent]);
    thisHist->SetMarkerColor (colors[iCent]);
    thisHist->SetMarkerStyle (kFullCircle);
    thisHist->SetMarkerSize (0.5);

    integral->SetLineColor (colors[iCent]);
    integral->SetMarkerColor (colors[iCent]);
    integral->SetMarkerStyle (kOpenCircle);
    integral->SetMarkerSize (0.5);

    integral->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");

    if (iCent == 0)
      integral->Draw ("e1");
    else
      integral->Draw ("same e1");
    thisHist->Draw ("same e1");
  }
  gPad->RedrawAxis ();

  myMarkerText (0.75, 0.88, kBlack, kFullCircle, "n_{Z}(#it{p}_{T})", 1.25, 0.04);
  myMarkerText (0.75, 0.80, kBlack, kOpenCircle, "N_{Z}(#it{p}_{T}) =  #int_{#it{p}_{T}}^{#infty}n(x)dx", 1.25, 0.04);
  myText (0.54, 0.88, colors[0], "#it{pp}", 0.04);
  for (short iCent = 1; iCent < numCentBins; iCent++) {
    myText (0.54, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
  }

  c->SaveAs (Form ("%s/z_pt_spectrum.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z mass spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotZMassSpectra () {
  TCanvas* c = new TCanvas ("ZMassCanvas", "", 800, 600);
  c->cd ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
  //for (short iSpc = 0; iSpc < 1; iSpc++) {
    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
    //for (short iCent = 0; iCent < 2; iCent++) {
      TH1D* thisHist = ZMYields[1][iCent][iSpc];
      if (thisHist->GetMaximum () > 0)
        thisHist->Scale (1. / (thisHist->GetMaximum ()));

      thisHist->SetFillColorAlpha (fillColors[iCent], fillAlpha);
      thisHist->SetLineColor (kBlack);
      thisHist->SetMarkerSize (0);
      thisHist->SetLineWidth (0);
      thisHist->GetYaxis ()->SetRangeUser (0, 1.3);

      thisHist->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
      thisHist->GetYaxis ()->SetTitle ("A.U.");

      thisHist->DrawCopy (iCent == numCentBins-1 ? "bar" : "bar same");
      thisHist->SetLineWidth (1);
      thisHist->Draw ("hist same");
    }
    gPad->RedrawAxis ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TH1D* thisHist = ZMYields[0][iCent][iSpc];
      if (thisHist->GetMaximum () > 0)
        thisHist->Scale (1. / (thisHist->GetMaximum ()));

      TGraphAsymmErrors* thisGraph = make_graph (thisHist);
      ResetXErrors (thisGraph);
      deltaize (thisGraph, 0.1*(-1.5+iCent));

      const int markerStyle = kFullCircle;
      //const int markerStyle = iSpc == 0 ? kOpenCircle : kFullCircle;
      thisGraph->SetMarkerStyle (markerStyle);
      thisGraph->SetMarkerSize (1);
      thisGraph->SetLineWidth (1);
      thisGraph->SetLineColor (colors[iCent]);
      thisGraph->SetMarkerColor (colors[iCent]);
      thisGraph->GetYaxis ()->SetRangeUser (0, 1.3);

      thisGraph->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
      thisGraph->GetYaxis ()->SetTitle ("A.U.");
      thisGraph->Draw ("P");

      const char* spc = iSpc == 0 ? "Z#rightarrowee Events":"Z#rightarrow#mu#mu Events";
      if (iCent == 0) {
        myText (0.66, 0.88, kBlack, spc, 0.04);
        myMarkerText (0.25, 0.88, colors[0], markerStyle, Form ("#it{pp}"), 1.25, 0.04);
        //myMarkerText (0.2, 0.88-0.07*iSpc, colors[0], markerStyle, Form ("#it{pp} %s", spc.c_str()), 1.25, 0.04);
      }
      else
        myMarkerText (0.25, 0.88-0.07*iCent, colors[iCent], markerStyle, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 1.25, 0.04);
        //myMarkerText (0.2, 0.74-0.07*iSpc, colors[iCent], markerStyle, Form ("Pb+Pb %s", spc.c_str()), 1.25, 0.04);
    }
    c->SaveAs (Form ("%s/z%s_mass_spectrum.pdf", plotPath.Data (), iSpc == 0 ? "ee":"mumu"));
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotZPhiYield () {
  TCanvas* c = new TCanvas ("ZPhiYieldCanvas", "", 800, 600);
  c->cd ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TH1D* thisHist = ZPhiYields[0][iCent];
    thisHist->Rebin (8);
    thisHist->Scale (1. / thisHist->Integral (), "width");

    TF1* fit = new TF1 ("fit", "[0]+[1]*cos(x)+[2]*cos(2*x)", -pi/2, 3*pi/2);
    thisHist->Fit (fit, "RN0");
    delete fit;
    fit = nullptr;

    TGraphAsymmErrors* thisGraph = make_graph (thisHist);
    deltaize (thisGraph, (1.5-iCent)*0.02, false);

    thisGraph->GetXaxis ()->SetTitle ("2|#phi_{Z} - #Psi_{2}|");
    thisGraph->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/d#Delta#phi");

    thisGraph->GetYaxis ()->SetRangeUser (0.1, 0.7);

    thisGraph->SetLineColor (colors[iCent]);
    thisGraph->SetMarkerColor (colors[iCent]);
    thisGraph->SetMarkerSize (0.75);
    if (iCent == 0)
      thisGraph->Draw ("ap");
    else
      thisGraph->Draw ("p");

  } // end loop over cents

  myText (0.66, 0.88, colors[0], "#it{pp}", 0.04);
  for (short iCent = 1; iCent < numCentBins; iCent++) {
    myText (0.66, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
  }

  c->SaveAs (Form ("%s/ZPhiYields.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotZMissingPt () {
  TCanvas* c = new TCanvas ("ZMissingPtCanvas", "", 600*numPhiTrkBins, 500);
  c->cd ();

  c->Divide (2, 1);

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        c->cd (iPhi+1);
        gPad->SetLogx ();
        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TH1D* thisHist = new TH1D (Form ("ZMissingPtAvg_%s_iPtZ%i_iPhi%i_iCent%i", spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          TH1D* integral = new TH1D (Form ("ZMissingPtInt_%s_iPtZ%i_iPhi%i_iCent%i", spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          thisHist->Sumw2 ();
          integral->Sumw2 ();
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            TH1D* px = ZMissingPt[iSpc][iPtZ][0][iPhi][iCent]->ProjectionX ("_px", iPtTrk+1, iPtTrk+1);
            TF1* fit = new TF1 ("fit", "gaus(0)", zMissingPtBins[0], zMissingPtBins[numZMissingPtBins]);
            px->Fit (fit, "RN0Q");
            thisHist->SetBinContent (iPtTrk+1, fit->GetParameter (1));
            thisHist->SetBinError (iPtTrk+1, fit->GetParError (1));
            if (px) delete px;

            px = ZMissingPt[iSpc][iPtZ][0][iPhi][iCent]->ProjectionX ("_px", 0, iPtTrk+1);
            px->Fit (fit, "RN0Q");
            integral->SetBinContent (iPtTrk+1, fit->GetParameter (1));
            integral->SetBinError (iPtTrk+1, fit->GetParError (1));
            if (px) delete px;

            if (fit) delete fit;
          }
          ZMissingPtAvgs[iSpc][iPtZ][0][iPhi][iCent] = thisHist;
          ZMissingPtInts[iSpc][iPtZ][0][iPhi][iCent] = integral;

          integral->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch}  [GeV]");
          integral->GetYaxis ()->SetTitle ("<#it{p}_{T}^{ ||}>  [GeV]");
          integral->GetXaxis ()->SetTitleSize (0.06);
          integral->GetYaxis ()->SetTitleSize (0.06);
          integral->GetXaxis ()->SetLabelSize (0.06);
          integral->GetYaxis ()->SetLabelSize (0.06);
          integral->GetXaxis ()->SetTitleOffset (1.2);
          integral->GetYaxis ()->SetTitleOffset (1.2);

          integral->GetYaxis ()->SetRangeUser (-32, 16);

          integral->SetFillColorAlpha (fillColors[iCent], fillAlpha);
          integral->SetLineColor (kBlack);
          integral->SetMarkerSize (0);
          integral->SetLineWidth (0);
          integral->DrawCopy (iCent == numCentBins-1 ? "b" : "b same");
          integral->SetLineWidth (1);
          integral->Draw ("hist same");
        }

        TLine* l = new TLine (ptTrkBins[0], 0, ptTrkBins[nPtTrkBins], 0);
        l->Draw ("same");

        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TGraphAsymmErrors* thisGraph = make_graph (ZMissingPtAvgs[iSpc][iPtZ][0][iPhi][iCent]);
          RecenterGraph (thisGraph);
          ResetXErrors (thisGraph);
          deltaize (thisGraph, (1.5-iCent)*0.04, false);

          thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch}  [GeV]");
          thisGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{ ||}>  [GeV]");

          thisGraph->SetLineColor (colors[iCent]);
          thisGraph->SetMarkerColor (colors[iCent]);
          thisGraph->SetMarkerSize (0.75);
          thisGraph->Draw ("P");
        } // end loop over cents
        if (iPhi == 0) {
          myText (0.25, 0.88, kBlack, "0 < #Delta#phi < #pi/8 or 7#pi/8 < #Delta#phi < #pi", 0.05);
          myText (0.25, 0.81, kBlack, Form ("%g < #it{p}_{T}^{ Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
        }
        else if (iPhi == 1)
          myText (0.25, 0.88, kBlack, "0 < #Delta#phi < #pi", 0.05);
      } // end loop over directions

      myText (0.66, 0.88, colors[0], "#it{pp}", 0.05);
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        myText (0.66, 0.88-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
      }

      c->SaveAs (Form ("%s/ZMissingPt.pdf", plotPath.Data ()));
    } // end loop over pT^Z bins
  } // end loop over species
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTtrk binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotTrkYield () {
  TCanvas* c = new TCanvas ("TrkYieldCanvas", "", 860*numCentBins, 1000);
  c->cd ();

  const double padRatio = 1.1; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;

  TPad** topPads = Get1DArray <TPad*> (numCentBins);
  TPad** bottomPads = Get1DArray <TPad*> (numCentBins);
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    TPad* topPad = new TPad ("topPad", "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
    TPad* bottomPad = new TPad ("bottomPad", "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);
    topPads[iCent] = topPad;
    bottomPads[iCent] = bottomPad;

    topPad->SetTopMargin (0.04);
    topPad->SetBottomMargin (0);
    topPad->SetLeftMargin (0.17);
    topPad->SetRightMargin (0.06);
    bottomPad->SetTopMargin (0);
    bottomPad->SetBottomMargin (0.20);
    bottomPad->SetLeftMargin (0.17);
    bottomPad->SetRightMargin (0.06);
    topPad->Draw ();
    bottomPad->Draw ();
  }

  double***** yieldNormFactor      = Get5DArray <double> (3, nZPtBins, 2, numCentBins, numPhiBins);
  double***** yieldNormFactorError = Get5DArray <double> (3, nZPtBins, 2, numCentBins, numPhiBins);
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iData = 0; iData < 2; iData++) {
        for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            double f = 0, e = 0;
            TH1D* h = ZPtSpecs[iData][iCent][iSpc];
            if (iPhi < numPhiBins-1)
              f = h->IntegralAndError (h->FindBin (zPtBins[iPtZ]), h->FindBin (zPtBins[iPtZ+1])-1, e);
            else
              f = h->IntegralAndError (1, 3, e);
            yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi]      = f * (phiHighBins[iPhi] - phiLowBins[iPhi]);
            yieldNormFactorError[iSpc][iPtZ][iData][iCent][iPhi] = e * (phiHighBins[iPhi] - phiLowBins[iPhi]);
          } // end loop over cents
        } // end loop over phi
      } // end loop over data types
    } // end loop over pT^Z bins
  } // end loop over species

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        TPad* topPad = topPads[iCent];
        TPad* bottomPad = bottomPads[iCent];

        for (short iData = 0; iData < 2; iData++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {

            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];

            //RescaleWithError (thisHist, yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi], yieldNormFactorError[iSpc][iPtZ][iData][iCent][iPhi]);
            if (yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi] > 0)
              thisHist->Scale (1. / yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi]);
          } // end loop over phi
        } // end loop over data types
        

        double min = 1e30;
        double max = 0;
        for (short iData = 0; iData < 2; iData++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];
            if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
            if (thisHist->GetMaximum () > max)  max = thisHist->GetMaximum ();
          } // end loop over phi
        } // end loop over data types

        topPad->cd ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        for (int iPhi = numPhiBins-2; iPhi >= 0; iPhi--) {
          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][1][iPhi][iCent];

          thisHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
          thisHist->SetMarkerSize (0);
          thisHist->SetLineColor (kBlack);
          thisHist->SetLineWidth (0);
          thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          thisHist->GetYaxis ()->SetTitle ("Per-trigger yield Y (#it{p}_{T})");
          thisHist->GetXaxis ()->SetTitleSize (0.03 / uPadY);
          thisHist->GetYaxis ()->SetTitleSize (0.03 / uPadY);
          thisHist->GetXaxis ()->SetLabelSize (0.03 / uPadY);
          thisHist->GetYaxis ()->SetLabelSize (0.03 / uPadY);
          thisHist->GetYaxis ()->SetTitleOffset (2.6 * uPadY);

          thisHist->GetYaxis ()->SetRangeUser (0.5*min, 2*max);
          if (iPhi == numPhiBins-2)
            thisHist->DrawCopy ("bar");
          else
            thisHist->DrawCopy ("bar same");
          thisHist->SetLineWidth (1);
          thisHist->Draw ("hist same");
        }
        gPad->RedrawAxis ();

        for (int iPhi = 0; iPhi < numPhiBins-1; iPhi++) {
          const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][0][iPhi][iCent];
          
          TGraphAsymmErrors* thisGraph = make_graph (thisHist);
          RecenterGraph (thisGraph);
          ResetXErrors (thisGraph);

          thisGraph->SetMarkerStyle (markerStyle);
          thisGraph->SetMarkerColor (colors[iPhi]);
          thisGraph->SetLineColor (colors[iPhi]);
          thisGraph->SetMarkerSize (1);
          thisGraph->SetLineWidth (2);
  
          thisGraph->GetYaxis ()->SetRangeUser (0.5*min, 2*max);
          thisGraph->Draw ("P");

          if (iCent == 0)
            myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.034/uPadY);
          else if (iCent == numCentBins-1) {
            if (iPhi == 0) {
              myText (0.43, 0.91, kBlack, "Data", 0.03/uPadY);
              myText (0.574, 0.91, kBlack, "MC", 0.03/uPadY);
              myText (0.683, 0.91, kBlack, "#Delta#phi", 0.03/uPadY);
            }

            TVirtualPad* cPad = gPad; // store current pad
            const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
            const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
            TBox* b = TBoxNDC (0.61-0.024, 0.85-0.06*iPhi-0.016, 0.61+0.024, 0.85-0.06*iPhi+0.016);
            b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            b->Draw ("l");
            myMarkerText (0.512, 0.852-0.06*iPhi, colors[iPhi], markerStyle, "", 1.4, 0.03/uPadY);
            myText (0.68, 0.85-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.03/uPadY);
            cPad->cd ();
          }

        } // end loop over phi

        bottomPad->cd ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        for (short iBkg = 0; iBkg < 2; iBkg++) {
          for (short iData = 0; iData < 2; iData++) {
            const char* data = iData == 0 ? "data" : "mc";
            for (int iPhi = numPhiBins-2; iPhi >= 0; iPhi--) {
              TH1D* subYield = new TH1D (Form ("subYieldPt_%s_%s_iPtZ%i_iBkg%i_iPhi%i_iCent%i", spc, data, iPtZ, iBkg, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
              subYield->Sumw2 ();
              //subYield->Rebin (8);

              subYield->Add (ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]);
              if (iBkg == 0)
                subYield->Add (ZTracksPt[iSpc][iPtZ][iData][0][iCent], -1);
              else
                subYield->Add (ZTracksPt[iSpc][iPtZ][iData][numPhiBins-1][iCent], -1);

              ZTracksSubYields[iSpc][iPtZ][iData][iBkg][iPhi][iCent] = subYield;
            } // end loop over phi
          } // end loop over data types
        } // end loop over bkgs
        
        min = 1e30;
        max = 0;
        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (short iData = 0; iData < 2; iData++) {
            for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
              TH1D* thisHist = ZTracksSubYields[iSpc][iPtZ][iData][iBkg][iPhi][iCent];
              if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
              if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
            } // end loop over phi
          } // end loop over data types
        } // end loop over bkgs

        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = numPhiBins-2; iPhi >= 1; iPhi--) {
            TH1D* subYield = ZTracksSubYields[iSpc][iPtZ][1][iBkg][iPhi][iCent];

            float delta = log10 (max) - log10 (min);
            subYield->GetYaxis ()->SetRangeUser (pow (10, log10 (min) - 0.1*delta), pow (10, log10 (max) + 0.1*delta));

            subYield->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            subYield->SetLineColor (kBlack);
            subYield->SetMarkerSize (0);
            subYield->SetLineWidth (0);
            subYield->SetMarkerStyle (kFullCircle);
            subYield->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            subYield->GetYaxis ()->SetTitle ("Signal Yield");
            subYield->GetXaxis ()->SetTitleSize (0.03 / uPadY);
            subYield->GetYaxis ()->SetTitleSize (0.03 / uPadY);
            subYield->GetXaxis ()->SetLabelSize (0.03 / uPadY);
            subYield->GetYaxis ()->SetLabelSize (0.03 / uPadY);
            subYield->GetXaxis ()->SetTitleOffset (1.5);
            subYield->GetYaxis ()->SetTitleOffset (2.8 * uPadY);

            if (iBkg == 0 && iPhi == numPhiBins-2)
              subYield->DrawCopy ("bar");
            else
              subYield->DrawCopy ("bar same");
            subYield->SetLineWidth (1);
            subYield->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();

          for (int iPhi = numPhiBins-2; iPhi >= 1; iPhi--) {
            const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
            TGraphAsymmErrors* subYield = make_graph (ZTracksSubYields[iSpc][iPtZ][0][iBkg][iPhi][iCent]);
            RecenterGraph (subYield);
            ResetXErrors (subYield);

            float delta = log10 (max) - log10 (min);
            subYield->GetYaxis ()->SetRangeUser (pow (10, log10 (min) - 0.1*delta), pow (10, log10 (max) + 0.1*delta));
            subYield->SetMarkerStyle (markerStyle);
            subYield->SetLineColor (colors[iPhi]);
            subYield->SetMarkerColor (colors[iPhi]);
            subYield->SetMarkerSize (1);
            subYield->SetLineWidth (2);
            subYield->SetMarkerStyle (kFullCircle);

            subYield->Draw ("P");
          } // end loop over phi
        } // end loop over bkgs
      } // end loop over cents
      

      topPads[0]->cd ();
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        topPads[iCent]->cd ();
        if (iCent == 0)
          myText (0.22, 0.06, kBlack, "#it{pp}", 0.03 / uPadY);// or d#it{p}_{T}^{1,2} > 5 GeV", 0.03 / uPadY);
        else
          myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.03 / uPadY);
      }

      c->SaveAs (Form ("%s/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));

    } // end loop over pT^Z bins
  } // end loop over species
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotIAARatios () {
  TCanvas* c = new TCanvas ("IAACanvas", "", 800*(numCentBins-1), 500);
  c->cd ();

  c->Divide (numCentBins-1, 1);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iData = 0; iData < 2; iData++) {
        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TH1D* ppHist = ZTracksSubYields[iSpc][iPtZ][0][iBkg][iPhi][0];

            for (short iCent = 1; iCent < numCentBins; iCent++) {
              if (!ZTracksIAARatios[iSpc][iPtZ][iData][iBkg][iPhi][iCent]) {
                TH1D* PbPbHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iData][iBkg][iPhi][iCent]->Clone ());
                PbPbHist->Divide (ppHist);
                ZTracksIAARatios[iSpc][iPtZ][iData][iBkg][iPhi][iCent] = PbPbHist;
              } 
            } // end loop over cents
          } // end loop over phi
        } // end loop over bkgs
      } // end loop over data types
    } // end loop over pT^Z bins
  } // end loop over species

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            c->cd (iCent);
            gPad->SetLogx ();
            
            TH1D* ratioHist = ZTracksIAARatios[iSpc][iPtZ][1][iBkg][iPhi][iCent];

            //ratioHist->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
            ratioHist->GetYaxis ()->SetRangeUser (0, 3);
            ratioHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            ratioHist->SetMarkerSize (0);
            ratioHist->SetLineColor (kBlack);
            ratioHist->SetLineWidth (0);
            ratioHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            ratioHist->GetYaxis ()->SetTitle ("I_{AA}");
            ratioHist->GetXaxis ()->SetTitleSize (0.06);
            ratioHist->GetYaxis ()->SetTitleSize (0.06);
            ratioHist->GetXaxis ()->SetLabelSize (0.06);
            ratioHist->GetYaxis ()->SetLabelSize (0.06);
            ratioHist->GetXaxis ()->SetTitleOffset (1.2);
            ratioHist->GetYaxis ()->SetTitleOffset (1.2);

            if (iBkg == 0 && iPhi == 1)
              ratioHist->DrawCopy ("bar");
            else
              ratioHist->DrawCopy ("bar same");
            ratioHist->SetLineWidth (1);
            ratioHist->Draw ("hist same");
          } // end loop over phi
        } // end loop over bkgs
        gPad->RedrawAxis ();

        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            c->cd (iCent);
            gPad->SetLogx ();
            //gPad->SetLogy ();
            
            TGraphAsymmErrors* ratioGraph = make_graph (ZTracksIAARatios[iSpc][iPtZ][0][iBkg][iPhi][iCent]);
            RecenterGraph (ratioGraph);
            ResetXErrors (ratioGraph);
            deltaize (ratioGraph, 1+(0.5*(numPhiBins-2)+0.5-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-2)+0.5

            //ratioGraph->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
            ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
            ratioGraph->SetLineColor (colors[iPhi]);
            ratioGraph->SetMarkerColor (colors[iPhi]);
            ratioGraph->SetMarkerSize (1.2);
            ratioGraph->SetLineWidth (2);
            ratioGraph->SetMarkerStyle (kFullCircle);

            ratioGraph->Draw ("P");

            if (iCent == 1)
              myText (0.49, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
            else if (iCent == numCentBins-1) {
              if (iPhi == 1) {
                myText (0.44, 0.91, kBlack, "Data", 0.05);
                myText (0.577, 0.91, kBlack, "MC", 0.05);
                myText (0.685, 0.91, kBlack, "#Delta#phi", 0.05);
              }
              TVirtualPad* cPad = gPad; // store current pad
              const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
              const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
              TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
              b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
              b->Draw ("l");
              myMarkerText (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
              myText (0.68, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
              //myText (0.56, 0.90-0.07*(iPhi-1), colors[iPhi], Form ("(%s, %s)", lo, hi), 0.06);
              cPad->cd ();
            }
          } // end loop over phi
        } // end loop over bkgs
      } // end loop over cents
      c->cd (1);

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
      } // end loop over cents

      c->SaveAs (Form ("%s/iaa_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species

}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotICPRatios () {
  TCanvas* c = new TCanvas ("ICPCanvas", "", 600*(numCentBins-2), 500);
  c->cd ();

  c->Divide (numCentBins-2, 1);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iData = 0; iData < 2; iData++) {
        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TH1D* periphHist = ZTracksSubYields[iSpc][iPtZ][iData][iBkg][iPhi][1];

            for (short iCent = 2; iCent < numCentBins; iCent++) {
              if (!ZTracksICPRatios[iSpc][iPtZ][iData][iBkg][iPhi][iCent]) {
                TH1D* centHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iData][iBkg][iPhi][iCent]->Clone ());
                centHist->Divide (periphHist);
                ZTracksICPRatios[iSpc][iPtZ][iData][iBkg][iPhi][iCent] = centHist;
              } 
            } // end loop over cents
          } // end loop over phi
        } // end loop over bkgs
      } // end loop over data types
    } // end loop over pT^Z bins
  } // end loop over species

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            c->cd (iCent-1);
            gPad->SetLogx ();
            
            TH1D* ratioHist = ZTracksICPRatios[iSpc][iPtZ][1][iBkg][iPhi][iCent];

            //ratioHist->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
            ratioHist->GetYaxis ()->SetRangeUser (0, 3);
            ratioHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            ratioHist->SetMarkerSize (0);
            ratioHist->SetLineColor (kBlack);
            ratioHist->SetLineWidth (0);
            ratioHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
            ratioHist->GetYaxis ()->SetTitle ("I_{CP}");
            ratioHist->GetXaxis ()->SetTitleSize (0.06);
            ratioHist->GetYaxis ()->SetTitleSize (0.06);
            ratioHist->GetXaxis ()->SetLabelSize (0.06);
            ratioHist->GetYaxis ()->SetLabelSize (0.06);
            ratioHist->GetXaxis ()->SetTitleOffset (1.2);
            ratioHist->GetYaxis ()->SetTitleOffset (1.2);

            if (iBkg == 0 && iPhi == 1)
              ratioHist->DrawCopy ("bar");
            else
              ratioHist->DrawCopy ("bar same");
            ratioHist->SetLineWidth (1);
            ratioHist->Draw ("hist same");
          } // end loop over phi
        } // end loop over bkgs
        gPad->RedrawAxis ();

        //for (short iBkg = 0; iBkg < 2; iBkg++) {
        for (short iBkg = 0; iBkg < 1; iBkg++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            c->cd (iCent-1);
            
            TGraphAsymmErrors* ratioGraph = make_graph (ZTracksICPRatios[iSpc][iPtZ][0][iBkg][iPhi][iCent]);
            RecenterGraph (ratioGraph);
            ResetXErrors (ratioGraph);
            deltaize (ratioGraph, 1+(0.5*(numPhiBins-2)+0.5-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-2)+0.5

            //ratioGraph->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
            ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
            ratioGraph->SetLineColor (colors[iPhi]);
            ratioGraph->SetMarkerColor (colors[iPhi]);
            ratioGraph->SetMarkerSize (1.2);
            ratioGraph->SetLineWidth (2);
            ratioGraph->SetMarkerStyle (kFullCircle);

            ratioGraph->Draw ("P");

            if (iCent == 2)
              myText (0.59, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
            else if (iCent == numCentBins-1) {
              if (iPhi == 1) {
                myText (0.555, 0.91, kBlack, "Data", 0.05);
                myText (0.682, 0.91, kBlack, "MC", 0.05);
                myText (0.785, 0.91, kBlack, "#Delta#phi", 0.05);
              }
              TVirtualPad* cPad = gPad; // store current pad
              const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
              const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
              TBox* b = TBoxNDC (0.71-0.024, 0.91-0.06*iPhi-0.016, 0.71+0.024, 0.91-0.06*iPhi+0.016);
              b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
              b->Draw ("l");
              myMarkerText (0.612, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
              myText (0.78, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
              //myText (0.56, 0.90-0.07*(iPhi-1), colors[iPhi], Form ("(%s, %s)", lo, hi), 0.06);
              cPad->cd ();
            }
          } // end loop over phi
        } // end loop over bkgs
      } // end loop over cents

      c->cd (1);
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);
      } // end loop over cents

      c->SaveAs (Form ("%s/icp_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void ZTrackHistMaker () {
  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TFile* fcalWeightsFile = new TFile ("FCalWeightsFile.root", "read");
  TH1D* fcalWeights = (TH1D*)fcalWeightsFile->Get ("FCalWeights");

  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    FCalSpec[iData] = new TH1D (Form ("FCalSpec_%s", data), "", 300, 0, 6000); 
    FCalSpec[iData]->Sumw2 ();
    FCalQ2Corr[iData] = new TH2D (Form ("FCalQ2Corr_%s", data), "", 300, 0, 6000, 150, 0, 300);
    FCalQ2Corr[iData]->Sumw2 ();
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      ZPhiYields[iData][iCent] = new TH1D (Form ("ZPhiYield_%s_iCent%i", data, iCent), "", 80, 0, pi);
      ZPhiYields[iData][iCent]->Sumw2 ();
      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : "mumu");
        ZTrackPtPhi[iData][iCent][iSpc] = new TH2D (Form ("ZTrackPtPhi_%s_%s_iCent%i", data, spc, iCent), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
        ZPtSpecs[iData][iCent][iSpc] = new TH1D (Form ("ZPtSpec_%s_%s_iCent%i", data, spc, iCent), "", 300, 0, 300);
        ZMYields[iData][iCent][iSpc] = new TH1D (Form ("ZMSpec_%s_%s_iCent%i", data, spc, iCent), "", 40, 76, 106);
        TrackSpec[iData][iCent][iSpc] = new TH1D (Form ("TrackSpec_%s_%s_iCent%i", data, spc, iCent), "", 100, 0, 100);
        //DRDists[iData][iCent][iSpc] = new TH2D (Form ("DRDist_%s_%s_iCent%i", data, spc, iCent), "", 100, 0, 0.1, 80, 0, 0.8);
        
        ZTrackPtPhi[iData][iCent][iSpc]->Sumw2 ();
        ZPtSpecs[iData][iCent][iSpc]->Sumw2 ();
        ZMYields[iData][iCent][iSpc]->Sumw2 ();
        TrackSpec[iData][iCent][iSpc]->Sumw2 ();
        //DRDists[iData][iCent][iSpc]->Sumw2 ();
      }
      ElectronSpec[iData][iCent] = new TH1D (Form ("ElectronSpec_%s_iCent%i", data, iCent), "", 250, 0, 250);
      MuonSpec[iData][iCent] = new TH1D (Form ("MuonSpec_%s_iCent%i", data, iCent), "", 250, 0, 250);

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent] = new TH2D (Form ("ZMissingPt_%s_%s_iPtZ%i_iPhi%i_iCent%i", data, spc, iPtZ, iPhi, iCent), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent] = new TH1D (Form ("ZTracksPt_%s_%s_iPtZ%i_iPhi%i_iCent%i", data, spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          }
        }
      }
    }
  }

  bool isEE = false, isMC = false;
  float fcal_et = 0, q2 = 0, psi2 = 0, z_pt = 0, z_eta = 0, z_phi = 0, z_m = 0, event_weight = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;//, *l_trk_pt = nullptr, *l_trk_eta = nullptr, *l_trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("isMC",         &isMC);
  PbPbTree->SetBranchAddress ("event_weight", &event_weight);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_eta",        &z_eta);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("z_m",          &z_m);
  PbPbTree->SetBranchAddress ("l1_pt",        &l1_pt);
  PbPbTree->SetBranchAddress ("l1_eta",       &l1_eta);
  PbPbTree->SetBranchAddress ("l1_phi",       &l1_phi);
  PbPbTree->SetBranchAddress ("l1_charge",    &l1_charge);
  PbPbTree->SetBranchAddress ("l2_pt",        &l2_pt);
  PbPbTree->SetBranchAddress ("l2_eta",       &l2_eta);
  PbPbTree->SetBranchAddress ("l2_phi",       &l2_phi);
  PbPbTree->SetBranchAddress ("l2_charge",    &l2_charge);
  PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
  PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
  PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);
  //PbPbTree->SetBranchAddress ("l_trk_pt",     &l_trk_pt);
  //PbPbTree->SetBranchAddress ("l_trk_eta",    &l_trk_eta);
  //PbPbTree->SetBranchAddress ("l_trk_phi",    &l_trk_phi);

  int nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    PbPbTree->GetEntry (iEvt);

    const short iData = isMC ? 1 : 0; // 0 for not MC (data), 1 for MC
    const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

    int iCent = 0;
    while (iCent < numCentBins) {
      if (fcal_et < centBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numCentBins-1)
      continue;

    int iPtZ = 0; // find z-pt bin
    while (iPtZ < nZPtBins) {
      if (z_pt < zPtBins[iPtZ+1])
        break;
      else
        iPtZ++;
    }

    FCalSpec[iData]->Fill (fcal_et, event_weight);
    FCalQ2Corr[iData]->Fill (fcal_et, q2, event_weight);
    if (isMC)
      event_weight *= fcalWeights->GetBinContent (fcalWeights->FindBin (fcal_et));

    ZPtSpecs[iData][iCent][iSpc]->Fill (z_pt, event_weight);
    ZMYields[iData][iCent][iSpc]->Fill (z_m, event_weight);
    if (isEE) {
      ElectronSpec[iData][iCent]->Fill (l1_pt, event_weight);
      ElectronSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }
    else {
      MuonSpec[iData][iCent]->Fill (l1_pt, event_weight);
      MuonSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    ZPhiYields[iData][iCent]->Fill (2*dphi, event_weight);

    for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
      for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }

    for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
      const float trkpt = trk_pt->at (iTrk);

      if (trkpt < 1)
        continue;

      TrackSpec[iData][iCent][iSpc]->Fill (trkpt, event_weight);
      //float minDR = 2;
      //int minLTrk = -1;
      //for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
      //  float dR = DeltaR (l_trk_eta->at (iLTrk), trk_eta->at (iTrk), l_trk_phi->at (iLTrk), trk_phi->at (iTrk));
      //  if (dR < minDR) {
      //    minDR = dR;
      //    minLTrk = iLTrk;
      //  }
      //}
      //if (isEE && minLTrk != -1) 
      //  DRDists[iData][iCent][0]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);
      //else if (minLTrk != -1) 
      //  DRDists[iData][iCent][1]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);

      // Add to missing pT (requires dphi in +/-pi/2 to +/-pi)
      dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
      bool awaySide = false;
      if (dphi > pi/2) {
        dphi = pi-dphi;
        awaySide = true;
      }

      short iPtTrk = 0;
      while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
        iPtTrk++;
      // start at the 1st phi bin and integrate outwards until the track is no longer contained 
      // e.g. so 7pi/8->pi is a subset of pi/2->pi
      short iPhi = 0;
      while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
        if (awaySide)
          trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi);
        else
          trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi);
        iPhi++;
      }

      // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
      dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
      short idPhi = 0;
      if (z_pt < 3)
        idPhi = numPhiBins-1;
      else if (5 < z_pt) {
        while (idPhi < numPhiBins-1) {
          if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
            break;
          else
            idPhi++;
        }
      }

      if (idPhi >= 0 && idPhi < numPhiBins)
        ZTracksPt[iSpc][iPtZ][iData][idPhi][iCent]->Fill (trkpt, event_weight);

      if (z_pt > 20) {
        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        ZTrackPtPhi[iData][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
      }
    } // end loop over tracks

    for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }
  } // end loop over Pb+Pb tree


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ppTree->SetBranchAddress ("isEE",         &isEE);
  ppTree->SetBranchAddress ("isMC",         &isMC);
  ppTree->SetBranchAddress ("event_weight", &event_weight);
  ppTree->SetBranchAddress ("psi2",         &psi2);
  ppTree->SetBranchAddress ("z_pt",         &z_pt);
  ppTree->SetBranchAddress ("z_eta",        &z_eta);
  ppTree->SetBranchAddress ("z_phi",        &z_phi);
  ppTree->SetBranchAddress ("z_m",          &z_m);
  ppTree->SetBranchAddress ("l1_pt",        &l1_pt);
  ppTree->SetBranchAddress ("l1_eta",       &l1_eta);
  ppTree->SetBranchAddress ("l1_phi",       &l1_phi);
  ppTree->SetBranchAddress ("l1_charge",    &l1_charge);
  ppTree->SetBranchAddress ("l2_pt",        &l2_pt);
  ppTree->SetBranchAddress ("l2_eta",       &l2_eta);
  ppTree->SetBranchAddress ("l2_phi",       &l2_phi);
  ppTree->SetBranchAddress ("l2_charge",    &l2_charge);
  ppTree->SetBranchAddress ("trk_pt",       &trk_pt);
  ppTree->SetBranchAddress ("trk_eta",      &trk_eta);
  ppTree->SetBranchAddress ("trk_phi",      &trk_phi);
  //ppTree->SetBranchAddress ("l_trk_pt",     &l_trk_pt);
  //ppTree->SetBranchAddress ("l_trk_eta",    &l_trk_eta);
  //ppTree->SetBranchAddress ("l_trk_phi",    &l_trk_phi);

  nEvts = ppTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    ppTree->GetEntry (iEvt);

    const short iData = isMC ? 1 : 0; // 0 for not MC (data), 1 for MC
    const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
    const short iCent = 0; // iCent = 0 for pp

    int iPtZ = 0; // find z-pt bin
    while (iPtZ < nZPtBins) {
      if (z_pt < zPtBins[iPtZ+1])
        break;
      else
        iPtZ++;
    }

    ZPtSpecs[iData][iCent][iSpc]->Fill (z_pt, event_weight);
    ZMYields[iData][iCent][iSpc]->Fill (z_m, event_weight);
    if (isEE) {
      ElectronSpec[iData][iCent]->Fill (l1_pt, event_weight);
      ElectronSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }
    else {
      MuonSpec[iData][iCent]->Fill (l1_pt, event_weight);
      MuonSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    ZPhiYields[iData][iCent]->Fill (2*dphi, event_weight);

    for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
      const float trkpt = trk_pt->at (iTrk);

      if (trkpt < 1)
        continue;

      TrackSpec[iData][iCent][iSpc]->Fill (trkpt, event_weight);
      //float minDR = 2;
      //int minLTrk = -1;
      //if (!isMC) {
      //  for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
      //    float dR = DeltaR (l_trk_eta->at (iLTrk), trk_eta->at (iTrk), l_trk_phi->at (iLTrk), trk_phi->at (iTrk));
      //    if (dR < minDR) {
      //      minDR = dR;
      //      minLTrk = iLTrk;
      //    }
      //  }
      //}
      //if (isEE && minLTrk != -1)
      //  DRDists[iData][iCent][0]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);
      //else if (minLTrk != -1)
      //  DRDists[iData][iCent][1]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);

      // Add to missing pT (requires dphi in -pi/2 to pi/2)
      dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
      bool awaySide = false;
      if (dphi > pi/2) {
        dphi = pi-dphi;
        awaySide = true;
      }

      short iPtTrk = 0;
      while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
        iPtTrk++;
      // start at the 1st phi bin and integrate outwards until the track is no longer contained 
      // e.g. so 7pi/8->pi is a subset of pi/2->pi
      short iPhi = 0;
      while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
        if (awaySide)
          trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi);
        else
          trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi);
        iPhi++;
      }
      
      // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
      dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
      short idPhi = 0;
      if (z_pt < 3)
        idPhi = numPhiBins-1;
      else if (5 < z_pt) {
        while (idPhi < numPhiBins-1) {
          if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
            break;
          else
            idPhi++;
        }
      }

      if (idPhi >= 0 && idPhi < numPhiBins)
        ZTracksPt[iSpc][iPtZ][iData][idPhi][iCent]->Fill (trkpt, event_weight);

      if (z_pt > 20) {
        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        ZTrackPtPhi[iData][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
      }
    } // end loop over tracks

    for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }
  } // end loop over pp tree


  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      ZTrackPtPhi[iData][iCent][2] = new TH2D (Form ("ZTrackPtPhi_%s_comb_iCent%i", data, iCent), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
      ZPtSpecs[iData][iCent][2] = new TH1D (Form ("ZPtSpec_%s_comb_iCent%i", data, iCent), "", 300, 0, 300);
      ZMYields[iData][iCent][2] = new TH1D (Form ("ZMSpec_%s_comb_iCent%i", data, iCent), "", 40, 76, 106);
      TrackSpec[iData][iCent][2] = new TH1D (Form ("TrackSpec_%s_comb_iCent%i", data, iCent), "", 100, 0, 100);
      //DRDists[iData][iCent][2] = new TH2D (Form ("DRDist_%s_comb_iCent%i", data, iCent), "", 100, 0, 0.1, 80, 0, 0.8);

      ZTrackPtPhi[iData][iCent][2]->Sumw2 ();
      ZPtSpecs[iData][iCent][2]->Sumw2 ();
      ZMYields[iData][iCent][2]->Sumw2 ();
      TrackSpec[iData][iCent][2]->Sumw2 ();
      //DRDists[iData][iCent][2]->Sumw2 ();

      ZTrackPtPhi[iData][iCent][2]->Add (ZTrackPtPhi[iData][iCent][0]);
      ZTrackPtPhi[iData][iCent][2]->Add (ZTrackPtPhi[iData][iCent][1]);
      ZPtSpecs[iData][iCent][2]->Add (ZPtSpecs[iData][iCent][0]);
      ZPtSpecs[iData][iCent][2]->Add (ZPtSpecs[iData][iCent][1]);
      ZMYields[iData][iCent][2]->Add (ZMYields[iData][iCent][0]);
      ZMYields[iData][iCent][2]->Add (ZMYields[iData][iCent][1]);
      TrackSpec[iData][iCent][2]->Add (TrackSpec[iData][iCent][0]);
      TrackSpec[iData][iCent][2]->Add (TrackSpec[iData][iCent][1]);
      //DRDists[iData][iCent][2]->Add (DRDists[iData][iCent][0]);
      //DRDists[iData][iCent][2]->Add (DRDists[iData][iCent][1]);

      for (short iPtZ = 0; iPtZ < nZPtBins; iPtZ++) {
        for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          ZMissingPt[2][iPtZ][iData][iPhi][iCent] = new TH2D (Form ("ZMissingPt_%s_comb_iPtZ%i_iPhi%i_iCent%i", data, iPtZ, iPhi, iCent), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
          ZMissingPt[2][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          ZMissingPt[2][iPtZ][iData][iPhi][iCent]->Add (ZMissingPt[0][iPtZ][iData][iPhi][iCent]);
          ZMissingPt[2][iPtZ][iData][iPhi][iCent]->Add (ZMissingPt[1][iPtZ][iData][iPhi][iCent]);
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          ZTracksPt[2][iPtZ][iData][iPhi][iCent] = new TH1D (Form ("ZTracksPt_%s_comb_iPtZ%i_iPhi%i_iCent%i", data, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          ZTracksPt[2][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          ZTracksPt[2][iPtZ][iData][iPhi][iCent]->Add (ZTracksPt[0][iPtZ][iData][iPhi][iCent]);
          ZTracksPt[2][iPtZ][iData][iPhi][iCent]->Add (ZTracksPt[1][iPtZ][iData][iPhi][iCent]);
        }
      }
    }
  }


  SaveHists ();
  LoadHists ();

  //PlotZPtSpectra ();
  //PlotZMassSpectra ();
  //PlotTrkYield ();

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  fcalWeightsFile->Close ();
  if (fcalWeightsFile) { delete fcalWeightsFile; fcalWeightsFile = nullptr; }

  Delete1DArray (trkPtProj, numPhiBins);
}
