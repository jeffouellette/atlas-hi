#ifndef __Analysis_h__
#define __Analysis_h__

#include "Params.h"
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <iostream>

using namespace atlashi;
using namespace std;

class Analysis {

  private:
  // TFile with histograms
  TFile* histFile = nullptr;
  bool histsLoaded = false;

  vector<TH1*> drawnHists;
  vector<TGraphAsymmErrors*> drawnGraphs;

  protected:
  const char* name;
  const char* directory;
  bool backgroundSubtracted = false;

  public:
  bool plotFill = false; // plot as filled (bar) graph or points w/ errors
  bool useAltMarker = false; // plot as open markers (instead of closed)

  // Event info distributions (for reweighting)
  TH3D* PbPbEventReweights = nullptr;
  TH1D* ppEventReweights   = nullptr;

  // Analysis checks
  TH1D* FCalSpec           = nullptr;
  TH2D* FCalQ2Corr         = nullptr;
  TH1D** ZPhiYields        = nullptr;
  TH1D*** ZPtSpecs         = nullptr;
  TH1D*** ZMYields         = nullptr;
  TH1D** ElectronSpec      = nullptr;
  TH1D** MuonSpec          = nullptr;
  TH1D*** TrackSpec        = nullptr;
  TH2D*** DRDists          = nullptr;
  TH2D*** ZTrackPtPhi      = nullptr;
  TH1D**** ZTrackPhiPtBins = nullptr;
  
  // Physics plots
  TH2D***** ZMissingPt     = nullptr;
  TH1D***** ZMissingPtAvgs = nullptr;
  TH1D***** ZMissingPtInts = nullptr;
  TH1D***** ZTracksPt      = nullptr;
  TH1D****  ZCounts        = nullptr;

  TH1D***** ZTracksSubYields = nullptr;
  TH1D***** ZTracksIAARatios = nullptr;
  TH1D***** ZTracksICPRatios = nullptr;

  Analysis () {
    name = "";
    directory = "";
    backgroundSubtracted = false;

    plotFill = false;
    useAltMarker = false;

    // Reweighting histograms
    PbPbEventReweights = nullptr;
    ppEventReweights   = nullptr;

    // Analysis checks
    FCalSpec        = nullptr;
    FCalQ2Corr      = nullptr;
    ZPhiYields      = Get1DArray <TH1D*> (numCentBins);                // iCent
    ZPtSpecs        = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    ZMYields        = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    ElectronSpec    = Get1DArray <TH1D*> (numCentBins);                // iCent
    MuonSpec        = Get1DArray <TH1D*> (numCentBins);                // iCent
    TrackSpec       = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
    DRDists         = Get2DArray <TH2D*> (numCentBins, 3);             // iCent, iSpc
    ZTrackPtPhi     = Get2DArray <TH2D*> (numCentBins, 3);             // iCent, iSpc (0=ee, 1=mumu, 2=combined)
    ZTrackPhiPtBins = Get3DArray <TH1D*> (nPtTrkBins, numCentBins, 3); // iPtTrk, iCent, iSpc
    
    // Physics plots
    ZMissingPt      = Get4DArray <TH2D*> (3, nPtZBins, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
    ZMissingPtAvgs  = Get4DArray <TH1D*> (3, nPtZBins, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
    ZMissingPtInts  = Get4DArray <TH1D*> (3, nPtZBins, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
    ZTracksPt       = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);    // iSpc, iPtZ, iPhi, iCent
    ZCounts         = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);    // iSpc, iPtZ, iCent

    ZTracksSubYields = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
    ZTracksIAARatios = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
    ZTracksICPRatios = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  }

  Analysis (Analysis* a) {
    // Reweighting histograms
    PbPbEventReweights = a->PbPbEventReweights;
    ppEventReweights   = a->ppEventReweights;

    // Analysis checks
    FCalSpec        = a->FCalSpec;
    FCalQ2Corr      = a->FCalQ2Corr;
    ZPhiYields      = a->ZPhiYields;
    ZPtSpecs        = a->ZPtSpecs;
    ZMYields        = a->ZMYields;
    ElectronSpec    = a->ElectronSpec;
    MuonSpec        = a->MuonSpec;
    TrackSpec       = a->TrackSpec;
    DRDists         = a->DRDists;
    ZTrackPtPhi     = a->ZTrackPtPhi;
    ZTrackPhiPtBins = a->ZTrackPhiPtBins;
    
    // Physics plots
    ZMissingPt      = a->ZMissingPt;
    ZMissingPtAvgs  = a->ZMissingPtAvgs;
    ZMissingPtInts  = a->ZMissingPtInts;
    ZTracksPt       = a->ZTracksPt;
    ZCounts         = a->ZCounts;
  }

  protected:
  void LabeldPhiPtTrk (const short iCent);
  void LabelZMassSpectra (const short iSpc, const short iCent);
  void LabelTrkYield (const short iCent, const short iPhi);

  void GetDrawnObjects ();
  void GetMinAndMax (double &min, double &max, const bool log = false);
  void SetMinAndMax (double min, double max);

  virtual void SubtractBackground ();

  public:
  const char* Name () { return name; }

  void CreateHists ();
  void LoadHists ();
  void SaveHists ();

  void Plot3DDist ();
  void PlotFCalDists ();
  void PlotdPhiPtTrk (const short pSpc = 2);
  void PlotLeptonPtSpectra ();
  void PlotLeptonTrackPtSpectra ();
  void PlotZPtSpectra ();
  void PlotZMassSpectra ();
  void PlotZPhiYield ();
  void PlotZMissingPt ();

  void CalculateIAA ();
  void CalculateICP ();

  void PlotTrkYield (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotIAARatios (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotICPRatios (const short pSpc = 2, const short pPtZ = nPtZBins-1);

};


void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}


void Analysis::GetDrawnObjects () {
  TList* primitives = gPad->GetListOfPrimitives ();
  drawnHists.clear ();
  for (int i = 0; i < primitives->GetSize (); i++) {
    TObject *obj = primitives->At (i);
    if (obj->IsA()->InheritsFrom (TH1::Class ())) {
      drawnHists.push_back ((TH1*)obj);
    }
    else if (obj->IsA()->InheritsFrom (TGraphAsymmErrors::Class ())) {
      drawnGraphs.push_back ((TGraphAsymmErrors*)obj);
    }
  }
}


void Analysis::GetMinAndMax (double &min, double &max, const bool log) {
  for (TH1* h : drawnHists) {
    const double _max = log ? h->GetMaximum (0) : h->GetMaximum ();
    const double _min = log ? h->GetMinimum (0) : h->GetMinimum ();

    if (_max > max) max = _max;
    if (_min < min) min = _min;
  }
  for (TGraphAsymmErrors* g : drawnGraphs) {
    double* ys = g->GetY ();
    for (int n = 0; n < g->GetN (); n++) {
      if (log && ys[n] <= 0)
        continue;
      if (ys[n] > max) max = ys[n];
      if (ys[n] < min) min = ys[n];
    }
  }
  return;
}


void Analysis::SetMinAndMax (double min, double max) {
  for (TH1* h : drawnHists) {
    h->GetYaxis ()->SetRangeUser (min, max);
  }
  for (TGraphAsymmErrors* g : drawnGraphs) {
    g->GetYaxis ()->SetRangeUser (min, max);
  }
  gPad->Update ();
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::CreateHists () {
  FCalSpec = new TH1D (Form ("FCalSpec_%s", name), "", 300, 0, 6000); 
  FCalSpec->Sumw2 ();
  FCalQ2Corr = new TH2D (Form ("FCalQ2Corr_%s", name), "", 300, 0, 6000, 150, 0, 300);
  FCalQ2Corr->Sumw2 ();
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    ZPhiYields[iCent] = new TH1D (Form ("ZPhiYield_iCent%i_%s", iCent, name), "", 80, 0, pi);
    ZPhiYields[iCent]->Sumw2 ();
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      ZTrackPtPhi[iCent][iSpc] = new TH2D (Form ("ZTrackPtPhi_%s_iCent%i_%s", spc, iCent, name), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
      ZPtSpecs[iCent][iSpc] = new TH1D (Form ("ZPtSpec_%s_iCent%i_%s", spc, iCent, name), "", 300, 0, 300);
      ZMYields[iCent][iSpc] = new TH1D (Form ("ZMSpec_%s_iCent%i_%s", spc, iCent, name), "", 40, 76, 106);
      TrackSpec[iCent][iSpc] = new TH1D (Form ("TrackSpec_%s_iCent%i_%s", spc, iCent, name), "", 100, 0, 100);
      
      ZTrackPtPhi[iCent][iSpc]->Sumw2 ();
      ZPtSpecs[iCent][iSpc]->Sumw2 ();
      ZMYields[iCent][iSpc]->Sumw2 ();
      TrackSpec[iCent][iSpc]->Sumw2 ();
    }
    ElectronSpec[iCent] = new TH1D (Form ("ElectronSpec_iCent%i_%s", iCent, name), "", 250, 0, 250);
    MuonSpec[iCent] = new TH1D (Form ("MuonSpec_iCent%i_%s", iCent, name), "", 250, 0, 250);

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          ZMissingPt[iSpc][iPtZ][iPhi][iCent] = new TH2D (Form ("ZMissingPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
          ZMissingPt[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          ZTracksPt[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("ZTracksPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name), "", nPtTrkBins, ptTrkBins);
          ZTracksPt[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        }
        ZCounts[iSpc][iPtZ][iCent] = new TH1D (Form ("ZCounts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name), "", 1, 0, 1);
        ZCounts[iSpc][iPtZ][iCent]->Sumw2 ();
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadHists () {
  SetupDirectories (directory, "ZTrackAnalysis/");
  if (histsLoaded)
    return;
  histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      ZTrackPtPhi[iCent][iSpc] = (TH2D*)histFile->Get (Form ("ZTrackPtPhi_%s_iCent%i_%s", spc, iCent, name));
      ZPtSpecs[iCent][iSpc] = (TH1D*)histFile->Get (Form ("ZPtSpec_%s_iCent%i_%s", spc, iCent, name));
      ZMYields[iCent][iSpc] = (TH1D*)histFile->Get (Form ("ZMSpec_%s_iCent%i_%s", spc, iCent, name));
      TrackSpec[iCent][iSpc] = (TH1D*)histFile->Get (Form ("TrackSpec_%s_iCent%i_%s", spc, iCent, name));
      DRDists[iCent][iSpc] = (TH2D*)histFile->Get (Form ("DRDist_%s_iCent%i_%s", spc, iCent, name));
    }
    ZPhiYields[iCent] = (TH1D*)histFile->Get (Form ("ZPhiYield_iCent%i_%s", iCent, name));
    ElectronSpec[iCent] = (TH1D*)histFile->Get (Form ("ElectronSpec_iCent%i_%s", iCent, name));
    MuonSpec[iCent] = (TH1D*)histFile->Get (Form ("MuonSpec_iCent%i_%s", iCent, name));
    
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          ZMissingPt[iSpc][iPtZ][iPhi][iCent] = (TH2D*)histFile->Get (Form ("ZMissingPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name));
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          ZTracksPt[iSpc][iPtZ][iPhi][iCent] = (TH1D*)histFile->Get (Form ("ZTracksPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name));
        }
        ZCounts[iSpc][iPtZ][iCent] = (TH1D*)histFile->Get (Form ("ZCounts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name));
      }
    }
  }
  FCalSpec = (TH1D*)histFile->Get (Form ("FCalSpec_%s", name));
  FCalQ2Corr = (TH2D*)histFile->Get (Form ("FCalQ2Corr_%s", name));
  
  histsLoaded = true;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::SaveHists () {
  SetupDirectories (directory, "ZTrackAnalysis/");
  histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "recreate");
  histFile->cd ();
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      SafeWrite (ZTrackPtPhi[iCent][iSpc]);
      SafeWrite (ZPtSpecs[iCent][iSpc]);
      SafeWrite (ZMYields[iCent][iSpc]);
      SafeWrite (TrackSpec[iCent][iSpc]);
      SafeWrite (DRDists[iCent][iSpc]);
    }
    SafeWrite (ZPhiYields[iCent]);
    SafeWrite (ElectronSpec[iCent]);
    SafeWrite (MuonSpec[iCent]);
    
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          SafeWrite (ZMissingPt[iSpc][iPtZ][iPhi][iCent]);
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SafeWrite (ZTracksPt[iSpc][iPtZ][iPhi][iCent]);
        }
        SafeWrite (ZCounts[iSpc][iPtZ][iCent]);
      }
    }
  }
  SafeWrite (FCalSpec);
  SafeWrite (FCalQ2Corr);
  
  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - xZTrk 3d distribution
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::Plot3DDist () {
  TCanvas* c = new TCanvas ("3dCanvas", "", 800, 600);
  c->cd ();

  c->SetLogy ();

  ZTrackPtPhi[numCentBins-1][2]->GetXaxis ()->SetTitle ("#phi_{Z} - #phi_{Trk}");
  ZTrackPtPhi[numCentBins-1][2]->GetYaxis ()->SetTitle ("#it{p}_{T}^{ ch} / #it{p}_{T}^{ Z}");

  //ZTrackPtPhi[0][numCentBins-1][2]->RebinY (2);
  ZTrackPtPhi[numCentBins-1][2]->Draw ("lego2");

  c->SaveAs (Form ("%s/ZTrackCorr.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotFCalDists () {
  TCanvas* c = new TCanvas ("FCalCanvas", "", 800, 600);
  c->cd ();

  c->SetLogy ();

  FCalSpec->Scale (1./FCalSpec->Integral ());

  FCalSpec->SetLineColor (kBlack);

  FCalSpec->GetYaxis ()->SetRangeUser (1.5e-4, 1e-1);

  FCalSpec->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
  FCalSpec->GetYaxis ()->SetTitle ("A.U.");

  FCalSpec->Draw ("hist");

  //myText (0.65, 0.88, kBlack, "2018 Pb+Pb", 0.04);
  //myText (0.65, 0.81, kBlue, "Pythia8 + Hijing", 0.04);

  c->SaveAs (Form ("%s/FCalDist.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - pTTrk 2d projections
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotdPhiPtTrk (const short pSpc) {
  TF1**** Fits = Get3DArray <TF1*> (nPtTrkBins, numCentBins, 3); // iPtTrk, iCent, iSpc
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    const char* canvasName = Form ("dPhiPtTrkCanvas_%s", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 1600, 1200);
      gDirectory->Add (c);
      c->cd ();
      c->Divide (2, 2);
    }
    c->cd ();
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      c->cd (iCent+1);
      GetDrawnObjects ();

      gPad->SetTopMargin (0.01);
      gPad->SetBottomMargin (0.12);
      gPad->SetRightMargin (0.01);
      gPad->SetLeftMargin (0.12);

      for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        ZTrackPhiPtBins[iPtTrk][iCent][iSpc] = ZTrackPtPhi[iCent][iSpc]->ProjectionX (Form ("ZTrackPhi_iPtTrk%i_iCent%i_%s", iPtTrk, iCent, name), iPtTrk+1, iPtTrk+1);
        TH1D* thisHist = ZTrackPhiPtBins[iPtTrk][iCent][iSpc];
        thisHist->Rebin (2);
        if (iPtTrk > 3)
          thisHist->Rebin (2);
        if (iCent != 0)
          thisHist->Rebin (2);
        if (ZPtSpecs[iCent][iSpc]->Integral () > 0)
          thisHist->Scale (1./ZPtSpecs[iCent][iSpc]->Integral (), "width");

        TF1* constFit = new TF1 (Form ("fit_iPtTrk%i_iCent%i_%s", iPtTrk, iCent, name), "[0]+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x)", -pi/2, 3*pi/2);
        thisHist->Fit (constFit, "RN0");
        //const float min = constFit->GetMinimum (-pi/2, 3*pi/2);
        const double min = thisHist->Integral (thisHist->GetXaxis ()->FindBin (-pi/3), thisHist->GetXaxis ()->FindBin (pi/3)) / (thisHist->GetXaxis      ()->FindBin (pi/3) - thisHist->GetXaxis ()->FindBin (-pi/3) + 1);

        for (int ix = 1; ix <= thisHist->GetNbinsX (); ix++)
          thisHist->SetBinContent (ix, thisHist->GetBinContent (ix) - min); 

        constFit->SetParameter (0, constFit->GetParameter (0) - min);
        Fits[iPtTrk][iCent][iSpc] = constFit;
      } // end loop over pT^trk bins
      
      double min = 1e30, max = 0;
      GetMinAndMax (min, max);
      for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
        TH1D* thisHist = ZTrackPhiPtBins[iPtTrk][iCent][iSpc];
        if (thisHist->GetMinimum () < min) min = thisHist->GetMinimum ();
        if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
      } // end loop over pT^trk bins
      SetMinAndMax (min, max);

      if (plotFill) {
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
          TH1D* thisHist = (TH1D*)ZTrackPhiPtBins[iPtTrk][iCent][iSpc]->Clone ();

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

          thisHist->DrawCopy (!canvasExists && iPtTrk == 0 ? "b" : "same b");
          thisHist->SetLineWidth (1);
          thisHist->Draw ("hist same");
        }
        gPad->RedrawAxis ();
      }
      else {
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
          TGraphAsymmErrors* thisGraph = make_graph (ZTrackPhiPtBins[iPtTrk][iCent][iSpc]);

          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
          thisGraph->GetYaxis ()->SetRangeUser (0, 1.3*max);

          thisGraph->SetMarkerStyle (markerStyle);
          thisGraph->SetLineColor (colors[iPtTrk]);
          thisGraph->SetMarkerColor (colors[iPtTrk]);

          thisGraph->Draw (!canvasExists && iPtTrk == 0 ? "AP" : "P");

          Fits[iPtTrk][iCent][iSpc]->SetLineColor (colors[iPtTrk]);
          Fits[iPtTrk][iCent][iSpc]->Draw ("same");
        } // end loop over pT^trk bins
      }

      LabeldPhiPtTrk (iCent);

    } // end loop over centrality

    c->SaveAs (Form ("%s/dPhi_pTtrk_%s.pdf", plotPath.Data (), spc));

  }
  Delete3DArray (Fits, nPtTrkBins, numCentBins, 3);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track dPhi distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::LabeldPhiPtTrk (const short iCent) {
  for (int iPtTrk = 0; iPtTrk < nPtTrkBins-1; iPtTrk++) {
    if (iCent == 0) {
      myText (0.2, 0.9, kBlack, "#it{pp}", 0.06);
      const float pt_lo = ptTrkBins[iPtTrk];
      const float pt_hi = ptTrkBins[iPtTrk+1];
      myText (0.2, 0.80-0.075*iPtTrk, colors[iPtTrk], Form ("%.1f < #it{p}_{T}^{ ch} < %.1f GeV", pt_lo, pt_hi), 0.06);
    }
    else {
      myText (0.2, 0.9, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot lepton Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotLeptonPtSpectra () {
  const char* canvasName = "LeptonPtCanvas";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1000, 600);
    gDirectory->Add (c);
  }
  c->cd ();
  c->SetLogy ();
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TH1D* thisHist = (TH1D*)ElectronSpec[iCent]->Clone ();
    thisHist->Rebin (5);
    thisHist->Scale (1./ZPtSpecs[iCent][0]->Integral (), "width");

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
    TH1D* thisHist = (TH1D*)MuonSpec[iCent]->Clone ();
    thisHist->Scale (1./ZPtSpecs[iCent][1]->Integral (), "width");
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
void Analysis::PlotLeptonTrackPtSpectra () {
  const char* canvasName = "LeptonTrackPtCanvas";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
  }
  c->cd ();
  c->SetLogy ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TH1D* thisHist = (TH1D*)TrackSpec[iCent][iSpc]->Clone (Form ("leptonTrackPt_iSpc%i_iCent%i", iSpc, iCent));
      thisHist->Rebin (5);
      thisHist->Scale (1./ZPtSpecs[iCent][iSpc]->Integral (ZPtSpecs[iCent][iSpc]->GetXaxis ()->FindBin (5), ZPtSpecs[iCent][iSpc]->GetNbinsX ()), "width");
      cout << iSpc << ", " << iCent << ", " << ZPtSpecs[iCent][iSpc]->Integral (ZPtSpecs[iCent][iSpc]->GetXaxis ()->FindBin (5), ZPtSpecs[iCent][iSpc]->GetNbinsX ()) << endl;

      thisHist->GetYaxis ()->SetRangeUser (6e-6, 450);

      thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
      if (iSpc == 0) {
        thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrowee} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }
      else if (iSpc == 1) {
        thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrow#mu#mu} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
      }
      else {
        thisHist->GetYaxis ()->SetTitle ("1/N_{Z#rightarrowll} dN_{ch}/d#it{p}_{T} [GeV^{-1}]");
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

    myText (0.66, 0.88, kBlack, iSpc == 0 ? "Z#rightarrowee Events" : (iSpc == 1 ?"Z#rightarrow#mu#mu Events" : "Z#rightarrowll Events"), 0.04);
    c->SaveAs (Form ("%s/%sTrackPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Comb")));
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotZPtSpectra () {
  const char* canvasName = "ZPtCanvas";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
  }
  c->cd ();
  c->SetLogy ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    //TH1D* thisHist = (TH1D*)ZPtSpecs[0][iCent][0]->Clone ();
    TH1D* thisHist = (TH1D*)ZPtSpecs[iCent][1]->Clone ();
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
void Analysis::PlotZMassSpectra () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    const char* canvasName = Form ("ZMassCanvas_%s", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 800, 600);
      gDirectory->Add (c);
    }
    c->cd ();

    if (plotFill) {
      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        TH1D* thisHist = ZMYields[iCent][iSpc];
        if (thisHist->GetMaximum () > 0)
          thisHist->Scale (1. / (thisHist->GetMaximum ()));

        thisHist->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        thisHist->SetLineColor (kBlack);
        thisHist->SetMarkerSize (0);
        thisHist->SetLineWidth (0);
        thisHist->GetYaxis ()->SetRangeUser (0, 1.3);

        thisHist->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
        thisHist->GetYaxis ()->SetTitle ("A.U.");

        thisHist->DrawCopy (!canvasExists && iCent == numCentBins-1 ? "bar" : "bar same");
        thisHist->SetLineWidth (1);
        thisHist->Draw ("hist same");

        LabelZMassSpectra (iSpc, iCent);
      }
      gPad->RedrawAxis ();
    }
    else {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        TH1D* thisHist = ZMYields[iCent][iSpc];
        if (thisHist->GetMaximum () > 0)
          thisHist->Scale (1. / (thisHist->GetMaximum ()));

        TGraphAsymmErrors* thisGraph = make_graph (thisHist);
        ResetXErrors (thisGraph);
        deltaize (thisGraph, 0.1*(-1.5+iCent));

        const int markerStyle = kFullCircle;
        thisGraph->SetMarkerStyle (markerStyle);
        thisGraph->SetMarkerSize (1);
        thisGraph->SetLineWidth (1);
        thisGraph->SetLineColor (colors[iCent]);
        thisGraph->SetMarkerColor (colors[iCent]);
        thisGraph->GetYaxis ()->SetRangeUser (0, 1.3);

        thisGraph->GetXaxis ()->SetTitle ("m_{Z} [GeV]");
        thisGraph->GetYaxis ()->SetTitle ("A.U.");
        thisGraph->Draw (!canvasExists && iCent == numCentBins-1 ? "AP" : "P");

        LabelZMassSpectra (iSpc, iCent);
      }
    }
    c->SaveAs (Form ("%s/z%s_mass_spectrum.pdf", plotPath.Data (), spc));
  }

}


void Analysis::LabelZMassSpectra (const short iSpc, const short iCent) {
  const char* spc = iSpc == 0 ? "Z#rightarrowee Events":"Z#rightarrow#mu#mu Events";
  if (iCent == 0) {
    myText (0.66, 0.88, kBlack, spc, 0.04);
    myMarkerText (0.25, 0.88, colors[0], kFullCircle, Form ("#it{pp}"), 1.25, 0.04);
  }
  else
    myMarkerText (0.25, 0.88-0.07*iCent, colors[iCent], kFullCircle, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 1.25, 0.04);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotZPhiYield () {
  const char* canvasName = "ZPhiYieldCanvas";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
  }
  c->cd ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TH1D* thisHist = ZPhiYields[iCent];
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
void Analysis::PlotZMissingPt () {
  const char* canvasName = "ZMissingPtCanvas";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600*numPhiTrkBins, 500);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (2, 1);
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        c->cd (iPhi+1);
        gPad->SetLogx ();
        for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
          TH1D* thisHist = new TH1D (Form ("ZMissingPtAvg_%s_iPtZ%i_iPhi%i_iCent%i", spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          TH1D* integral = new TH1D (Form ("ZMissingPtInt_%s_iPtZ%i_iPhi%i_iCent%i", spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          thisHist->Sumw2 ();
          integral->Sumw2 ();
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            TH1D* px = ZMissingPt[iSpc][iPtZ][iPhi][iCent]->ProjectionX ("_px", iPtTrk+1, iPtTrk+1);
            TF1* fit = new TF1 ("fit", "gaus(0)", zMissingPtBins[0], zMissingPtBins[numZMissingPtBins]);
            px->Fit (fit, "RN0Q");
            thisHist->SetBinContent (iPtTrk+1, fit->GetParameter (1));
            thisHist->SetBinError (iPtTrk+1, fit->GetParError (1));
            if (px) delete px;

            px = ZMissingPt[iSpc][iPtZ][iPhi][iCent]->ProjectionX ("_px", 0, iPtTrk+1);
            px->Fit (fit, "RN0Q");
            integral->SetBinContent (iPtTrk+1, fit->GetParameter (1));
            integral->SetBinError (iPtTrk+1, fit->GetParError (1));
            if (px) delete px;

            if (fit) delete fit;
          }
          ZMissingPtAvgs[iSpc][iPtZ][iPhi][iCent] = thisHist;
          ZMissingPtInts[iSpc][iPtZ][iPhi][iCent] = integral;

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
          TGraphAsymmErrors* thisGraph = make_graph (ZMissingPtAvgs[iSpc][iPtZ][iPhi][iCent]);
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
// Subtracts default (HM) background from track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::SubtractBackground () {
  if (backgroundSubtracted)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* subYield = new TH1D (Form ("defaultSubYield_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name), "", nPtTrkBins, ptTrkBins);
          subYield->Sumw2 ();

          subYield->Add (ZTracksPt[iSpc][iPtZ][iPhi][iCent]);
          subYield->Add (ZTracksPt[iSpc][iPtZ][0][iCent], -1);

          ZTracksSubYields[iSpc][iPtZ][iPhi][iCent] = subYield;
        } // end loop over phi
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  backgroundSubtracted = true;
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTtrk binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotTrkYield (const short pSpc, const short pPtZ) {
  const double padRatio = 1.1; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0 / (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("TrkYieldCanvas_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 860*numCentBins, 1000);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* topPadName = Form ("topPad_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
        const char* bottomPadName = Form ("bottomPad_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);

        TPad* topPad = nullptr, *bottomPad = nullptr;
        if (!canvasExists) {
          topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
          bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

          gDirectory->Add (topPad);
          gDirectory->Add (bottomPad);

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
        else {
          topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
          bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
        }

        topPad->cd ();
        GetDrawnObjects ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        double min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
          if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
          if (thisHist->GetMaximum () > max)  max = thisHist->GetMaximum ();
        } // end loop over phi
        min = 0.5*min;
        max = 2*max;
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];

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

            thisHist->GetYaxis ()->SetRangeUser (min, max);
            thisHist->DrawCopy (!canvasExists && iPhi == numPhiBins-1 ? "bar" : "bar same");
            thisHist->SetLineWidth (1);
            thisHist->Draw ("hist same");
          }
          gPad->RedrawAxis ();
        }
        else {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
            
            TGraphAsymmErrors* thisGraph = make_graph (thisHist);
            RecenterGraph (thisGraph);
            ResetXErrors (thisGraph);

            thisGraph->SetMarkerStyle (markerStyle);
            thisGraph->SetMarkerColor (colors[iPhi]);
            thisGraph->SetLineColor (colors[iPhi]);
            thisGraph->SetMarkerSize (1);
            thisGraph->SetLineWidth (2);
  
            thisGraph->GetYaxis ()->SetRangeUser (min, max);
            thisGraph->Draw (!canvasExists && iPhi == 0 ? "AP" : "P");
          } // end loop over phi
        }

        if (!canvasExists)
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++)
            LabelTrkYield (iCent, iPhi);

        bottomPad->cd ();
        GetDrawnObjects ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        if (!backgroundSubtracted)
          SubtractBackground ();
        
        min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* thisHist = ZTracksSubYields[iSpc][iPtZ][iPhi][iCent];
          if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
          if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
        } // end loop over phi
        float delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* subYield = ZTracksSubYields[iSpc][iPtZ][iPhi][iCent];

            subYield->GetYaxis ()->SetRangeUser (min, max);

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

            subYield->DrawCopy (!canvasExists && iPhi == numPhiBins-1 ? "bar" : "bar same");
            subYield->SetLineWidth (1);
            subYield->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();
        }
        else {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGraphAsymmErrors* subYield = make_graph (ZTracksSubYields[iSpc][iPtZ][iPhi][iCent]);
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

            if (!canvasExists && iPhi == numPhiBins-1)
              subYield->Draw ("AP");
            else
              subYield->Draw ("P");
          } // end loop over phi
        }
      } // end loop over cents
      
      c->SaveAs (Form ("%s/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track pT distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::LabelTrkYield (const short iCent, const short iPhi) {
  const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);

  if (iCent == 0)
    myText (0.22, 0.06, kBlack, "#it{pp}", 0.06);
  else
    myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

  if (iCent == 0)
    myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
  else if (iCent == numCentBins-1) {
    if (iPhi == 0) {
      myText (0.43, 0.91, kBlack, "Data", 0.06);
      myText (0.574, 0.91, kBlack, "MC", 0.06);
      myText (0.683, 0.91, kBlack, "#Delta#phi", 0.06);
    }

    TVirtualPad* cPad = gPad; // store current pad
    const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
    const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
    TBox* b = TBoxNDC (0.61-0.024, 0.85-0.06*iPhi-0.016, 0.61+0.024, 0.85-0.06*iPhi+0.016);
    b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    b->Draw ("l");
    myMarkerText (0.512, 0.852-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.06);
    myText (0.68, 0.85-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.06);
    cPad->cd ();
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::CalculateIAA () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
        TH1D* ppHist = ZTracksSubYields[iSpc][iPtZ][iPhi][0];

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          if (!ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent]) {
            TH1D* PbPbHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iPhi][iCent]->Clone ());
            PbPbHist->Divide (ppHist);
            ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
          } 
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotIAARatios (const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("IAACanvas_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 860*numCentBins, 1000);
        c->Divide (numCentBins-1, 1);
        gDirectory->Add (c);
      }

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TH1D* ratioHist = ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent];

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

            ratioHist->DrawCopy (!canvasExists && iPhi == 1 ? "bar" : "bar same");
            ratioHist->SetLineWidth (1);
            ratioHist->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();
        }
        else {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TGraphAsymmErrors* ratioGraph = make_graph (ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (ratioGraph);
            ResetXErrors (ratioGraph);
            deltaize (ratioGraph, 1+(0.5*(numPhiBins-1)-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-1)

            ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
            ratioGraph->SetLineColor (colors[iPhi]);
            ratioGraph->SetMarkerColor (colors[iPhi]);
            ratioGraph->SetMarkerSize (1.2);
            ratioGraph->SetLineWidth (2);
            ratioGraph->SetMarkerStyle (kFullCircle);

            ratioGraph->Draw (!canvasExists && iPhi == 1 ? "AP" : "P");

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
              cPad->cd ();
            }
          } // end loop over phi
        }
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
// Calculates subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::CalculateICP () {
  if (!backgroundSubtracted)
    SubtractBackground ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
        TH1D* periphHist = ZTracksSubYields[iSpc][iPtZ][iPhi][1];

        for (short iCent = 2; iCent < numCentBins; iCent++) {
          if (!ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent]) {
            TH1D* centHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iPhi][iCent]->Clone ());
            centHist->Divide (periphHist);
            ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent] = centHist;
          }
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotICPRatios (const short pSpc, const short pPtZ) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("ICPCanvas_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 600*(numCentBins-1), 500);
        c->Divide (numCentBins-2, 1);
        gDirectory->Add (c);
      }

      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        gPad->SetLogx ();

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TH1D* ratioHist = ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent];

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

            ratioHist->DrawCopy (!canvasExists && iPhi == 1 ? "bar" : "bar same");
            ratioHist->SetLineWidth (1);
            ratioHist->Draw ("hist same");
          } // end loop over phi
          gPad->RedrawAxis ();
        }
        else {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            c->cd (iCent-1);
            
            TGraphAsymmErrors* ratioGraph = make_graph (ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (ratioGraph);
            ResetXErrors (ratioGraph);
            deltaize (ratioGraph, 1+(0.5*(numPhiBins-1)-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-1)

            //ratioGraph->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
            ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
            ratioGraph->SetLineColor (colors[iPhi]);
            ratioGraph->SetMarkerColor (colors[iPhi]);
            ratioGraph->SetMarkerSize (1.2);
            ratioGraph->SetLineWidth (2);
            ratioGraph->SetMarkerStyle (kFullCircle);

            ratioGraph->Draw (!canvasExists && iPhi == 1 ? "AP" : "P");

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
        }
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

#endif
