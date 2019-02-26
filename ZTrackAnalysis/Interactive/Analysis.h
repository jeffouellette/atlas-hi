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

  public:

  // Event info distributions (for reweighting)
  TH3D** PbPbEventInfoDist  = nullptr;
  TH3D** PbPbEventReweights = nullptr;
  TH1D** ppEventInfoDist    = nullptr;
  TH1D** ppEventReweights   = nullptr;

  // Analysis checks
  TH1D** FCalSpec           = nullptr;
  TH2D** FCalQ2Corr         = nullptr;
  TH1D*** ZPhiYields        = nullptr;
  TH1D**** ZPtSpecs         = nullptr;
  TH1D**** ZMYields         = nullptr;
  TH1D*** ElectronSpec      = nullptr;
  TH1D*** MuonSpec          = nullptr;
  TH1D**** TrackSpec        = nullptr;
  TH2D**** DRDists          = nullptr;
  TH2D**** ZTrackPtPhi      = nullptr;
  TH1D***** ZTrackPhiPtBins = nullptr;
  
  // Physics plots
  TH2D****** ZMissingPt     = nullptr;
  TH1D****** ZMissingPtAvgs = nullptr;
  TH1D****** ZMissingPtInts = nullptr;
  TH1D****** ZTracksPt      = nullptr;
  TH1D****** ZTracksDefSig  = nullptr;
  TH1D***** ZCounts         = nullptr;

  Analysis () {
    // Reweighting histograms
    PbPbEventInfoDist  = Get1DArray <TH3D*> (2); // iData
    PbPbEventReweights = Get1DArray <TH3D*> (2); // iData
    ppEventInfoDist    = Get1DArray <TH1D*> (2); // iData
    ppEventReweights   = Get1DArray <TH1D*> (2); // iData

    // Analysis checks
    FCalSpec        = Get1DArray <TH1D*> (2);                             // iData
    FCalQ2Corr      = Get1DArray <TH2D*> (2);                             // iData
    ZPhiYields      = Get2DArray <TH1D*> (2, numCentBins);                // iData, iCent
    ZPtSpecs        = Get3DArray <TH1D*> (2, numCentBins, 3);             // iData, iCent, iSpc
    ZMYields        = Get3DArray <TH1D*> (2, numCentBins, 3);             // iData, iCent, iSpc
    ElectronSpec    = Get2DArray <TH1D*> (2, numCentBins);                // iData, iCent
    MuonSpec        = Get2DArray <TH1D*> (2, numCentBins);                // iData, iCent
    TrackSpec       = Get3DArray <TH1D*> (2, numCentBins, 3);             // iData, iCent, iSpc
    DRDists         = Get3DArray <TH2D*> (2, numCentBins, 3);             // iData, iCent, iSpc
    ZTrackPtPhi     = Get3DArray <TH2D*> (2, numCentBins, 3);             // iData, iCent, iSpc (0=ee, 1=mumu, 2=combined)
    ZTrackPhiPtBins = Get4DArray <TH1D*> (2, nPtTrkBins, numCentBins, 3); // iData, iPtTrk, iCent, iSpc
    
    // Physics plots
    ZMissingPt      = Get5DArray <TH2D*> (3, nPtZBins, 2, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iData, iPhi, iCent
    ZMissingPtAvgs  = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iData, iPhi, iCent
    ZMissingPtInts  = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiTrkBins, numCentBins); // iSpc, iPtZ, iData, iPhi, iCent
    ZTracksPt       = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiBins, numCentBins);    // iSpc, iPtZ, iData, iPhi, iCent
    ZTracksDefSig   = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiBins, numCentBins);    // iSpc, iPtZ, iData, iPhi, iCent
    ZCounts         = Get4DArray <TH1D*> (3, nPtZBins, 2, numCentBins);    // iSpc, iPtZ, iData, iCent
  }

  Analysis (Analysis* a) {
    // Reweighting histograms
    PbPbEventInfoDist  = a->PbPbEventInfoDist;
    PbPbEventReweights = a->PbPbEventReweights;
    ppEventInfoDist    = a->ppEventInfoDist;
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
    // DRDists        = a->DRDists;
    ZTrackPtPhi     = a->ZTrackPtPhi;
    ZTrackPhiPtBins = a->ZTrackPhiPtBins;
    
    // Physics plots
    ZMissingPt      = a->ZMissingPt;
    ZMissingPtAvgs  = a->ZMissingPtAvgs;
    ZMissingPtInts  = a->ZMissingPtInts;
    ZTracksPt       = a->ZTracksPt;
  }

  public:
  void LoadHists (const char* name);
  void SaveHists ();

  void Plot3DDist ();
  void PlotFCalDists ();
  void PlotdPhiPtTrk ();
  void PlotLeptonPtSpectra ();
  void PlotLeptonTrackPtSpectra ();
  void PlotZPtSpectra ();
  void PlotZMassSpectra ();
  void PlotZPhiYield ();
  void PlotZMissingPt ();
  void PlotTrkYield ();
  void PlotIAARatios ();
  void PlotICPRatios ();

};


void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadHists (const char* name) {
  if (histsLoaded)
    return;
  histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");

  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        ZTrackPtPhi[iData][iCent][iSpc] = (TH2D*)histFile->Get (Form ("ZTrackPtPhi_%s_%s_iCent%i_%s", data, spc, iCent, name));
        ZPtSpecs[iData][iCent][iSpc] = (TH1D*)histFile->Get (Form ("ZPtSpec_%s_%s_iCent%i_%s", data, spc, iCent, name));
        ZMYields[iData][iCent][iSpc] = (TH1D*)histFile->Get (Form ("ZMSpec_%s_%s_iCent%i_%s", data, spc, iCent, name));
        TrackSpec[iData][iCent][iSpc] = (TH1D*)histFile->Get (Form ("TrackSpec_%s_%s_iCent%i_%s", data, spc, iCent, name));
        DRDists[iData][iCent][iSpc] = (TH2D*)histFile->Get (Form ("DRDist_%s_%s_iCent%i_%s", data, spc, iCent, name));
      }
      ZPhiYields[iData][iCent] = (TH1D*)histFile->Get (Form ("ZPhiYield_%s_iCent%i_%s", data, iCent, name));
      ElectronSpec[iData][iCent] = (TH1D*)histFile->Get (Form ("ElectronSpec_%s_iCent%i_%s", data, iCent, name));
      MuonSpec[iData][iCent] = (TH1D*)histFile->Get (Form ("MuonSpec_%s_iCent%i_%s", data, iCent, name));
      
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent] = (TH2D*)histFile->Get (Form ("ZMissingPt_%s_%s_iPtZ%i_iPhi%i_iCent%i_%s", data, spc, iPtZ, iPhi, iCent, name));
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent] = (TH1D*)histFile->Get (Form ("ZTracksPt_%s_%s_iPtZ%i_iPhi%i_iCent%i_%s", data, spc, iPtZ, iPhi, iCent, name));
          }
          ZCounts[iSpc][iPtZ][iData][iCent] = (TH1D*)histFile->Get (Form ("ZCounts_%s_%s_iPtZ%i_iCent%i_%s", data, spc, iPtZ, iCent, name));
        }
      }
    }
    PbPbEventInfoDist[iData] = (TH3D*)histFile->Get (Form ("PbPbEventInfoDist_%s_%s", data, name));
    ppEventInfoDist[iData] = (TH1D*)histFile->Get (Form ("ppEventInfoDist_%s_%s", data, name));

    FCalSpec[iData] = (TH1D*)histFile->Get (Form ("FCalSpec_%s_%s", data, name));
    FCalQ2Corr[iData] = (TH2D*)histFile->Get (Form ("FCalQ2Corr_%s_%s", data, name));
  }
  histsLoaded = true;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::SaveHists () {
  histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "recreate");
  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        SafeWrite (ZTrackPtPhi[iData][iCent][iSpc]);
        SafeWrite (ZPtSpecs[iData][iCent][iSpc]);
        SafeWrite (ZMYields[iData][iCent][iSpc]);
        SafeWrite (TrackSpec[iData][iCent][iSpc]);
        SafeWrite (DRDists[iData][iCent][iSpc]);
      }
      SafeWrite (ZPhiYields[iData][iCent]);
      SafeWrite (ElectronSpec[iData][iCent]);
      SafeWrite (MuonSpec[iData][iCent]);
      
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            SafeWrite (ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]);
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            SafeWrite (ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]);
          }
          SafeWrite (ZCounts[iSpc][iPtZ][iData][iCent]);
        }
      }
    }
    SafeWrite (PbPbEventInfoDist[iData]);
    SafeWrite (ppEventInfoDist[iData]);
    SafeWrite (FCalSpec[iData]);
    SafeWrite (FCalQ2Corr[iData]);
  }
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

  ZTrackPtPhi[0][numCentBins-1][2]->GetXaxis ()->SetTitle ("#phi_{Z} - #phi_{Trk}");
  ZTrackPtPhi[0][numCentBins-1][2]->GetYaxis ()->SetTitle ("#it{p}_{T}^{ ch} / #it{p}_{T}^{ Z}");

  //ZTrackPtPhi[0][numCentBins-1][2]->RebinY (2);
  ZTrackPtPhi[0][numCentBins-1][2]->Draw ("lego2");

  c->SaveAs (Form ("%s/ZTrackCorr.pdf", plotPath.Data ()));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::PlotFCalDists () {
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
void Analysis::PlotdPhiPtTrk () {
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
void Analysis::PlotLeptonPtSpectra () {
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
void Analysis::PlotLeptonTrackPtSpectra () {
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
void Analysis::PlotZPtSpectra () {
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
void Analysis::PlotZMassSpectra () {
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
void Analysis::PlotZPhiYield () {
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
void Analysis::PlotZMissingPt () {
  TCanvas* c = new TCanvas ("ZMissingPtCanvas", "", 600*numPhiTrkBins, 500);
  c->cd ();

  c->Divide (2, 1);

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
void Analysis::PlotTrkYield () {
  TCanvas* c = nullptr;
  if (gDirectory->Get ("TrkYieldCanvas") != nullptr)
    c = dynamic_cast<TCanvas*>(gDirectory->Get ("TrkYieldCanvas"));
  else
    c = new TCanvas ("TrkYieldCanvas", "", 860*numCentBins, 1000);
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

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        TPad* topPad = topPads[iCent];
        TPad* bottomPad = bottomPads[iCent];

        for (short iData = 0; iData < 2; iData++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {

            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];
            double yieldNormFactor = ZCounts[iSpc][iPtZ][iData][iCent]->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
            double yieldNormFactorError = ZCounts[iSpc][iPtZ][iData][iCent]->GetBinError (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);

            //RescaleWithError (thisHist, yieldNormFactor, yieldNormFactorError);
            if (yieldNormFactor > 0)
              thisHist->Scale (1. / yieldNormFactor);
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

        for (short iData = 0; iData < 2; iData++) {
          const char* data = iData == 0 ? "data" : "mc";
          for (int iPhi = numPhiBins-2; iPhi >= 0; iPhi--) {
            TH1D* subYield = new TH1D (Form ("defSubYield_%s_%s_iPtZ%i_iPhi%i_iCent%i", spc, data, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
            subYield->Sumw2 ();
            //subYield->Rebin (8);

            subYield->Add (ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]);
            subYield->Add (ZTracksPt[iSpc][iPtZ][iData][0][iCent], -1);

            ZTracksDefSig[iSpc][iPtZ][iData][iPhi][iCent] = subYield;
          } // end loop over phi
        } // end loop over data types
        
        min = 1e30;
        max = 0;
        for (short iData = 0; iData < 2; iData++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TH1D* thisHist = ZTracksDefSig[iSpc][iPtZ][iData][iPhi][iCent];
            if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
            if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
          } // end loop over phi
        } // end loop over data types

        for (int iPhi = numPhiBins-2; iPhi >= 1; iPhi--) {
          TH1D* subYield = ZTracksDefSig[iSpc][iPtZ][1][iPhi][iCent];

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

          if (iPhi == numPhiBins-2)
            subYield->DrawCopy ("bar");
          else
            subYield->DrawCopy ("bar same");
          subYield->SetLineWidth (1);
          subYield->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();

        for (int iPhi = numPhiBins-2; iPhi >= 1; iPhi--) {
          const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
          TGraphAsymmErrors* subYield = make_graph (ZTracksDefSig[iSpc][iPtZ][0][iPhi][iCent]);
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

#endif
