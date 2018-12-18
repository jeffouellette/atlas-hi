#include "PhotonAnalysisHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace atlashi;

namespace JetCalibration {

void PhotonAnalysisHist () {

  SetAtlasStyle ();

  SetupDirectories ("PhotonAnalysis/", "JetCalibration/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* inTree = (TTree*)inFile->Get ("jeffsphotons");
  float ppt = 0, peta = 0, pphi = 0;
  double evtWeight = 0;
  bool isPeriodA = false, isMC = false;
  inTree->SetBranchAddress ("photon_pt", &ppt);
  inTree->SetBranchAddress ("photon_eta", &peta);
  inTree->SetBranchAddress ("photon_phi", &pphi);
  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);

  const int nPhotons = inTree->GetEntries ();

  TH2D**** photonSpectrum = Get3DArray <TH2D*> (3, 2, 2); // iPer, iData, iWgt
  TH2D**** photonEtaPhi = Get3DArray <TH2D*> (3, 2, 2); // iPer, iData, iWgt

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* per = iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB");

    for (short iData = 0; iData < 2; iData++) {
      const char* data = iData == 0 ? "data" : "mc";

      for (short iWgt = 0; iWgt < 2; iWgt++) {
        const char* wgt = iWgt == 0 ? "weighted" : "unweighted";

        photonSpectrum[iPer][iData][iWgt] = new TH2D (Form ("photonSpectrum_%s_%s_%s", per, data, wgt), "", 50, logspace (20, 500, 50), 240, 0, 2.40);
        photonSpectrum[iPer][iData][iWgt]->Sumw2 ();

        photonEtaPhi[iPer][iData][iWgt] = new TH2D (Form ("photonEtaPhi_%s_%s_%s", per, data, wgt), "", 240, -2.4, 2.4, numphibins, phibins);
        photonEtaPhi[iPer][iData][iWgt]->Sumw2 ();
      }
    }
  }


  for (int iPhoton = 0; iPhoton < nPhotons; iPhoton++) {
    inTree->GetEntry (iPhoton);

    const short iPer = isPeriodA ? 0 : 1;
    const short iMC = isMC ? 1 : 0;

    photonSpectrum[iPer][iMC][0]->Fill (ppt, peta, evtWeight); // period A or B, data or MC, weighted
    photonSpectrum[iPer][iMC][1]->Fill (ppt, peta); // etc.
    photonSpectrum[2][iMC][0]->Fill (ppt, peta, evtWeight);
    photonSpectrum[2][iMC][1]->Fill (ppt, peta);

    photonEtaPhi[iPer][iMC][0]->Fill (peta, pphi, evtWeight); // period A or B, data or MC, weighted
    photonEtaPhi[iPer][iMC][1]->Fill (peta, pphi); // etc.
    photonEtaPhi[2][iMC][0]->Fill (peta, pphi, evtWeight);
    photonEtaPhi[2][iMC][1]->Fill (peta, pphi);
  }


  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  TCanvas* etaPhiCanvas = new TCanvas ("etaPhiCanvas", "", 800, 600);
  FormatTH2Canvas (etaPhiCanvas);

  for (short iPer = 0; iPer < 3; iPer++) {
    canvas->cd ();
    gPad->SetLogx ();
    gPad->SetLogy ();

    for (short iData = 1; iData >= 0; iData--) {
      for (short iEta = 0; iEta < 2; iEta++) {
        TH1D* thisHist = photonSpectrum[iPer][iData][0]->ProjectionX ("_py", iEta==0?1:153,iEta==0?137:240);
        thisHist->Scale (1.0/(thisHist->Integral () * (iEta==0?1.37:0.88)), "width");
        //thisHist->GetYaxis ()->SetRangeUser (1e-2, 1e8);

        TGraphAsymmErrors* thisGraph = make_graph (thisHist);
        deltaize (thisGraph, 1+0.005*(iEta-numetabins/2), true);

        thisGraph->GetYaxis ()->SetRangeUser (1e-8, 5e-1);
        thisGraph->SetMarkerStyle (iData == 0 ? kFullCircle : kOpenCircle);
        thisGraph->SetMarkerSize (0.5);
        thisGraph->SetLineColor (iEta == 0 ? kBlack : kBlue);
        thisGraph->SetMarkerColor (iEta == 0 ? kBlack : kBlue);

        thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
        thisGraph->GetYaxis ()->SetTitle ("1 / N^{#gamma} d^{2}N^{#gamma} / d#it{p}_{T} d#eta #left[GeV^{-1}#right]");

        if (iEta == 0 && iData == 1)
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("p");

        if (thisHist) { delete thisHist; thisHist = NULL; }
        if (thisGraph) { delete thisGraph; thisGraph = NULL; }
      }
    }

    myText (0.25, 0.48, kBlack, "Barrel, #left|#eta#right| < 1.37", 0.04);
    myText (0.25, 0.42, kBlue, "Endcaps, 1.52 < #left|#eta#right| < 2.40", 0.04);

    canvas->SaveAs (Form ("%s/%s/photonSpectrum.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

    etaPhiCanvas->cd ();
    gPad->SetLogx (false);
    gPad->SetLogy (false);
    gStyle->SetPalette (kRainBow);
    photonEtaPhi[iPer][0][0]->Draw ("colz");
    photonEtaPhi[iPer][0][0]->GetXaxis ()->SetTitle ("Photon #eta");
    photonEtaPhi[iPer][0][0]->GetYaxis ()->SetTitle ("Photon #phi");
    photonEtaPhi[iPer][0][0]->GetXaxis ()->SetTitleOffset (1.1);
    photonEtaPhi[iPer][0][0]->GetYaxis ()->SetTitleOffset (1.1);
    etaPhiCanvas->SaveAs (Form ("%s/%s/photonEtaPhi_data.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

    etaPhiCanvas->cd ();
    gPad->SetLogx (false);
    gPad->SetLogy (false);
    gStyle->SetPalette (kRainBow);
    photonEtaPhi[iPer][1][0]->Draw ("colz");
    photonEtaPhi[iPer][1][0]->GetXaxis ()->SetTitle ("Photon #eta");
    photonEtaPhi[iPer][1][0]->GetYaxis ()->SetTitle ("Photon #phi");
    photonEtaPhi[iPer][1][0]->GetXaxis ()->SetTitleOffset (1.1);
    photonEtaPhi[iPer][1][0]->GetYaxis ()->SetTitleOffset (1.1);
    etaPhiCanvas->SaveAs (Form ("%s/%s/photonEtaPhi_mc.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

  }


  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");

  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iData = 0; iData < 2; iData++) {
      for (short iWgt = 0; iWgt < 2; iWgt++) {
        photonSpectrum[iPer][iData][iWgt]->Write ();
        photonEtaPhi[iPer][iData][iWgt]->Write ();
      }
    }
  }

  Delete3DArray (photonSpectrum, 3, 2, 2);
  Delete3DArray (photonEtaPhi, 3, 2, 2);

  outFile->Close ();

  return;
}

} // end namespace
