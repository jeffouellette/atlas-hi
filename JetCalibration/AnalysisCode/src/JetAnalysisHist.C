#include "JetAnalysisHist.h"
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

void JetAnalysisHist () {

  SetAtlasStyle ();

  SetupDirectories ("JetAnalysis/", "JetCalibration/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* inTree = (TTree*)inFile->Get ("jeffsjets");
  float jpt = 0, jeta = 0, jphi = 0, je = 0, ppt = 0, peta = 0, pphi = 0;
  double evtWeight = 0;
  bool isPeriodA = false, isMC = false;
  inTree->SetBranchAddress ("jet_pt", &jpt);
  inTree->SetBranchAddress ("jet_eta", &jeta);
  inTree->SetBranchAddress ("jet_phi", &jphi);
  inTree->SetBranchAddress ("jet_e", &je);
  //inTree->SetBranchAddress ("photon_pt", &ppt);
  //inTree->SetBranchAddress ("photon_eta", &peta);
  //inTree->SetBranchAddress ("photon_phi", &pphi);
  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);

  const long long nJets = inTree->GetEntries ();

  TH2D**** jetSpectrum = Get3DArray <TH2D*> (3, 2, 2); // iPer, iData, iWgt
  TH2D**** jetEtaPhi = Get3DArray <TH2D*> (3, 2, 2); // iPer, iData, iWgt

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* per = iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB");

    for (short iData = 0; iData < 2; iData++) {
      const char* data = iData == 0 ? "data" : "mc";

      for (short iWgt = 0; iWgt < 2; iWgt++) {
        const char* wgt = iWgt == 0 ? "weighted" : "unweighted";

        jetSpectrum[iPer][iData][iWgt] = new TH2D (Form ("jetSpectrum_%s_%s_%s", per, data, wgt), "", 40, logspace (20, 500, 40), numetabins, etabins);
        jetSpectrum[iPer][iData][iWgt]->Sumw2 ();

        jetEtaPhi[iPer][iData][iWgt] = new TH2D (Form ("jetEtaPhi_%s_%s_%s", per, data, wgt), "", 450, -4.5, 4.5, numphibins, phibins);
        jetEtaPhi[iPer][iData][iWgt]->Sumw2 ();
      }
    }
  }


  for (long long iJet = 0; iJet < nJets; iJet++) {
    inTree->GetEntry (iJet);

    //if (ppt < 60)
    //  continue;

    jetSpectrum[(short)isPeriodA][(short)isMC][0]->Fill (jpt, jeta, evtWeight); // period A or B, data or MC, weighted
    jetSpectrum[(short)isPeriodA][(short)isMC][1]->Fill (jpt, jeta); // etc.
    jetSpectrum[2][(short)isMC][0]->Fill (jpt, jeta, evtWeight);
    jetSpectrum[2][(short)isMC][1]->Fill (jpt, jeta);

    jetEtaPhi[(short)isPeriodA][(short)isMC][0]->Fill (jeta, jphi, evtWeight); // period A or B, data or MC, weighted
    jetEtaPhi[(short)isPeriodA][(short)isMC][1]->Fill (jeta, jphi); // etc.
    jetEtaPhi[2][(short)isMC][0]->Fill (jeta, jphi, evtWeight);
    jetEtaPhi[2][(short)isMC][1]->Fill (jeta, jphi);
  }


  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  TCanvas* etaPhiCanvas = new TCanvas ("etaPhiCanvas", "", 800, 600);
  FormatTH2Canvas (etaPhiCanvas);
  const Color_t colors[9] = {kMagenta, kRed, kBlue, 8, kOrange+1, kCyan+2};
  const float scales[numetabins] = {-3, -1, 1, 2, 0, -2};

  for (short iPer = 0; iPer < 3; iPer++) {
    canvas->cd ();
    gPad->SetLogx ();
    gPad->SetLogy ();

    for (short iData = 1; iData >= 0; iData--) {
      for (short iEta = 0; iEta < numetabins; iEta++) {
        TH1D* thisHist = jetSpectrum[iPer][iData][0]->ProjectionX ("_py", iEta+1, iEta+1);
        thisHist->Scale (pow (10, scales[iEta])/(thisHist->Integral () * (etabins[iEta+1]-etabins[iEta])), "width");

        TGraphAsymmErrors* thisGraph = make_graph (thisHist);
        deltaize (thisGraph, 1+0.005*(iEta-numetabins/2), true);

        thisGraph->GetXaxis ()->SetRangeUser (20, 500);
        thisGraph->GetYaxis ()->SetRangeUser (1e-11, 1e1);
        thisGraph->SetMarkerStyle (iData == 0 ? kFullCircle : kOpenCircle);
        thisGraph->SetMarkerSize (0.5);

        thisGraph->SetLineColor (colors[iEta]);
        thisGraph->SetMarkerColor (colors[iEta]);

        thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet} #left[GeV#right]");
        thisGraph->GetYaxis ()->SetTitle ("1 / N^{jet} d^{2}N^{jet} / d#it{p}_{T} d#eta #left[GeV^{-1}#right]");

        if (iData == 1 && iEta == 0)
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("p");

        myText (0.20, 0.2+0.06*iEta, colors[iEta], Form ("%g < #eta < %g (x 10^{%g})", etabins[iEta], etabins[iEta+1], scales[iEta]), 0.04);

        if (thisHist) { delete thisHist; thisHist = NULL; }
        if (thisGraph) { delete thisGraph; thisGraph = NULL; }
      }
    }

    canvas->SaveAs (Form ("%s/%s/jetSpectrum.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

    etaPhiCanvas->cd ();
    gPad->SetLogx (false);
    gPad->SetLogy (false);
    gStyle->SetPalette (kRainBow);
    jetEtaPhi[iPer][0][0]->GetXaxis ()->SetTitle ("Jet #eta");
    jetEtaPhi[iPer][0][0]->GetYaxis ()->SetTitle ("Jet #phi");
    jetEtaPhi[iPer][0][0]->GetXaxis ()->SetTitleOffset (1.1);
    jetEtaPhi[iPer][0][0]->GetYaxis ()->SetTitleOffset (1.1);
    jetEtaPhi[iPer][0][0]->Draw ("colz");
    etaPhiCanvas->SaveAs (Form ("%s/%s/jetEtaPhi_data.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

    etaPhiCanvas->cd ();
    gPad->SetLogx (false);
    gPad->SetLogy (false);
    gStyle->SetPalette (kRainBow);
    jetEtaPhi[iPer][1][0]->GetXaxis ()->SetTitle ("Jet #eta");
    jetEtaPhi[iPer][1][0]->GetYaxis ()->SetTitle ("Jet #phi");
    jetEtaPhi[iPer][1][0]->GetXaxis ()->SetTitleOffset (1.1);
    jetEtaPhi[iPer][1][0]->GetYaxis ()->SetTitleOffset (1.1);
    jetEtaPhi[iPer][1][0]->Draw ("colz");
    etaPhiCanvas->SaveAs (Form ("%s/%s/jetEtaPhi_mc.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

  }


  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");

  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iData = 0; iData < 2; iData++) {
      for (short iWgt = 0; iWgt < 2; iWgt++) {
        jetSpectrum[iPer][iData][iWgt]->Write ();
        jetEtaPhi[iPer][iData][iWgt]->Write ();
      }
    }
  }

  Delete3DArray (jetSpectrum, 3, 2, 2);
  Delete3DArray (jetEtaPhi, 3, 2, 2);

  outFile->Close ();

  return;
}

} // end namespace
