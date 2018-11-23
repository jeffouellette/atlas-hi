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

  TTree* photonTree = (TTree*)inFile->Get ("jeffsphotons");
  float ppt = 0, peta = 0, pphi = 0, evtWeight = 0;
  bool isPeriodA = false, isMC = false;
  photonTree->SetBranchAddress ("photon_pt", &ppt);
  photonTree->SetBranchAddress ("photon_eta", &peta);
  photonTree->SetBranchAddress ("photon_phi", &pphi);
  photonTree->SetBranchAddress ("evt_weight", &evtWeight);
  photonTree->SetBranchAddress ("isMC", &isMC);
  photonTree->SetBranchAddress ("isPeriodA", &isPeriodA);

  const int nPhotons = photonTree->GetEntries ();

  TH2D**** photonSpectrum = Get3DArray <TH2D*> (3, 2, 2); // iPer, iData, iWgt
  TH2D**** photonEtaPhi = Get3DArray <TH2D*> (3, 2, 2); // iPer, iData, iWgt

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* per = iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB");

    for (short iData = 0; iData < 2; iData++) {
      const char* data = iData == 0 ? "data" : "mc";

      for (short iWgt = 0; iWgt < 2; iWgt++) {
        const char* wgt = iWgt == 0 ? "weighted" : "unweighted";

        photonSpectrum[iPer][iData][iWgt] = new TH2D (Form ("photonSpectrum_%s_%s_%s", per, data, wgt), "", numpbins, pbins, 240, 0, 2.40);
        photonSpectrum[iPer][iData][iWgt]->Sumw2 ();

        photonEtaPhi[iPer][iData][iWgt] = new TH2D (Form ("photonEtaPhi_%s_%s_%s", per, data, wgt), "", 240, -2.4, 2.4, numphibins, phibins);
        photonEtaPhi[iPer][iData][iWgt]->Sumw2 ();
      }
    }
  }


  for (int iPhoton = 0; iPhoton < nPhotons; iPhoton++) {
    photonTree->GetEntry (iPhoton);

    photonSpectrum[(short)isPeriodA][(short)isMC][0]->Fill (ppt, peta, evtWeight); // period A or B, data or MC, weighted
    photonSpectrum[(short)isPeriodA][(short)isMC][1]->Fill (ppt, peta); // etc.
    photonSpectrum[2][(short)isMC][0]->Fill (ppt, peta, evtWeight);
    photonSpectrum[2][(short)isMC][1]->Fill (ppt, peta);

    photonEtaPhi[(short)isPeriodA][(short)isMC][0]->Fill (peta, pphi, evtWeight); // period A or B, data or MC, weighted
    photonEtaPhi[(short)isPeriodA][(short)isMC][1]->Fill (peta, pphi); // etc.
    photonEtaPhi[2][(short)isMC][0]->Fill (peta, pphi, evtWeight);
    photonEtaPhi[2][(short)isMC][1]->Fill (peta, pphi);
  }


  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);

  for (short iPer = 0; iPer < 3; iPer++) {
    gPad->SetLogx ();
    gPad->SetLogy ();

    for (short iEta = 0; iEta < 2; iEta++) {
      TH1D* thisHist = photonSpectrum[iPer][0][0]->ProjectionX ("_py", iEta==0?1:153,iEta==0?137:240);
      thisHist->Scale (1.0/(iEta==0?1.37:0.88), "width");

      TGraphAsymmErrors* thisGraph = make_graph (thisHist);
      deltaize (thisGraph, 1+0.002*(iEta-numetabins/2), true);

      thisGraph->SetLineColor (iEta==0? kBlack:kBlue);
      thisGraph->SetMarkerColor (iEta==0? kBlack:kBlue);

      thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle ("Counts");

      ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw (iEta == 0 ? "ap":"p");

      if (thisHist) { delete thisHist; thisHist = NULL; }
      if (thisGraph) { delete thisGraph; thisGraph = NULL; }
    }

    myText (0.65, 0.88, kBlack, "Barrel, #left|#eta#right| < 1.37", 0.04);
    myText (0.65, 0.81, kBlue, "Endcaps, 1.52 < #left|#eta#right| < 2.40", 0.04);

    canvas->SaveAs (Form ("%s/%s/photonSpectrum.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

    gPad->SetLogx (false);
    gPad->SetLogy (false);
    gStyle->SetPalette (kRainBow);
    photonEtaPhi[iPer][0][0]->Draw ("colz");
    canvas->SaveAs (Form ("%s/%s/photonEtaPhi_data.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

    gPad->SetLogx (false);
    gPad->SetLogy (false);
    gStyle->SetPalette (kRainBow);
    photonEtaPhi[iPer][1][0]->Draw ("colz");
    canvas->SaveAs (Form ("%s/%s/photonEtaPhi_mc.pdf", plotPath.Data (), iPer == 0 ? "PeriodA" : (iPer == 1 ? "PeriodB" : "PeriodAB")));

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
