#include "ElectronContaminationStudyHist.h"

#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace pPb8TeV2016JetCalibration {

void ElectronContaminationStudyHist () {

  // Setup trigger vectors
  SetupDirectories("ElectronContamination/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<TString> zeeJetSampleIds(0);
  zeeJetSampleIds.push_back("Pbp_ZeeJet_Overlay");
  zeeJetSampleIds.push_back("pPb_ZeeJet_Overlay");

  TH2D* fakePhotonSpectrum = new TH2D ("fakePhotonSpectrum", ";#it{p}_{T}^{e} #left[GeV#right];#eta;#sigma_{Z#rightarrow ee} #times #SigmaN_{fake #gamma} / #SigmaN_{evt} #left[#mub#right]", numpebins, pebins, numeetabins, eetabins);
  fakePhotonSpectrum->Sumw2();
  TH2D* fakePhotonCounts = new TH2D ("fakePhotonCounts", ";#it{p}_{T}^{e} #left[GeV#right];#eta;# fake photons", numpebins, pebins, numeetabins, eetabins);
  fakePhotonCounts->Sumw2();
  TH2D* allElectronCounts = new TH2D ("allElectronCounts", ";#it{p}_{T}^{e} #left[GeV#right];#eta;# truth electrons", numpebins, pebins, numeetabins, eetabins);
  allElectronCounts->Sumw2();
  TH2D* allGammaCounts = new TH2D ("allGammaCounts", ";#it{p}_{T}^{e} #left[GeV#right];#eta;# total photons", numpebins, pebins, numeetabins, eetabins);
  allGammaCounts->Sumw2();
  TH2D* electronSpectrum = new TH2D ("electronSpectrum", ";#it{p}_{T}^{e} #left[GeV#right];#eta;#sigma_{Z#rightarrow ee} #times #SigmaN_{e} / #SigmaN_{evt} #left[#mub#right]", numpebins, pebins, numeetabins, eetabins);
  electronSpectrum->Sumw2();
  TH2D* truthElectronRecoPhotonCounts = new TH2D ("truthElectronRecoPhotonCounts", "", numpebins, pebins, numeetabins, eetabins);
  truthElectronRecoPhotonCounts->Sumw2();

  TH2D* truthElectronRecoElectronCounts = new TH2D ("truthElectronRecoElectronCounts", "", numpebins, pebins, numeetabins, eetabins);
  truthElectronRecoElectronCounts->Sumw2();
  

  {
   TSystemDirectory dir(rootPath.Data(), rootPath.Data());
   TList* sysfiles = dir.GetListOfFiles();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile *sysfile;
   TString fname;
   TString histName;
   TIter next(sysfiles);
   int numFiles = 0;
   while ((sysfile=(TSystemFile*)next())) {
    fname = sysfile->GetName();
    if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {

     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains(zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");

       fakePhotonSpectrum->Add((TH2D*)thisFile->Get(Form("fakePhotonSpectrum_dataSet%s", zeeJetSampleId.Data())));
       fakePhotonCounts->Add((TH2D*)thisFile->Get(Form("fakePhotonCounts_dataSet%s", zeeJetSampleId.Data())));
       allElectronCounts->Add((TH2D*)thisFile->Get(Form("allElectronCounts_dataSet%s", zeeJetSampleId.Data())));
       allGammaCounts->Add((TH2D*)thisFile->Get(Form("allGammaCounts_dataSet%s", zeeJetSampleId.Data())));
       electronSpectrum->Add((TH2D*)thisFile->Get(Form("electronSpectrum_dataSet%s", zeeJetSampleId.Data())));
       truthElectronRecoPhotonCounts->Add((TH2D*)thisFile->Get(Form("truthElectronRecoPhotonCounts_dataSet%s", zeeJetSampleId.Data())));
       truthElectronRecoElectronCounts->Add((TH2D*)thisFile->Get(Form("truthElectronRecoElectronCounts_dataSet%s", zeeJetSampleId.Data())));

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
    }
   }
   cout << numFiles << " files read in." << endl;
  }
  /**** End loop over input files ****/


  /**** Divide by number of events ****/
  fakePhotonSpectrum->Scale(1.0/(302313 + 167761));
  electronSpectrum->Scale(1.0/(302313 + 167761));
  truthElectronRecoPhotonCounts->Scale(1.0/(302313 + 167761));
  truthElectronRecoElectronCounts->Scale(1.0/(302313 + 167761));


  /**** Canvas definitions ****/

  TStyle* myStyle = AtlasStyle();
  myStyle->SetPalette(55);
  myStyle->SetPadRightMargin(0.18);
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
  FormatTH2Canvas (canvas);

  gPad->SetLogx(true);
  gPad->SetLogz(true);

  fakePhotonSpectrum->GetYaxis()->SetTitleOffset(1.1);
  fakePhotonSpectrum->Draw("colz");
  fakePhotonSpectrum->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs (Form ("%s/fakePhotonSpectrum.pdf", plotPath.Data()));

  gPad->SetLogz(true); 
  fakePhotonCounts->Draw("colz");
  fakePhotonCounts->GetYaxis()->SetTitleOffset(1.1);
  fakePhotonCounts->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs (Form ("%s/fakePhotonCounts.pdf", plotPath.Data()));

  gPad->SetLogz(false);
  TH2D* fakePhotonRate = new TH2D("fakePhotonRate", ";#it{p}_{T}^{#gamma} #left[GeV#right];#eta^{#gamma};# fake photons / # truth electron", numpebins, pebins, numeetabins, eetabins);
  for (short iP = 1; iP <= numpebins; iP++) {
   for (short iEta = 1; iEta <= numeetabins; iEta++) {
    if (allElectronCounts->GetBinContent (iP, iEta) != 0) {
     fakePhotonRate->SetBinContent (iP, iEta, fakePhotonCounts->GetBinContent (iP, iEta) / allElectronCounts->GetBinContent (iP, iEta));
    }
   }
  }
  fakePhotonRate->GetYaxis()->SetTitleOffset(1.1);
  fakePhotonRate->GetZaxis()->SetTitleOffset(1.3);
  fakePhotonRate->Draw("colz");
  fakePhotonRate->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs(Form ("%s/fakePhotonRate.pdf", plotPath.Data()));

  gPad->SetLogz(true);
  electronSpectrum->GetYaxis()->SetTitleOffset(1.1);
  electronSpectrum->GetZaxis()->SetTitleOffset(1.3);
  electronSpectrum->Draw("colz");
  electronSpectrum->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs(Form("%s/electronSpectrum.pdf", plotPath.Data()));


  TCanvas* truthElectronCanvas = new TCanvas ("truthElectronCanvas", "", 800, 600);
  truthElectronCanvas->Divide (3, 2, -1);  
//  FormatTH2Canvas (truthElectronCanvas, false);

  truthElectronRecoElectronCounts->GetXaxis()->SetTitle("Truth #it{E}_{T} #left[GeV#right]");
  truthElectronRecoElectronCounts->GetYaxis()->SetTitle("Reco. counts / #SigmaN_{evt}");
  truthElectronRecoPhotonCounts->GetXaxis()->SetTitle("Truth #it{E}_{T} #left[GeV#right]");
  truthElectronRecoPhotonCounts->GetYaxis()->SetTitle("Reco. counts / #SigmaN_{evt}");

  int padnum = 1;
  for (int iEta = 0; iEta < numeetabins; iEta++) {
   if (iEta == 1 || iEta == 3) continue;

   truthElectronCanvas->cd(padnum);
   gPad->SetLogy(true);
   gPad->SetLogx(true);

   TH1D* recoElectronSlice = truthElectronRecoElectronCounts->ProjectionX(Form("electron_%i",iEta), iEta, iEta+1);
   TH1D* recoPhotonSlice = truthElectronRecoPhotonCounts->ProjectionX(Form("photon_%i",iEta), iEta, iEta+1);

   recoElectronSlice->GetYaxis()->SetTitle("Counts / Event / GeV");
   recoElectronSlice->GetYaxis()->SetTitleOffset(1.1);
   recoPhotonSlice->GetYaxis()->SetTitle("Counts / Event / GeV");
   recoPhotonSlice->GetYaxis()->SetTitleOffset(1.1);

   recoElectronSlice->Scale(1, "width");
   recoPhotonSlice->Scale(1, "width");

   recoElectronSlice->SetAxisRange(1e-7, 1, "Y");
   recoPhotonSlice->SetAxisRange(1e-7, 1, "Y");

   recoElectronSlice->SetLineColor(kBlack);
   recoElectronSlice->SetMarkerColor(kBlack);

   recoPhotonSlice->SetLineColor(kRed);
   recoPhotonSlice->SetMarkerColor(kRed);

   recoElectronSlice->Draw("e1");
   recoPhotonSlice->Draw("e1 same");

   myText (0.55, 0.82, kBlack, Form ("%g < #eta < %g", eetabins[iEta], eetabins[iEta+1]));
   //canvas->SaveAs(Form("%s/truthElectronRecoCounts_etabin%i.pdf", plotPath.Data(), iEta));

   truthElectronCanvas->cd(padnum+3);
   gPad->SetLogy(false);
   gPad->SetLogx(true);

   TH1D* ratio = truthElectronRecoPhotonCounts->ProjectionX(Form("ratio_%i", iEta), iEta, iEta+1);
   ratio->Scale(1, "width"); // necessary since denominator has width divided
   ratio->Divide(recoElectronSlice);
   ratio->GetYaxis()->SetTitle("N^{e#rightarrow#gamma} / N^{e#rightarrowe}");
   ratio->GetYaxis()->SetTitleOffset(1.1);
   ratio->SetAxisRange(0, 0.12, "Y");

   ratio->Draw("e1");
   myText (0.55, 0.82, kBlack, Form ("%g < #eta < %g", eetabins[iEta], eetabins[iEta+1]));
   //canvas->SaveAs(Form("%s/electronMisIdRate_etabin%i.pdf", plotPath.Data(), iEta));

   padnum++;

  }
  truthElectronCanvas->SaveAs(Form("%s/truthElectronReco.pdf", plotPath.Data()));
  

  TFile* outFile = new TFile(TString(rootPath) + "electronContaminationStudy.root", "recreate");

  fakePhotonCounts->Write();
  if (fakePhotonCounts) delete fakePhotonCounts;
  fakePhotonSpectrum->Write();
  if (fakePhotonSpectrum) delete fakePhotonSpectrum;
  allElectronCounts->Write();
  if (allElectronCounts) delete allElectronCounts;
  allGammaCounts->Write();
  if (allGammaCounts) delete allGammaCounts;
  fakePhotonRate->Write();
  if (fakePhotonRate) delete fakePhotonRate;
  electronSpectrum->Write();
  if (electronSpectrum) delete electronSpectrum;
  truthElectronRecoElectronCounts->Write();
  if (truthElectronRecoElectronCounts) delete truthElectronRecoElectronCounts;
  truthElectronRecoPhotonCounts->Write();
  if (truthElectronRecoPhotonCounts) delete truthElectronRecoPhotonCounts;
  truthElectronCanvas->Write();
  if (truthElectronCanvas) delete truthElectronCanvas;

  outFile->Write();
  if (outFile) delete outFile;
}

} // end namespace
