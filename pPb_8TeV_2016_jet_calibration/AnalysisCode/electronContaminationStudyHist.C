#include "../Params.C"
#include "../../Initialization.C"

static const double pebins[18] = {20, 25, 35, 45, 55, 65, 75, 85, 105, 125, 150, 175, 200, 250, 300, 350, 400, 550};
static const short numpebins = sizeof(pebins)/sizeof(pebins[0]) - 1;

static const double eetabins[6] = { -2.37, -1.56, -1.37, 1.37, 1.56, 2.37};
static const short numeetabins = sizeof(eetabins)/sizeof(eetabins[0]) - 1;


void electronContaminationStudyHist () {

  // Setup trigger vectors
  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

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
   const TString contpath = rootPath + "/electronContamination/";
   TSystemDirectory dir(contpath.Data(), contpath.Data());
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
       cout << "Reading in " << contpath+fname << endl;
       TFile* thisFile = new TFile(contpath + fname, "READ");

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
  for (int pebin = 1; pebin <= numpebins; pebin++) {
   for (int eetabin = 1; eetabin <= numeetabins; eetabin++) {
    if (allElectronCounts->GetBinContent (pebin, eetabin) != 0) {
     fakePhotonRate->SetBinContent (pebin, eetabin, fakePhotonCounts->GetBinContent (pebin, eetabin) / allElectronCounts->GetBinContent (pebin, eetabin));
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

  
  FormatTH2Canvas (canvas, false);
  gPad->SetLogy(true);

  truthElectronRecoElectronCounts->GetXaxis()->SetTitle("Truth #it{E}_{T} #left[GeV#right]");
  truthElectronRecoElectronCounts->GetYaxis()->SetTitle("Reco. counts / #SigmaN_{evt}");
  truthElectronRecoPhotonCounts->GetXaxis()->SetTitle("Truth #it{E}_{T} #left[GeV#right]");
  truthElectronRecoPhotonCounts->GetYaxis()->SetTitle("Reco. counts / #SigmaN_{evt}");

  for (int eetabin = 0; eetabin < numeetabins; eetabin++) {
   if (eetabin == 1 || eetabin == 3) continue;
   TH1D* recoElectronSlice = truthElectronRecoElectronCounts->ProjectionX(Form("electron_%i",eetabin), eetabin, eetabin+1);
   TH1D* recoPhotonSlice = truthElectronRecoPhotonCounts->ProjectionX(Form("photon_%i",eetabin), eetabin, eetabin+1);

   recoElectronSlice->GetYaxis()->SetTitle("Counts / #SigmaN_{evt} / GeV");
   recoElectronSlice->GetYaxis()->SetTitleOffset(1.1);
   recoPhotonSlice->GetYaxis()->SetTitle("Counts / #SigmaN_{evt} / GeV");
   recoPhotonSlice->GetYaxis()->SetTitleOffset(1.1);

   recoElectronSlice->Scale(1, "width");
   recoPhotonSlice->Scale(1, "width");

   recoElectronSlice->SetLineColor(kBlack);
   recoElectronSlice->SetMarkerColor(kBlack);

   recoPhotonSlice->SetLineColor(kRed);
   recoPhotonSlice->SetMarkerColor(kRed);

   recoElectronSlice->Draw("e1");
   recoPhotonSlice->Draw("e1 same");

   myText (0.55, 0.82, kBlack, Form ("%g < #eta < %g", eetabins[eetabin], eetabins[eetabin+1]));
   canvas->SaveAs(Form("%s/truthElectronRecoCounts_etabin%i.pdf", plotPath.Data(), eetabin));
  }
  

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

  outFile->Write();
  if (outFile) delete outFile;
}

