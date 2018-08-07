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

  TH2D* fakePhotonSpectrum = new TH2D ("fakePhotonSpectrum", ";#it{p}_{T}^{e} #left[GeV#right];#eta;#sigma_{Z#rightarrow ee} #times N_{fake #gamma} / N_{evt} #left[#mub#right]", numpebins, pebins, numeetabins, eetabins);
  fakePhotonSpectrum->Sumw2();
  TH2D* fakePhotonCounts = new TH2D ("fakePhotonCounts", ";#it{p}_{T}^{e} #left[GeV#right];#eta;# fake photons", numpebins, pebins, numeetabins, eetabins);
  fakePhotonCounts->Sumw2();
  TH2D* allElectronCounts = new TH2D ("allElectronCounts", ";#it{p}_{T}^{e} #left[GeV#right];#eta;# truth electrons", numpebins, pebins, numeetabins, eetabins);
  allElectronCounts->Sumw2();
  TH2D* allGammaCounts = new TH2D ("allGammaCounts", ";#it{p}_{T}^{e} #left[GeV#right];#eta;# total photons", numpebins, pebins, numeetabins, eetabins);
  allGammaCounts->Sumw2();
  TH2D* electronSpectrum = new TH2D ("electronSpectrum", ";#it{p}_{T}^{e} #left[GeV#right];#eta;#sigma_{Z#rightarrow ee} #times N_{e} / N_{evt} #left[#mub#right]", numpebins, pebins, numeetabins, eetabins);
  electronSpectrum->Sumw2();
  

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


  /**** Canvas definitions ****/

  TStyle* myStyle = AtlasStyle();
  myStyle->SetPalette(55);
  myStyle->SetPadRightMargin(0.18);
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
  gPad->SetLogx(true);

  gPad->SetLogz(true);

  fakePhotonSpectrum->Draw("colz");
  fakePhotonSpectrum->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs (Form ("%s/fakePhotonSpectrum.pdf", plotPath.Data()));

  gPad->SetLogz(true); 
  fakePhotonCounts->Draw("colz");
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
  fakePhotonRate->Draw("colz");
  fakePhotonRate->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs(Form ("%s/fakePhotonRate.pdf", plotPath.Data()));

  gPad->SetLogz(true);
  electronSpectrum->Draw("colz");
  electronSpectrum->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs(Form("%s/electronSpectrum.pdf", plotPath.Data()));
  

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

  outFile->Write();
  if (outFile) delete outFile;
}

