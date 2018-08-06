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

  TH2D* electronContamination = new TH2D ("electronContamination", ";#it{p}_{T}^{e} #left[GeV#right];#eta;Misidentified photon count", numpebins, pebins, numeetabins, eetabins);
  electronContamination->Sumw2();
  TH2D* electronSpectrum = new TH2D ("electronSpectrum", ";#it{p}_{T}^{e} #left[GeV#right];#eta;#sigma_{Z#rightarrow ee} #times N_{e} / N_{evt} #left[mb#right]", numpebins, pebins, numeetabins, eetabins);
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

       electronContamination->Add((TH2D*)thisFile->Get(Form("electronContamination_dataSet%s", zeeJetSampleId.Data())));
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
  TPad* pad = new TPad("pad", "", 0, 0, 1, 1);

  pad->SetLogz(true);

  electronContamination->Draw("colz");
  electronContamination->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs(Form("%s/electronContamination.pdf", plotPath.Data()));

  electronSpectrum->Draw("colz");
  electronSpectrum->GetZaxis()->SetTitleOffset(1.3);
  canvas->SaveAs(Form("%s/electronSpectrum.pdf", plotPath.Data()));
  

  TFile* outFile = new TFile(TString(rootPath) + "electronContaminationStudy.root", "recreate");

  electronContamination->Write();
  if (electronContamination) delete electronContamination;
  electronSpectrum->Write();
  if (electronSpectrum) delete electronSpectrum;

  outFile->Write();
  if (outFile) delete outFile;
}

