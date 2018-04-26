#include "../Params.C"
#include "../../Initialization.C"

void ZGammaJetCrossCheckHist () {

  // Setup trigger vectors
  setupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  vector<int> dataSets(0);
  for (int i = 0; i < sizeof(full_run_list)/sizeof(full_run_list[0]); i++) dataSets.push_back(full_run_list[i]);

  const int nhists = 6;
  TH2F* zjethists[nhists];
  for (int nhist = 0; nhist < nhists; nhist++) {
    zjethists[nhist] = new TH2F(Form("Zjet_hist%i", nhist), "", numpbins, pbins, 40, 0, 1.6);
    zjethists[nhist]->Sumw2();
  }

  {
    TSystemDirectory dir(rootPath.c_str(), rootPath.c_str());
    TList* sysfiles = dir.GetListOfFiles();
    if (sysfiles) {
      TSystemFile *sysfile;
      TString fname;
      TString histName;
      TIter next(sysfiles);
      TVectorD* infoVec;
      int numFiles = 0;
      while ((sysfile=(TSystemFile*)next())) {
        fname = sysfile->GetName();
        if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
          if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl;
          for (int dataSet : dataSets) {
            if (fname.Contains(to_string(dataSet))) {
              numFiles++;
              TFile* thisFile = new TFile(rootPath + fname, "READ");
              infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", dataSet));
              for (int nhist = 0; nhist < nhists; nhist++) {
                zjethists[nhist]->Add((TH2F*)thisFile->Get(Form("Zjet_run%i_hist%i", dataSet, nhist)));
              }
              thisFile->Close();
              delete thisFile;
              break;
            }
          }
        }
      }
      cout << numFiles << " files read in." << endl;
    }
  }  
 
  TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
  gPad->SetLogx();
  for (int nhist = 0; nhist < nhists; nhist++) {
    TProfile* zjethist = zjethists[nhist]->ProfileX();
    zjethist->SetAxisRange(0, 1.6, "Y");
    zjethist->Draw("");
    canvas->SaveAs(Form("%s/zjet%i.pdf", plotPath.c_str(), nhist));
  }
}
