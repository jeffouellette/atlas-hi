#include "ZGammaJetCrossCheckHist.C"

void GammaMCValidOverlayComp () {

  setupDirectories("", "pPb_8TeV_2016_jet_calibration/");
  vector<TString> gammaJetOverlaySampleIds(0);
  vector<TString> gammaJetValidSampleIds(0);
  for (short i = 0; i < 6; i++) {
    gammaJetValidSampleIds.push_back(TString("Pbp_Valid_GammaJet_Slice") + to_string(i+1));
    gammaJetValidSampleIds.push_back(TString("pPb_Valid_GammaJet_Slice") + to_string(i+1));
    gammaJetOverlaySampleIds.push_back(TString("Pbp_Overlay_GammaJet_Slice") + to_string(i+1));
    gammaJetOverlaySampleIds.push_back(TString("pPb_Overlay_GammaJet_Slice") + to_string(i+1));
  }
  for (short i = 0; i < 12; i++) {
    cout << gammaJetValidSampleIds[i] << endl;
    cout << gammaJetOverlaySampleIds[i] << endl;
  }

  int g_n[3][2][numetabins+1] = {{}, {}};
  TH2F* gJetHists[3][numetabins][2];
  for (short etabin = 0; etabin < numetabins; etabin++) {

    for (short periodType = 0; periodType < 3; periodType++) {
      string period = "periodA";
      if (periodType == 1) period = "periodB";
      else if (periodType == 2) period = "periodAB";

      for (short dType = 0; dType < 2; dType++) { // overlay = 1, valid = 0
        gJetHists[periodType][etabin][dType] = new TH2F(Form("gJetPtRatio_hist%i_%s_%s", etabin, period.c_str(), (dType==0 ? "valid":"overlay")), ";#it{p}_{T}^{#gamma} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{#gamma}#cos#Delta#phi", numpbins, pbins, numxjrefbins, xjrefbins);
        gJetHists[periodType][etabin][dType]->Sumw2();
      }
    }
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
   
          for (TString gammaJetValidSampleId : gammaJetValidSampleIds) { // check for gamma jet MC
            if (fname.Contains(gammaJetValidSampleId)) {
              numFiles++;
              cout << "Reading in " << rootPath+fname << endl;
              TFile* thisFile = new TFile(rootPath + fname, "READ");
              infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", gammaJetValidSampleId.Data()));

              for (short etabin = 0; etabin < numetabins; etabin++) {
                const short periodType = (gammaJetValidSampleId.Contains("pPb") ? 0 : 1);
                const short act_etabin = (gammaJetValidSampleId.Contains("pPb") ? numetabins - etabin - 1 : etabin); // period A condition

                g_n[periodType][0][etabin] += (*infoVec)[2+2*numetabins+etabin];
                g_n[2][0][act_etabin] += (*infoVec)[2+2*numetabins+etabin];

                // Only add the statistical error plots for MC (don't need to consider systematics)
                gJetHists[periodType][etabin][0]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetValidSampleId.Data(), etabin)));
                gJetHists[2][etabin][0]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetValidSampleId.Data(), act_etabin)));
              }

              thisFile->Close();
              if (thisFile) delete thisFile;
              break;
            }
          }

          for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
            if (fname.Contains(gammaJetOverlaySampleId)) {
              numFiles++;
              cout << "Reading in " << rootPath+fname << endl;
              TFile* thisFile = new TFile(rootPath + fname, "READ");
              infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", gammaJetOverlaySampleId.Data()));

              for (short etabin = 0; etabin < numetabins; etabin++) {
                const short periodType = (gammaJetOverlaySampleId.Contains("pPb") ? 0 : 1);
                const short act_etabin = (gammaJetOverlaySampleId.Contains("pPb") ? numetabins - etabin - 1 : etabin); // period A condition

                g_n[periodType][1][etabin] += (*infoVec)[2+2*numetabins+etabin]; // don't flip for period A & period B separately
                g_n[2][1][act_etabin] += (*infoVec)[2+2*numetabins+etabin]; // but do flip when combining period A & B

                // Only add the statistical error plots for MC (don't need to consider systematics)
                gJetHists[periodType][etabin][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetOverlaySampleId.Data(), etabin)));
                gJetHists[2][etabin][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetOverlaySampleId.Data(), act_etabin)));
              }

              thisFile->Close();
              if (thisFile) delete thisFile;
              break;
            }
          }
        }
      }
      cout << numFiles << " files read in." << endl;
    }
  }

  for (short periodType = 0; periodType < 3; periodType++) {
    for (short dType = 0; dType < 2; dType++) {
      g_n[periodType][dType][numetabins] = 0;
      for (short etabin = 0; etabin < numetabins; etabin++) {
        g_n[periodType][dType][numetabins] += g_n[periodType][dType][etabin];
      }
    }
  }

  TLine* lines[8] = {};
  for (short i = 0; i < 8; i++) {
    lines[i] = new TLine(pbins[0], 0.95+0.025*i, pbins[numpbins], 0.95+0.025*i);
    if (0.95+0.025*i == 1) lines[i]->SetLineStyle(1);
    else lines[i]->SetLineStyle(3);
  }

  TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
  const double padRatio = 1.5; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double lPadY = 1.0/(padRatio+1.0);
  const double uPadY = 1.0 - lPadY;
  TPad* upperPad = new TPad("upperPad", "", 0, lPadY, 1, 1);
  TPad* lowerPad = new TPad("lowerPad", "", 0, 0, 1, lPadY);
  upperPad->SetBottomMargin(0);
  upperPad->SetLeftMargin(-0.20);
  lowerPad->SetTopMargin(0);
  lowerPad->SetBottomMargin(0.30);
  lowerPad->SetLeftMargin(-0.20);
  upperPad->Draw();
  lowerPad->Draw();

  TProfile *gJetValidHist, *gJetOverlayHist;
  TH1D *gJetHist_rat;
  for (short periodType = 0; periodType < 3; periodType++) {
    string period = "Period A";
    if (periodType == 1) period = "Period B";
    else if (periodType == 2) period = "Period A+B";
    for (short etabin = 0; etabin < numetabins; etabin++) {

      /**** Plots GammaJet info ****/
      upperPad->cd();
      gPad->SetLogx();
      gJetValidHist = gJetHists[periodType][etabin][0]->ProfileX();
      gJetValidHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
      gJetValidHist->SetAxisRange(0.65, 1.45, "Y");
      gJetValidHist->SetMarkerColor(data_color);
      gJetValidHist->SetLineColor(data_color);
      gJetValidHist->GetXaxis()->SetLabelSize(0.04/uPadY);
      gJetValidHist->GetYaxis()->SetLabelSize(0.04/uPadY);
      gJetValidHist->GetYaxis()->SetTitleSize(0.04/uPadY);
      gJetValidHist->GetYaxis()->SetTitleOffset(uPadY);

      gJetOverlayHist = gJetHists[periodType][etabin][1]->ProfileX();
      gJetOverlayHist->SetMarkerColor(mc_color);
      gJetOverlayHist->SetLineColor(mc_color);

      gJetValidHist->DrawCopy("E1 X0");
      gJetOverlayHist->DrawCopy("SAME E1 X0");
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("Validation Sample (Signal Only) (%i events)", g_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Data Overlay Sample (%i events)", g_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
      if (periodType == 0) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{A} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
      else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{B} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
      myText(0.155, 0.28, kBlack, "#gamma + Jet", 0.04/uPadY);
      myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

      lowerPad->cd();
      gPad->SetLogx();
      gJetHist_rat = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hist%i", etabin)), gJetHists[periodType][etabin][0], gJetHists[periodType][etabin][1], numpbins, pbins, numxjrefbins, xjrefbins);

      gJetHist_rat->SetYTitle("Valid. / Overlay");
      gJetHist_rat->SetAxisRange(0.95, 1.05, "Y");
      gJetHist_rat->GetYaxis()->SetNdivisions(405);
      gJetHist_rat->GetXaxis()->SetTitleSize(0.04/lPadY);
      gJetHist_rat->GetYaxis()->SetTitleSize(0.04/lPadY);
      gJetHist_rat->GetXaxis()->SetTitleOffset(1);
      gJetHist_rat->GetYaxis()->SetTitleOffset(lPadY);
      gJetHist_rat->GetXaxis()->SetLabelSize(0.04/lPadY);
      gJetHist_rat->GetYaxis()->SetLabelSize(0.04/lPadY);
      gJetHist_rat->GetXaxis()->SetTickLength(0.08);

      gJetHist_rat->Draw("e1"); 
      for (TLine* line : lines) line->Draw("SAME");
      canvas->SaveAs(Form("%s/ValidCompare/gamma_jet%i_period%s.pdf", plotPath.c_str(), etabin, (string(periodType%2==0 ? "A":"") + string(periodType>0 ? "B":"")).c_str()));
      if (gJetHist_rat) delete gJetHist_rat;
    }
  }

}
