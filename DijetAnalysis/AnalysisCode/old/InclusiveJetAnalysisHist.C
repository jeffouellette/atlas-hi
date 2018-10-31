#include "../Util.C"
#include "../TreeVariables.C"

void InclusiveJetAnalysisHist() {

  Initialize(0, true);

  //if (debugStatements) {
  // cout << "Status: In IdealPtAnalysisHist.C (breakpoint A): Building trigger pt histograms with " << numDataSets << " runs being used" << endl;
  // cout << "Status: In IdealPtAnalysisHist.C (breakpoint B): Numtrigs = " << numtrigs << endl;
  // cout << "Status: In IdealPtAnalysisHist.C (breakpoint C): Numetabins = " << numetabins << endl;
  // cout << "Status: In IdealPtAnalysisHist.C (breakpoint D): Numpbins = " << numpbins << endl;
  // cout << "Status: In IdealPtAnalysisHist.C (breakpoint E): ptPath = " << ptPath << endl;
  //}

  const bool cutEtaPhiPlot = false;
  const double ymin = 5e-8;
  const double ymax = 1e7;
//    const Style_t mkstyles[2] = {kFullCircle, kOpenCircle};
  const Style_t mkstyles[2] = {kFullDotMedium, kFullDotMedium};
  const Color_t mkcolors[8] = {kOrange+4, kOrange-3, kRed, kViolet, kBlue, kCyan-2, kGreen+3, kTeal+9};

  double hscale, deta;
  TH1F* thisHist;

  if (debugStatements) cout << "Status: In IdealPtAnalysisHist.C (breakpoint F): Initialized histArr histograms..." << endl;
  TH1F* jetPtHistArr[numetabins][2][3]; // for pt spectra
  TH1F* fcalSumEt[2];
  TH2F* jetEtaPhiHist = new TH2F("jetEtaPhiHist", ";#eta;#phi;Counts / d#eta d#phi", 98, -4.9, 4.9, 80, 0, 2*pi);
  jetEtaPhiHist->Sumw2();
  TH2F* subJetEtaPhiHist = new TH2F("subJetEtaPhiHist", ";#eta;#phi;Counts / d#eta d#phi", 98, -4.9, 4.9, 80, 0, 2*pi);
  subJetEtaPhiHist->Sumw2();
  for (short dir = 0; dir < 2; dir++) {
   fcalSumEt[dir] = new TH1F(Form("fcal%sSumEt", (dir==0?"A":"C")), ";#Sigma #it{E}_{T} #left[GeV#right];Counts", 100, -50, 250);
   fcalSumEt[dir]->Sumw2();
  }
  TH1F* jetPtResponseCalib = new TH1F("jetPtResponseCalib", ";Jet #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 50, 0, 1.4);
  jetPtResponseCalib->Sumw2();
  TH1F* jetPtResponseReco = new TH1F("jetPtResponseReco", ";Jet #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 50, 0, 1.4);
  jetPtResponseReco->Sumw2();
  TH1F* jetEnergyResponseCalib = new TH1F("jetEnergyResponseCalib", ";Jet #it{E}^{reco} / #it{E}^{truth};", 50, 0, 1.4);
  jetEnergyResponseCalib->Sumw2();
  TH1F* jetEnergyResponseReco = new TH1F("jetEnergyResponseReco", ";Jet #it{E}^{reco} / #it{E}^{truth};", 50, 0, 1.4);
  jetEnergyResponseReco->Sumw2();

  for (short pType = 0; pType < 3; pType++) {
   const char* periodStr = GetPeriodStr(pType); 
   for (short dType = 0; dType < 2; dType++) {
    const char* dataStr = GetDataStr(dType);
    for (int etabin = 0; etabin < numetabins; etabin++) {
     const TString name = Form("jetpt_etabin%i_dType%s_%s", etabin, dataStr, periodStr);
     //jetPtHistArr[etabin][dType][pType] = new TH1F(name, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}d#eta #left[nb GeV^{-1}#right]", numpbins, pbins);
     jetPtHistArr[etabin][dType][pType] = new TH1F(name, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}N/Ad#it{p}_{T}d#eta #left[nb GeV^{-1}#right]", numpbins, pbins);
     jetPtHistArr[etabin][dType][pType]->Sumw2();
     //numeratorHistArr[etabin][dType] = new TH1F(Form("numerator_etabin%i_%s", etabin, dataStr.c_str()), ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma_{A}/d#sigma_{B}", numpbins, pbins);
     //numeratorHistArr[etabin][dType]->Sumw2();
     //denominatorHistArr[etabin][dType] = new TH1F(Form("denominator_etabin%i_s", etabin, dataStr.c_str()), ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma_{A}/d#sigma_{B}", numpbins, pbins);
     //denominatorHistArr[etabin][dType]->Sumw2();
    }
   }
  }
  if (debugStatements) {
   cout << "Status: In IdealPtAnalysisHist.C (breakpoint G): histArr histograms initialized." << endl;
   cout << "Status: In IdealPtAnalysisHist.C (breakpoint H): Starting loop over triggers..." << endl;
  }

  double totalLuminosity = 0;

  /**** Fill eta-phi correlation plot with results from event loops ****/
  {
   TSystemDirectory dir(ptPath.c_str(), ptPath.c_str());
   TList* sysfiles = dir.GetListOfFiles();
   if (sysfiles) {
    TSystemFile *sysfile;
    TString fname;
    TString histName;
    TIter next(sysfiles);
    TVectorD* infoVec;

    while ((sysfile=(TSystemFile*)next())) {
     fname = sysfile->GetName();
     if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
      if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl; 
      for (int runNumber : runNumbers) {
       const bool pA = runNumber < 313500;
       if (SkipRun(runNumber)) continue;
       if (fname.Contains(to_string(runNumber))) {
        TFile* thisFile = new TFile(ptPath + fname, "READ");
        //totalLuminosity += (TVectorD*)thisFile->Get(Form("lum_vec_%i", dataSet))[0];

        // quickly check the parameters stored in this root file
        infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i_per%s", runNumber, (pA?"A":"B")));

        histName = Form("etaPhiHist_dSet%i_%s", runNumber, (pA?"pPb":"Pbp"));
        jetEtaPhiHist->Add((TH2F*)thisFile->Get(histName));
        histName = Form("subJetEtaPhiHist_dSet%i_%s", runNumber, (pA?"pPb":"Pbp"));
        subJetEtaPhiHist->Add((TH2F*)thisFile->Get(histName));

        thisFile->Close();
        delete thisFile;
        break;
       }
      }
      for (TString mcSample : mcSamples) {
       const bool pA = mcSample.Contains("pPb");
       if (SkipMC(pA)) continue; 
       if (fname.Contains(mcSample)) {
        TFile* thisFile = new TFile(ptPath + fname, "READ");
        //totalLuminosity += (TVectorD*)thisFile->Get(Form("lum_vec_%i", dataSet))[0];

        // quickly check the parameters stored in this root file
        infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", mcSample.Data()));

        histName = Form("etaPhiHist_dSet%s", mcSample.Data());
        jetEtaPhiHist->Add((TH2F*)thisFile->Get(histName));
        histName = Form("subJetEtaPhiHist_dSet%s", mcSample.Data());
        subJetEtaPhiHist->Add((TH2F*)thisFile->Get(histName));

        thisFile->Close();
        delete thisFile;
        break;
       }
      }
     }
    }
   }
  }
  /**** End fill eta-phi correlation plot ****/


  /**** Fill in HEC gap region with average ****/
  //int nbins_x = jetEtaPhiHist->GetNbinsX();
  //int nbins_y = jetEtaPhiHist->GetNbinsY();
  //if (cutEtaPhiPlot) {
  // for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
  //  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
  //   double x = jetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
  //   double y = jetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
  //   if (lowerEtaCut < x && x < upperEtaCut && lowerPhiCut < y && y < upperPhiCut) {
  //    jetEtaPhiHist->SetBinContent(bin_x+1, bin_y+1, 0);
  //   }
  //  }
  // }
  //}
  //else {
  // for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
  //  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {  
  //   double x = jetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
  //   double y = jetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
  //   if ((lowerPhiCut < y && y < upperPhiCut) && (lowerEtaCut < x && x < upperEtaCut)) { // if true, we are in the disabled HEC
  //    double integral_dy = 0;
  //    int numbins_dy = 0; 
  //    for (int bin_y_prime = 0; bin_y_prime < nbins_y; bin_y_prime++) {
  //     double y_prime = jetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
  //     if (!(lowerPhiCut < y_prime && y_prime < upperPhiCut)) { // if true we are outside of the disabled HEC
  //      integral_dy += jetEtaPhiHist->GetBinContent(bin_x+1, bin_y_prime+1);
  //      numbins_dy++; 
  //     }
  //    }
  //    if (numbins_dy != 0) integral_dy = integral_dy / (double)numbins_dy; // take the average
  //    jetEtaPhiHist->SetBinContent(bin_x+1, bin_y+1, integral_dy);
  //   }
  //  }
  // }
  //}
  //// average out HEC region on just subleading jets plot
  //if (cutEtaPhiPlot) {
  // for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
  //  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
  //   double x = subJetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
  //   double y = subJetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
  //   if (lowerEtaCut < x && x < upperEtaCut && lowerPhiCut < y && y < upperPhiCut) {
  //    subJetEtaPhiHist->SetBinContent(bin_x+1, bin_y+1, 0);
  //   }
  //  }
  // }
  //}
  //else {
  // for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
  //  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {  
  //   double x = subJetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
  //   double y = subJetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
  //   if ((lowerPhiCut < y && y < upperPhiCut) && (lowerEtaCut < x && x < upperEtaCut)) { // if true, we are in the disabled HEC
  //    double integral_dy = 0;
  //    int numbins_dy = 0; 
  //    for (int bin_y_prime = 0; bin_y_prime < nbins_y; bin_y_prime++) {
  //     double y_prime = subJetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
  //     if (!(lowerPhiCut < y_prime && y_prime < upperPhiCut)) { // if true we are outside of the disabled HEC
  //      integral_dy += subJetEtaPhiHist->GetBinContent(bin_x+1, bin_y_prime+1);
  //      numbins_dy++; 
  //     }
  //    }
  //    if (numbins_dy != 0) integral_dy = integral_dy / (double)numbins_dy; // take the average
  //    subJetEtaPhiHist->SetBinContent(bin_x+1, bin_y+1, integral_dy);
  //   }
  //  }
  // }
  //}
  /**** End fill in HEC region ****/


  /**** Fill summed histograms with results from event loops ****/
  {
   TSystemDirectory dir(ptPath.c_str(), ptPath.c_str());
   TList* sysfiles = dir.GetListOfFiles();
   if (sysfiles) {
    TSystemFile *sysfile;
    TString fname;
    TString histName;
    TIter next(sysfiles);
    TVectorD* infoVec;

    while ((sysfile=(TSystemFile*)next())) {
     fname = sysfile->GetName();
     if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
      if (debugStatements) cout << "Status: In IdealPtAnalysisHist.C (breakpoint I): Found " << fname.Data() << endl; 

      for (int runNumber : runNumbers) {
       const bool pA = runNumber < 313500;
       if (SkipRun(pA)) continue;

       if (fname.Contains(to_string(runNumber))) {
        TFile* thisFile = new TFile(ptPath + fname, "READ");
        //totalLuminosity += (TVectorD*)thisFile->Get(Form("lum_vec_%i", dataSet))[0];

        // quickly check the parameters stored in this root file
        infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", runNumber));
        totalLuminosity += (*infoVec)[4];

        int actetabin; // used to flip period A pseudorapidities
        for (int etabin = 0; etabin < numetabins; etabin++) {
         if (pA) actetabin = numetabins - etabin - 1;
         else actetabin = etabin;
         histName = Form("jetPtSpectrum_dSet%i_%s_etabin%i", runNumber, (pA?"pPb":"Pbp"), actetabin);
         thisHist = (TH1F*)thisFile->Get(histName);

         //if (scaleBackHEC) {
         // double scaleHEC = 0;
         // double numerator = 0;
         // double denominator = 0;  
         // for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
         //  double x = jetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
         //  if (x < etabins[actetabin] || etabins[actetabin+1] < x) continue;
         //  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
         //   double y = jetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
         //   double content = jetEtaPhiHist->GetBinContent(bin_x+1, bin_y+1);
         //   numerator += content;
         //   if (!InDisabledHEC(x, y))
         //    denominator += content; // if true we are in the HEC, add to "total" events
         //  }
         // }
         // if (denominator != 0) scaleHEC = numerator/denominator;
      
         // thisHist->Scale(scaleHEC);
         //}

         const short pType = (pA ? 0:1); // 0 for period A
         jetPtHistArr[etabin][0][pType]->Add(thisHist);
         jetPtHistArr[etabin][0][2]->Add(thisHist);

        }
        for (short dir = 0; dir < 2; dir++) {
         histName = Form("fcal%sSumEt_dSet%i_%s", (dir==0?"A":"C"), runNumber, (pA?"pPb":"Pbp"));
         thisHist = (TH1F*)thisFile->Get(histName);
         if (pA)
          fcalSumEt[(dir+1)%2]->Add(thisHist);
         else
          fcalSumEt[dir]->Add(thisHist); 
        }
        thisFile->Close();
        delete thisFile;
        break;
       }
      }

      for (TString mcSample : mcSamples) {
       const bool pA = mcSample.Contains("pPb");
       if (SkipMC(pA)) continue;

       if (fname.Contains(mcSample)) {
        TFile* thisFile = new TFile(ptPath + fname, "READ");

        // quickly check the parameters stored in this root file
        infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", mcSample.Data()));
        //totalLuminosity += (*infoVec)[4];

        jetEnergyResponseCalib->Add((TH1F*)thisFile->Get(Form("jetEnergyResponseCalib_dSet%s", mcSample.Data())));
        jetEnergyResponseReco->Add((TH1F*)thisFile->Get(Form("jetEnergyResponseReco_dSet%s", mcSample.Data())));
        jetPtResponseCalib->Add((TH1F*)thisFile->Get(Form("jetPtResponseCalib_dSet%s", mcSample.Data())));
        jetPtResponseReco->Add((TH1F*)thisFile->Get(Form("jetPtResponseReco_dSet%s", mcSample.Data())));

        int actetabin; // used to flip period A pseudorapidities
        for (int etabin = 0; etabin < numetabins; etabin++) {
         if (pA) actetabin = numetabins - etabin - 1;
         else actetabin = etabin;
         histName = Form("jetPtSpectrum_dSet%s_etabin%i", mcSample.Data(), actetabin);
         thisHist = (TH1F*)thisFile->Get(histName);

         //if (scaleBackHEC) {
         // double scaleHEC = 0;
         // double numerator = 0;
         // double denominator = 0;  
         // for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
         //  double x = jetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
         //  if (x < etabins[actetabin] || etabins[actetabin+1] < x) continue;
         //  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
         //   double y = jetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
         //   double content = jetEtaPhiHist->GetBinContent(bin_x+1, bin_y+1);
         //   numerator += content;
         //   if (!InDisabledHEC(x, y))
         //    denominator += content; // if true we are in the HEC, add to "total" events
         //  }
         // }
         // if (denominator != 0) scaleHEC = numerator/denominator;
      
         // thisHist->Scale(scaleHEC);
         //}

         const short pType = (pA ? 0:1); // 0 for period A
         jetPtHistArr[etabin][1][pType]->Add(thisHist);
         jetPtHistArr[etabin][1][2]->Add(thisHist);
        }
        for (short dir = 0; dir < 2; dir++) {
         histName = Form("fcal%sSumEt_dSet%s", (dir==0?"A":"C"), mcSample.Data());
         thisHist = (TH1F*)thisFile->Get(histName);
         if (pA)
          fcalSumEt[2-dir-1]->Add(thisHist);
         else
          fcalSumEt[dir]->Add(thisHist); 
        }
        thisFile->Close();
        delete thisFile;
        break;
       }
      }
     }
    }
   }
  }
  /**** End fill summed histograms ****/


  /** Plotting routines **/

  // Plot best-selected pt spectra
  double* jetPtHistArrScales;
  if (scaleAnalyses) jetPtHistArrScales = linspace(-1.5, 1.5, numetabins/2 - 1);
  else jetPtHistArrScales = linspace(0, 0, numetabins/2 - 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts
  TCanvas* canvas = new TCanvas("trigCanvas", "", 800, 600);
  canvas->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTicks();
  gStyle->SetErrorX(0);
  canvas->Draw();
  bool plotted = false;

  for (short dType = 0; dType < 2; dType++) {
   if (!runData && dType == 0) continue;
   if (!runMC && dType == 1) continue;

   for (int etabin = 0; etabin < numetabins; etabin++) {
    thisHist = jetPtHistArr[etabin][dType][2];
    for (int pbin = 0; pbin < numpbins; pbin++) {
     if (kinematicLumiVec[pbin][etabin] != 0) {
      thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / kinematicLumiVec[pbin][etabin]);
      thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / kinematicLumiVec[pbin][etabin]);
     }
     else if (debugStatements) {
      cout << "Warning: In IdealPtAnalysisHist.C (breakpoint J): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
     }
    }
    deta = etabins[etabin+1] - etabins[etabin];
    thisHist->Scale(1e3*TMath::Power(10, jetPtHistArrScales[(int)(numetabins/2 - 0.5 -TMath::Abs(etabin - numetabins/2 + 0.5))])/deta, "width"); // separate different etabins
//      thisHist->SetMarkerStyle(kDot);
    const Style_t kStyle = mkstyles[etabin%2];
    const Color_t kColor = mkcolors[etabin%8];
    thisHist->SetMarkerStyle(kStyle);
    thisHist->SetMarkerColor(kColor);
    thisHist->SetLineColor(kColor);
    thisHist->SetMinimum(ymin);
    thisHist->SetMaximum(ymax);
    thisHist->GetYaxis()->SetTitleOffset(1.35);
    thisHist->GetXaxis()->SetTickLength(0.02);
    thisHist->GetYaxis()->SetTickLength(0.02);

    if (plotted) thisHist->Draw("same e1");
    else {
     thisHist->Draw("e1");
     plotted = true;
    }
    
    const float textx = 0.47 + (etabin>=(numetabins/2))*0.25;
    const float texty = 0.91 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
    const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[etabin], etabins[etabin+1], jetPtHistArrScales[(int)((0.5*(numetabins-1))-TMath::Abs(etabin-(0.5*(numetabins-1))))]);
    myMarkerText (textx, texty, kColor, kFullCircle, text);
   }

   myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity));
   myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

   string histName;
   if (numetabins > 1) {
    histName = "ptSpectra_allTrigs";
   }
   else histName = "ptSpectra_allTrigs_integratedEta";

   if (runPeriodA && !runPeriodB) {
    myText (0.19, 0.33, kBlack, "Period A");
    histName = histName + "_periodA";
   }
   else if (!runPeriodA && runPeriodB) {
    myText (0.19, 0.33, kBlack, "Period B");
    histName = histName + "_periodB";
   }
   else {
    myText (0.19, 0.33, kBlack, "Period A+B");
   }
   if (highPtJetsOnly) histName += "_highPtJetsOnly";
   if (runPeriodA || runPeriodB) canvas->SaveAs((plotPath + "ptSpectra/" + histName + ".pdf").c_str());
  }


  /**** Plot fCal Sum Et distributions ****/
  {
   Color_t fcalColors[2] = {kBlack, kBlue};
   gPad->SetLogx(false);
   gPad->SetLogy(true);
   for (short dir = 0; dir < 2; dir++) {
    fcalSumEt[dir]->SetMarkerColor(fcalColors[dir]);
    fcalSumEt[dir]->SetLineColor(fcalColors[dir]);
    fcalSumEt[dir]->GetYaxis()->SetRangeUser(1, 1200);
    //fcalSumEt[dir]->SetMarkerStyle(kFullDotMedium);
    if (dir == 0) fcalSumEt[dir]->Draw("e1 x0");
    else fcalSumEt[dir]->Draw("same e1 x0");
    const char* text = (dir==0?"Proton-going direction":"Lead-going direction");
    myMarkerText(0.6, 0.85 - dir*0.08, fcalColors[dir], kFullCircle, text);
   }
   canvas->SaveAs((plotPath + "mc_fcal_distributions.pdf").c_str());
  }


  /**** Plot jet energy response ****/
  {
   gPad->SetLogx(false);
   gPad->SetLogy(true);

   float ncounts = jetEnergyResponseCalib->Integral();
   jetEnergyResponseCalib->SetLineColor(kBlue);
   jetEnergyResponseCalib->SetMarkerColor(kBlue);
   jetEnergyResponseCalib->Scale(1./ncounts, "width");

   TF1* calibFit = new TF1("calibFit", "gaus(0)", 0, 2.0);
   jetEnergyResponseCalib->Fit(calibFit, "RLN");
   float m = calibFit->GetParameter(1);
   float s = calibFit->GetParameter(2);
   if (calibFit) delete calibFit;

   calibFit = new TF1("calibFit2", "gaus(0)", m-1.5*s, m+1.5*s);
   jetEnergyResponseCalib->Fit(calibFit, "RLN");

   jetEnergyResponseCalib->GetYaxis()->SetTitle("Fractional counts / 0.01 GeV");
   jetEnergyResponseCalib->GetYaxis()->SetRangeUser(1e-4, 20);
   jetEnergyResponseCalib->Draw("e1 x0");
   calibFit->SetLineColor(kBlue);
   calibFit->Draw("same");

   ncounts = jetEnergyResponseReco->Integral();
   jetEnergyResponseReco->SetLineColor(kBlack);
   jetEnergyResponseReco->SetMarkerColor(kBlack);
   jetEnergyResponseReco->Scale(1./ncounts, "width");

   TF1* recoFit = new TF1("recoFit", "gaus(0)", 0, 2.0);
   jetEnergyResponseReco->Fit(recoFit, "RLN");
   m = recoFit->GetParameter(1);
   s = recoFit->GetParameter(2);
   if (recoFit) delete recoFit;

   recoFit = new TF1("recoFit2", "gaus(0)", m-1.5*s, m+1.5*s);
   jetEnergyResponseReco->Fit(recoFit, "RLN");

   jetEnergyResponseReco->GetYaxis()->SetTitle("Fractional counts / 0.01 GeV");
   jetEnergyResponseReco->GetYaxis()->SetRangeUser(1e-4, 20);
   jetEnergyResponseReco->Draw("same e1 x0");
   recoFit->SetLineColor(kBlack);
   recoFit->Draw("same");

   myText(0.19, 0.89, kBlack, "Jet energy response");
   myText(0.19, 0.83, kBlack, Form("%i jets", (int)ncounts));
   myText(0.19, 0.77, kBlack, Form("JES = %.2f #pm %.2f", calibFit->GetParameter(1), calibFit->GetParError(1)));
   myText(0.19, 0.71, kBlack, Form("JER = %.2f #pm %.2f", calibFit->GetParameter(2), calibFit->GetParError(2)));
   canvas->SaveAs(Form("%s/jetEnergyResponse.pdf", plotPath.c_str()));
   if (calibFit) delete calibFit;
   if (recoFit) delete recoFit;
  }


  /**** Plots jet momentum response ****/
  {
   gPad->SetLogx(false);
   gPad->SetLogy(true);

   float ncounts = jetPtResponseCalib->Integral();
   jetPtResponseCalib->SetLineColor(kBlue);
   jetPtResponseCalib->SetMarkerColor(kBlue);
   jetPtResponseCalib->Scale(1./ncounts, "width");

   TF1* calibFit = new TF1("calibFit", "gaus(0)", 0, 2.0);
   jetPtResponseCalib->Fit(calibFit, "RLN");
   float m = calibFit->GetParameter(1);
   float s = calibFit->GetParameter(2);
   if (calibFit) delete calibFit;

   calibFit = new TF1("calibFit2", "gaus(0)", m-1.5*s, m+1.5*s);
   jetPtResponseCalib->Fit(calibFit, "RLN");
   
   jetPtResponseCalib->GetYaxis()->SetTitle("Fractional counts / 0.01 GeV");
   jetPtResponseCalib->GetYaxis()->SetRangeUser(1e-4, 20);
   jetPtResponseCalib->Draw("e1 x0");
   calibFit->SetLineColor(kBlue);
   calibFit->Draw("same");

   ncounts = jetPtResponseReco->Integral();
   jetPtResponseReco->SetLineColor(kBlack);
   jetPtResponseReco->SetMarkerColor(kBlack);
   jetPtResponseReco->Scale(1./ncounts, "width");

   TF1* recoFit = new TF1("recoFit", "gaus(0)", 0, 2.0);
   jetPtResponseReco->Fit(recoFit, "RLN");
   m = recoFit->GetParameter(1);
   s = recoFit->GetParameter(2);
   if (recoFit) delete recoFit;

   recoFit = new TF1("recoFit2", "gaus(0)", m-1.5*s, m+1.5*s);
   jetPtResponseReco->Fit(recoFit, "RLN");

   jetPtResponseReco->GetYaxis()->SetTitle("Fractional counts / 0.01 GeV");
   jetPtResponseReco->GetYaxis()->SetRangeUser(1e-4, 20);
   jetPtResponseReco->Draw("same e1 x0");
   recoFit->SetLineColor(kBlack);
   recoFit->Draw("same");

   myText(0.19, 0.89, kBlack, "Jet momentum response");
   myText(0.19, 0.83, kBlack, Form("%i jets", (int)ncounts));
   myText(0.19, 0.77, kBlack, Form("JES = %.2f #pm %.2f", calibFit->GetParameter(1), calibFit->GetParError(1)));
   myText(0.19, 0.71, kBlack, Form("JER = %.2f #pm %.2f", calibFit->GetParameter(2), calibFit->GetParError(2)));
   canvas->SaveAs(Form("%s/jetPtResponse.pdf", plotPath.c_str()));
  }

  /**** Plot pt spectrum ratios ****/
  if (runPeriodA && runPeriodB) {
   delete canvas;
   delete[] jetPtHistArrScales;
   jetPtHistArrScales = linspace(-2, 1, numetabins/2 - 1);
   canvas = new TCanvas("ratioCanvas", "", 800, 600);
   canvas->cd();
   //gPad->SetLogy();
   gPad->SetLogx();
   gPad->SetTicks();
   gStyle->SetErrorX(0);
   canvas->Draw();

   double*** periodALuminosities = new double**[numpbins];
   double*** periodBLuminosities = new double**[numpbins];
   for (short pbin = 0; pbin < numpbins; pbin++) {
    periodALuminosities[pbin] = new double*[numetabins];
    periodBLuminosities[pbin] = new double*[numetabins];
    for (short etabin = 0; etabin < numetabins; etabin++) { 
     periodALuminosities[pbin][etabin] = new double[2];
     periodBLuminosities[pbin][etabin] = new double[2];
     for (short dType = 0; dType < 2; dType++) {
      periodALuminosities[pbin][etabin][dType] = 0;
      periodBLuminosities[pbin][etabin][dType] = 0;
     }
    }
   }

   // Data luminosities 
   Trigger* bestTrigger;
   double bestLumi;
   int thisRunNumber;
   double** allLuminosities = GetTriggerLuminosities();
   for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
    thisRunNumber = runNumbers[rnIndex];
    const bool pA = thisRunNumber < 313500;
    /**** Need to get kinematic trigger vec for this run ****/
    SetBestTriggers(rnIndex);

    if (pA) {
     for (int pbin = 0; pbin < numpbins; pbin++) {
      for (int etabin = 0; etabin < numetabins; etabin++) {
       bestTrigger = kinematicTriggerVec[pbin][etabin];
       if (bestTrigger == NULL) continue;
       const int actetabin = numetabins - etabin - 1;
       bestLumi = allLuminosities[bestTrigger->index][rnIndex];
       periodALuminosities[pbin][actetabin][0] += bestLumi;
      }
     }
    }
    else if (!pA) {
     for (int pbin = 0; pbin < numpbins; pbin++) {
      for (int etabin = 0; etabin < numetabins; etabin++) {
       bestTrigger = kinematicTriggerVec[pbin][etabin];
       if (bestTrigger == NULL) continue;
       bestLumi = allLuminosities[bestTrigger->index][rnIndex];
       periodBLuminosities[pbin][etabin][0] += bestLumi;
      }
     }
    }
   }
   delete[] allLuminosities;
 
   // TODO: MC luminosities - set to 0 for now
   for (int pbin = 0; pbin < numpbins; pbin++) {
    for (int etabin = 0; etabin < numetabins; etabin++) { 
     periodALuminosities[pbin][etabin][1] = 0;
     periodBLuminosities[pbin][etabin][1] = 0;
    }
   }

   TLine* lineDrawer = new TLine();
   lineDrawer->SetLineStyle(6);
   lineDrawer->SetLineWidth(2);

   for (short dType = 0; dType < 2; dType++) {
    for (int etabin = 0; etabin < numetabins; etabin++) {
     // Scale numerator by luminosity to get a cross section
     thisHist = jetPtHistArr[etabin][dType][0];
     for (int pbin = 0; pbin < numpbins; pbin++) {
      double thisLumi = periodALuminosities[pbin][etabin][dType];
      if (thisLumi != 0) {
       thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / thisLumi);
       thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / thisLumi);
      }
      else if (debugStatements) cout << "Warning: In IdealPtAnalysisHist.C (breakpoint J): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
     }
     // Scale denominator by luminosity to get a cross section
     thisHist = jetPtHistArr[etabin][dType][1];
     for (int pbin = 0; pbin < numpbins; pbin++) {
      double thisLumi = periodBLuminosities[pbin][etabin][dType];
      if (thisLumi != 0) {
       thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / thisLumi);
       thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / thisLumi);
      }
      else if (debugStatements) cout << "Warning: In IdealPtAnalysisHist.C (breakpoint J): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
     }
     // Now divide the histograms to get a ratio of the cross sections
     jetPtHistArr[etabin][dType][0]->Divide(jetPtHistArr[etabin][dType][1]);

     // Now go through the plotting routines
     thisHist = jetPtHistArr[etabin][dType][0];
     double scale = jetPtHistArrScales[(int)(numetabins/2 - 0.5 - TMath::Abs(etabin - numetabins/2 + 0.5))];
     for (int pbin = 0; pbin < numpbins; pbin++) {
      if (thisHist->GetBinContent(pbin+1) == 0.) thisHist->SetBinContent(pbin+1, -5);
      else thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) + scale);
     }

//   thisHist->SetMarkerStyle(kDot);
     const Style_t kStyle = mkstyles[etabin%2];
     const Color_t kColor = mkcolors[etabin%8];
     thisHist->SetMarkerStyle(kStyle);
     thisHist->SetMarkerColor(kColor);
     thisHist->SetLineColor(kColor);
     thisHist->SetMinimum(-2);
     thisHist->SetMaximum(3.5);
     thisHist->GetYaxis()->SetTitleOffset(1.35);
     thisHist->GetXaxis()->SetTickLength(0.02);
     thisHist->GetYaxis()->SetTickLength(0.02);
     if (etabin == 0) thisHist->Draw("e1");
     else thisHist->Draw("same e1");

     lineDrawer->SetLineColor(kColor);
     if (etabins[etabin] < 0) lineDrawer->DrawLine(pbins[0], scale+1, pbins[numpbins], scale+1);
     else lineDrawer->DrawLine(0.5*(pbins[0]+pbins[1]), scale+1, 0.5*(pbins[numpbins-1]+pbins[numpbins]), scale+1);

     const float textx = 0.47 + (etabin>=(numetabins/2))*0.25;
     const float texty = 0.91 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
     const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[etabin], etabins[etabin+1], scale);
     myMarkerText (textx, texty, kColor, kFullCircle, text);
    }
   } 

   myText (0.19, 0.91, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity));
   myText (0.19, 0.85, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

   string histName = "ptSpectra_ratio_allTrigs";
   if (highPtJetsOnly) histName += "_highPtJetsOnly";
   canvas->SaveAs((plotPath + "ptSpectra/" + histName + ".pdf").c_str());
   delete lineDrawer;
  }
  /**** End plot ratios ****/


  /**** Plot eta-phi jet phase space diagram ****/
  delete canvas;
  canvas = new TCanvas("etaPhiCanvas", "", 800, 600);
  canvas->cd();

  const Int_t Number = 3;
  Double_t Red[Number]    = { 0.00, 1.00, 1.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 0.00, 0.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=100;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  canvas->SetRightMargin(0.18);
  canvas->SetTopMargin(-0.02);
  jetEtaPhiHist->GetXaxis()->SetTitleOffset(0.8);
  jetEtaPhiHist->GetYaxis()->SetTitleOffset(0.8);
  jetEtaPhiHist->GetZaxis()->SetTitleOffset(1.2);

  jetEtaPhiHist->Scale(1., "width");
  jetEtaPhiHist->Draw("colz");
  string histName = "inclusive_jets_eta_phi_correlation";
  myText (0.19, 0.39, kBlack, "Inclusive jets");
  if (runPeriodA && !runPeriodB) {
   myText (0.19, 0.33, kBlack, "Period A");
   histName = histName + "_periodA";
  }
  else if (!runPeriodA && runPeriodB) {
   myText (0.19, 0.33, kBlack, "Period B");
   histName = histName + "_periodB";
  }
  else {
   myText (0.19, 0.33, kBlack, "Period A & B");
  }
  if (cutEtaPhiPlot) {
   histName += "_withoutDroppedBins";
  }

  TLine* lineDrawer = new TLine();
  lineDrawer->SetLineStyle(7);
  lineDrawer->SetLineColor(kBlack);
  lineDrawer->DrawLine(lowerEtaCut+0.02, lowerPhiCut+0.02, upperEtaCut-0.02, lowerPhiCut+0.02);
  lineDrawer->DrawLine(lowerEtaCut+0.02, upperPhiCut-0.02, upperEtaCut-0.02, upperPhiCut-0.02);
  lineDrawer->DrawLine(lowerEtaCut+0.02, lowerPhiCut+0.02, lowerEtaCut+0.02, upperPhiCut-0.02);
  lineDrawer->DrawLine(upperEtaCut-0.02, lowerPhiCut+0.02, upperEtaCut-0.02, upperPhiCut-0.02);

  if (highPtJetsOnly) histName += "_highPtJetsOnly";
      
  canvas->SaveAs((plotPath + histName + ".pdf").c_str());

  histName = "jetEtaPhiHist";
  if (highPtJetsOnly) histName += "_highPtJetsOnly";
  TFile* output = new TFile((rootPath + histName + ".root").c_str(), "RECREATE");
  jetEtaPhiHist->Write(); // save the eta phi distribution for further use
  subJetEtaPhiHist->Write();
  output->Close();

  delete output;
  delete jetEtaPhiHist;
  /**** End plot eta-phi diagram ****/


  /**** Plot eta-phi jet phase space diagram ****/
  delete canvas;
  canvas = new TCanvas("subleadingEtaPhiCanvas", "", 800, 600);
  canvas->cd();

  canvas->SetRightMargin(0.18);
  canvas->SetTopMargin(-0.02);
  subJetEtaPhiHist->GetXaxis()->SetTitleOffset(0.8);
  subJetEtaPhiHist->GetYaxis()->SetTitleOffset(0.8);
  subJetEtaPhiHist->GetZaxis()->SetTitleOffset(1.2);

  subJetEtaPhiHist->Scale(1., "width");
  subJetEtaPhiHist->Draw("colz");
  histName = "subleading_jets_eta_phi_correlation";
  myText (0.19, 0.39, kBlack, "Subleading jets");
  if (runPeriodA && !runPeriodB) {
   myText (0.19, 0.33, kBlack, "Period A");
   histName = histName + "_periodA";
  }
  else if (!runPeriodA && runPeriodB) {
   myText (0.19, 0.33, kBlack, "Period B");
   histName = histName + "_periodB";
  }
  else {
   myText (0.19, 0.33, kBlack, "Period A & B");
  }
  if (cutEtaPhiPlot) {
   histName += "_withoutDroppedBins";
  }

  lineDrawer->SetLineStyle(7);
  lineDrawer->SetLineColor(kBlack);
  lineDrawer->DrawLine(lowerEtaCut+0.02, lowerPhiCut+0.02, upperEtaCut-0.02, lowerPhiCut+0.02);
  lineDrawer->DrawLine(lowerEtaCut+0.02, upperPhiCut-0.02, upperEtaCut-0.02, upperPhiCut-0.02);
  lineDrawer->DrawLine(lowerEtaCut+0.02, lowerPhiCut+0.02, lowerEtaCut+0.02, upperPhiCut-0.02);
  lineDrawer->DrawLine(upperEtaCut-0.02, lowerPhiCut+0.02, upperEtaCut-0.02, upperPhiCut-0.02);

  if (highPtJetsOnly) histName += "_highPtJetsOnly";
      
  canvas->SaveAs((plotPath + histName + ".pdf").c_str());

  delete subJetEtaPhiHist;
  /**** End plot eta-phi diagram ****/


  /**** Free memory ****/
  for (short pType = 0; pType < 2; pType++) {
   for (short dType = 0; dType < 2; dType++) {
    for (int etabin = 0; etabin < numetabins; etabin++) {
     delete jetPtHistArr[etabin][dType][pType];
    }
   }
  }
  delete lineDrawer;
  delete canvas;
  delete[] jetPtHistArrScales;
  /**** End free memory ****/


  if (debugStatements) cout << "Status: In IdealPtAnalysisHist.C (breakpoint K): Finished plotting pt spectrum" << endl;
  return;
}
