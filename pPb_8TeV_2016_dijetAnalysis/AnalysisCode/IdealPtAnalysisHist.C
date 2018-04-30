#include "../Util.C"

void IdealPtAnalysisHist() {

  initialize(0, true);
  const bool isMC = useDataVersion == 0;

  vector<int>* dataSets = NULL;
  if (!isMC) dataSets = getRunNumbers();
  else dataSets = getMCSamples();
  
  //vector<int>* thisRunNumbers = getRunNumbers();

  const int numDataSets = (*dataSets).size();
  if (debugStatements) {
    cout << "Status: In IdealPtAnalysisHist.C (breakpoint A): Building trigger pt histograms with " << numDataSets << " runs being used" << endl;
    cout << "Status: In IdealPtAnalysisHist.C (breakpoint B): Numtrigs = " << numtrigs << endl;
    cout << "Status: In IdealPtAnalysisHist.C (breakpoint C): Numetabins = " << numetabins << endl;
    cout << "Status: In IdealPtAnalysisHist.C (breakpoint D): Numpbins = " << numpbins << endl;
    cout << "Status: In IdealPtAnalysisHist.C (breakpoint E): ptPath = " << ptPath << endl;
  }

  const bool cutEtaPhiPlot = false;
  const double ymin = 5e-8;
  const double ymax = 1e7;
//    const Style_t mkstyles[2] = {kFullCircle, kOpenCircle};
  const Style_t mkstyles[2] = {kFullDotMedium, kFullDotMedium};
  const Color_t mkcolors[8] = {kOrange+4, kOrange-3, kRed, kViolet, kBlue, kCyan-2, kGreen+3, kTeal+9};

  double hscale, deta;
  TH1D* thisHist;

  if (debugStatements) cout << "Status: In IdealPtAnalysisHist.C (breakpoint F): Initialized histArr histograms..." << endl;
  TH1D* histArr[numetabins]; // for pt spectra
  TH1D* numeratorHistArr[numetabins]; // for period A - period B pt spectrum ratio
  TH1D* denominatorHistArr[numetabins];
  TH2D* etaPhiHist = new TH2D("etaPhiHist", ";#eta;#phi;Counts / d#eta d#phi", 98, -4.9, 4.9, 100, 0, 2*pi);
  etaPhiHist->Sumw2();
  TH2D* subleadingEtaPhiHist = new TH2D("subleadingEtaPhiHist", ";#eta;#phi;Counts / d#eta d#phi", 98, -4.9, 4.9, 100, 0, 2*pi);
  subleadingEtaPhiHist->Sumw2();
  
  for (int etabin = 0; etabin < numetabins; etabin++) {
    histArr[etabin] = new TH1D(Form("best_statistics_etabin%i", etabin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}d#eta #left[nb GeV^{-1}#right]", numpbins, pbins);
    histArr[etabin]->Sumw2();
    numeratorHistArr[etabin] = new TH1D(Form("best_statistics_numerator_etabin%i", etabin), ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma_{A}/d#sigma_{B}", numpbins, pbins);
    numeratorHistArr[etabin]->Sumw2();
    denominatorHistArr[etabin] = new TH1D(Form("best_statistics_denominator_etabin%i", etabin), ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma_{A}/d#sigma_{B}", numpbins, pbins);
    denominatorHistArr[etabin]->Sumw2();
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
          for (int dataSet : *dataSets) {
            if (!isMC && skipRun(dataSet)) continue;
            else if (isMC && skipMC(dataSet)) continue; 
            if (fname.Contains(to_string(dataSet))) {
              TFile* thisFile = new TFile(ptPath + fname, "READ");
              //totalLuminosity += (TVectorD*)thisFile->Get(Form("lum_vec_%i", dataSet))[0];

              // quickly check the parameters stored in this root file
              infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", dataSet));
              assert ((int)(*infoVec)[0] == dataSet);
              assert ((int)(*infoVec)[1] == numetabins);
              assert ((int)(*infoVec)[2] == numtrigs);
              assert ((int)(*infoVec)[3] == numpbins);

              histName = Form("etaPhiHist_dataset%i", dataSet);
              etaPhiHist->Add((TH2D*)thisFile->Get(histName));
              histName = Form("subleadingEtaPhiHist_dataset%i", dataSet);
              subleadingEtaPhiHist->Add((TH2D*)thisFile->Get(histName));

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
  int nbins_x = etaPhiHist->GetNbinsX();
  int nbins_y = etaPhiHist->GetNbinsY();
  if (cutEtaPhiPlot) {
    for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
      for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
        double x = etaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
        double y = etaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
        if (lowerEtaCut < x && x < upperEtaCut && lowerPhiCut < y && y < upperPhiCut) {
          etaPhiHist->SetBinContent(bin_x+1, bin_y+1, 0);
        }
      }
    }
  }
  else {
    for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
      for (int bin_y = 0; bin_y < nbins_y; bin_y++) {  
        double x = etaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
        double y = etaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
        if ((lowerPhiCut < y && y < upperPhiCut) && (lowerEtaCut < x && x < upperEtaCut)) { // if true, we are in the disabled HEC
          double integral_dy = 0;
          int numbins_dy = 0; 
          for (int bin_y_prime = 0; bin_y_prime < nbins_y; bin_y_prime++) {
            double y_prime = etaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
            if (!(lowerPhiCut < y_prime && y_prime < upperPhiCut)) { // if true we are outside of the disabled HEC
              integral_dy += etaPhiHist->GetBinContent(bin_x+1, bin_y_prime+1);
              numbins_dy++; 
            }
          }
          if (numbins_dy != 0) integral_dy = integral_dy / (double)numbins_dy; // take the average
          etaPhiHist->SetBinContent(bin_x+1, bin_y+1, integral_dy);
        }
      }
    }
  }
  // average out HEC region on just subleading jets plot
  if (cutEtaPhiPlot) {
    for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
      for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
        double x = subleadingEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
        double y = subleadingEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
        if (lowerEtaCut < x && x < upperEtaCut && lowerPhiCut < y && y < upperPhiCut) {
          subleadingEtaPhiHist->SetBinContent(bin_x+1, bin_y+1, 0);
        }
      }
    }
  }
  else {
    for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
      for (int bin_y = 0; bin_y < nbins_y; bin_y++) {  
        double x = subleadingEtaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
        double y = subleadingEtaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
        if ((lowerPhiCut < y && y < upperPhiCut) && (lowerEtaCut < x && x < upperEtaCut)) { // if true, we are in the disabled HEC
          double integral_dy = 0;
          int numbins_dy = 0; 
          for (int bin_y_prime = 0; bin_y_prime < nbins_y; bin_y_prime++) {
            double y_prime = subleadingEtaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
            if (!(lowerPhiCut < y_prime && y_prime < upperPhiCut)) { // if true we are outside of the disabled HEC
              integral_dy += subleadingEtaPhiHist->GetBinContent(bin_x+1, bin_y_prime+1);
              numbins_dy++; 
            }
          }
          if (numbins_dy != 0) integral_dy = integral_dy / (double)numbins_dy; // take the average
          subleadingEtaPhiHist->SetBinContent(bin_x+1, bin_y+1, integral_dy);
        }
      }
    }
  }
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
          if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl; 
          for (int dataSet : *dataSets) {
            if (!isMC && skipRun(dataSet)) continue;
            else if (isMC && skipMC(dataSet)) continue;
            if (fname.Contains(to_string(dataSet))) {
              TFile* thisFile = new TFile(ptPath + fname, "READ");
              //totalLuminosity += (TVectorD*)thisFile->Get(Form("lum_vec_%i", dataSet))[0];

              // quickly check the parameters stored in this root file
              infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", dataSet));
              assert ((int)(*infoVec)[0] == dataSet);
              assert ((int)(*infoVec)[1] == numetabins);
              assert ((int)(*infoVec)[2] == numtrigs);
              assert ((int)(*infoVec)[3] == numpbins);
              totalLuminosity += (*infoVec)[4];

              int actetabin; // used to flip period A pseudorapidities
              for (int etabin = 0; etabin < numetabins; etabin++) {
                if (isPeriodA(dataSet)) actetabin = numetabins - etabin - 1;
                else actetabin = etabin;
                histName = Form("trig_pt_counts_dataset%i_etabin%i", dataSet, actetabin);
                thisHist = (TH1D*)thisFile->Get(histName);

                double scaleHEC = 0;
                double numerator = 0;
                double denominator = 0;  
                for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
                  double x = etaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
                  if (x < etabins[actetabin] || etabins[actetabin+1] < x) continue;
                  for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
                    double y = etaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
                    double content = etaPhiHist->GetBinContent(bin_x+1, bin_y+1);
                    numerator += content;
                    if (!(lowerEtaCut < x && x < upperEtaCut && lowerPhiCut < y && y < upperPhiCut)) denominator += content; // if true we are in the HEC, add to "total" events
                  }
                }
                if (denominator != 0) scaleHEC = numerator/denominator;
          
                thisHist->Scale(scaleHEC);

                histArr[etabin]->Add(thisHist);
                if (!runPeriodA || !runPeriodB) continue;
                if (isPeriodA(dataSet)) numeratorHistArr[etabin]->Add(thisHist);
                else denominatorHistArr[etabin]->Add(thisHist);
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
  double* histArrScales;
  if (scaleAnalyses) histArrScales = linspace(-1.5, 1.5, numetabins/2 - 1);
  else histArrScales = linspace(0, 0, numetabins/2 - 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts
  TCanvas* canvas = new TCanvas("trigCanvas", "", 800, 600);
  canvas->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTicks();
  gStyle->SetErrorX(0);
  canvas->Draw();
  for (int etabin = 0; etabin < numetabins; etabin++) {
    thisHist = histArr[etabin];
    for (int pbin = 0; pbin < numpbins; pbin++) {
      if (kinematicLumiVec[pbin + etabin*numpbins] != 0) {
        thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / (kinematicLumiVec[pbin + etabin*numpbins]));
        thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / (kinematicLumiVec[pbin + etabin*numpbins]));
      }
      else if (debugStatements) cout << "Warning: In IdealPtAnalysisHist.C (breakpoint J): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
    }
    deta = etabins[etabin+1] - etabins[etabin];
    thisHist->Scale(1e3*TMath::Power(10, histArrScales[(int)(numetabins/2 - 0.5 -TMath::Abs(etabin - numetabins/2 + 0.5))])/deta, "width"); // separate different etabins
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
    if (etabin == 0) thisHist->Draw("e1");
    else thisHist->Draw("same e1");
    
    const float textx = 0.47 + (etabin>=(numetabins/2))*0.25;
    const float texty = 0.91 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
    const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[etabin], etabins[etabin+1], histArrScales[(int)((0.5*(numetabins-1))-TMath::Abs(etabin-(0.5*(numetabins-1))))]);
    myMarkerText (textx, texty, kColor, kFullCircle, text);
  }

  myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity));
  myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

  string histName;
  if (numetabins > 1) {
    histName = "ptSpectra_combinedTriggers_etabinned";
  }
  else histName = "ptSpectra_combinedTriggers_0eta490";

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
  if (highPtJetsOnly) histName += "_highPtJetsOnly";
  if (runPeriodA || runPeriodB) canvas->SaveAs((plotPath + "ptSpectra/" + histName + ".pdf").c_str());


  /**** Plot pt spectrum ratios ****/
  if (runPeriodA && runPeriodB) {
    delete canvas;
    delete[] histArrScales;
    histArrScales = linspace(-2, 1, numetabins/2 - 1);
    canvas = new TCanvas("ratioCanvas", "", 800, 600);
    canvas->cd();
    //gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicks();
    gStyle->SetErrorX(0);
    canvas->Draw();

    double* periodALuminosities = new double[numpbins*numetabins];
    double* periodBLuminosities = new double[numpbins*numetabins];

    if (isMC) {
      for (int i = 0; i < numpbins*numetabins; i++) {
        periodALuminosities[i] = 0;
        periodBLuminosities[i] = 0;
      }
    } else {
      for (int i = 0; i < numpbins*numetabins; i++) {
        periodALuminosities[i] = 0;
        periodBLuminosities[i] = 0;
      }
      const double* allLuminosities = getTriggerLuminosities();
      
      Trigger* bestTrigger;
      double bestLumi;
      int thisRunNumber;
      for (int rnIndex = 0; rnIndex < numDataSets; rnIndex++) {
        thisRunNumber = (*dataSets)[rnIndex];
        /**** Need to get kinematic trigger vec for this run ****/
        setBestTriggers(rnIndex);

        if (isPeriodA(thisRunNumber)) {
          for (int pbin = 0; pbin < numpbins; pbin++) {
            for (int etabin = 0; etabin < numetabins; etabin++) {
              bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
              if (bestTrigger == NULL) continue;
              const int actetabin = numetabins - etabin - 1;
              bestLumi = allLuminosities[rnIndex + (bestTrigger->index)*numDataSets];
              periodALuminosities[pbin + actetabin*numpbins] += bestLumi;
            }
          }
        } else if (!isPeriodA(thisRunNumber)) {
          for (int pbin = 0; pbin < numpbins; pbin++) {
            for (int etabin = 0; etabin < numetabins; etabin++) {
              bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
              if (bestTrigger == NULL) continue;
              bestLumi = allLuminosities[rnIndex + (bestTrigger->index)*numDataSets];
              periodBLuminosities[pbin + etabin*numpbins] += bestLumi;
            }
          }
        }
      }
      delete[] allLuminosities;
    }

    TLine* lineDrawer = new TLine();
    lineDrawer->SetLineStyle(6);
    lineDrawer->SetLineWidth(2);
    
    for (int etabin = 0; etabin < numetabins; etabin++) {
      // Scale numerator by luminosity to get a cross section
      thisHist = numeratorHistArr[etabin];
      for (int pbin = 0; pbin < numpbins; pbin++) {
        double thisLumi = periodALuminosities[pbin + etabin*numpbins];
        if (thisLumi != 0) {
          thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / thisLumi);
          thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / thisLumi);
        }
        else if (debugStatements) cout << "Warning: In IdealPtAnalysisHist.C (breakpoint J): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
      }
      // Scale denominator by luminosity to get a cross section
      thisHist = denominatorHistArr[etabin];
      for (int pbin = 0; pbin < numpbins; pbin++) {
        double thisLumi = periodBLuminosities[pbin + etabin*numpbins];
        if (thisLumi != 0) {
          thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / thisLumi);
          thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / thisLumi);
        }
        else if (debugStatements) cout << "Warning: In IdealPtAnalysisHist.C (breakpoint J): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
      }
      // Now divide the histograms to get a ratio of the cross sections
      numeratorHistArr[etabin]->Divide(denominatorHistArr[etabin]);

      // Now go through the plotting routines
      thisHist = numeratorHistArr[etabin];
      double scale = histArrScales[(int)(numetabins/2 - 0.5 - TMath::Abs(etabin - numetabins/2 + 0.5))];
      for (int pbin = 0; pbin < numpbins; pbin++) {
        if (thisHist->GetBinContent(pbin+1) == 0.) thisHist->SetBinContent(pbin+1, -5);
        else thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) + scale);
      }

//    thisHist->SetMarkerStyle(kDot);
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
  
    myText (0.19, 0.91, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity));
    myText (0.19, 0.85, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

    histName = "ptSpectra_ratio_combinedTriggers_etabinned";
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
  etaPhiHist->GetXaxis()->SetTitleOffset(0.8);
  etaPhiHist->GetYaxis()->SetTitleOffset(0.8);
  etaPhiHist->GetZaxis()->SetTitleOffset(1.2);

  etaPhiHist->Scale(1., "width");
  etaPhiHist->Draw("colz");
  histName = "inclusive_jets_eta_phi_correlation";
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

  histName = "etaPhiHist";
  if (highPtJetsOnly) histName += "_highPtJetsOnly";
  TFile* output = new TFile((rootPath + histName + ".root").c_str(), "RECREATE");
  etaPhiHist->Write(); // save the eta phi distribution for further use
  subleadingEtaPhiHist->Write();
  output->Close();

  delete output;
  delete etaPhiHist;
  /**** End plot eta-phi diagram ****/


  /**** Plot eta-phi jet phase space diagram ****/
  delete canvas;
  canvas = new TCanvas("subleadingEtaPhiCanvas", "", 800, 600);
  canvas->cd();

  canvas->SetRightMargin(0.18);
  canvas->SetTopMargin(-0.02);
  subleadingEtaPhiHist->GetXaxis()->SetTitleOffset(0.8);
  subleadingEtaPhiHist->GetYaxis()->SetTitleOffset(0.8);
  subleadingEtaPhiHist->GetZaxis()->SetTitleOffset(1.2);

  subleadingEtaPhiHist->Scale(1., "width");
  subleadingEtaPhiHist->Draw("colz");
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

  delete subleadingEtaPhiHist;
  /**** End plot eta-phi diagram ****/


  /**** Free memory ****/
  for (int etabin = 0; etabin < numetabins; etabin++) {
    delete histArr[etabin];
    delete numeratorHistArr[etabin];
    delete denominatorHistArr[etabin];
  }
  delete lineDrawer;
  delete canvas;
  delete[] histArrScales;
  delete dataSets;
  /**** End free memory ****/


  if (debugStatements) cout << "Status: In IdealPtAnalysisHist.C (breakpoint K): Finished plotting pt spectrum" << endl;
  return;
}
