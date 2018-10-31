#include "TriggerPtAnalysisHist.h"
#include "Util.h"
#include "Params.h"
#include <GlobalParams.h>

using namespace atlashi;

namespace JetAnalysis {

void TriggerPtAnalysisHist(int thisRunNumber) {

  if (skipRun(thisRunNumber)) return;

  initialize(thisRunNumber, false);

  /**** Generate list of physics triggers ****/
  vector<Trigger*> triggerSubList(0);
  for (Trigger* trig : triggerVec) {
    if (trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber && trig->name != minbiasTriggerName) triggerSubList.push_back(trig);
  }
  if (debugStatements) {
    cout << "Status: In TriggerPtAnalysisHist.C (breakpoint A): Processing run " << thisRunNumber << " with triggers:" << endl;
    for (Trigger* trig : triggerSubList) {
      cout << "\t" << trig->name << endl;
    }
  }
  const int numtrigs_sublist = triggerSubList.size();
  const double* trigbins = linspace(0, numtrigs_sublist+1, numtrigs_sublist+1);
  double* triggerLuminosities = getTriggerLuminosities();
  /**** End generate list of physics triggers ****/

  double minval = 1e0;
  double maxval;
  const double textangle = 30.;

  TFile* thisFile = new TFile(Form("%srun_%i.root", trigPath.c_str(), thisRunNumber), "READ");
  TVectorD* lum_vec = (TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber));
  const double luminosity = (*lum_vec)[0];

  int rnIndex = 0;
  int numruns;
  {
    vector<int>* runNumbers = getRunNumbers();
    numruns = runNumbers->size();
    while ((*runNumbers)[rnIndex] != thisRunNumber) rnIndex++;
    delete runNumbers;
  }

  TH1D* histArr[numtrigs_sublist];
  for (int t = 0; t < numtrigs_sublist; t++) {
    Trigger* trig = triggerSubList[t];
    TString histName = Form("integratedCounts_%s_run%i", trig->name.c_str(), thisRunNumber);
    histArr[t] = (TH1D*)thisFile->Get(histName);
    double thisLumi = 0.;
    for (int pbin = 0; pbin < numpbins && thisLumi == 0.; pbin++) {
      for (int etabin = 0; etabin < numetabins && thisLumi == 0.; etabin++) {
        thisLumi = triggerLuminosities[rnIndex + (trig->index + (pbin + etabin*numpbins)*numtrigs)*numruns];
      }
    }
    if (thisLumi != 0.) histArr[t]->SetBinContent(t+1, 1e3*histArr[t]->GetBinContent(t+1)/thisLumi);
  }
  
  TCanvas* canvas = new TCanvas(Form("canvas_%i", thisRunNumber), "", 800, 600);
  canvas->cd();

  TLatex* description = new TLatex();
  description->SetTextAlign(22);
  description->SetTextFont(42);
  description->SetTextSize(0.024);

  gPad->SetLogy();

  const Style_t mkstyles[8] = {kDot, kDot, kDot, kDot, kDot, kDot, kDot, kDot};
  const Color_t mkcolors[20] = {30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49};

  double maxbincontent = 0;
  {
    double thisbincontent;
    for (int i = 0; i < numtrigs_sublist; i++) {
      thisbincontent = histArr[i]->GetBinContent(i+1);
      if (maxbincontent < thisbincontent) maxbincontent = thisbincontent;
    }
    maxval = 6 * maxbincontent;
  }

  int numticks = 0;
  for (int t = 0; t < numtrigs_sublist; t++) {
    numticks += histArr[t]->Integral();
    histArr[t]->SetMarkerStyle(mkstyles[t%8]);
    histArr[t]->SetMarkerColor(mkcolors[t%20]);
    histArr[t]->SetLineColor(mkcolors[t%20]);
    histArr[t]->SetFillColor(mkcolors[t%20]);
    histArr[t]->SetMinimum(minval);
    histArr[t]->SetMaximum(maxval);
    histArr[t]->GetXaxis()->SetLabelOffset(999);
    histArr[t]->GetXaxis()->SetLabelSize(0);
    histArr[t]->GetXaxis()->SetTickLength(0);

    histArr[t]->GetYaxis()->SetTitleOffset(1.3);
    histArr[t]->GetYaxis()->SetTickLength(0.0075);

    if (t == 0) histArr[t]->Draw("BAR");
    else histArr[t]->Draw("BAR, SAME");

    TLatex* text = description->DrawLatex(t+0.5, TMath::Power(10, TMath::Log10(minval) + 0.01*(TMath::Log10(maxval)-TMath::Log10(minval))), triggerSubList[t]->name.c_str());
    text->SetTextAngle(90);
    text->SetTextAlign(12);
  }

  canvas->Draw();

  myText (0.46, 0.9, kBlack, Form("Run %i, %.1f nb^{-1}, #sqrt{#it{s}} = 8.16 TeV", thisRunNumber, luminosity));
  myText (0.7, 0.84, kBlack, Form("N^{total}_{counts} = %i", numticks)); 

  canvas->SaveAs(Form("%scounts/run_trig_%i.pdf", plotPath.c_str(), thisRunNumber));
  if (debugStatements) cout << Form("Status: In TriggerPtAnalysisHist.C (breakpoint B): Triggers for run number %i finished", thisRunNumber) << endl;

  thisFile->Close();
  delete [] triggerLuminosities;
  return;
}

} // end namespace
