#include "InclusiveRpPbAnalysisHist.h"
#include "InclusiveRpPbAnalysis.h"
#include "Util.h"
#include "Params.h"
#include "TreeVariables.h"
#include <GlobalParams.h>

using namespace atlashi;

namespace JetAnalysis {

void InclusiveRpPbAnalysisHist() {

  TH1D** ppHistArr = setupPPConfiguration(); // creates the pp histograms and setups the relevant environment variables

  initialize(0, true);
  std::vector<int>* thisRunNumbers = getRunNumbers();

  const int numruns = (*thisRunNumbers).size();
  const int numhists = numtrigs * numruns * numetabins;
  if (debugStatements) {
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint A): Building trigger pt histograms with " << numruns << " runs being used" << endl;
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint B): Numtrigs = " << numtrigs << endl;
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint C): Numetabins = " << numppYbins << endl;
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint D): Numpbins = " << numpbins << endl;
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint E): ptPath = " << ptPath << endl;
  }

  const bool shiftBins = false;
  double ymin = 0;
  double ymax = 7.5;
  if (shiftBins) ymax = 15;
  const Style_t mkstyles[2] = {kFullCircle, kOpenCircle};
  const Color_t mkcolors[8] = {kOrange-3, kGreen, kCyan, kBlue, kViolet, kRed, kMagenta, kYellow};

  TH1D* thisHist;

  if (debugStatements) cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint F): Initializing pPbHistArr histograms..." << endl;
  TH1D* pPbHistArr[numpPbYbins];
  for (int pPbYbin = 0; pPbYbin < numpPbYbins; pPbYbin++) {
    int ppYbin_equiv = getppYbin(0.5*(pPbYbins[pPbYbin]+pPbYbins[pPbYbin+1]) - etaCoM);
    pPbHistArr[pPbYbin] = new TH1D(Form("pPb_spectrum_pPbYbin%i", pPbYbin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T} #left[nb GeV^{-1}#right]", numppPtbins, &(ppPtbins[ppYbin_equiv*(numppPtbins+1)]));
    pPbHistArr[pPbYbin]->Sumw2();
  }
  TH1D* RpPbHistArr[numppYbins];
  for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
    RpPbHistArr[ppYbin] = new TH1D(Form("R_pPb_ppYbin%i", ppYbin), ";#it{p}_{T}^{jet} [GeV];R_{#it{pPb}}", numppPtbins, &(ppPtbins[ppYbin*(numppPtbins+1)]));
    RpPbHistArr[ppYbin]->Sumw2();
  }

  if (debugStatements) {
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint G): pPbHistArr histograms initialized." << endl;
    cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint H): Starting loop over triggers..." << endl;
  }

  vector<int>* runNumbers = getRunNumbers();
  double totalpPbLuminosity = 0; // nb^-1
  double totalppLuminosity = 20.2; // fb^-1

  /**** Fill summed histograms with results from event loops ****/
  {
    TSystemDirectory dir(ptPath.c_str(), RpPbPath.c_str());
    TList* sysfiles = dir.GetListOfFiles();
    if (sysfiles) {
      TSystemFile *sysfile;
      TString fname;
      TString histName;
      TIter next(sysfiles);
      TVectorD* run_vec;
      TVectorD* lum_vec;

      while ((sysfile=(TSystemFile*)next())) {
        fname = sysfile->GetName();
        if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
          if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl; 
          for (int thisRunNumber : *runNumbers) {
            if (skipRun(thisRunNumber)) continue;
            if (fname.Contains(to_string(thisRunNumber))) {
              TFile* thisFile = new TFile(RpPbPath + fname, "READ");

              // quickly check the parameters stored in this root file
              run_vec = (TVectorD*)thisFile->Get(Form("run_vec_%i", thisRunNumber));
              assert ((*run_vec)[0] == thisRunNumber);
              assert ((*run_vec)[1] == numppYbins);
              assert ((*run_vec)[2] == numpPbYbins);
              assert ((*run_vec)[3] == numtrigs);
              assert ((*run_vec)[4] == numppPtbins);
              totalpPbLuminosity += (*run_vec)[5];

              for (int pPbYbin = 0; pPbYbin < numpPbYbins; pPbYbin++) {
                histName = Form("pPb_spectrum_run%i_pPbYbin%i", thisRunNumber, pPbYbin);
                pPbHistArr[pPbYbin]->Add((TH1D*)thisFile->Get(histName));
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


  /**** Add histograms to match pp binning in eta, then divide the pPb spectrum by the pp spectrum ****/
  for (int pPbYbin = 0; pPbYbin < numpPbYbins; pPbYbin++) {
    const int ppYbin = getppYbin(0.5*(pPbYbins[pPbYbin] + pPbYbins[pPbYbin+1]) - etaCoM); // takes the bin center, then shifts "out" of the CoM frame so it looks like we're in the lab frame, then finds the corresponding pp y bin
    RpPbHistArr[ppYbin]->Add(pPbHistArr[pPbYbin]);
    //cout << "Adding pPbYbin " << pPbYbins[pPbYbin] << ", " << pPbYbins[pPbYbin+1] << " to ppYbin " << ppYbins[ppYbin] << ", " << ppYbins[ppYbin+1] << endl;
  }
  for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
    const double deta = ppYbins[ppYbin+1] - ppYbins[ppYbin]; // always = 0.5
    RpPbHistArr[ppYbin]->Scale(1e3/deta, "width"); // convert nb to pb
    RpPbHistArr[ppYbin]->Divide(ppHistArr[ppYbin]); // now divide by the pp spectrum
  }
  /**** End add histograms ****/


  /**** Plot imported pp spectrum ****/
  TCanvas* canvas = new TCanvas("ppCanvas", "", 800, 600);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetTicks();
  for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
    thisHist = ppHistArr[ppYbin];
    const Color_t kColor = mkcolors[TMath::Abs(numppYbins/2 - (ppYbin+(ppYbin>= numppYbins/2)))%(numppYbins/2)];
    const Style_t kStyle = mkstyles[ppYbin >= numppYbins/2];
    thisHist->SetMarkerStyle(kStyle);
    thisHist->SetMarkerColor(kColor);
    thisHist->SetLineColor(kColor);
    thisHist->SetMinimum(1e-8);
    thisHist->SetMaximum(1e5);
    //thisHist->GetXaxis()->SetRangeUser(60, 2550);
    if (ppYbin == 0) thisHist->Draw("e1");
    else thisHist->Draw("same e1");
    thisHist->GetXaxis()->SetRangeUser(60, 2550);
  }
  string histName = "pp_spectrum";
  if (runPeriodA || runPeriodB) canvas->SaveAs((plotPath + "R_pPb/" + histName + ".pdf").c_str());
  /**** End plot pp spectrum ****/


  /**** Plot calculated R_pPb histograms ****/
  const double histArrShifts[numppYbins] = {};
  //const double histArrShifts[numppYbins] = {0, 1.5, 3, 4.5, 6, 7.5};
  //const double histArrShifts[numppYbins] = {7.5, 6, 4.5, 3, 1.5, 0};
  const int canvasIndices[numppYbins] = {1, 5, 9, 10, 6, 2, 3, 7, 11, 12, 8, 4};
  delete canvas;
  canvas = new TCanvas("RpPbCanvas", "", 2000, 1500);
  canvas->Divide(4, 3, 0.0, 0.0);
  TLine* lineDrawer = new TLine();  
  gPad->SetTicks();
  canvas->Draw();
  for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
    canvas->cd(canvasIndices[ppYbin]);
    //canvas->cd((4*ppYbin)%numppYbins+(4*ppYbin)/numppYbins+1);
    gPad->SetLogx();
    thisHist = RpPbHistArr[ppYbin];
    const Style_t kStyle = mkstyles[0];
    const Color_t kColor = mkcolors[TMath::Abs(numppYbins/2 - ppYbin - 1) - (ppYbin>=numppYbins/2)];

    if (shiftBins) {
      for (int ppPtbin = 0; ppPtbin < numppPtbins; ppPtbin++) {
        RpPbHistArr[ppYbin]->SetBinContent(ppPtbin+1, RpPbHistArr[ppYbin]->GetBinContent(ppPtbin+1) + histArrShifts[ppYbin]); // shifts each histogram appropriately
      }
      lineDrawer->SetLineColor(kColor);
      lineDrawer->SetLineStyle(9);
      //lineDrawer->DrawLine(ppPtbins[0], histArrShifts[ppYbin]+1, ppPtbins[numppPtbins], histArrShifts[ppYbin]+1);
    }

    thisHist->SetMarkerStyle(kStyle);
    thisHist->SetMarkerColor(kColor);
    thisHist->SetLineColor(kColor);
    thisHist->SetMinimum(ymin);
    thisHist->SetMaximum(ymax);
    thisHist->GetYaxis()->SetLabelSize(0.07);
    thisHist->GetYaxis()->SetTitleSize(0.07);
    thisHist->GetXaxis()->SetLabelSize(0.07);
    thisHist->GetXaxis()->SetTitleSize(0.07);

    thisHist->GetYaxis()->SetTitleOffset(1.1);
    thisHist->GetXaxis()->SetTitleOffset(1.1);
    thisHist->GetXaxis()->SetTickLength(0.02);
    thisHist->GetYaxis()->SetTickLength(0.02);
    thisHist->SetAxisRange(60, 2550, "X");
    if (ppYbin == 0) thisHist->Draw("e1");
    else thisHist->Draw("same e1");

    if (shiftBins) lineDrawer->DrawLine(ppPtbins[0], histArrShifts[ppYbin]+1, ppPtbins[numppPtbins], histArrShifts[ppYbin]+1);

    const float textx = 0.26; //- 0.04*shiftBins;
    const float texty = 0.85; // - (ppYbin)*0.05;
    char* text;
    if (shiftBins) text = Form("%g < #it{y}_{CoM} < %g (+%g)", ppYbins[ppYbin], ppYbins[ppYbin+1], histArrShifts[ppYbin]); /*histArrScales[(int)((0.5*(numppYbins-1))-TMath::Abs(ppYbin-(0.5*(numppYbins-1))))]);*/
    else text = Form("%g < #it{y}_{CoM} < %g", ppYbins[ppYbin], ppYbins[ppYbin+1]);
    myText (textx, texty, kBlack, text,  0.1);
  }

  canvas->cd(12);
  myText (0.04, 0.325, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{#it{s}} = 8.16 TeV", totalpPbLuminosity), 0.061);
  myText (0.04, 0.225, kBlack, Form("2012 #it{pp}, %.1f fb^{-1}, #sqrt{#it{s}} = 8 TeV", totalppLuminosity), 0.061);

  histName = "R_pPb_combinedTriggers";

  if (shiftBins) histName = histName + "_shifted";

  if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
  else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

  if (runPeriodA || runPeriodB) canvas->SaveAs((plotPath + "R_pPb/" + histName + ".pdf").c_str());
  /**** End plot R_pPb ****/


  /**** Free memory and quit ****/
  for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
    delete pPbHistArr[ppYbin];
    delete RpPbHistArr[ppYbin];
  }
  delete canvas;
//    delete[] histArrScales;
  delete runNumbers;

  if (debugStatements) cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint J): Finished plotting pt spectrum" << endl;
  return;
}

} // end namespace
