#ifndef __AnalyzeTrackMomentumResolution_C__
#define __AnalyzeTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <AtlasUtils.h>
#include <Utilities.h>

#include "Params.h"

using namespace atlashi;

typedef TGraphAsymmErrors TGAE;

const int numFinerEtachBins = 40;
const double* finerEtachBins = linspace (-2.5, 2.5, numFinerEtachBins);

const double etachBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int numEtachBins = sizeof (etachBins) / sizeof (etachBins[0]) - 1;

const double pTchBinsPP[30] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60};
const int numPtchBinsPP = sizeof (pTchBinsPP) / sizeof (pTchBinsPP[0]) - 1;
const double pTchBinsPbPb[28] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 15, 20, 30, 40, 60};
const int numPtchBinsPbPb = sizeof (pTchBinsPbPb) / sizeof (pTchBinsPbPb[0]) - 1;

const bool isPbPb = false;

TH1D**** h_tms = nullptr;
TH2D** h2_avg_tms = nullptr;
TH2D** h2_avg_tmr = nullptr;
TH1D*** h_avg_tms = nullptr;
TH1D*** h_avg_tmr = nullptr;


void AnalyzeTrackMomentumResolution () {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/TrackingMomentumResolution/Nominal/trackingMomentumResolution.root", rootPath.Data ()), "read");

  h_tms = new TH1D***[numCentBins];
  for (int iCent = 0; iCent < 1; iCent++) {
    h_tms[iCent] = new TH1D**[(isPbPb ? numPtchBinsPP : numPtchBinsPbPb)];
    for (int iPtch = 0; iPtch < (isPbPb ? numPtchBinsPP : numPtchBinsPbPb); iPtch++) {
      h_tms[iCent][iPtch] = new TH1D*[numFinerEtachBins];
      for (int iEta = 0; iEta < numFinerEtachBins; iEta++) {
        h_tms[iCent][iPtch][iEta] = (TH1D*) inFile->Get (Form ("h_tms_iCent%i_iPtch%i_iEta%i", iCent, iPtch, iEta));
      } // end loop over iEta
    } // end loop over iPtch
  } // end loop over iCent



  TFile* outFile = new TFile (Form ("%s/TrackingMomentumResolution/Nominal/trackingMomentumResolutionFactors.root", rootPath.Data ()), "recreate");

  h2_avg_tms = new TH2D*[numCentBins];
  h2_avg_tmr = new TH2D*[numCentBins];
  h_avg_tms = new TH1D**[numCentBins];
  h_avg_tmr = new TH1D**[numCentBins];

  for (int iCent = 0; iCent < 1; iCent++) {
    h2_avg_tms[iCent] = new TH2D (Form ("h2_avg_tms_iCent%i", iCent), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMS [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb), numFinerEtachBins, finerEtachBins);
    h2_avg_tmr[iCent] = new TH2D (Form ("h2_avg_tmr_iCent%i", iCent), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMR [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb), numFinerEtachBins, finerEtachBins);

    h_avg_tms[iCent] = new TH1D*[numEtachBins+1];
    h_avg_tmr[iCent] = new TH1D*[numEtachBins+1];
    for (int iEta = 0; iEta <= numEtachBins; iEta++) {
      h_avg_tms[iCent][iEta] = new TH1D (Form ("h_avg_tms_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T}^{truth} [GeV];TMS [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb));
      h_avg_tmr[iCent][iEta] = new TH1D (Form ("h_avg_tmr_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T}^{truth} [GeV];TMR [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb));
    } // end loop over iEta
  } // end loop over iCent



  for (int iCent = 0; iCent < 1; iCent++) {

    for (int iPtch = 0; iPtch < (isPbPb ? numPtchBinsPP : numPtchBinsPbPb); iPtch++) {

      TH1D** h_tms_integratedEta = new TH1D*[numEtachBins+1];
      for (int iEta = 0; iEta <= numEtachBins; iEta++) {
        h_tms_integratedEta[iEta] = new TH1D (Form ("h_tms_integratedEta_iEta%i", iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 200, 0.5, 1.5);
        h_tms_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < numFinerEtachBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtachBins[iFinerEta] + finerEtachBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < numEtachBins) {
          if (etachBins[iEta] < binCenter && binCenter < etachBins[iEta+1])
            break;
          iEta++;
        }

        h_tms_integratedEta[iEta]->Add (h_tms[iCent][iPtch][iFinerEta]);
        h_tms_integratedEta[numEtachBins]->Add (h_tms[iCent][iPtch][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.5, 1.5);
        fit->SetParameter (0, h_tms[iCent][iPtch][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_tms[iCent][iPtch][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", mean-4*sigma, mean+4*sigma);

        fit->SetParameter (0, h_tms[iCent][iPtch][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_tms[iCent][iPtch][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double tms = mean;
        const double tms_err = mean_err;
        
        const double tmr = sigma / mean;
        const double tmr_err = fabs (tmr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_avg_tms[iCent]->SetBinContent (iPtch+1, iFinerEta+1, tms * 100);
        h2_avg_tms[iCent]->SetBinError (iPtch+1, iFinerEta+1, tms_err * 100);
        h2_avg_tmr[iCent]->SetBinContent (iPtch+1, iFinerEta+1, tmr * 100);
        h2_avg_tmr[iCent]->SetBinError (iPtch+1, iFinerEta+1, tmr_err * 100);
      }

      for (int iEta = 0; iEta <= numEtachBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.5, 1.5);
        fit->SetParameter (0, h_tms_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_tms_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", mean-4*sigma, mean+4*sigma);

        fit->SetParameter (0, h_tms_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_tms_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double tms = mean;
        const double tms_err = mean_err;

        const double tmr = sigma / mean;
        const double tmr_err = fabs (tmr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_avg_tms[iCent][iEta]->SetBinContent (iPtch+1, tms * 100);
        h_avg_tms[iCent][iEta]->SetBinError (iPtch+1, tms_err * 100);
        h_avg_tmr[iCent][iEta]->SetBinContent (iPtch+1, tmr * 100);
        h_avg_tmr[iCent][iEta]->SetBinError (iPtch+1, tmr_err * 100);

        delete h_tms_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_tms_integratedEta;
    } // end loop over iPtch
  } // end loop over iCent


  for (int iCent = 0; iCent < 1; iCent++) {
    h2_avg_tms[iCent]->Write ();
    h2_avg_tmr[iCent]->Write ();

    for (int iEta = 0; iEta <= numEtachBins; iEta++) {
      h_avg_tms[iCent][iEta]->Write ();
      h_avg_tmr[iCent][iEta]->Write ();
    } // end loop over iEta
  } // end loop over iCent

  outFile->Close ();
  
}

#endif
