#include "ZMassCalcHist.h"
#include "Params.h"
#include "Utils.h"

#include <GlobalParams.h>
#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace pPb8TeV2016JetCalibration {

void ZMassCalcHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("ZMassCalc/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back (TString ("Pbp") + (runValidation ? "_Signal":"_Overlay") + "_GammaJet_Slice" + to_string (i+1));
   gammaJetSampleIds.push_back (TString ("pPb") + (runValidation ? "_Signal":"_Overlay") + "_GammaJet_Slice" + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds (0);
  dijetSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  TH1D**** zMassSpectra = Get3DArray <TH1D*> (2, 2, numzetabins+1);

  for (short iEta = 0; iEta <= numzetabins; iEta++) {
   for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
    const TString dataType = (iData == 0 ? "data":"mc");
    for (short iSpc = 0; iSpc < 2; iSpc++) {
     const TString species = (iSpc == 0 ? "mumu":"ee");

     zMassSpectra[iSpc][iData][iEta] = new TH1D (Form ("z%sMassSpectrum_%s_iEta%i", species.Data (), dataType.Data (), iEta), "", 50, 60, 110);
     zMassSpectra[iSpc][iData][iEta]->Sumw2 ();
    }
   }
  }

  int*** nZeeMass = Get3DArray <int> (3, 2, numzetabins+1);
  int*** nZmumuMass = Get3DArray <int> (3, 2, numzetabins+1);

  {
   TSystemDirectory dir (rootPath.Data (), rootPath.Data ());
   TList* sysfiles = dir.GetListOfFiles ();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile *sysfile;
   TString fname;
   TString histName;
   TIter next (sysfiles);
   TVectorD *nZeeMassVec, *nZmumuMassVec; 
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);
       //infoVec = (TVectorD*)thisFile->Get (Form ("infoVec_%i", runNumber));
       nZeeMassVec = (TVectorD*)thisFile->Get (Form ("nZeeMassVec_%i", runNumber));
       nZmumuMassVec = (TVectorD*)thisFile->Get (Form ("nZmumuMassVec_%i", runNumber));

       for (short iEta = 0; iEta <= numzetabins; iEta++) {
        //const bool flipEta = runNumber < 313500 && iEta < numzetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? (numzetabins - iEta - 1) : iEta);

        nZeeMass[iPer][0][iEta] += (*nZeeMassVec)[iEta];
        nZeeMass[2][0][iEta] += (*nZeeMassVec)[act_iEta];

        nZmumuMass[iPer][0][iEta] += (*nZmumuMassVec)[iEta];
        nZmumuMass[2][0][iEta] += (*nZmumuMassVec)[act_iEta];
       }

       for (short iSpc = 0; iSpc < 2; iSpc++) {
        const TString species = (iSpc == 0 ? "mumu":"ee");

        for (short iEta = 0; iEta <= numzetabins; iEta++) {
         zMassSpectra[iSpc][0][iEta]->Add ( (TH1D*)thisFile->Get (Form ("z%sMassSpectrum_dataSet%i_data_iEta%i", species.Data (), runNumber, iEta)));
        }
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (zeeJetSampleId.Contains ("pPb") ? 0 : 1);
       nZeeMassVec = (TVectorD*)thisFile->Get (Form ("nZeeMassVec_%s", zeeJetSampleId.Data ()));

       for (short iEta = 0; iEta <= numzetabins; iEta++) {
        //const bool flipEta = zeeJetSampleId.Contains ("pPb") && iEta < numzetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? numzetabins - iEta - 1 : iEta); // period A condition

        nZeeMass[iPer][1][iEta] += (*nZeeMassVec)[iEta];
        nZeeMass[2][1][iEta] += (*nZeeMassVec)[act_iEta];
       }

       for (short iEta = 0; iEta <= numzetabins; iEta++) {
        zMassSpectra[1][1][iEta]->Add ( (TH1D*)thisFile->Get (Form ("zeeMassSpectrum_dataSet%s_mc_iEta%i", zeeJetSampleId.Data (), iEta)));
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if Z->mumu sample
     for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
      if (fname.Contains (zmumuJetSampleId)) { // if Z->mumu sample do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (zmumuJetSampleId.Contains ("pPb") ? 0 : 1);
       nZmumuMassVec = (TVectorD*)thisFile->Get (Form ("nZmumuMassVec_%s", zmumuJetSampleId.Data ()));

       for (short iEta = 0; iEta <= numzetabins; iEta++) {
        //const bool flipEta = zmumuJetSampleId.Contains ("pPb") && iEta < numzetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? numzetabins - iEta - 1 : iEta); // period A condition

        nZmumuMass[iPer][1][iEta] += (*nZmumuMassVec)[iEta];
        nZmumuMass[2][1][iEta] += (*nZmumuMassVec)[act_iEta];
       }

       for (short iEta = 0; iEta <= numzetabins; iEta++) {
        zMassSpectra[0][1][iEta]->Add ( (TH1D*)thisFile->Get (Form ("zmumuMassSpectrum_dataSet%s_mc_iEta%i", zmumuJetSampleId.Data (), iEta)));
       }

       thisFile->Close ();
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
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  const double padRatio = 1.5; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetBottomMargin (0);
  topPad->SetLeftMargin (-0.20);
  bottomPad->SetTopMargin (0);
  bottomPad->SetBottomMargin (0.30);
  bottomPad->SetLeftMargin (-0.20);
  topPad->Draw ();
  bottomPad->Draw ();


  /**** create new lines for Z mass spectra ****/
  TLine* lines[5] = {};
  for (short i = 0; i < 5; i++) {
   lines[i] = new TLine (60, 0.6+0.2*i, 110, 0.6+0.2*i);
   if (0.6+0.2*i == 1) lines[i]->SetLineStyle (1);
   else lines[i]->SetLineStyle (3);
  }


  /**** Plot dilepton mass spectra ****/
  for (short iEta = 0; iEta <= numzetabins; iEta++) {
   for (short iSpc = 0; iSpc < 2; iSpc++) {
    topPad->cd ();
    topPad->SetLogx (0);
    double mean[2] = {};
    double mean_err[2] = {};
    double sigma[2] = {};
    double sigma_err[2] = {};
    TF1* fits[2];
    //RooZfit* fits[2];
    for (short iData = 0; iData < 2; iData++) {
     TH1D* thisHist = zMassSpectra[iSpc][iData][iEta];
     Color_t color = (iData==0 ? dataColor : mcOverlayColor);
     thisHist->GetXaxis ()->SetTitle ("#font[12]{ll} Invariant Mass #left[GeV#right]");
     thisHist->GetYaxis ()->SetTitle ("Counts / Area under Curve");
     thisHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     thisHist->GetYaxis ()->SetTitleOffset (1.1*uPadY);
     thisHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     thisHist->SetLineColor (color);
     thisHist->SetMarkerColor (color);
//     thisHist->GetXaxis ()->SetNdivisions (50802, false);
     thisHist->Scale (1.0, "width");

     //TF1* gausFit = new TF1 (Form ("fit_guess_iSpc%i_iData%i", iSpc, iData), "gaus (0)", Z_mass - 10, Z_mass + 10);
     //thisHist->Fit (gausFit, "R", "L");
     //double m = gausFit->GetParameter (1);
     //double s = gausFit->GetParameter (2);
     //const double scale = 1.0 / thisHist->Integral (thisHist->FindBin (m-Z_mass_fitNsigma*s), thisHist->FindBin (m+Z_mass_fitNsigma*s));
     //thisHist->Scale (scale);
     //
     //RooZfit* fit = new RooZfit (thisHist, color);
     //fits[iData] = fit;

     //mean[iData] = fit->mean->getVal ();
     //mean_err[iData] = fit->mean->getError ();
     //sigma[iData] = fit->sigma->getVal ();
     //sigma_err[iData] = fit->sigma->getError ();

     TF1* fit = new TF1 (Form ("fit_guess_iSpc%i_iData%i", iSpc, iData), "gaus (0)", Z_mass - 10, Z_mass + 10);
     thisHist->Fit (fit, "R", "L");
     double m = fit->GetParameter (1);
     double s = fit->GetParameter (2);

     TF1* fit_better = new TF1 (Form ("fit_better_iSpc%i_iData%i", iSpc, iData), "gaus (0)", m-Z_mass_fitNsigma*s, m+Z_mass_fitNsigma*s);
     thisHist->Fit (fit_better, "R", "L");
     m = fit_better->GetParameter (1);
     s = fit_better->GetParameter (2);
     mean[iData] = m;
     mean_err[iData] = fit_better->GetParError (1);
     sigma[iData] = s;
     sigma_err[iData] = fit_better->GetParError (2);
     const double scale = 1.0 / thisHist->Integral (thisHist->FindBin (m-Z_mass_fitNsigma*s), thisHist->FindBin (m+Z_mass_fitNsigma*s));
     thisHist->Scale (scale);

     fits[iData] = fit_better;
     fit_better->SetLineColor (color);
     fit_better->SetParameter (0, scale*fit_better->GetParameter (0));

     thisHist->SetAxisRange (0, 0.20, "Y");
     thisHist->GetYaxis ()->ChangeLabel (1, -1, -1, -1, -1, -1, " ");
    }
    for (short iData = 0; iData < 2; iData++) {
     TH1D* thisHist = zMassSpectra[iSpc][iData][iEta];
     if (iData == 0) thisHist->Draw ("p");
     else thisHist->Draw ("hist same");
     //fits[iData]->plot->Draw ("same");
     fits[iData]->Draw ("same");
     if (iEta < numzetabins) {
      if (iSpc == 0)
       myText (0.175, 0.15, kBlack, Form ("%g < #eta_{Lab}^{Z} < %g", zetabins[iEta], zetabins[iEta+1]), 0.04/uPadY);
      else if (iSpc == 1)
       myText (0.175, 0.15, kBlack, Form ("%g < #eta_{Lab}^{Z} < %g", zetabins[iEta], zetabins[iEta+1]), 0.04/uPadY);
     }
    }
    if (iSpc == 0) {
     myText (0.175, 0.88, kBlack, "Z (#mu#mu) events", 0.04/uPadY);
     myMarkerText (0.175, 0.80, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV (%i events)", nZmumuMass[2][0][iEta]), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.55, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZmumuMass[2][1][iEta]), 1.25, 0.04/uPadY);
    }
    else if (iSpc == 1) {
     myText (0.175, 0.88, kBlack, "Z (ee) events", 0.04/uPadY);
     myMarkerText (0.175, 0.80, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV (%i events)", nZeeMass[2][0][iEta]), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.55, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZeeMass[2][1][iEta]), 1.25, 0.04/uPadY);
    }
    myText (0.175, 0.72, kBlack, Form ("m_{Z}^{data} = %.2f #pm %.2f GeV", mean[0], mean_err[0]), 0.04/uPadY);
    myText (0.175, 0.64, kBlack, Form ("#sigma_{Z}^{data} = %.2f #pm %.2f GeV", sigma[0], sigma_err[0]), 0.04/uPadY);
    myText (0.175, 0.47, kBlack, Form ("m_{Z}^{mc} = %.2f #pm %.2f GeV", mean[1], mean_err[1]), 0.04/uPadY);
    myText (0.175, 0.39, kBlack, Form ("#sigma_{Z}^{mc} = %.2f #pm %.2f GeV", sigma[1], sigma_err[1]), 0.04/uPadY);

    bottomPad->cd ();
    bottomPad->SetLogx (0);
    TH1D* thisHist = (TH1D*)zMassSpectra[iSpc][0][iEta]->Clone (Form ("invMass_iSpc%i_clone", iSpc));
    thisHist->Divide (zMassSpectra[iSpc][1][iEta]);
    thisHist->GetXaxis ()->SetTitle (Form ("#font[12]{%s} Invariant Mass #left[GeV#right]", (iSpc == 1 ? "ee":"#mu#mu")));
    thisHist->GetYaxis ()->SetTitle ("Data / MC");
    thisHist->GetXaxis ()->SetTitleSize (0.04/dPadY);
    thisHist->GetYaxis ()->SetTitleSize (0.04/dPadY);
    thisHist->GetXaxis ()->SetTitleOffset (1);
    thisHist->GetYaxis ()->SetTitleOffset (1.1*dPadY);
    thisHist->GetYaxis ()->CenterTitle (true);
    thisHist->GetXaxis ()->SetLabelSize (0.04/dPadY);
    thisHist->GetYaxis ()->SetLabelSize (0.04/dPadY);
    thisHist->SetLineColor (kBlack);
    thisHist->SetMarkerColor (kBlack);
    //thisHist->SetAxisRange (0.6, 1.6, "Y");
    thisHist->SetAxisRange (0.2, 2.2, "Y");

    thisHist->GetYaxis ()->ChangeLabel (-1, -1, -1, -1, -1, -1, " ");
    thisHist->GetXaxis ()->SetTickLength (0.08);
//    thisHist->GetXaxis ()->SetNdivisions (50802, false);
    thisHist->GetYaxis ()->SetNdivisions (405, false);
    thisHist->Draw ("p");
    for (TLine* line : lines) line->Draw ("same");

    if (iEta != numzetabins) {
     if (iSpc == 0) canvas->SaveAs (Form ("%s/zmumu_mass_comparison_iEta%i.pdf", plotPath.Data (), iEta));
     else if (iSpc == 1) canvas->SaveAs (Form ("%s/zee_mass_comparison_iEta%i.pdf", plotPath.Data (), iEta));
    }
    else {
     if (iSpc == 0) canvas->SaveAs (Form ("%s/zmumu_mass_comparison.pdf", plotPath.Data ()));
     else if (iSpc == 1) canvas->SaveAs (Form ("%s/zee_mass_comparison.pdf", plotPath.Data ()));
    }

    for (short fit = 0; fit < 2; fit++) if (fits[fit]) delete fits[fit];
   }
  }

  return;
}

} // end namespace
