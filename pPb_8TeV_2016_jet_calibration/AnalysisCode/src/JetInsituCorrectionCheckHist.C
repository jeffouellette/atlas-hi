#include "JetInsituCorrectionCheckHist.h"
#include "Params.h"

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
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

void JetInsituCorrectionCheckHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("JetInsituCorrectionCheck/", "pPb_8TeV_2016_jet_calibration/");

  // Setup lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);

  TH2D* jetPreInsituSpectrum = new TH2D ("jetPreInsituSpectrum", ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPreInsituSpectrum->Sumw2 ();
  TH2D* jetPostInsituSpectrum = new TH2D ("jetPostInsituSpectrum", ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPostInsituSpectrum->Sumw2 ();
  TH2D* jetPostInsituSpectrumSysHi = new TH2D ("jetPostInsituSpectrumSysHi", ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPostInsituSpectrumSysHi->Sumw2 ();
  TH2D* jetPostInsituSpectrumSysLo = new TH2D ("jetPostInsituSpectrumSysLo", ";#it{p}_{T}^{Pre-Insitu} #left[GeV#right];#eta;", 66, 20, 350, numetabins, etabins);
  jetPostInsituSpectrumSysLo->Sumw2 ();

  TH2D* jetInsituResponse[numpbins]; 
  TH2D* jetInsituResponseSysHi[numpbins];
  TH2D* jetInsituResponseSysLo[numpbins];

  for (short iP = 0; iP < numpbins; iP++) {
   jetInsituResponse[iP] = new TH2D (Form ("jetInsituResponse_iP%i", iP), ";#it{p}_{T}^{Insitu} / #it{p}_{T}^{EtaJES};#eta;", 200, 0.9, 1.1, numetabins, etabins);
   jetInsituResponse[iP]->Sumw2 ();
   jetInsituResponseSysHi[iP] = new TH2D (Form ("jetInsituResponseSysHi_iP%i", iP), ";#it{p}_{T}^{Insitu} / #it{p}_{T}^{EtaJES};#eta;", 200, 0.9, 1.1, numetabins, etabins);
   jetInsituResponseSysHi[iP]->Sumw2 ();
   jetInsituResponseSysLo[iP] = new TH2D (Form ("jetInsituResponseSysLo_iP%i", iP), ";#it{p}_{T}^{Insitu} / #it{p}_{T}^{EtaJES};#eta;", 200, 0.9, 1.1, numetabins, etabins);
   jetInsituResponseSysLo[iP]->Sumw2 ();
  }
 

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
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data () << endl;

     // do this if is a data sample
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       jetPreInsituSpectrum->Add ( (TH2D*)thisFile->Get (Form ("jetPreInsituSpectrum_dataSet%i", runNumber)));
       jetPostInsituSpectrum->Add ( (TH2D*)thisFile->Get (Form ("jetPostInsituSpectrum_dataSet%i", runNumber)));
       jetPostInsituSpectrumSysHi->Add ( (TH2D*)thisFile->Get (Form ("jetPostInsituSpectrumSysHi_dataSet%i", runNumber)));
       jetPostInsituSpectrumSysLo->Add ( (TH2D*)thisFile->Get (Form ("jetPostInsituSpectrumSysLo_dataSet%i", runNumber)));

       for (short iP = 0; iP < numpbins; iP++) {
        jetInsituResponse[iP]->Add ( (TH2D*)thisFile->Get (Form ("jetInsituResponse_dataSet%i_iP%i", runNumber, iP)));
        jetInsituResponseSysHi[iP]->Add ( (TH2D*)thisFile->Get (Form ("jetInsituResponseSysHi_dataSet%i_iP%i", runNumber, iP)));
        jetInsituResponseSysLo[iP]->Add ( (TH2D*)thisFile->Get (Form ("jetInsituResponseSysLo_dataSet%i_iP%i", runNumber, iP)));
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


  /**** Define local histograms, graphs, etc. ****/
  TFile* outFile = new TFile (rootPath + "JetInsituCorrectionCheck.root", "recreate");


  /**** Canvas definitions ****/
  TCanvas* jetInsituSpectrumCanvas = new TCanvas ("jetInsituSpectrumCanvas", "", 530 * (numetabins/2), 1200);
  jetInsituSpectrumCanvas->Draw ();

  const double padYwidth = 0.5;
  const double padRatio = 1.8; // ratio of size of upper pad to lower pad
  const int nx = numetabins/2;
  //const int ny = 4;
  double* lPadXbounds = new double[nx];
  double* uPadXbounds = new double[nx];
  for (short i = 0; i < nx-1; i++) {
   lPadXbounds[i+1] = 0.04 + (0.955/nx) * (i+1); 
   uPadXbounds[i] = 0.04 + (0.955/nx) * (i+1);
  }
  lPadXbounds[0] = 0;
  uPadXbounds[nx-1] = 1.;

  double* padLeftMargins = new double[nx];
  double* padRightMargins = new double[nx];
  padLeftMargins[0] = 0.04/lPadXbounds[1];
  padRightMargins[nx-1] = 0.005/ (1.-lPadXbounds[nx-1]);

//  const double lPadXbounds[3] = {0, 0.39, 0.69};
//  const double uPadXbounds[3] = {0.39, 0.69, 1.0};
  const double lPadYbounds[4] = {1-padYwidth/padRatio, 1-padYwidth, padYwidth* (1-1/padRatio), 0};
  const double uPadYbounds[4] = {1, 1-padYwidth/padRatio, padYwidth, padYwidth* (1-1/padRatio)};
  const double padTopMargins[4] = {0.02, 0, 0.02, 0};
  const double padBottomMargins[4] = {0, 0.32, 0, 0.32};

  TPad*** pads = new TPad**[nx]; //[nx][4];
  for (int ix = 0; ix < nx; ix++) {
   pads[ix] = new TPad*[4];
   const double lPadX = lPadXbounds[ix];
   const double uPadX = uPadXbounds[ix];
   for (int iy = 0; iy < 4; iy++) {
    const double lPadY = lPadYbounds[iy];
    const double uPadY = uPadYbounds[iy];
    pads[ix][iy] = new TPad (Form ("pad_%i_%i", ix, iy), "", lPadX, lPadY, uPadX, uPadY);
    pads[ix][iy]->Draw ();

    pads[ix][iy]->SetLeftMargin (padLeftMargins[ix]);
    pads[ix][iy]->SetRightMargin (padRightMargins[ix]);
    pads[ix][iy]->SetBottomMargin (padBottomMargins[iy]);
    pads[ix][iy]->SetTopMargin (padTopMargins[iy]);

    //if (iy%2 == 0) {
    // pads[ix][iy]->SetBottomMargin (0);
    // pads[ix][iy]->SetTopMargin (0.08);
    //}
    //else {
    // pads[ix][iy]->SetTopMargin (0);
    // pads[ix][iy]->SetBottomMargin (0.08);
    //}
    //if (ix != 0) pads[ix][iy]->SetLeftMargin (0);
    //else pads[ix][iy]->SetLeftMargin (0.07/0.31);
    //if (ix != 2) pads[ix][iy]->SetRightMargin (0);
    //else pads[ix][iy]->SetRightMargin (0.01/0.31);
   }
  }

  for (int iEta = 0; iEta < numetabins; iEta++) {
   const int ix = (iEta < nx ? (nx-iEta-1) : (iEta-nx));
   const int iy = 2* (iEta/nx);

   pads[ix][iy]->cd ();
   pads[ix][iy]->SetLogy (true);
   pads[ix][iy]->SetLogx (true);

   TH1D* preSpectrum = jetPreInsituSpectrum->ProjectionX (Form ("preSpectrumProjX_%i", iEta), iEta+1, iEta+1);
   TH1D* postSpectrum = jetPostInsituSpectrum->ProjectionX (Form ("postSpectrumProjX_%i", iEta), iEta+1, iEta+1);
   TH1D* postSpectrumSysHi = jetPostInsituSpectrumSysHi->ProjectionX (Form ("postSpectrumSysHiProjX_%i", iEta), iEta+1, iEta+1);
   TH1D* postSpectrumSysLo = jetPostInsituSpectrumSysLo->ProjectionX (Form ("postSpectrumSysLoProjX_%i", iEta), iEta+1, iEta+1);
   preSpectrum->Rebin (2);
   postSpectrum->Rebin (2);
   postSpectrumSysHi->Rebin (2);
   postSpectrumSysLo->Rebin (2);

   preSpectrum->Scale (1., "width");
   postSpectrum->Scale (1., "width");
   postSpectrumSysHi->Scale (1., "width");
   postSpectrumSysLo->Scale (1., "width");

   preSpectrum->SetLineColor (kBlack);
   preSpectrum->SetMarkerColor (kBlack);
   preSpectrum->SetAxisRange (2e0, 2e5, "Y");
   preSpectrum->GetYaxis ()->SetTitle ("Counts / GeV");
   preSpectrum->GetYaxis ()->SetTitleOffset (0.9);
   preSpectrum->GetXaxis ()->SetTitleSize (0.03/ (padYwidth/padRatio));
   preSpectrum->GetYaxis ()->SetTitleSize (0.03/ (padYwidth/padRatio));
   preSpectrum->GetXaxis ()->SetLabelSize (0.03/ (padYwidth/padRatio));
   preSpectrum->GetYaxis ()->SetLabelSize (0.03/ (padYwidth/padRatio));

   postSpectrum->SetLineColor (38);
   postSpectrum->SetMarkerColor (38);
   postSpectrum->SetMarkerStyle (7);

   TGraphAsymmErrors* postSys = new TGraphAsymmErrors (postSpectrum->GetNbinsX ());
   CalcSystematics (postSys, postSpectrum, postSpectrumSysHi, postSpectrumSysLo);
   postSys->SetFillColor (38);
   postSys->SetFillStyle (1001);

   preSpectrum->Draw ("hist");
   postSpectrum->DrawCopy ("same e1 x0");
   postSys->Draw ("3");

   const double textx = padLeftMargins[ix] + 0.5* (1.-padLeftMargins[ix]-padRightMargins[ix]);
   myText (textx, 0.85, kBlack, Form ("%g < #eta_{Lab} < %g", etabins[iEta], etabins[iEta+1]), 0.03*3);
   
   pads[ix][iy+1]->cd ();
   pads[ix][iy+1]->SetLogx (true);

   postSpectrum->Divide (preSpectrum);
   postSpectrumSysHi->Divide (preSpectrum);
   postSpectrumSysLo->Divide (preSpectrum);

   TGraphAsymmErrors* postSysRatio = new TGraphAsymmErrors (postSpectrum);
   CalcSystematics (postSysRatio, postSpectrum, postSpectrumSysHi, postSpectrumSysLo);

   postSpectrum->SetAxisRange (0.75, 1.25, "Y");
   postSpectrum->GetYaxis ()->SetTitle ("Insitu / EtaJES");
   postSpectrum->GetYaxis ()->SetNdivisions (403);
   postSpectrum->GetXaxis ()->SetTitleSize (0.03/ (padYwidth* (1-1/padRatio)));
   postSpectrum->GetYaxis ()->SetTitleSize (0.03/ (padYwidth* (1-1/padRatio)));
   postSpectrum->GetYaxis ()->CenterTitle (true);
   postSpectrum->GetXaxis ()->SetLabelSize (0.03/ (padYwidth* (1-1/padRatio)));
   postSpectrum->GetYaxis ()->SetLabelSize (0.03/ (padYwidth* (1-1/padRatio)));
   postSpectrum->GetXaxis ()->SetTickLength (0.08);
   postSpectrum->GetXaxis ()->SetTitleOffset (1.1);

   postSysRatio->SetFillColor (38);
   postSysRatio->SetFillStyle (1001);

   postSpectrum->DrawCopy ("e1 x0");
   postSysRatio->Draw ("3");

   postSys->Write ();
   postSysRatio->Write ();
  }
  jetInsituSpectrumCanvas->SaveAs (Form ("%s/JetInsituSpectrum.pdf", plotPath.Data ()));
  jetPreInsituSpectrum->Write ();
  jetPostInsituSpectrum->Write ();
  jetPostInsituSpectrumSysHi->Write ();
  jetPostInsituSpectrumSysLo->Write ();
  jetInsituSpectrumCanvas->Write ();


  TCanvas* jetInsituResponseCanvas = new TCanvas ("jetInsituResponseCanvas", "", 1600, 1200);
  jetInsituResponseCanvas->Divide (4, 3);
  jetInsituResponseCanvas->Draw ();
  //const double padYwidth = 0.48;
  //const double padRatio = 1.5; // ratio of size of upper pad to lower pad
  //const double lPadXbounds[3] = {0, 0.38, 0.69};
  //const double uPadXbounds[3] = {0.38, 0.69, 1.0};
  //const double lPadYbounds[4] = {1-padYwidth/padRatio, 1-padYwidth, padYwidth* (1-1/padRatio), 0};
  //const double uPadYbounds[4] = {1, 1-padYwidth/padRatio, padYwidth, padYwidth* (1-1/padRatio)};

  for (short iP = 0; iP < numpbins; iP++) {
   jetInsituResponse[iP]->Write ();
   jetInsituResponseSysHi[iP]->Write ();
   jetInsituResponseSysLo[iP]->Write ();
  }

  for (short iEta = 0; iEta < numetabins; iEta++) {

   for (short iP = 0; iP < numpbins; iP++) {
    if (iP <= 11) jetInsituResponseCanvas->cd (iP+1);
    TPad* thisPad = (TPad*)jetInsituResponseCanvas->GetPad (iP);
    thisPad->SetTopMargin (0.1);
    //gPad->SetLogy (true);

    TH1D* thisHist = jetInsituResponse[iP]->ProjectionX (Form ("response_iP%i_iEta%i", iP, iEta), iEta+1, iEta+1);
    TH1D* sysHi = jetInsituResponseSysHi[iP]->ProjectionX (Form ("response_syshi_iP%i_iEta%i", iP, iEta), iEta+1, iEta+1);
    TH1D* sysLo = jetInsituResponseSysLo[iP]->ProjectionX (Form ("response_syslo_iP%i_iEta%i", iP, iEta), iEta+1, iEta+1);

    if (thisHist->Integral () != 0) {
     thisHist->Scale (1./thisHist->Integral ());
     sysHi->Scale (1./sysHi->Integral ());
     sysLo->Scale (1./sysLo->Integral ());
    }

    //thisHist->Rebin (2);
    //sysHi->Rebin (2);
    //sysLo->Rebin (2);

    thisHist->SetLineColor (kBlack);
    thisHist->SetMarkerColor (kBlack);

    TGraphAsymmErrors* sysGraph = new TGraphAsymmErrors (thisHist);
    CalcSystematics (sysGraph, thisHist, sysHi, sysLo);

    sysHi->SetLineColor (kRed);
    sysHi->SetMarkerColor (kRed);
    sysLo->SetLineColor (kBlue);
    sysLo->SetMarkerColor (kBlue);

    thisHist->GetXaxis ()->SetTitleSize (0.08);
    thisHist->GetXaxis ()->SetLabelSize (0.08);
    thisHist->GetXaxis ()->SetTitleOffset (0.9);
    thisHist->GetXaxis ()->CenterTitle (true);
    thisHist->GetYaxis ()->SetTitleSize (0.08);
    thisHist->GetYaxis ()->SetLabelSize (0.08);
    thisHist->GetXaxis ()->SetNdivisions (503);

    if (iP <= 11) {
     thisHist->Draw ("e1");
     //sysHi->Draw ("same p");
     //sysLo->Draw ("same p");

     //sysGraph->SetFillStyle (1001);
     //sysGraph->SetFillColor (14);
     //sysGraph->Draw ("2");

     myText (0.18, 0.95, kBlack, Form ("%gGeV < #it{p}_{T} < %gGeV", pbins[iP], pbins[iP+1]), 0.08);
    }

   }
   jetInsituResponseCanvas->cd (1);
   myText (0.2, 0.8, kBlack, Form ("%g < #eta_{Lab} < %g", etabins[iEta], etabins[iEta+1]), 0.03*3);
   jetInsituResponseCanvas->SaveAs (Form ("%s/JetInsituResponse_%i.pdf", plotPath.Data (), iEta));
  }  

  outFile->Write ();
  if (outFile) delete outFile;

  return;
}

} // end namespace
