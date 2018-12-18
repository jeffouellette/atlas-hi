#include "ZMassCalcHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <GlobalParams.h>
#include <ArrayTemplates.h>

#include <TTree.h>
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
#include <TLorentzVector.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace JetCalibration {

void ZMassCalcHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("ZMassCalc/", "JetCalibration/");

  TH3D*** zMassSpectra = Get2DArray <TH3D*> (2, 2);
  TH3D*** zMassCounts = Get2DArray <TH3D*> (2, 2);

  const int numzmassbins = 50;
  const double* zmassbins = linspace (60, 110, numzmassbins);

  for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
    const char* data = (iData == 0 ? "data":"mc");
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      const char* spc = (iSpc == 0 ? "mumu":"ee");

      zMassSpectra[iSpc][iData] = new TH3D (Form ("z%sMassSpectrum_%s", spc, data), "", numzmassbins, zmassbins, numzetabins, zetabins, numphibins, phibins);
      zMassSpectra[iSpc][iData]->Sumw2 ();

      zMassCounts[iSpc][iData] = new TH3D (Form ("z%sMassCounts_%s", spc, data), "", numzmassbins, zmassbins, numzetabins, zetabins, numphibins, phibins);
      zMassCounts[iSpc][iData]->Sumw2 ();
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Load analyzed TTrees
  //////////////////////////////////////////////////////////////////////////////
  float zpt = 0, zeta = 0, zphi = 0, zm = 0, l1eta = 0, l2eta = 0;
  double evtWeight = 0;
  bool isMC = false, isPeriodA = false, isZee = false;

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* inTree = (TTree*)inFile->Get ("jeffsztree");

  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("isZee", &isZee);
  inTree->SetBranchAddress ("Z_pt", &zpt);
  inTree->SetBranchAddress ("Z_eta", &zeta);
  inTree->SetBranchAddress ("Z_phi", &zphi);
  inTree->SetBranchAddress ("Z_m", &zm);
  inTree->SetBranchAddress ("l1_eta", &l1eta);
  inTree->SetBranchAddress ("l2_eta", &l2eta);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);

  //////////////////////////////////////////////////////////////////////////////
  // Fill desired histograms
  //////////////////////////////////////////////////////////////////////////////
  const int nZJets = inTree->GetEntries ();
  TLorentzVector Z;
  for (int zJet = 0; zJet < nZJets; zJet++) {
    inTree->GetEntry (zJet);

    if (1.37 < fabs (l1eta) || 1.37 < fabs (l2eta))
      continue;

    const short iSpc = isZee ? 1 : 0; // 0 = mumu, 1 = ee
    const short iData = isMC ? 1 : 0;
    //const short iPer = isPeriodA ? 0 : 1;

    Z.SetPtEtaPhiM (zpt, zeta, zphi, zm);

    zMassSpectra[iSpc][iData]->Fill (zm, Z.Rapidity (), zphi, evtWeight);
    zMassCounts[iSpc][iData]->Fill (zm, Z.Rapidity (), zphi);
  }
  inTree = NULL;
  inFile->Close ();
  if (inFile) { delete inFile; inFile = NULL; }


  //////////////////////////////////////////////////////////////////////////////
  // Save histograms for interactive access
  //////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iData = 0; iData < 2; iData++) {
      zMassSpectra[iSpc][iData]->Write ();
      zMassCounts[iSpc][iData]->Write ();
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Canvas definitions
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 800);
  const double padRatio = 1.2; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetTopMargin (0.04);
  topPad->SetBottomMargin (0);
  topPad->SetLeftMargin (0.15);
  topPad->SetRightMargin (0.04);
  bottomPad->SetTopMargin (0);
  bottomPad->SetBottomMargin (0.20);
  bottomPad->SetLeftMargin (0.15);
  bottomPad->SetRightMargin (0.04);
  topPad->Draw ();
  bottomPad->Draw ();

  //////////////////////////////////////////////////////////////////////////////
  // Plotting elements 
  //////////////////////////////////////////////////////////////////////////////
  TLine* lines[5] = {};
  for (short i = 0; i < 5; i++) {
    const float dzm = 0.2;
    lines[i] = new TLine (60, 0.6+dzm*i, 110, 0.6+dzm*i);
    if (0.6+dzm*i == 1) lines[i]->SetLineStyle (1);
    else lines[i]->SetLineStyle (3);
  }


  //////////////////////////////////////////////////////////////////////////////
  // Plot Z mass spectra
  //////////////////////////////////////////////////////////////////////////////
  for (short iEta = 0; iEta <= numzetabins; iEta++) {
    const short eta_lo = (iEta != numzetabins ? iEta+1 : 1);
    const short eta_hi = (iEta != numzetabins ? iEta+1 : numzetabins);

    for (short iSpc = 0; iSpc < 2; iSpc++) {
      const char* spc = (iSpc == 0 ? "mumu" : "ee");
      const Color_t dataColor = kBlack;

      topPad->cd ();
      topPad->SetLogx (0);

      double mean[2] = {};
      double mean_err[2] = {};
      double sigma[2] = {};
      double sigma_err[2] = {};

      TH1D** hists = Get1DArray <TH1D*> (2);
      TF1** fits = Get1DArray <TF1*> (2);
      //RooZfit* fits[2];
      for (short iData = 1; iData >= 0; iData--) {
        const char* data = (iData == 0 ? "data" : "mc");
        const Color_t color = (iData == 0 ? dataColor : mcOverlayColor);

        TH2D* proj2d = Project2D ("proj2d", zMassSpectra[iSpc][iData], "x", "z", eta_lo, eta_hi, exclusive && iEta == numzetabins);
        TH1D* thisHist = proj2d->ProjectionX (Form ("zMassSpectra_%s_%s_iEta%i", spc, data, iEta));
        if (proj2d) { delete proj2d; proj2d = NULL; }

        //TH1D* thisHist = zMassSpectra[iSpc][iData]->ProjectionX (Form ("zMassSpectra_%s_%s_iEta%i", spc, data, iEta), eta_lo, eta_hi);

        hists[iData] = thisHist;
        thisHist->Scale (1.0, "width");

        thisHist->SetLineColor (color);
        thisHist->SetMarkerColor (color);

        thisHist->GetXaxis ()->SetTitle (Form ("#font[12]{%s} Invariant Mass #left[GeV#right]", (iSpc == 1 ? "ee":"#mu#mu")));
        thisHist->GetYaxis ()->SetTitle ("Counts / Area under Curve");
        thisHist->GetYaxis ()->SetTitleSize (0.035/uPadY);
        thisHist->GetYaxis ()->SetTitleOffset (1.5*uPadY);
        thisHist->GetYaxis ()->SetLabelSize (0.035/uPadY);
//        thisHist->GetXaxis ()->SetNdivisions (50802, false);

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

        thisHist->GetYaxis ()->SetRangeUser (0, 0.20);
        thisHist->GetYaxis ()->ChangeLabel (1, -1, -1, -1, -1, -1, " ");
      }

      for (short iData = 1; iData >= 0; iData--) {
        TH1D* thisHist = hists[iData];
        if (iData == 1) thisHist->Draw ("hist");
        else thisHist->Draw ("p same");
        //fits[iData]->plot->Draw ("same");
        fits[iData]->Draw ("same");
        if (eta_lo != 1 || eta_hi != numzetabins) {
          if (!exclusive || iEta != numzetabins)
            myText (0.205, 0.15, kBlack, Form ("%g < y_{det}^{Z} < %g", zetabins[eta_lo-1], zetabins[eta_hi]), 0.035/uPadY);
          else if (fabs (zetabins[eta_lo-1]) != fabs (zetabins[eta_hi]))
            myText (0.205, 0.15, kBlack, Form ("y_{det}^{Z} < %g #intersection %g < y_{det}^{Z}", zetabins[eta_lo-1], zetabins[eta_hi]), 0.035/uPadY);
          else
            myText (0.205, 0.15, kBlack, Form ("#left|y_{det}^{Z}#right| > %g", fabs (zetabins[eta_hi])), 0.035/uPadY);
        }
      }

      int countsData = 0, countsMC = 0;
      if (exclusive && iEta == numzetabins) {
        countsData = zMassCounts[iSpc][0]->Integral () - zMassCounts[iSpc][0]->Integral (1, numzmassbins, eta_lo, eta_hi, 1, numphibins);
        countsMC = zMassCounts[iSpc][1]->Integral () - zMassCounts[iSpc][1]->Integral (1, numzmassbins, eta_lo, eta_hi, 1, numphibins);
      }
      else {
        countsData = zMassCounts[iSpc][0]->Integral (1, numzmassbins, eta_lo, eta_hi, 1, numphibins);
        countsMC = zMassCounts[iSpc][1]->Integral (1, numzmassbins, eta_lo, eta_hi, 1, numphibins);
      }
      if (iSpc == 0)
        myText (0.205, 0.88, kBlack, "Z#rightarrow#mu#mu", 0.035/uPadY);
      else if (iSpc == 1)
        myText (0.205, 0.88, kBlack, "Z#rightarrowee", 0.035/uPadY);
      myMarkerText (0.225, 0.80, dataColor, kFullCircle, Form ("2016 Data (%i events)", countsData), 1.25, 0.035/uPadY);
      myMarkerText (0.225, 0.55, mcOverlayColor, kFullCircle, Form ("Pythia8 MC + Overlay (%i events)", countsMC), 1.25, 0.035/uPadY);
      myText (0.205, 0.72, kBlack, Form ("m_{Z}^{data} = %.2f #pm %.2f GeV", mean[0], mean_err[0]), 0.035/uPadY);
      myText (0.205, 0.64, kBlack, Form ("#sigma_{Z}^{data} = %.2f #pm %.2f GeV", sigma[0], sigma_err[0]), 0.035/uPadY);
      myText (0.205, 0.47, kBlack, Form ("m_{Z}^{mc} = %.2f #pm %.2f GeV", mean[1], mean_err[1]), 0.035/uPadY);
      myText (0.205, 0.39, kBlack, Form ("#sigma_{Z}^{mc} = %.2f #pm %.2f GeV", sigma[1], sigma_err[1]), 0.035/uPadY);

      bottomPad->cd ();
      bottomPad->SetLogx (0);

      TH1D* thisHist = (TH1D*)(hists[0]->Clone (Form ("dataOverMC_%s_iEta%i", spc, iEta)));
      thisHist->Divide (hists[1]);

      thisHist->SetLineColor (dataColor);
      thisHist->SetMarkerColor (dataColor);

      //thisHist->GetYaxis ()->SetRangeUser (0.6, 1.6);
      thisHist->GetYaxis ()->SetRangeUser (0.2, 2.2);
      thisHist->GetXaxis ()->SetTitle (Form ("#font[12]{%s} Invariant Mass #left[GeV#right]", (iSpc == 1 ? "ee":"#mu#mu")));
      thisHist->GetYaxis ()->SetTitle ("Data / MC");
      thisHist->GetXaxis ()->SetTitleSize (0.035/dPadY);
      thisHist->GetYaxis ()->SetTitleSize (0.035/dPadY);
      thisHist->GetXaxis ()->SetTitleOffset (1);
      thisHist->GetYaxis ()->SetTitleOffset (1.5*dPadY);
      thisHist->GetYaxis ()->CenterTitle (true);
      thisHist->GetXaxis ()->SetLabelSize (0.035/dPadY);
      thisHist->GetYaxis ()->SetLabelSize (0.035/dPadY);

      thisHist->GetYaxis ()->ChangeLabel (-1, -1, -1, -1, -1, -1, " ");
      thisHist->GetXaxis ()->SetTickLength (0.08);
//      thisHist->GetXaxis ()->SetNdivisions (50802, false);
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

      if (thisHist) { delete thisHist; thisHist = NULL; }
      Delete1DArray (hists, 2);
      Delete1DArray (fits, 2);
    }
  }

  outFile->Close ();
  if (outFile) delete outFile;

  return;
}

} // end namespace
