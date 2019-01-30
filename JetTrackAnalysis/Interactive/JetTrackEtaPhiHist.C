#include "Params.h"
#include "TemplateFitting.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <TLine.h>
#include <TH2D.h>
#include <TFile.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <iostream>

using namespace atlashi;
using namespace std;
using namespace TemplateFitting;

using namespace JetTrackAnalysis;

void JetTrackEtaPhiHist () {

  SetAtlasStyle ();

  SetupDirectories ("JetTrackEtaPhi/", "JetTrackAnalysis/");

  TFile* inDataSignalFile = new TFile (Form ("%s/SignalEvents/data_combined.root", rootPath.Data ()), "read");
  TFile* inDataMixedFile = new TFile (Form ("%s/MixedEvents/data_combined.root", rootPath.Data ()), "read");
  TFile* inMCSignalFile = new TFile (Form ("%s/SignalEvents/mc_combined.root", rootPath.Data ()), "read");
  TFile* inMCMixedFile = new TFile (Form ("%s/MixedEvents/mc_combined.root", rootPath.Data ()), "read");

  TH1D****** JetPtDists = Get5DArray <TH1D*> (2, 2, numCentBins, numPtBins, numVertZBins+1);
  TH2D******* Correlations = Get6DArray <TH2D*> (2, 2, 2, numCentBins, numPtBins, numVertZBins+1); // iData (2), iPer (2), iSignal (2), iCent (numCentBins), iPt (numPtBins), iVertZ (numVertZBins+1)
  //TH2D* JetPtSpectrum[numCentBins] = {};
  //TH2D* JetCounts[numCentBins] = {};
  //TH2D* TrackPtSpectrum[numCentBins] = {};
  //TH2D* TrackCounts[numCentBins] = {};

  //TH1D* TrackSelectionHist[numCentBins] = {};
  //TH1D* JetSelectionHist[numCentBins] = {};

  //TH2D* DijetDeltaEtaDeltaPhi[numCentBins] = {};
  //TH2D* DijetEta1Eta2[numCentBins] = {};
  //TH1D* DijetAj[numCentBins] = {};

  //double ajmax = 0;
  //double ajmin = 1e100;
  for (short iData = 0; iData < 2; iData++) {
    for (short iSignal = 0; iSignal < 2; iSignal++) { // 0 = signal, 1 = mixed
      TFile* f = NULL;
      if      (iData == 0 && iSignal == 0)
        f = inDataSignalFile;
      else if (iData == 0 && iSignal == 1)
        f = inDataMixedFile;
      else if (iData == 1 && iSignal == 0)
        f = inMCSignalFile;
      else
        f = inMCMixedFile;

      for (short iPer = 0; iPer < 2; iPer++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iPt = 0; iPt < numPtBins; iPt++) {
            for (short iVertZ = 0; iVertZ < numVertZBins; iVertZ++) {

              const char* histName = Form ("Correlation_iVertZ%i_iCent%i_iPt%i_%s_p%s_%s", iVertZ, iCent, iPt, iData == 0 ? "data":"mc", iPer == 0 ? "A":"B", iSignal == 0 ? "SE":"ME");
              if (!f->GetListOfKeys ()->Contains (histName)) {
                //cout << "Cannot find " << histName << endl;
                continue;
              }

              Correlations[iData][iPer][iSignal][iCent][iPt][iVertZ] = (TH2D*)f->Get (histName);
            }
          }
          //JetPtSpectrum[iCent] = (TH2D*)inFile->Get (Form ("JetPtSpectrum_cent%i", iCent));
          //JetPtSpectrum[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

          //JetCounts[iCent] = (TH2D*)inFile->Get (Form ("JetCounts_cent%i", iCent)); // divide by centrality bin width
          //JetCounts[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

          //TrackPtSpectrum[iCent] = (TH2D*)inFile->Get (Form ("TrackPtSpectrum_cent%i", iCent));
          //TrackPtSpectrum[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

          //TrackCounts[iCent] = (TH2D*)inFile->Get (Form ("TrackCounts_cent%i", iCent));
          //TrackCounts[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

          //TrackSelectionHist[iCent] = (TH1D*)inFile->Get (Form ("TrackSelectionHist_cent%i", iCent));
          //JetSelectionHist[iCent] = (TH1D*)inFile->Get (Form ("JetSelectionHist_cent%i", iCent));

          //DijetDeltaEtaDeltaPhi[iCent] = (TH2D*)inFile->Get (Form ("DijetDeltaEtaDeltaPhi_cent%i", iCent));
          //DijetDeltaEtaDeltaPhi[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

          //DijetEta1Eta2[iCent] = (TH2D*)inFile->Get (Form ("DijetEta1Eta2_cent%i", iCent));
          //DijetEta1Eta2[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

          //DijetAj[iCent] = (TH1D*)inFile->Get (Form ("DijetAj_cent%i", iCent));
          //DijetAj[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));
          //if (DijetAj[iCent]->GetMaximum () > ajmax)
          //  ajmax = DijetAj[iCent]->GetMaximum ();
          //if (DijetAj[iCent]->GetMinimum () > ajmin)
          //  ajmin = DijetAj[iCent]->GetMinimum ();
        }
      }
    }
  }
  for (short iData = 0; iData < 2; iData++) {
    TFile* f = iData == 0 ? inDataSignalFile : inMCSignalFile;
    for (short iPer = 0; iPer < 2; iPer++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (short iPt = 0; iPt < numPtBins; iPt++) {
          for (short iVertZ = 0; iVertZ < numVertZBins; iVertZ++) {
            const char* histName = Form ("JetPtDist_iVertZ%i_iCent%i_iPt%i_%s_p%s", iVertZ, iCent, iPt, iData == 0 ? "data":"mc", iPer == 0 ? "A":"B");
            JetPtDists[iData][iPer][iCent][iPt][iVertZ] = (TH1D*)f->Get (histName);
          }
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Add up histograms with vertices between -5 and 5 along z
  //////////////////////////////////////////////////////////////////////////////
  for (short iData = 0; iData < 2; iData++) {
    for (short iSignal = 0; iSignal < 2; iSignal++) { // 0 = signal, 1 = mixed
      for (short iPer = 0; iPer < 2; iPer++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iPt = 0; iPt < numPtBins; iPt++) {
            if (!Correlations[iData][iPer][iSignal][iCent][iPt][0])
              continue;

            const char* histName = Form ("Correlation_iCent%i_iPt%i_%s_p%s_%s", iCent, iPt, iData == 0 ? "data":"mc", iPer == 0 ? "A":"B", iSignal == 0 ? "SE":"ME");
            Correlations[iData][iPer][iSignal][iCent][iPt][numVertZBins] = new TH2D (histName, "", 160, -8, 8, numphibins, phibins);
            for (short iVertZ = 10; iVertZ < 20; iVertZ++) {
              Correlations[iData][iPer][iSignal][iCent][iPt][numVertZBins]->Add (Correlations[iData][iPer][iSignal][iCent][iPt][iVertZ]);
            }
          }
        }
      }
    }
  }
  for (short iData = 0; iData < 2; iData++) {
    for (short iPer = 0; iPer < 2; iPer++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (short iPt = 0; iPt < numPtBins; iPt++) {
          if (!JetPtDists[iData][iPer][iCent][iPt][0])
            continue;

          const char* histName = Form ("JetPtDist_iCent%i_iPt%i_%s_p%s", iCent, iPt, iData == 0 ? "data":"mc", iPer == 0 ? "A":"B");
          JetPtDists[iData][iPer][iCent][iPt][numVertZBins] = new TH1D (histName, "", numPtBins, ptBins);
          for (short iVertZ = 10; iVertZ < 20; iVertZ++) {
            JetPtDists[iData][iPer][iCent][iPt][numVertZBins]->Add (JetPtDists[iData][iPer][iCent][iPt][iVertZ]);
          }
        }
      }
    }
  }

//  for (short iData = 0; iData < 2; iData++) {
//    for (short iCent = 0; iCent < numCentBins; iCent++) {
//      for (short iPt = 0; iPt < numPtBins; iPt++) {
//        Correlations[iData][2][iCent][iPt] = new TH2D (Form ("Correlations_%s_cent%i_pt%i_C", iData == 0 ? "data":"MC", iCent, iPt), "", 80, 0, 8, numphibins, phibins);
//          Correlations[iData][iSignal][iCent][iPt] = (TH2D*)f->Get (Form ("Correlations_cent%i_pt%i_%s", iCent, iPt, iSignal == 0 ? "SE":"ME"));

  const Color_t colors[9] = {kBlack, kRed, kBlue, kMagenta, 8, kCyan+1, kOrange, kViolet, kGray};

  //TCanvas* c = new TCanvas ("c", "", 900, 800);
  //gStyle->SetPalette (kRainBow);
  ////gPad->SetLogz();
  //Correlations[0]->RebinX (2);
  //Correlations[0]->RebinY (2);
  //Correlations[0]->Draw ("lego2");
  //c->SaveAs (Form ("%s/JetTrackEtaPhiCorr.pdf", plotPath.Data ()));

  //TCanvas* DijetDeltaEtaDeltaPhiCanvas = new TCanvas ("DijetDeltaEtaDeltaPhiCanvas", "", 800, 600);
  //FormatTH2Canvas (DijetDeltaEtaDeltaPhiCanvas);

  //TCanvas* DijetDeltaPhiCanvas = new TCanvas ("DijetDeltaPhiCanvas", "", 800, 600);

  //TCanvas* DijetEta1Eta2Canvas = new TCanvas ("DijetEta1Eta2Canvas", "", 800, 600);
  //FormatTH2Canvas (DijetEta1Eta2Canvas);

  //TCanvas* DijetAjCanvas = new TCanvas ("DijetAjCanvas", "", 800, 600);
  //for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //  DijetDeltaEtaDeltaPhiCanvas->cd ();
  //  gPad->SetLogz ();

  //  DijetDeltaEtaDeltaPhi[iCent]->GetYaxis ()->SetTitleOffset (1);
  //  DijetDeltaEtaDeltaPhi[iCent]->GetZaxis ()->SetTitle ("Relative counts");
  //  DijetDeltaEtaDeltaPhi[iCent]->Draw ("colz");

  //  myText (0.7, 0.88, kBlack, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 0.045);
  //  DijetDeltaEtaDeltaPhiCanvas->SaveAs (Form ("%s/DijetDeltaEtaDeltaPhi_cent%i.pdf", plotPath.Data (), iCent));

  //  DijetDeltaPhiCanvas->cd ();
  //  gPad->SetLogy ();
  //  TH1D* th = DijetDeltaEtaDeltaPhi[iCent]->ProjectionY ();
  //  th->GetYaxis ()->SetRangeUser (1e6, 1e9);
  //  th->SetLineColor (colors[iCent]);
  //  th->SetMarkerColor (colors[iCent]);
  //  th->GetYaxis ()->SetTitle ("Relative counts");
  //  th->GetYaxis ()->SetTitleOffset (1);
  //  th->GetXaxis ()->SetTitleOffset (1);
  //
  //  if (iCent == numCentBins-1) {
  //    th->DrawCopy ("e1");
  //    TLine* line1 = new TLine (3*pi/4, 1e6, 3*pi/4, 1e9);
  //    line1->SetLineStyle (7);
  //    line1->Draw ();
  //    TLine* line2 = new TLine (5*pi/4, 1e6, 5*pi/4, 1e9);
  //    line2->SetLineStyle (7);
  //    line2->Draw ();
  //  }
  //  else
  //    th->DrawCopy ("same e1");
  //  if (th) { delete th; th = NULL; }
  //  myMarkerText (0.25, 0.9 - iCent*0.05, colors[iCent], kFullCircle, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.045/gPad->GetHNDC ());

  //  DijetEta1Eta2Canvas->cd ();
  //  gPad->SetLogz ();
  //  DijetEta1Eta2[iCent]->GetXaxis ()->SetTitle ("#eta^{lead}_{det}");
  //  DijetEta1Eta2[iCent]->GetXaxis ()->SetTitleOffset (1);
  //  DijetEta1Eta2[iCent]->GetYaxis ()->SetTitle ("#eta^{sublead}_{det}");
  //  DijetEta1Eta2[iCent]->GetYaxis ()->SetTitleOffset (1);
  //  DijetEta1Eta2[iCent]->Draw ("colz");
  //  
  //  myText (0.7, 0.88, kBlack, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 0.045);
  //  DijetEta1Eta2Canvas->SaveAs (Form ("%s/DijetEta1Eta2corr_cent%i.pdf", plotPath.Data (), iCent));

  //  DijetAjCanvas->cd ();
  //  gPad->SetLogy ();
  //  DijetAj[iCent]->GetYaxis ()->SetRangeUser (0.5*ajmin, 2*ajmax);
  //  DijetAj[iCent]->SetMarkerColor (colors[iCent]);
  //  DijetAj[iCent]->SetLineColor (colors[iCent]);
  //  if (iCent == numCentBins-1)
  //    DijetAj[iCent]->Draw ("e1");
  //  else
  //    DijetAj[iCent]->Draw ("same e1");
  //}
  //DijetDeltaPhiCanvas->SaveAs (Form ("%s/DijetDeltaPhi.pdf", plotPath.Data ()));
  //DijetAjCanvas->SaveAs (Form ("%s/DijetAjCorr.pdf", plotPath.Data ()));


  TCanvas* JetTrackPhiCanvas = new TCanvas ("JetTrackPhiCanvas", "", 2000, 600);
  TPad* leftPad = new TPad ("leftPad", "", 0, 0, 0.333, 1);
  TPad* middlePad = new TPad ("middlePad", "", 0.333, 0, 1-0.333, 1);
  TPad* rightPad = new TPad ("rightPad", "", 1-0.333, 0, 1, 1);
  leftPad->SetLeftMargin (0.20);
  leftPad->SetRightMargin (0.02);
  middlePad->SetLeftMargin (0.20);
  middlePad->SetRightMargin (0.02);
  rightPad->SetRightMargin (0.02);
  rightPad->SetLeftMargin (0.20);
  JetTrackPhiCanvas->Draw ();
  leftPad->Draw ();
  middlePad->Draw ();
  rightPad->Draw ();

  TFile* outFile = new TFile (Form ("%s/perTriggerYields.root", rootPath.Data ()), "recreate");
  TH1D****** perTriggerYields = Get5DArray <TH1D*> (2, 2, 3, numCentBins, numPtBins);

  for (short iData = 0; iData < 2; iData++) {
    for (short iPer = 0; iPer < 2; iPer++) {
      for (short iPt = 0; iPt < numPtBins; iPt++) {
        bool skip = false;
        for (short iEta = 0; iEta < 3; iEta++) {
          const int eta_lo = (iEta == 0 ? -6 : (iEta == 1 ? -2 : 2));
          const int eta_hi = (iEta == 0 ? -2 : (iEta == 1 ? 2  : 6));

          if (iEta == 0)
            leftPad->cd ();
          else if (iEta == 1)
            middlePad->cd ();
          else
            rightPad->cd ();
          gPad->SetLogy (false);

          for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
            TH2D* signal2d = Correlations[iData][iPer][0][iCent][iPt][numVertZBins];
            TH2D* mixed2d  = Correlations[iData][iPer][1][iCent][iPt][numVertZBins];

            if (!signal2d || !mixed2d) {
              skip = true;
              break;
            }

            float ntrig = JetPtDists[iData][iPer][iCent][iPt][numVertZBins]->Integral ();

            TH1D* signal = NULL, *mixed = NULL;

            //TH1D* signal_left  = signal2d->ProjectionY ("signal_left",  signal2d->GetXaxis ()->FindBin (-eta_hi), signal2d->GetXaxis ()->FindBin (-eta_lo));
            //TH1D* mixed_left   = mixed2d ->ProjectionY ("mixed_left",   mixed2d ->GetXaxis ()->FindBin (-eta_hi), mixed2d ->GetXaxis ()->FindBin (-eta_lo));
            //TH1D* signal_right = signal2d->ProjectionY ("signal_right", signal2d->GetXaxis ()->FindBin (eta_lo),  signal2d->GetXaxis ()->FindBin (eta_hi));
            //TH1D* mixed_right  = mixed2d ->ProjectionY ("mixed_right",  mixed2d ->GetXaxis ()->FindBin (eta_lo),  mixed2d ->GetXaxis ()->FindBin (eta_hi));

            //signal = new TH1D ("signal", "", numphibins, phibins);
            //mixed  = new TH1D ("mixed",  "", numphibins, phibins); 

            //signal->Add (signal_left, signal_right);
            //mixed ->Add (mixed_left,  mixed_right);

            //if (signal_left)  { delete signal_left;  signal_left = NULL; }
            //if (mixed_left)   { delete mixed_left;   mixed_left = NULL; }
            //if (signal_right) { delete signal_right; signal_right = NULL; }
            //if (mixed_right)  { delete mixed_right;  mixed_right = NULL; }

            if (iPer == 0) {
              signal = signal2d->ProjectionY ("signal",  signal2d->GetXaxis ()->FindBin (-eta_hi), signal2d->GetXaxis ()->FindBin (-eta_lo));
              mixed  = mixed2d ->ProjectionY ("mixed",   mixed2d ->GetXaxis ()->FindBin (-eta_hi), mixed2d ->GetXaxis ()->FindBin (-eta_lo));
            }
            else {
              signal = signal2d->ProjectionY ("signal", signal2d->GetXaxis ()->FindBin (eta_lo),  signal2d->GetXaxis ()->FindBin (eta_hi));
              mixed  = mixed2d ->ProjectionY ("mixed",  mixed2d ->GetXaxis ()->FindBin (eta_lo),  mixed2d ->GetXaxis ()->FindBin (eta_hi));
            }

            signal->Rebin (8);
            mixed ->Rebin (8);

            TH1D* thisHist = (TH1D*)signal->Clone (Form ("perTriggerYield_iCent%i_iPt%i_%s_p%s_%s", iCent, iPt, iData == 0 ? "data":"mc", iPer == 0 ? "A":"B", iEta == 0 ? "Pbdir": (iEta == 1 ? "src":"pdir")));
            thisHist->Divide (signal, mixed); // construct C = S/B

            // construct Y = 1/2pi integral(B)/Ntrig * C
            if (ntrig > 0) thisHist->Scale (1./ntrig);
            thisHist->Scale (mixed->Integral () / (2*pi), "width");

            cout << "Integral of Y: " << thisHist->Integral ("width") << endl;

            if (ZYAM) {
              const float min = thisHist->GetMinimum ();
              for (int ix = 1; ix <= thisHist->GetNbinsX (); ix++) {
                thisHist->SetBinContent (ix, thisHist->GetBinContent (ix) - min); 
              }
            }

            if (signal) { delete signal; signal = NULL; }
            if (mixed)  { delete mixed;  mixed = NULL; }

            perTriggerYields[iData][iPer][iEta][iCent][iPt] = thisHist;
          }

          if (skip)
            continue;

          double max = 0, min=2e32;
          for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
            TH1D* h = perTriggerYields[iData][iPer][iEta][iCent][iPt];
            if (h->GetMaximum () > max)
              max = h->GetMaximum ();
            if (h->GetMinimum () < min)
              min = h->GetMinimum ();
          }

          for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
            TH1D* thisHist = perTriggerYields[iData][iPer][iEta][iCent][iPt];

            if (ZYAM) {
              thisHist->GetYaxis ()->SetRangeUser (0, 1.6* max);
              thisHist->GetYaxis ()->SetTitle ("dY (#Delta#phi) / d#Delta#phi, ZYAM");
            }
            else {
              thisHist->GetYaxis ()->SetRangeUser (0.9* min, 1.6* max);
              thisHist->GetYaxis ()->SetTitle ("dY (#Delta#phi) / d#Delta#phi");
            }
            thisHist->GetXaxis ()->SetTitle ("#phi_{det}^{Jet} - #phi_{det}^{Trk}");
            thisHist->GetXaxis ()->SetTitleOffset (1);
            thisHist->GetYaxis ()->SetTitleOffset (1.6);
            thisHist->SetLineColor (colors[iCent]);
            thisHist->SetMarkerColor (colors[iCent]);
            if (iCent == numCentBins-1)
              thisHist->DrawCopy ("e1 x0");
            else
              thisHist->DrawCopy ("same e1 x0");

            if (iEta == 0) {
              myMarkerText (0.76, 0.70 - iCent*0.05, colors[iCent], kFullCircle, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.045/gPad->GetHNDC ());
            }
            thisHist->Write ();
          }

          if (iEta == 0) {
            myText (0.65, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/gPad->GetHNDC ());
            myText (0.65, 0.82, kBlack, Form ("%g < #Delta#eta_{proton} < %g", (float)eta_lo, (float)eta_hi), 0.045/gPad->GetHNDC ()); 
            myText (0.65, 0.77, kBlack, iData == 0 ? "Data":"MC", 0.045/gPad->GetHNDC ());
          }
          else if (iEta == 1) {

          }
          else {
            myText (0.25, 0.89, kBlack, "#left|#Delta#eta_{Trk}^{Sublead.}#right| > 2", 0.045/gPad->GetHNDC ());
            myText (0.25, 0.82, kBlack, Form ("%g < #Delta#eta_{proton} < %g", (float)eta_lo, (float)eta_hi), 0.045/gPad->GetHNDC ()); 
            myText (0.25, 0.76, kBlack, Form ("%g GeV < #it{p}_{T}^{lead} < %g GeV", ptBins[iPt], ptBins[iPt+1]), 0.045/gPad->GetHNDC ());
            myText (0.25, 0.70, kBlack, "10 GeV < #it{p}_{T}^{sublead}", 0.045/gPad->GetHNDC ());
          }
        } // end loop over eta bins
        if (!skip)
          JetTrackPhiCanvas->SaveAs (Form ("%s/%s/JetTrackPhiCorr_%s_p%s_iPt%i.pdf", plotPath.Data (), iData == 0 ? "Data":"MC", iData == 0 ? "data":"mc", iPer == 0 ? "A":"B", iPt));
        JetTrackPhiCanvas->Draw ();

        cout << "Press enter to continue..." << flush;
        //cin.get ();
        gPad->WaitPrimitive();
        cout << endl;
      } // end loop over pT cuts
    } // end loop over periods
  } // end loop over data sets

  cout << "Calculating example template fit..." << endl;
  Fitting* result;
  double v22, v22err;
  result = TemplateFit(perTriggerYields[0][0][1][0][0], perTriggerYields[0][0][1][numCentBins-1][0]);
	result->GetVnnAndError(v22, v22err, 0);
	result->c1->SaveAs (Form ("%s/test.pdf", plotPath.Data ()));

  if (result) { delete result; result = NULL; }


  Delete5DArray (perTriggerYields, 2, 2, 3, numCentBins, numPtBins);

  outFile->Close ();
}
