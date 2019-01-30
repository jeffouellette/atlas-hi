#include "PhotonAnalysisHist.h"
#include "Params.h"

#include <Utils.h>
#include <GlobalParams.h>
#include <ArrayTemplates.h>
#include <Trigger.h>

#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TText.h>

#include <iostream>
#include <string>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace std;
using namespace atlashi;

namespace offlineAnalyses {

float lumi_int = 0;

const char* GetShowerShape (const int iShowerShape) {
  switch (iShowerShape) {
    case 0:  return "Rhad1";
    case 1:  return "Rhad";
    case 2:  return "e277";
    case 3:  return "Reta";
    case 4:  return "Rphi";
    case 5:  return "weta1";
    case 6:  return "weta2";
    case 7:  return "wtots1";
    case 8:  return "f1";
    case 9:  return "f3";
    case 10: return "fracs1"; // fside
    case 11: return "DeltaE";
    case 12: return "Eratio";
    default: return "";
  }
}


void PhotonAnalysisHist () {

  SetAtlasStyle ();
  SetupDirectories ("PhotonAnalysis/", "OfflineAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TH1D**** ptSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 2);

  TH1D*** xJgNgamma = Get2DArray <TH1D*> (photonTrigN, 3);
  TH1D*** xJgDist = Get2DArray <TH1D*> (photonTrigN, 3);

  TH1D** photonCounts = Get1DArray <TH1D*> (photonTrigN);
  TH2D** etaPhiMap = Get1DArray <TH2D*> (photonTrigN);

  TH1D***** etconeDists = Get4DArray <TH1D*> (3, photonTrigN, 2, 2);
  TH2D**** etconeFCalDists = Get3DArray <TH2D*> (3, photonTrigN, 2);

  TH1D***** etconeEventPlaneDists = Get4DArray <TH1D*> (3, photonTrigN, 2, 3);

  TH1D*** photonEventPlaneDists = Get2DArray <TH1D*> (photonTrigN, 2);

  TH1D*** fcal_et = Get2DArray <TH1D*> (photonTrigN, 2);

  TH1D**** sidebandSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 4);
  TH1D**** showerShapeDists = Get3DArray <TH1D*> (photonTrigN+1, 2, 13);

  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
    photonCounts[iTrig] = (TH1D*)inFile->Get (Form ("photonCounts_%s", photonTrigNames[iTrig].c_str ()));
    etaPhiMap[iTrig] = (TH2D*)inFile->Get (Form ("etaPhiMap_%s", photonTrigNames[iTrig].c_str()));
    for (short iEta = 0; iEta < 2; iEta++) {
      const char* eta = iEta==0?"barrel":"endcap";

      for (short iAngle = 0; iAngle < 3; iAngle++) {
        const char* angle = iAngle==0?"0_pi6":(iAngle==1?"pi6_pi3":"pi3_pi2");
        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          const char etcone = iEtcone==0?'2':(iEtcone==1?'3':'4');
          etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle] = (TH1D*)inFile->Get (Form ("etcone%c0_%s_%s_%s", etcone, photonTrigNames[iTrig].c_str (), eta, angle));
        }
      }

      for (short iCent = 0; iCent < 2; iCent++) {
        const char* cent = iCent==0?"cent":"periph";

        ptSpectrum[iTrig][iEta][iCent] = (TH1D*)inFile->Get (Form ("ptSpectrum_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent));

        for (short iAngle = 0; iAngle < 3; iAngle++) {
          const char* angle = iAngle==0?"0_pi6":(iAngle==1?"pi6_pi3":"pi3_pi2");
          if (iEta == 0 && iCent == 0) {
            xJgDist[iTrig][iAngle] = (TH1D*)inFile->Get (Form ("xJgDist_%s_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent, angle));
            xJgNgamma[iTrig][iAngle] = (TH1D*)inFile->Get (Form ("xJgNgamma_%s_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent, angle));
          }
          else {
            xJgDist[iTrig][iAngle]->Add ((TH1D*)inFile->Get (Form ("xJgDist_%s_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent, angle)));
            xJgNgamma[iTrig][iAngle]->Add ((TH1D*)inFile->Get (Form ("xJgNgamma_%s_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent, angle)));
          }
        }

        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          const char etcone = iEtcone==0?'2':(iEtcone==1?'3':'4');
          etconeDists[iEtcone][iTrig][iEta][iCent] = (TH1D*)inFile->Get (Form ("etcone%c0_%s_%s_%s", etcone, photonTrigNames[iTrig].c_str(), eta, cent));
        }
      }

      for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
        const char etcone = iEtcone==0?'2':(iEtcone==1?'3':'4');
        etconeFCalDists[iEtcone][iTrig][iEta] = (TH2D*)inFile->Get (Form ("etcone%c0FCal_%s_%s", etcone, photonTrigNames[iTrig].c_str(), eta));
      }

      photonEventPlaneDists[iTrig][iEta] = (TH1D*)inFile->Get (Form ("photonEventPlaneDist_%s_%s", photonTrigNames[iTrig].c_str (), eta));

      fcal_et[iTrig][iEta] = (TH1D*)inFile->Get (Form ("fcal_et_%s_%s", photonTrigNames[iTrig].c_str(), eta));

      for (short iSide = 0; iSide < 4; iSide++) {
        char side = iSide==0?'a':(iSide==1?'b':(iSide==2?'c':'d'));
        sidebandSpectrum[iTrig][iEta][iSide] = (TH1D*)inFile->Get (Form ("sidebandSpectrum_%s_%s_%c", photonTrigNames[iTrig].c_str(), eta, side));
      }

      for (short iShower = 0; iShower < 13; iShower++) {
        showerShapeDists[iTrig][iEta][iShower] = (TH1D*)inFile->Get (Form ("showerShapeDist_%s_%s_%s", photonTrigNames[iTrig].c_str(), eta, GetShowerShape (iShower)));
      }
    }
  }

  TFile* inFile_mc = new TFile (Form ("%s/outfile_MC_GammaJet_SC_DP50_70.lst_Tight.root", rootPath.Data ()), "read");
  for (short iEta = 0; iEta < 2; iEta++) {
    showerShapeDists[photonTrigN][iEta][0] = (TH1D*)inFile_mc->Get (Form ("hph_Rhad_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][1] = (TH1D*)inFile_mc->Get (Form ("hph_Rhad1_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][2] = (TH1D*)inFile_mc->Get (Form ("hph_e277_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][3] = (TH1D*)inFile_mc->Get (Form ("hph_Reta_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][4] = (TH1D*)inFile_mc->Get (Form ("hph_Rphi_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][5] = (TH1D*)inFile_mc->Get (Form ("hph_weta1_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][6] = (TH1D*)inFile_mc->Get (Form ("hph_weta2_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][7] = (TH1D*)inFile_mc->Get (Form ("hph_wtots1_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][8] = (TH1D*)inFile_mc->Get (Form ("hph_f1_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][9] = (TH1D*)inFile_mc->Get (Form ("hph_f3_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][10] = (TH1D*)inFile_mc->Get (Form ("hph_fracs1_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][11] = (TH1D*)inFile_mc->Get (Form ("hph_DeltaE_eta%i", iEta));
    showerShapeDists[photonTrigN][iEta][12] = (TH1D*)inFile_mc->Get (Form ("hph_Eratio_eta%i", iEta));
    for (short iShower = 0; iShower < 13; iShower++) {
      TH1D* h = showerShapeDists[photonTrigN][iEta][iShower];
      h->Rebin (showerShapeMCRebinFactors[iShower]);
      //const float integral = h->Integral (h->FindBin (showerShapeBinsLow[iShower]), h->FindBin (showerShapeBinsHigh[iShower]));
      TAxis *axis = h->GetXaxis();
      const double xmin = showerShapeBinsLow[iShower];
      const double xmax = showerShapeBinsHigh[iShower];

      const int bmin = axis->FindBin (xmin); //in your case xmin=-1.5
      const int bmax = axis->FindBin (xmax); //in your case xmax=0.8
      double integral = h->Integral (bmin,bmax);
      integral -= h->GetBinContent (bmin) * (xmin-axis->GetBinLowEdge (bmin)) / axis->GetBinWidth (bmin);
      integral -= h->GetBinContent (bmax) * (axis->GetBinUpEdge (bmax)-xmax) / axis->GetBinWidth (bmax);
      if (integral > 0)
        h->Scale (1. / integral, "width");
    }
  }


  const Color_t colors[4] = {kBlack, kBlue, kRed, 8};
  const TString evtPlaneStrs[3] = {"0 < #Delta#phi < #pi/6 (\"In plane\")", "#pi/6 < #Delta#phi < #pi/3", "#pi/3 < #Delta#phi < #pi/2 (\"Out of plane\")"};
  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
    TCanvas* canvas = new TCanvas (Form ("canvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 600);
    canvas->cd ();

    lumi_int = 0;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot photon counts
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gPad->SetLogy (false);
    photonCounts[iTrig]->GetYaxis ()->SetTitle ("Counts / #mub^{-1}");
    photonCounts[iTrig]->Draw ("e1 x0");
    for (short iRunNumber = 0; iRunNumber < numRunNumbers; iRunNumber++) {
      TText* t = new TText ();
      t->SetTextAngle (89);
      if (photonCounts[iTrig]->GetBinContent (iRunNumber+1) > 0) {
        t->DrawText (iRunNumber+0.25,  0.38*photonCounts[iTrig]->GetMaximum (), to_string(runNumbers[iRunNumber]).c_str ());
        lumi_int += lumis[iRunNumber];
      }
    }
    myText (0.63, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
    myText (0.63, 0.81, kBlack, "Tight photons > 60 GeV", 0.04);
    myText (0.63, 0.74, kBlack, "Etcone30 < 10 GeV", 0.04);
    canvas->SaveAs (Form ("%s/photonCounts/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot pT spectra of photons
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gPad->SetLogy ();
    double max = 0;
    TH1D** ptSpectra = Get1DArray <TH1D*> (2);
    //canvas->SetLeftMargin (-0.1);
    for (int iEta = 0; iEta < 2; iEta++) {
      ptSpectra[iEta] = new TH1D (Form ("ptSpectrum_%s", iEta==0?"barrel":"endcap"), "", 600, 50, 650);
      ptSpectra[iEta]->Sumw2 ();
      for (int iCent = 0; iCent < 2; iCent++) {
        TH1D* thisHist = ptSpectrum[iTrig][iEta][iCent];
        ptSpectra[iEta]->Add (thisHist);
      }
      ptSpectra[iEta]->Rebin(25);
      //ptSpectra[iEta]->Scale (1/pi);
      if (max < ptSpectra[iEta]->GetMaximum ())
        max = ptSpectra[iEta]->GetMaximum ();
    }
    TGraphAsymmErrors** ptSpectraGraphs = Get1DArray <TGraphAsymmErrors*> (2);
    for (int iEta = 0; iEta < 2; iEta++) {
      //for (int iCent = 0; iCent < 2; iCent++) {
      TGraphAsymmErrors* thisGraph = make_graph (ptSpectra[iEta], 0.5);
      ptSpectraGraphs[iEta] = thisGraph;

      thisGraph->SetMarkerColor (colors[iEta]);
      thisGraph->SetLineColor (colors[iEta]);
      thisGraph->SetLineWidth (2);
      //thisGraph->SetLineStyle (iCent==0?2:1);

      thisGraph->GetXaxis ()->SetTitle ("Photon #it{p}_{T} #left[GeV#right]");
      //thisGraph->GetYaxis ()->SetTitle ("Entries / bin");
      thisGraph->GetYaxis ()->SetTitle ("Counts");

      thisGraph->GetXaxis ()->SetRangeUser (50, 550);
      thisGraph->GetYaxis ()->SetRangeUser (0.3, 3*max);

      //thisGraph->GetYaxis ()->SetTitleOffset (0.2);
      //thisGraph->GetYaxis ()->SetTitleSize (0.04);
      //thisGraph->GetYaxis ()->SetLabelSize (0);

      if (iEta == 0)
        thisGraph->Draw ("ap");
      else
        thisGraph->Draw ("p");
    }
    //myText (0.56, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.045);
    myText (0.56, 0.88, kBlack, "#bf{#it{ATLAS}} Preliminary", 0.045);
    myText (0.56, 0.82, kBlack, Form ("#sqrt{s_{NN}} = 5.02 TeV, %.2f nb^{-1}", lumi_int*1e-3), 0.045);
    myMarkerText (0.56, 0.76, kBlack, kFullCircle, "0 < #left|#eta#right| < 1.37", 1.25, 0.045);
    myMarkerText (0.56, 0.70, kBlue, kFullCircle, "1.52 < #left|#eta#right| < 2.37", 1.25, 0.045);
    //myText (0.56, 0.665, kBlack, "Dashed: #Sigma#it{E}_{T}^{FCal} > 2 TeV", 0.045);
    //myText (0.56, 0.605, kBlack, "Solid: #Sigma#it{E}_{T}^{FCal} < 2 TeV", 0.045);
    canvas->SaveAs (Form ("%s/ptSpectrum/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));
    Delete1DArray (ptSpectra, 2);
    Delete1DArray (ptSpectraGraphs, 2);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot xJg distributions 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    max = 0;
    gPad->SetLogy (false);
    for (short iAngle = 0; iAngle < 3; iAngle++) {
      TH1D* thisHist = xJgDist[iTrig][iAngle];
      thisHist->Rebin (8);
      thisHist->Scale (1. / xJgNgamma[iTrig][iAngle]->Integral (), "width");
      if (max < thisHist->GetMaximum ())
        max = thisHist->GetMaximum ();
    }
    
    for (short iAngle = 0; iAngle < 3; iAngle++) {
      TH1D* thisHist = xJgDist[iTrig][iAngle];

      TGraphAsymmErrors* thisGraph = make_graph (thisHist);
      deltaize (thisGraph, -0.01+0.01*iAngle, false);

      thisGraph->GetYaxis ()->SetRangeUser (0, 1.1*max);
      //thisGraph->SetMarkerStyle (iCent==0?24:20);
      thisGraph->SetMarkerColor (colors[iAngle]);
      thisGraph->SetLineColor (colors[iAngle]);
      thisGraph->GetXaxis ()->SetTitle ("x_{J}^{#gamma}");
      thisGraph->GetYaxis ()->SetTitle ("1 / N_{#gamma}  dN_{pairs}^{#gamma+jet} / dx_{J}^{#gamma}");
      if (iAngle == 0)
        ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("ap");
      else
        ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("p");

      myText (0.6, 0.805-0.06*iAngle, colors[iAngle], evtPlaneStrs[iAngle].Data (), 0.04);

      if (thisGraph) delete thisGraph;
    }
    
    
    myText (0.6, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
    myText (0.6, 0.62, kBlack, "Tight, iso photons", 0.04);
    myText (0.6, 0.56, kBlack, "#it{p}_{T}^{#gamma} > 60 GeV, #it{p}_{T}^{Jet} > 30 GeV", 0.04);
    myText (0.6, 0.5, kBlack, "Anti-#it{k}_{T} R=0.4 HI Jets",  0.04);
    canvas->SaveAs (Form ("%s/xJgDist/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot ET cone distributions
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TString etconeStrs[3] = {"Etcone20", "Etcone30", "Etcone40"};
    gPad->SetLogy (false);
    for (int iEta = 0; iEta < 2; iEta++) {
      max = 0;
      for (short iCent = 0; iCent < 2; iCent++) {
        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          if (etconeDists[iEtcone][iTrig][iEta][iCent]->GetMaximum () > max)
            max = etconeDists[iEtcone][iTrig][iEta][iCent]->GetMaximum ();
        }
      }

      for (short iCent = 0; iCent < 2; iCent++) {
        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          TH1D* thisHist = etconeDists[iEtcone][iTrig][iEta][iCent];
          thisHist->GetYaxis ()->SetRangeUser (0.5, 1.1*max);
          thisHist->SetLineStyle (iCent==0?2:1);
          thisHist->SetMarkerColor (colors[iEtcone]);
          thisHist->SetLineColor (colors[iEtcone]);

          thisHist->GetXaxis ()->SetTitle ("Etcone Value #left[GeV#right]");
          thisHist->GetYaxis ()->SetTitle ("Counts");

          if (iEtcone == 0 && iCent == 0) {
            thisHist->Draw ("hist");
            myText (0.56, 0.54, kBlack, "Dashed: #Sigma#it{E}_{T}^{FCal} > 2 TeV", 0.04);
            myText (0.56, 0.485, kBlack, "Solid: #Sigma#it{E}_{T}^{FCal} < 2 TeV", 0.04);
          }
          else
            thisHist->Draw ("same hist");

          if (iCent == 0)
            myText (0.56, 0.69-0.05*iEtcone, colors[iEtcone], etconeStrs[iEtcone], 0.04);
        }
      }
      //myText (0.6, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
      myText (0.56, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.56, 0.82, kBlack, (string(iEta==0?"Barrel":"Endcaps") + ", tight photons > 20 GeV").c_str(), 0.04);
      myText (0.56, 0.76, kBlack, Form ("#sqrt{s_{NN}} = 5.02 TeV, %.1f #mub^{-1}", lumi_int), 0.04);
      canvas->SaveAs (Form ("%s/etcones/%s_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"));
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot ET cone distributions, binned by angle w.r.t. the event plane
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* etconeEventPlaneCanvas = new TCanvas (Form ("etconeEventPlaneCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 1000);
    etconeEventPlaneCanvas->Divide (1, 3, -1);
    etconeEventPlaneCanvas->cd ();
    gPad->SetLogy (false);

    for (int iEta = 0; iEta < 2; iEta++) {
      for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
        etconeEventPlaneCanvas->cd (iEtcone+1);
        gPad->SetRightMargin (0.01);
        gPad->SetLeftMargin (0.16);
        max = 0;
        for (short iAngle = 0; iAngle < 3; iAngle++) {
          if (etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle]->GetMaximum () > max)
            max = etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle]->GetMaximum ();
        }
        for (short iAngle = 0; iAngle < 3; iAngle++) {
       
          TH1D* thisHist = etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle];
          thisHist->GetYaxis ()->SetRangeUser (0, 1.1*max);
          thisHist->SetMarkerColor (colors[iAngle]);
          thisHist->SetLineColor (colors[iAngle]);
          thisHist->SetMarkerSize (0.5);

          thisHist->GetXaxis ()->SetTitle ("Isolation Etcone Value #left[GeV#right]");
          thisHist->GetYaxis ()->SetTitle ("Counts");

          thisHist->GetXaxis ()->SetTitleOffset (0.9);
          thisHist->GetYaxis ()->SetTitleOffset (1.7*gPad->GetHNDC ());
          thisHist->GetYaxis ()->CenterTitle (true);
          thisHist->GetXaxis ()->SetTitleSize (0.024/gPad->GetHNDC ());
          thisHist->GetYaxis ()->SetTitleSize (0.024/gPad->GetHNDC ());
          thisHist->GetXaxis ()->SetLabelSize (0.024/gPad->GetHNDC ());
          thisHist->GetYaxis ()->SetLabelSize (0.024/gPad->GetHNDC ());

          if (iAngle == 0) {
            thisHist->Draw ("hist");
          }
          else
            thisHist->Draw ("same hist");

          myText (0.19, 0.88-0.09*iAngle, colors[iAngle], Form ("RMS=%.2f#pm%.2f", thisHist->GetRMS (), thisHist->GetRMSError ()), 0.024/gPad->GetHNDC ());
          myText (0.19, 0.58-0.09*iAngle, colors[iAngle], Form ("#mu=%.2f#pm%.2f", thisHist->GetMean (), thisHist->GetMeanError ()), 0.024/gPad->GetHNDC ());

          if (iAngle == 0)
            myText (0.55, 0.88, kBlack, etconeStrs[iEtcone], 0.024/gPad->GetHNDC ());
          if (iEtcone == 1 && iAngle == 0) {
            myText (0.55, 0.76, kBlack, (TString(iEta==0?"Barrel":"Endcaps") + ", tight photons").Data (), 0.024/gPad->GetHNDC ());
            myText (0.55, 0.64, kBlack, "#it{p}_{T}^{#gamma} > 20 GeV", 0.024/gPad->GetHNDC ());
          }
          else if (iEtcone == 2 && iAngle == 0) {
            myText (0.55, 0.76, kBlack, photonTrigNames[iTrig].c_str (), 0.024/gPad->GetHNDC ());
            myText (0.55, 0.64, kBlack, Form ("%.1f #mub^{-1}", lumi_int), 0.024/gPad->GetHNDC ());
          }
          if (iEtcone == 0)
            myText (0.55, 0.76-0.12*iAngle, colors[iAngle], evtPlaneStrs[iAngle], 0.024/gPad->GetHNDC ());
        }
      }
      etconeEventPlaneCanvas->SaveAs (Form ("%s/etcones_eventplane/%s_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"));
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot mean, sigma ET cone as a function of FCal Sum ET
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* etconeFCalCanvas = new TCanvas ("etconeFCalCanvas", "", 800, 1000);
    etconeFCalCanvas->Divide (1, 3, -1);
    etconeFCalCanvas->cd ();
    gPad->SetLogy (false);

    for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
      etconeFCalCanvas->cd (iEtcone+1);
      gPad->SetRightMargin (0.01);
      gPad->SetLeftMargin (0.16);

      TGraphAsymmErrors*** graphArr = Get2DArray <TGraphAsymmErrors*> (2, 2);
      max = 0;
      for (short iEta = 0; iEta < 2; iEta++) {
        TH2D* thisHist = etconeFCalDists[iEtcone][iTrig][iEta];
        TGraphAsymmErrors* thisGraph1 = new TGraphAsymmErrors (); // mean
        TGraphAsymmErrors* thisGraph2 = new TGraphAsymmErrors (); // std. dev.

        for (int ix = 1; ix <= thisHist->GetNbinsX (); ix++) {
          TH1D* proj = thisHist->ProjectionY ("_px", ix, ix);

          TF1* fit = new TF1 ("fit", "gaus(0)", -1.5*proj->GetRMS(), 1.5*proj->GetRMS ());

          proj->Fit (fit, "R0Q");
          float m = fit->GetParameter (1);
          float merr = fit->GetParError (1);
          float s = fit->GetParameter (2);
          float serr = fit->GetParError (2);

          thisGraph1->SetPoint (ix-1, thisHist->GetXaxis ()->GetBinCenter (ix), m); 
          thisGraph1->SetPointError (ix-1, 0.5 * thisHist->GetXaxis ()->GetBinWidth (ix), 0.5 * thisHist->GetXaxis ()->GetBinWidth (ix), merr, merr);
          if (m + merr > max)
            max = m + merr;

          thisGraph2->SetPoint (ix-1, thisHist->GetXaxis ()->GetBinCenter (ix), s); 
          thisGraph2->SetPointError (ix-1, 0.5 * thisHist->GetXaxis ()->GetBinWidth (ix), 0.5 * thisHist->GetXaxis ()->GetBinWidth (ix), serr, serr);
          if (s + serr > max)
            max = s + serr; 

          graphArr[iEta][0] = thisGraph1;
          graphArr[iEta][1] = thisGraph2;

          if (fit) delete fit;
          if (proj) delete proj;
        }
      }

      for (short iEta = 0; iEta < 2; iEta++) {
        TGraphAsymmErrors* thisGraph1 = graphArr[iEta][0];
        TGraphAsymmErrors* thisGraph2 = graphArr[iEta][1];

        deltaize (thisGraph1, -0.06+0.12*iEta);
        deltaize (thisGraph2, -0.06+0.12*iEta);

        thisGraph1->GetYaxis ()->SetRangeUser (0, 1.1*max);
        thisGraph2->GetYaxis ()->SetRangeUser (0, 1.1*max);

        thisGraph1->SetMarkerStyle (20);
        thisGraph2->SetMarkerStyle (24);

        thisGraph1->SetLineColor (colors[iEta]);
        thisGraph1->SetMarkerColor (colors[iEta]);
        thisGraph2->SetLineColor (colors[iEta]);
        thisGraph2->SetMarkerColor (colors[iEta]);

        thisGraph1->SetMarkerSize (0.8);
        thisGraph2->SetMarkerSize (0.8);

        thisGraph1->GetXaxis ()->SetTitle ("FCal #Sigma#it{E}_{T} #left[TeV#right]");
        thisGraph1->GetYaxis ()->SetTitle ("#mu, #sigma #left[GeV#right]");

        thisGraph1->GetXaxis ()->SetTitleOffset (0.9);
        thisGraph1->GetYaxis ()->SetTitleOffset (1.7*gPad->GetHNDC ());
        thisGraph1->GetYaxis ()->CenterTitle (true);
        thisGraph1->GetXaxis ()->SetTitleSize (0.032/gPad->GetHNDC ());
        thisGraph1->GetYaxis ()->SetTitleSize (0.032/gPad->GetHNDC ());
        thisGraph1->GetXaxis ()->SetLabelSize (0.032/gPad->GetHNDC ());
        thisGraph1->GetYaxis ()->SetLabelSize (0.032/gPad->GetHNDC ());

        if (iEta == 0)
          ( (TGraphAsymmErrors*)thisGraph1->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)thisGraph1->Clone ())->Draw ("p");
        ( (TGraphAsymmErrors*)thisGraph2->Clone ())->Draw ("p");

        if (iEta == 0) {
          float texty = gPad->GetBottomMargin() + 0.9*(1-gPad->GetBottomMargin()-gPad->GetTopMargin());
          myText (0.2, texty, kBlack, Form ("Etcone%i0", iEtcone+2), 0.024/gPad->GetHNDC ());
        }

        if (iEtcone == 0) {
          float texty = gPad->GetBottomMargin() + (0.8-0.1*iEta)*(1-gPad->GetBottomMargin()-gPad->GetTopMargin());
          myText (0.2, texty, colors[iEta], iEta==0?"Barrel":"Endcaps", 0.024/gPad->GetHNDC ());
        }
        else if (iEtcone == 1) {
          float texty = gPad->GetBottomMargin() + 0.8*(1-gPad->GetBottomMargin()-gPad->GetTopMargin());
          myMarkerText (0.225, texty, kBlack, 20, "#mu = <etconeXX>", 1.25, 0.024/gPad->GetHNDC ());
          texty = gPad->GetBottomMargin() + 0.7*(1-gPad->GetBottomMargin()-gPad->GetTopMargin());
          myMarkerText (0.225, texty, kBlack, 24, "#sigma = #sigma (etconeXX)", 1.25, 0.024/gPad->GetHNDC ());
        }
        else if (iEtcone == 2) {
          float texty = gPad->GetBottomMargin() + 0.8*(1-gPad->GetBottomMargin()-gPad->GetTopMargin());
          myText (0.2, texty, kBlack, photonTrigNames[iTrig].c_str (), 0.024/gPad->GetHNDC ());
          texty = gPad->GetBottomMargin() + 0.7*(1-gPad->GetBottomMargin()-gPad->GetTopMargin());
          myText (0.2, texty, kBlack, "Tight photons > 20 GeV", 0.024/gPad->GetHNDC ());
        }

      }
      Delete2DArray (graphArr, 2, 2);
    }
    etconeFCalCanvas->SaveAs (Form ("%s/etcones_fcal/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot FCal sum Et distributions
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    canvas->cd ();
    max = 0;
    for (short iEta = 0; iEta < 2; iEta++) {
      fcal_et[iTrig][iEta]->Rebin(4);
      if (max < fcal_et[iTrig][iEta]->GetMaximum ())
        max = fcal_et[iTrig][iEta]->GetMaximum ();
    }    
    
    for (short iEta = 0; iEta < 2; iEta++) {
      TH1D* thisHist = fcal_et[iTrig][iEta];
      thisHist->GetYaxis ()->SetRangeUser (0.5, 1.3*max);
      thisHist->SetMarkerColor (colors[iEta]);
      thisHist->SetLineColor (colors[iEta]);

      thisHist->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} #left[TeV#right]");
      thisHist->GetYaxis ()->SetTitle ("Counts");

      if (iEta == 0)
        thisHist->Draw ("hist");
      else
        thisHist->Draw ("same hist");

      myText (0.2, 0.74-0.06*iEta, colors[iEta], iEta==0?"Barrel":"Endcaps", 0.04);
    }    
    //myText (0.2, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
    myText (0.65, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
    myText (0.2, 0.88, kBlack, "Tight photons, #it{p}_{T}^{#gamma} > 20 GeV", 0.04);
    myText (0.2, 0.82, kBlack, Form ("#sqrt{s_{NN}} = 5.02 TeV, %.1f #mub^{-1}", lumi_int), 0.04);
    canvas->SaveAs (Form ("%s/fcal_et/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot FCal sum Et distributions
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    canvas->cd ();
    max = 0;
    for (short iEta = 0; iEta < 2; iEta++) {
      photonEventPlaneDists[iTrig][iEta]->Rebin(10);
      if (max < photonEventPlaneDists[iTrig][iEta]->GetMaximum ())
        max = photonEventPlaneDists[iTrig][iEta]->GetMaximum ();
    }    
    
    for (short iEta = 0; iEta < 2; iEta++) {
      TH1D* thisHist = photonEventPlaneDists[iTrig][iEta];
      TGraphAsymmErrors* thisGraph = make_graph (thisHist);
      deltaize (thisGraph, (-0.5+iEta)*pi/90., false);

      thisGraph->SetMarkerColor (colors[iEta]);
      thisGraph->SetLineColor (colors[iEta]);

      thisGraph->GetXaxis ()->SetTitle ("#phi - #Psi_{2}");
      thisGraph->GetYaxis ()->SetTitle ("Counts");

      thisGraph->GetXaxis ()->SetRangeUser (0, pi/2);
      thisGraph->GetYaxis ()->SetRangeUser (0.5, 1.1*max);

      if (iEta == 0)
        thisGraph->Draw ("ap");
      else
        thisGraph->Draw ("p");

      myText (0.25, 0.28-0.06*iEta, colors[iEta], iEta==0?"Barrel":"Endcaps", 0.04);
    }    
    myText (0.6, 0.32, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
    myText (0.6, 0.26, kBlack, "Tight, isolated photons", 0.04);
    myText (0.6, 0.20, kBlack, "#it{p}_{T}^{#gamma} > 60 GeV", 0.04);
    canvas->SaveAs (Form ("%s/photonEventPlaneDists/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot eta phi map of photons
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* etaPhiMapCanvas = new TCanvas (Form ("etaPhiMapCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 600);
    etaPhiMapCanvas->cd ();
    etaPhiMapCanvas->SetRightMargin (0.18);
    etaPhiMapCanvas->SetLeftMargin (0.12);
    gStyle->SetPalette (kRainBow);
    etaPhiMap[iTrig]->GetXaxis ()->SetTitle ("#eta");
    etaPhiMap[iTrig]->GetYaxis ()->SetTitle ("#phi");
    etaPhiMap[iTrig]->GetZaxis ()->SetTitle ("Counts");
    etaPhiMap[iTrig]->GetXaxis ()->SetTitleOffset (1);
    etaPhiMap[iTrig]->GetYaxis ()->SetTitleOffset (1);
    etaPhiMap[iTrig]->Draw ("colz");
    //myText (0.6, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
    //myText (0.2, 0.32, kBlack, Form ("#sqrt{s_{NN}} = 5.02 TeV, %.1f #mub^{-1}", lumi_int), 0.04);
    //myText (0.2, 0.26, kBlack, "Tight photons, #it{p}_{T}^{#gamma} > 20 GeV", 0.04);

    etaPhiMapCanvas->SaveAs (Form ("%s/etaPhiMaps/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ())); 


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determine and plot purity factors
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* purityCanvas = new TCanvas (Form ("purityCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 600);
    purityCanvas->cd ();
    TGraphAsymmErrors** purityGraphs = Get1DArray<TGraphAsymmErrors*> (2);
    for (short iEta = 0; iEta < 2; iEta++) {
      TH1D* purityHist = new TH1D (Form ("purityHist_%s_%s", photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"), "", numpbins, pbins);
      for (int ix = 1; ix <= purityHist->GetNbinsX (); ix++) {
        const float A = sidebandSpectrum[iTrig][iEta][0]->GetBinContent (ix);
        const float B = sidebandSpectrum[iTrig][iEta][1]->GetBinContent (ix);
        const float C = sidebandSpectrum[iTrig][iEta][2]->GetBinContent (ix);
        const float D = sidebandSpectrum[iTrig][iEta][3]->GetBinContent (ix);
        const float eA = sidebandSpectrum[iTrig][iEta][0]->GetBinError (ix);
        const float eB = sidebandSpectrum[iTrig][iEta][1]->GetBinError (ix);
        const float eC = sidebandSpectrum[iTrig][iEta][2]->GetBinError (ix);
        const float eD = sidebandSpectrum[iTrig][iEta][3]->GetBinError (ix);
        if (A != 0. && D != 0.) {
          float purity = 1 - B*C / (A*D);
          float ePurity = sqrt (pow (purity * eA/A, 2) + pow (purity * eD/D, 2) + pow (C*eB / (A*D), 2) + pow (B*eC / (A*D), 2));
          purityHist->SetBinContent (ix, purity); 
          purityHist->SetBinError (ix, ePurity);
        }
      }

      purityGraphs[iEta] = make_graph (purityHist); 
      TGraphAsymmErrors* thisGraph = purityGraphs[iEta];
      deltaize (thisGraph, -1+(2*iEta), false);

      thisGraph->SetMarkerColor (colors[iEta]);
      thisGraph->SetLineColor (colors[iEta]);

      thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle ("1 - BC/AD");

      thisGraph->GetXaxis ()->SetRangeUser (0, 200);
      thisGraph->GetYaxis ()->SetRangeUser (0, 1.1);
      thisGraph->GetXaxis ()->SetTitleOffset (1.4);
      thisGraph->GetYaxis ()->SetTitleOffset (1.1);

      if (iEta == 0)
        thisGraph->Draw ("ap");
      else
        thisGraph->Draw ("p");

      if (iEta == 0) {
        myText (0.66, 0.24, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
        //myText (0.2, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
        myText (0.2, 0.88, kBlack, Form ("#sqrt{s_{NN}} = 5.02 TeV, %.1f #mub^{-1}", lumi_int), 0.04);
      }
      myText (0.2, 0.82-0.06*iEta, colors[iEta], iEta==0?"0 < |#eta| < 1.37":"1.52 < |#eta| < 2.37", 0.04);

    }
    purityCanvas->SaveAs (Form ("%s/purities/%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ())); 
    Delete1DArray (purityGraphs, 2);
    


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot sideband spectra and ratios
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* sidebandsCanvas = new TCanvas (Form ("sidebandsCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 800);
    sidebandsCanvas->cd ();
    const double padRatio = 1.2; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
    const double dPadY = 1.0/ (padRatio+1.0);
    const double uPadY = 1.0 - dPadY;
    TPad* topPad = new TPad (Form ("topPad_%s", photonTrigNames[iTrig].c_str ()), "", 0, dPadY, 1, 1);
    TPad* bottomPad = new TPad (Form ("bottomPad_%s", photonTrigNames[iTrig].c_str ()), "", 0, 0, 1, dPadY);
    topPad->SetBottomMargin (0);
    topPad->SetRightMargin (-0.20);
    bottomPad->SetTopMargin (0);
    bottomPad->SetBottomMargin (0.30);
    bottomPad->SetRightMargin (-0.20);
    topPad->Draw ();
    bottomPad->Draw ();
    for (short iEta = 0; iEta < 2; iEta++) {
      max = 0;
      for (short iSide = 0; iSide < 4; iSide++)
        max = std::max (sidebandSpectrum[iTrig][iEta][iSide]->GetMaximum (), max);

      topPad->cd ();
      gPad->SetLogy ();
      for (short iSide = 0; iSide < 4; iSide++) {
        TH1D* thisHist = sidebandSpectrum[iTrig][iEta][iSide];
        TGraphAsymmErrors* thisGraph = make_graph (thisHist);
        deltaize (thisGraph, 0.5*(-3+2*iSide), false);

        thisGraph->GetXaxis ()->SetRangeUser (50, 205);
        thisGraph->GetYaxis ()->SetRangeUser (0.5, 2*max);

        thisGraph->SetLineColor (colors[iSide]);
        thisGraph->SetMarkerColor (colors[iSide]);

        thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
        thisGraph->GetYaxis ()->SetTitle ("N");

        thisGraph->GetXaxis ()->SetTitleSize (0.04/uPadY);
        thisGraph->GetYaxis ()->SetTitleSize (0.04/uPadY);
        thisGraph->GetXaxis ()->SetLabelSize (0.04/uPadY);
        thisGraph->GetYaxis ()->SetLabelSize (0.04/uPadY);
        thisGraph->GetXaxis ()->SetTitleOffset (1.4);
        thisGraph->GetYaxis ()->SetTitleOffset (1.4*uPadY);

        if (iSide == 0)
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("p");

        TString sidebandStr = "";
        switch (iSide) {
          case 0: sidebandStr = sidebandStr + "A (tight, etcone 30 < 10 GeV)"; break;
          case 1: sidebandStr = sidebandStr + "B (tight, etcone 30 > 12 GeV)"; break;
          case 2: sidebandStr = sidebandStr + "C (non-tight, etcone 30 < 10 GeV)"; break;
          case 3: sidebandStr = sidebandStr + "D (non-tight, etcone30 > 12 GeV)"; break;
        }
        myText (0.2, 0.3-0.065*iSide, colors[iSide], sidebandStr, 0.03/uPadY);

        if (thisGraph) delete thisGraph;
      }
      myText (0.65, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.03/uPadY);
      myText (0.52, 0.76, kBlack, Form ("#sqrt{s_{NN}} = 5.02 TeV, %.1f #mub^{-1}", lumi_int), 0.03/uPadY);
      myText (0.2, 0.36, kBlack, /*(photonTrigNames[iTrig] + */(iEta==0?"Barrel, 0 < |#eta| < 1.37":"Endcaps, 1.52 < |#eta| < 2.37"), 0.03/uPadY);

      bottomPad->cd ();
      max = 0;
      for (short iSide = 1; iSide < 4; iSide++) {
        TH1D* thisHist = sidebandSpectrum[iTrig][iEta][iSide];
        thisHist->Divide (sidebandSpectrum[iTrig][iEta][0]);
        max = std::max (max, thisHist->GetMaximum ());
      }
      for (short iSide = 1; iSide < 4; iSide++) {
        TH1D* thisHist = sidebandSpectrum[iTrig][iEta][iSide];
        TGraphAsymmErrors* thisGraph = make_graph (thisHist);
        deltaize (thisGraph, 0.5*(-3+2*iSide), false);

        thisGraph->GetXaxis ()->SetRangeUser (50, 205);
        thisGraph->GetYaxis ()->SetRangeUser (0, 1);

        thisGraph->SetLineColor (colors[iSide]);
        thisGraph->SetMarkerColor (colors[iSide]);

        thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
        thisGraph->GetYaxis ()->SetTitle ("N / N_{A}");

        thisGraph->GetXaxis ()->SetTitleSize (0.04/dPadY);
        thisGraph->GetYaxis ()->SetTitleSize (0.04/dPadY);
        thisGraph->GetXaxis ()->SetLabelSize (0.04/dPadY);
        thisGraph->GetYaxis ()->SetLabelSize (0.04/dPadY);
        thisGraph->GetXaxis ()->SetTitleOffset (1.4);
        thisGraph->GetYaxis ()->SetTitleOffset (1.4*dPadY);

        if (iSide == 1)
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)thisGraph->Clone ())->Draw ("p");

        const TString sidebandStr = TString (iSide==1?"B":(iSide==2?"C":"D")) + " / A";
        myText (0.78, 0.95-0.07*iSide, colors[iSide], sidebandStr, 0.03/dPadY);

        if (thisGraph) delete thisGraph;
      }
      sidebandsCanvas->SaveAs (Form ("%s/sidebands/%s_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"));
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot shower shapes
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* showerShapeCanvas = new TCanvas (Form ("showerShapeCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 800);
    showerShapeCanvas->Divide (3, 3);
    showerShapeCanvas->cd ();
    for (short iEta = 0; iEta < 2; iEta++) {
      int iCanvas = 1;
      for (short iShower = 3; iShower < 13; iShower++) {
        if (iShower == 2 || iShower == 7 || iShower == 10) continue; // don't plot e277 or wtots1 for now!
 
        showerShapeCanvas->cd (iCanvas++);

        TH1D* thisHist = showerShapeDists[iTrig][iEta][iShower];
        thisHist->Rebin (20);
        const float integral = thisHist->Integral ();
        if (integral > 0)
          thisHist->Scale (1. / integral, "width");
        TH1D* thisHist_mc = showerShapeDists[photonTrigN][iEta][iShower];
        thisHist_mc->GetXaxis ()->SetRangeUser (showerShapeBinsLow[iShower], showerShapeBinsHigh[iShower]);
        max = std::max (thisHist->GetMaximum (), thisHist_mc->GetMaximum ());

        thisHist->GetYaxis ()->SetRangeUser (0, 1.01*max);
        thisHist_mc->GetYaxis ()->SetRangeUser (0, 1.01*max);

        thisHist_mc->GetXaxis ()->SetTitle (GetShowerShape (iShower));
        thisHist_mc->GetYaxis ()->SetTitle ("Counts");

        thisHist_mc->GetXaxis ()->SetTitleSize (0);
        if (iCanvas == 1)
          thisHist_mc->GetYaxis ()->SetTitleSize (0.12);
        else
          thisHist_mc->GetYaxis ()->SetTitleSize (0);

        thisHist_mc->GetXaxis ()->SetTitleOffset (1);
        thisHist_mc->GetYaxis ()->SetTitleOffset (0.4);

        thisHist_mc->GetXaxis ()->SetLabelSize (0.1);
        thisHist_mc->GetYaxis ()->SetLabelSize (0);

        thisHist_mc->GetXaxis ()->SetNdivisions (404);
        thisHist_mc->GetYaxis ()->SetNdivisions (404);

        thisHist_mc->SetLineWidth (2);
        thisHist_mc->SetLineColor (kRed);

        thisHist_mc->Draw ("hist");
        thisHist->Draw ("same hist");

        myText (0.2, 0.88, kBlack, GetShowerShape (iShower), 0.1);
      }
      showerShapeCanvas->cd (9);
      gPad->Clear ();
      //myText (0.1, 0.9, kBlack, photonTrigNames[iTrig].c_str (), 0.1);
      myText (0.1, 0.9, kBlack, "#bf{#it{ATLAS}} Internal", 0.1);
      myText (0.1, 0.75, kBlack, "Tight, iso photons", 0.1);
      myText (0.1, 0.6, kBlack, (iEta == 0 ? "0 < |#eta| < 1.37" : "1.37 < |#eta| < 2.37"), 0.1);
      myText (0.1, 0.35, kBlack, Form ("Data, %.0f #mub^{-1}", lumi_int), 0.1);
      myText (0.1, 0.2, kRed, "Pythia8 + HIJING", 0.1);

      showerShapeCanvas->SaveAs (Form ("%s/showerShapes/%s_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"));
    }

    if (canvas)
      delete canvas;
    if (etconeEventPlaneCanvas)
      delete etconeEventPlaneCanvas;
    if (etconeFCalCanvas)
      delete etconeFCalCanvas;
    if (etaPhiMapCanvas)
      delete etaPhiMapCanvas;
    if (purityCanvas)
      delete purityCanvas;
    if (sidebandsCanvas)
      delete sidebandsCanvas;
    if (showerShapeCanvas)
      delete showerShapeCanvas;
  }
}

} // end namespace
