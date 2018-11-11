#include "PhotonAnalysisHist.h"
#include "Params.h"

#include <Utils.h>
#include <GlobalParams.h>
#include <ArrayTemplates.h>
#include <Trigger.h>

#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

#include <iostream>
#include <string>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace std;
using namespace atlashi;

namespace offlineAnalyses {

const char* GetShowerShape (const int iShowerShape) {
  switch (iShowerShape) {
   case 0:  return "rhad1";
   case 1:  return "rhad";
   case 2:  return "e277";
   case 3:  return "reta";
   case 4:  return "rphi";
   case 5:  return "weta1";
   case 6:  return "weta2";
   case 7:  return "wtots1";
   case 8:  return "f1";
   case 9:  return "fracs1"; // fside
   case 10: return "deltae";
   case 11: return "eratio";
   default: return "";
  }
}


void PhotonAnalysisHist () {

  SetAtlasStyle ();
  SetupDirectories ("PhotonAnalysis/", "OfflineAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TH1D**** ptSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 2);
  TH1D*** xJgDist = Get2DArray <TH1D*> (photonTrigN, 2);

  TH2D** etaPhiMap = Get1DArray <TH2D*> (photonTrigN);

  TH1D**** etcone20 = Get3DArray <TH1D*> (photonTrigN, 2, 2);
  TH1D**** etcone30 = Get3DArray <TH1D*> (photonTrigN, 2, 2);
  TH1D**** etcone40 = Get3DArray <TH1D*> (photonTrigN, 2, 2);

  TH1D*** fcal_et = Get2DArray <TH1D*> (photonTrigN, 2);

  TH1D**** sidebandSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 4);
  TH1D**** showerShapeDists = Get3DArray <TH1D*> (photonTrigN+1, 2, 12);

  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
   etaPhiMap[iTrig] = (TH2D*)inFile->Get (Form ("etaPhiMap_%s", photonTrigNames[iTrig].c_str()));
   for (short iEta = 0; iEta < 2; iEta++) {
    const char* eta = iEta==0?"barrel":"endcap";

    for (short iCent = 0; iCent < 2; iCent++) {
     const char* cent = iCent==0?"cent":"periph";

     ptSpectrum[iTrig][iEta][iCent] = (TH1D*)inFile->Get (Form ("ptSpectrum_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent));
     if (iEta == 0)
      xJgDist[iTrig][iCent] = (TH1D*)inFile->Get (Form ("xJgDist_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent));
     else
      xJgDist[iTrig][iCent]->Add ((TH1D*)inFile->Get (Form ("xJgDist_%s_%s_%s", photonTrigNames[iTrig].c_str (), eta, cent)));

     etcone20[iTrig][iEta][iCent] = (TH1D*)inFile->Get (Form ("etcone20_%s_%s_%s", photonTrigNames[iTrig].c_str(), eta, cent));
     etcone30[iTrig][iEta][iCent] = (TH1D*)inFile->Get (Form ("etcone30_%s_%s_%s", photonTrigNames[iTrig].c_str(), eta, cent));
     etcone40[iTrig][iEta][iCent] = (TH1D*)inFile->Get (Form ("etcone40_%s_%s_%s", photonTrigNames[iTrig].c_str(), eta, cent));
    }

    fcal_et[iTrig][iEta] = (TH1D*)inFile->Get (Form ("fcal_et_%s_%s", photonTrigNames[iTrig].c_str(), eta));

    for (short iSide = 0; iSide < 4; iSide++) {
     char side = iSide==0?'a':(iSide==1?'b':(iSide==2?'c':'d'));
     sidebandSpectrum[iTrig][iEta][iSide] = (TH1D*)inFile->Get (Form ("sidebandSpectrum_%s_%s_%c", photonTrigNames[iTrig].c_str(), eta, side));
    }

    for (short iShower = 0; iShower < 12; iShower++) {
     showerShapeDists[iTrig][iEta][iShower] = (TH1D*)inFile->Get (Form ("showerShapeDist_%s_%s_%s", photonTrigNames[iTrig].c_str(), eta, GetShowerShape (iShower)));
    }
   }
  }

  TFile* inFile_mc = new TFile (Form ("%s/outfile_MC_GammaJet_SC_DP50_70.lst_NoCut.root", rootPath.Data ()), "read");
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
   showerShapeDists[photonTrigN][iEta][9] = (TH1D*)inFile_mc->Get (Form ("hph_fracs1_eta%i", iEta));
   showerShapeDists[photonTrigN][iEta][10] = (TH1D*)inFile_mc->Get (Form ("hph_DeltaE_eta%i", iEta));
   showerShapeDists[photonTrigN][iEta][11] = (TH1D*)inFile_mc->Get (Form ("hph_Eratio_eta%i", iEta));
   for (short iShower = 0; iShower < 12; iShower++) {
    const float integral = showerShapeDists[photonTrigN][iEta][iShower]->Integral ();
    if (integral > 0)
     showerShapeDists[photonTrigN][iEta][iShower]->Scale (1. / integral, "width");
   }
  }

  const Color_t colors[3] = {kBlack, kBlue, kRed};
  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
   TCanvas* canvas = new TCanvas (Form ("canvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 600);
   canvas->cd ();
   gPad->SetLogy ();

   float max = 0;
   for (int iEta = 0; iEta < 2; iEta++) {
    for (int iCent = 0; iCent < 2; iCent++) {
     TH1D* thisHist = ptSpectrum[iTrig][iEta][iCent];
     if (max < thisHist->GetMaximum ())
      max = thisHist->GetMaximum ();
    }
   }

   float xJgScales[2] = {0, 0};
   for (int iEta = 0; iEta < 2; iEta++) {
    for (int iCent = 0; iCent < 2; iCent++) {
     TH1D* thisHist = ptSpectrum[iTrig][iEta][iCent];

     xJgScales[iCent] += thisHist->Integral ();

     thisHist->GetYaxis ()->SetRangeUser (0.5, 1.1*max);
     thisHist->SetMarkerColor (colors[iEta]);
     thisHist->SetLineColor (colors[iEta]);
     thisHist->SetLineStyle (iCent==0?2:1);
     thisHist->GetXaxis ()->SetTitle ("Offline Tight #it{p}_{T}^{#gamma} #left[GeV#right]");
     thisHist->GetYaxis ()->SetTitle ("Counts");
     if (iEta == 0 && iCent == 0)
      thisHist->Draw ("hist");
     else
      thisHist->Draw ("same hist");
    }
   }
   myText (0.6, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
   myText (0.6, 0.805, kBlack, "Barrel", 0.04);
   myText (0.6, 0.745, kBlue, "Endcaps", 0.04);
   myText (0.6, 0.665, kBlack, "Dashed: #Sigma#it{E}_{T}^{FCal} > 2 TeV", 0.04);
   myText (0.6, 0.605, kBlack, "Solid: #Sigma#it{E}_{T}^{FCal} < 2 TeV", 0.04);
   canvas->SaveAs (Form ("%s/ptSpectrum_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


   max = 0;
   gPad->SetLogy (false);
   for (int iCent = 0; iCent < 2; iCent++) {
    TH1D* thisHist = xJgDist[iTrig][iCent];
    thisHist->Rebin (8);
    thisHist->Scale (1. / xJgScales[iCent], "width");
    if (max < thisHist->GetMaximum ())
     max = thisHist->GetMaximum ();
   }
   for (int iCent = 0; iCent < 2; iCent++) {
    TH1D* thisHist = xJgDist[iTrig][iCent];

    thisHist->GetYaxis ()->SetRangeUser (0, 1.1*max);
    thisHist->SetMarkerColor (colors[iCent]);
    thisHist->SetLineColor (colors[iCent]);
    thisHist->GetXaxis ()->SetTitle ("x_{J}^{#gamma}");
    thisHist->GetYaxis ()->SetTitle ("1 / N_{#gamma}  dN_{pairs}^{#gamma+jet} / dx_{J}^{#gamma}");
    if (iCent == 0)
     thisHist->Draw ("e1 x0");
    else
     thisHist->Draw ("same e1 x0");
   }
   
   myText (0.6, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
   myText (0.6, 0.805, colors[0], "#Sigma#it{E}_{T}^{FCal} > 2 TeV", 0.04);
   myText (0.6, 0.745, colors[1], "#Sigma#it{E}_{T}^{FCal} < 2 TeV", 0.04);
   canvas->SaveAs (Form ("%s/xJgDist_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));


   TH1D**** etconePlots[3] = {etcone20, etcone30, etcone40};
   TString etconeStrs[3] = {"Etcone20", "Etcone30", "Etcone40"};
   gPad->SetLogy (false);
   for (int iEta = 0; iEta < 2; iEta++) {
    max = 0;
    for (short iCent = 0; iCent < 2; iCent++) {
     for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
      if (etconePlots[iEtcone][iTrig][iEta][iCent]->GetMaximum () > max)
       max = etconePlots[iEtcone][iTrig][iEta][iCent]->GetMaximum ();
     }
    }

    for (short iCent = 0; iCent < 2; iCent++) {
     for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
      TH1D* thisHist = etconePlots[iEtcone][iTrig][iEta][iCent];
      thisHist->GetYaxis ()->SetRangeUser (0.5, 1.1*max);
      thisHist->SetLineStyle (iCent==0?2:1);
      thisHist->SetMarkerColor (colors[iEtcone]);
      thisHist->SetLineColor (colors[iEtcone]);

      thisHist->GetXaxis ()->SetTitle ("Isolation Etcone Value #left[GeV#right]");
      thisHist->GetYaxis ()->SetTitle ("Counts");

      if (iEtcone == 0 && iCent == 0) {
       thisHist->Draw ("hist");
       myText (0.6, 0.53, kBlack, "Dashed: #Sigma#it{E}_{T}^{FCal} > 2 TeV", 0.04);
       myText (0.6, 0.475, kBlack, "Solid: #Sigma#it{E}_{T}^{FCal} < 2 TeV", 0.04);
      }
      else
       thisHist->Draw ("same hist");

      if (iCent == 0)
       myText (0.6, 0.725-0.055*iEtcone, colors[iEtcone], etconeStrs[iEtcone], 0.04);
     }
    }
    myText (0.6, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
    myText (0.6, 0.81, kBlack, (string(iEta==0?"Barrel":"Endcap") + ", tight photons > 20 GeV").c_str(), 0.04);
    canvas->SaveAs (Form ("%s/etcones_%s_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"));
   }

   max = 0;
   for (short iEta = 0; iEta < 2; iEta++) {
    fcal_et[iTrig][iEta]->Rebin(2);
    if (max < fcal_et[iTrig][iEta]->GetMaximum ())
     max = fcal_et[iTrig][iEta]->GetMaximum ();
   }    
   
   for (short iEta = 0; iEta < 2; iEta++) {
    TH1D* thisHist = fcal_et[iTrig][iEta];
    thisHist->GetYaxis ()->SetRangeUser (0.5, 1.1*max);
    thisHist->SetMarkerColor (colors[iEta]);
    thisHist->SetLineColor (colors[iEta]);

    thisHist->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} #left[TeV#right]");
    thisHist->GetYaxis ()->SetTitle ("Counts");

    if (iEta == 0)
     thisHist->Draw ("hist");
    else
     thisHist->Draw ("same hist");

    myText (0.25, 0.88-0.07*iEta, colors[iEta], iEta==0?"Barrel":"Endcaps", 0.04);
   }    
   myText (0.6, 0.88, kBlack, photonTrigNames[iTrig].c_str (), 0.04);
   myText (0.6, 0.81, kBlack, "Tight photons, #it{p}_{T}^{#gamma} > 20 GeV", 0.04);
   canvas->SaveAs (Form ("%s/fcal_et_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ()));

   TCanvas* etaPhiMapCanvas = new TCanvas (Form ("etaPhiMapCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 600);
   etaPhiMapCanvas->cd ();
   FormatTH2Canvas (etaPhiMapCanvas);
   gStyle->SetPalette (kRainBow);
   etaPhiMap[iTrig]->GetXaxis ()->SetTitle ("#eta");
   etaPhiMap[iTrig]->GetYaxis ()->SetTitle ("#phi");
   etaPhiMap[iTrig]->Draw ("colz");
   etaPhiMapCanvas->SaveAs (Form ("%s/etaPhiMap_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str ())); 

   TCanvas* showerShapeCanvas = new TCanvas (Form ("showerShapeCanvas_%s", photonTrigNames[iTrig].c_str ()), "", 800, 600);
   showerShapeCanvas->Divide (4, 3);
   showerShapeCanvas->cd ();
   for (short iEta = 0; iEta < 2; iEta++) {
    for (short iShower = 0; iShower < 12; iShower++) {

     showerShapeCanvas->cd (iShower+1);
     TH1D* thisHist = showerShapeDists[iTrig][iEta][iShower];
     thisHist->Rebin (20);
     const float integral = thisHist->Integral ();
     if (integral > 0)
      thisHist->Scale (1. / integral, "width");

     max = std::max (thisHist->GetMaximum (), showerShapeDists[photonTrigN][iEta][iShower]->GetMaximum ());
     thisHist->GetYaxis ()->SetRangeUser (0, 1.1*max);

     thisHist->GetYaxis ()->SetTitle ("Counts / Total");

     thisHist->Draw ("hist");

     thisHist = showerShapeDists[photonTrigN][iEta][iShower];
     thisHist->SetLineColor (kBlue);
     thisHist->Draw ("same hist");

     myText (0.2, 0.9, kBlack, GetShowerShape (iShower), 0.12);
    }
    showerShapeCanvas->SaveAs (Form ("%s/showerShapes_%s_%s.pdf", plotPath.Data (), photonTrigNames[iTrig].c_str (), iEta==0?"barrel":"endcap"));
   }

   if (canvas)
    delete canvas;
   if (etaPhiMapCanvas)
    delete etaPhiMapCanvas;
   if (showerShapeCanvas)
    delete showerShapeCanvas;
  }
}

} // end namespace
