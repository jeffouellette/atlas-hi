#include "PhotonAnalysisHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TVirtualFitter.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace atlashi;

namespace JetCalibration {


void PhotonAnalysisHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("PhotonAnalysis/", "JetCalibration/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TH2D***** photonSpectrum = Get4DArray <TH2D*> (3, 3, 2, 2); // iPer, iData, iPCut, iWgt
  TH2D***** photonEtaPhi = Get4DArray <TH2D*> (3, 3, 2, 2); // iPer, iData, iPCut, iWgt

  for (short iPer = 0; iPer < 3; iPer++) {
   const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

   for (short iData = 0; iData < 3; iData++) {
    if (iData == 2 && skipSignalMC)
     continue; // ignore signal MC if not desired

    for (short iPCut = 0; iPCut < 2; iPCut++) {
     const TString pCut = (iPCut == 0 ? "tight" : (iPCut == 1 ? "lnt" : "signal")); // lnt = loose non-tight

     for (short iWgt = 0; iWgt < 2; iWgt++) {
      const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

      photonSpectrum[iPer][iData][iPCut][iWgt] = (TH2D*)inFile->Get (Form ("photonSpectrum_%s_%s_%s_%s", period.Data (), data.Data (), pCut.Data (), weight.Data ()));
      photonEtaPhi[iPer][iData][iPCut][iWgt] = (TH2D*)inFile->Get (Form ("photonEtaPhi_%s_%s_%s_%s", period.Data (), data.Data (), pCut.Data (), weight.Data ()));
     }
    }
   }
  }


  /**** Canvas definitions ****/
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);

  for (short iPer = 0; iPer < 3; iPer++) {
   for (short iPCut = 0; iPCut < 3; iPCut++) {
    
   }
  }

  return;
}

} // end namespace
