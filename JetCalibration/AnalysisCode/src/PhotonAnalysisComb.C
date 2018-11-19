#include "PhotonAnalysisComb.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TVirtualFitter.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <iostream>

using namespace atlashi;

namespace JetCalibration {


void PhotonAnalysisComb () {

  // Setup trigger vectors
  SetupDirectories ("PhotonAnalysis/", "JetCalibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetOverlaySampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetOverlaySampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetOverlaySampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> gammaJetSignalSampleIds (0);
  if (!skipSignalMC) {
   for (short i = 0; i < 6; i++) {
    gammaJetSignalSampleIds.push_back (TString ("pPb_Signal_GammaJet_Slice") + to_string (i+1));
   }
  }

  // initialize histograms
  TH2D***** photonSpectrum = Get4DArray <TH2D*> (3, 3, 2, 2); // iPer, iData, iErr, tight/lnt, weighted/non-weighted
  TH2D***** photonEtaPhi = Get4DArray <TH2D*> (3, 3, 2, 2); // iPer, iData, iErr, tight/lnt, weighted/non-weighted

  for (short iPer = 0; iPer < 3; iPer++) {
   const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

   for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
    if (iData == 2 && skipSignalMC)
     continue; // ignore signal MC if not desired
    const TString data = (iData == 0 ? "data" : (iData == 1 ? "mc_overlay" : "mc_signal"));

    for (short iPCut = 0; iPCut < 2; iPCut++) {
     const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

     for (short iWgt = 0; iWgt < 2; iWgt++) {
      const TString weight = (iWgt == 0 ? "unweighted":"weighted");

      photonSpectrum[iPer][iData][iPCut][iWgt] = new TH2D (Form ("photonSpectrum_%s_%s_%s_%s", period.Data (), data.Data (), pCut.Data (), weight.Data ()), "", numpbins, pbins, numetabins, etabins);
      photonSpectrum[iPer][iData][iPCut][iWgt]->Sumw2 ();

      photonEtaPhi[iPer][iData][iPCut][iWgt] = new TH2D (Form ("photonEtaPhi_%s_%s_%s_%s", period.Data (), data.Data (), pCut.Data (), weight.Data ()), "", 48, -2.4, 2.4, numphibins, phibins);
      photonEtaPhi[iPer][iData][iPCut][iWgt]->Sumw2 ();
     }
    }
    
   }
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
     if (debugStatements) cout << "Status: In PhotonAnalysisComb.C: Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);

       for (short iPCut = 0; iPCut < 2; iPCut++) {
        const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

        for (short iWgt = 0; iWgt < 2; iWgt++) {
         const TString weight = (iWgt == 0 ? "unweighted":"weighted");

         TH2D* thisHist = (TH2D*)thisFile->Get (Form ("photonSpectrum_dataSet%i_data_%s_%s", runNumber, pCut.Data (), weight.Data ()));
         photonSpectrum[iPer][iData][iPCut][iWgt]->Add (thisHist);
         photonSpectrum[2][iData][iPCut][iWgt]->Add (thisHist);

         thisHist = (TH2D*)thisFile->Get (Form ("photonEtaPhi_dataSet%i_data_%s_%s", runNumber, pCut.Data (), weight.Data ()));
         photonEtaPhi[iPer][iData][iPCut][iWgt]->Add (thisHist);
         photonEtaPhi[2][iData][iPCut][iWgt]->Add (thisHist);
        }
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet OVERLAY MC sample
     for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetOverlaySampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (gammaJetOverlaySampleId.Contains ("pPb") ? 0 : 1);

       for (short iPCut = 0; iPCut < 2; iPCut++) {
        const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

        for (short iWgt = 0; iWgt < 2; iWgt++) {
         const TString weight = (iWgt == 0 ? "unweighted":"weighted");

         TH2D* thisHist = (TH2D*)thisFile->Get (Form ("photonSpectrum_%s_mc_overlay_%s_%s", gammaJetOverlaySampleId.Data (), pCut.Data (), weight.Data ()));
         photonSpectrum[iPer][iData][iPCut][iWgt]->Add (thisHist);
         photonSpectrum[2][iData][iPCut][iWgt]->Add (thisHist);

         thisHist = (TH2D*)thisFile->Get (Form ("photonEtaPhi_%s_mc_overlay_%s_%s", gammaJetOverlaySampleId.Data (), pCut.Data (), weight.Data ()));
         photonEtaPhi[iPer][iData][iPCut][iWgt]->Add (thisHist);
         photonEtaPhi[2][iData][iPCut][iWgt]->Add (thisHist);
        }
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet SIGNAL MC sample
     for (TString gammaJetSignalSampleId : gammaJetSignalSampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetSignalSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (gammaJetSignalSampleId.Contains ("pPb") ? 0 : 1);

       for (short iPCut = 0; iPCut < 2; iPCut++) {
        const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

        for (short iWgt = 0; iWgt < 2; iWgt++) {
         const TString weight = (iWgt == 0 ? "unweighted":"weighted");

         TH2D* thisHist = (TH2D*)thisFile->Get (Form ("photonSpectrum_%s_mc_signal_%s_%s", gammaJetSignalSampleId.Data (), pCut.Data (), weight.Data ()));
         photonSpectrum[iPer][iData][iPCut][iWgt]->Add (thisHist);
         photonSpectrum[2][iData][iPCut][iWgt]->Add (thisHist);

         thisHist = (TH2D*)thisFile->Get (Form ("photonEtaPhi_%s_mc_signal_%s_%s", gammaJetSignalSampleId.Data (), pCut.Data (), weight.Data ()));
         photonEtaPhi[iPer][iData][iPCut][iWgt]->Add (thisHist);
         photonEtaPhi[2][iData][iPCut][iWgt]->Add (thisHist);
        }
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


  /**** Save output ****/
  TFile* outFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "recreate");
  outFile->cd ();

  for (short iPer = 0; iPer < 3; iPer++) {

   for (short iData = 0; iData < 3; iData++) {
    if (iData == 2 && skipSignalMC)
     continue; // ignore signal MC if not desired

    for (short iPCut = 0; iPCut < 2; iPCut++) {
     for (short iWgt = 0; iWgt < 2; iWgt++) {
      photonSpectrum[iPer][iData][iPCut][iWgt]->Write ();
      photonEtaPhi[iPer][iData][iPCut][iWgt]->Write ();
     }
    }
   }
  }

  Delete5DArray (photonSpectrum, 3, 3, 2, 2);
  Delete5DArray (photonEtaPhi, 3, 3, 2, 2);

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }
  /**** End save output ****/

  return;
}

} // end namespace
