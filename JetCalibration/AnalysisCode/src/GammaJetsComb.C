#include "GammaJetsComb.h"
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

/**
 * Returns a new TH3D equal to tight - lnt.
 */
void SubtractLooseNonTight (TH3D* signal, const TH3D* tight, const TH3D* lnt, const TH2D* tightCounts, const TH2D* lntCounts) {
  const int nx = signal->GetNbinsX ();
  const int ny = signal->GetNbinsY ();
  const int nz = signal->GetNbinsZ ();
  for (int ix = 1; ix <= nx; ix++) {
   for (int iy = 1; iy <= ny; iy++) {
    const double ntight = tightCounts->GetBinContent (ix, iy);
    const double ntighterr = tightCounts->GetBinError (ix, iy);
    const double nlnt = lntCounts->GetBinContent (ix, iy); 
    const double nlnterr = lntCounts->GetBinError (ix, iy);

    double w = 0, werr = 0;
    if (nlnt != 0) {
     w = ntight / nlnt;
     werr += pow (w * nlnterr / nlnt, 2);
    } 
    else
     w = 0;

    if (ntight != 0)
     werr += pow (w * ntighterr / ntight, 2);

    werr = sqrt (werr);

    for (int iz = 1; iz <= nz; iz++) {
     const double t = tight->GetBinContent (ix, iy, iz);
     const double terr = tight->GetBinError (ix, iy, iz);
     const double l = lnt->GetBinContent (ix, iy, iz);
     const double lerr = lnt->GetBinError (ix, iy, iz);

     signal->SetBinContent (ix, iy, iz, t - w * l);
     signal->SetBinError (ix, iy, iz, sqrt (pow (terr, 2) + pow (werr * l, 2) + pow (w * lerr, 2)));
    }
   }
  }
  return;
}


void GammaJetsComb () {

  // Setup trigger vectors
  SetupDirectories ("GammaJets/", "JetCalibration/");

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

  TH3D****** gJetHists = Get5DArray <TH3D*> (2, 3, 3, 3, 3); // iYear, iPer, iData, iErr, iPCut
  TH2D****** gJetCounts = Get5DArray <TH2D*> (2, 3, 3, 3, 2); // iYear, iPer, iData, iPCut, iWgt

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iPer = 0; iPer < 3; iPer++) {
    const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

    for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
     if (iData == 2 && skipSignalMC)
      continue; // ignore signal MC if not desired

     const TString dataType = (iData == 0 ? "data" : (iData == 1 ? "mc_overlay" : "mc_signal"));

     for (short iPCut = 0; iPCut < 3; iPCut++) {
      const TString pCut = (iPCut == 0 ? "tight" : (iPCut == 1 ? "lnt" : "signal")); // lnt = loose non-tight

      for (short iErr = 0; iErr < 3; iErr++) {
       if (iErr != 1 && iData != 0)
        continue; // ignore systematics for MC

       const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

       const TString name = Form ("gJetPtRatio_%s_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data (), pCut.Data ());
       gJetHists[iYear][iPer][iData][iErr][iPCut] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
       gJetHists[iYear][iPer][iData][iErr][iPCut]->Sumw2 ();
      }

      for (short iWgt = 0; iWgt < 2; iWgt++) {
       const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

       const TString name = Form ("gJetCounts_%s_%s_%s_%s_%s", year.Data (), dataType.Data (), period.Data (), pCut.Data (), weight.Data ());
       gJetCounts[iYear][iPer][iData][iPCut][iWgt] = new TH2D (name, "", numpbins, pbins, numetabins, etabins);
       gJetCounts[iYear][iPer][iData][iPCut][iWgt]->Sumw2 ();
      }
     }
    }
   }
  }

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0) continue;

   const TString path = (iYear == 0 ? rootPath + "/2015factors/" : rootPath);

   TSystemDirectory dir (path.Data (), path.Data ());
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
     if (debugStatements) cout << "Status: In GammaJetsComb.C: Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);

       for (short iPCut = 0; iPCut < 2; iPCut++) {
        const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

        for (short iErr = 0; iErr < 3; iErr++) {
         const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

         TH3D* temp = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%i_data_%s_%s", runNumber, error.Data (), pCut.Data ()));
         gJetHists[iYear][iPer][0][iErr][iPCut]->Add (temp);
         gJetHists[iYear][2][0][iErr][iPCut]->Add (temp);
        }

        for (short iWgt = 0; iWgt < 2; iWgt++) {
         const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

         TH2D* temp = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%i_data_%s_%s", runNumber, pCut.Data (), weight.Data ()));
         gJetCounts[iYear][iPer][0][iPCut][iWgt]->Add (temp);
         gJetCounts[iYear][2][0][iPCut][iWgt]->Add (temp);
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
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (gammaJetOverlaySampleId.Contains ("pPb") ? 0 : 1);

       for (short iPCut = 0; iPCut < 2; iPCut++) {
        const TString pCut = (iPCut == 0 ? "tight" : "lnt"); // lnt = loose non-tight

        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_overlay_stat_%s", gammaJetOverlaySampleId.Data (), pCut.Data ()));
        gJetHists[iYear][iPer][1][1][iPCut]->Add (temp3);
        gJetHists[iYear][2][1][1][iPCut]->Add (temp3);

        for (short iWgt = 0; iWgt < 2; iWgt++) {
         const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

         TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_overlay_%s_%s", gammaJetOverlaySampleId.Data (), pCut.Data (), weight.Data ()));
         gJetCounts[iYear][iPer][1][iPCut][iWgt]->Add (temp2);
         gJetCounts[iYear][2][1][iPCut][iWgt]->Add (temp2);
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
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (gammaJetSignalSampleId.Contains ("pPb") ? 0 : 1);

       for (short iPCut = 0; iPCut < 3; iPCut++) {
        const TString pCut = (iPCut == 0 ? "tight" : (iPCut == 1 ? "lnt" : "signal")); // lnt = loose non-tight

        TH3D* temp = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_signal_stat_%s", gammaJetSignalSampleId.Data (), pCut.Data ()));
        gJetHists[iYear][iPer][2][1][iPCut]->Add (temp);
        gJetHists[iYear][2][2][1][iPCut]->Add (temp);

        for (short iWgt = 0; iWgt < 2; iWgt++) {
         const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

         TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_signal_%s_%s", gammaJetSignalSampleId.Data (), pCut.Data (), weight.Data ()));
         gJetCounts[iYear][iPer][2][iPCut][iWgt]->Add (temp2);
         gJetCounts[iYear][2][2][iPCut][iWgt]->Add (temp2);
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


  /**** Subtract hadronic background from tight distributions ****/
  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0) continue;

   for (short iPer = 0; iPer < 3; iPer++) {
    for (short iErr = 0; iErr < 3; iErr++) { // subtract for data
     if (subtractBackground) {
      SubtractLooseNonTight (gJetHists[iYear][iPer][0][iErr][2], // signal
                             gJetHists[iYear][iPer][0][iErr][0], // tight
                             gJetHists[iYear][iPer][0][iErr][1], // loose non tight
                             gJetCounts[iYear][iPer][0][0][1], // tight weighted counts
                             gJetCounts[iYear][iPer][0][1][1]); // loose non tight weighted counts
     }
     else {
      delete gJetHists[iYear][iPer][0][iErr][2];
      gJetHists[iYear][iPer][0][iErr][2] = gJetHists[iYear][iPer][0][iErr][0];
     }
    }
    // subtract for overlay MC
    if (subtractBackground) {
     SubtractLooseNonTight (gJetHists[iYear][iPer][1][1][2],
                            gJetHists[iYear][iPer][1][1][0],
                            gJetHists[iYear][iPer][1][1][1],
                            gJetCounts[iYear][iPer][1][0][1],
                            gJetCounts[iYear][iPer][1][1][1]);
    }
    else {
     delete gJetHists[iYear][iPer][1][1][2];
     gJetHists[iYear][iPer][1][1][2] = gJetHists[iYear][iPer][1][1][0];
    }

    if (!skipSignalMC) { // subtract for signal MC
     if (subtractBackground) {
      SubtractLooseNonTight (gJetHists[iYear][iPer][2][1][2],
                             gJetHists[iYear][iPer][2][1][0],
                             gJetHists[iYear][iPer][2][1][1],
                             gJetCounts[iYear][iPer][2][0][1],
                             gJetCounts[iYear][iPer][2][1][1]);
     }
     else {
      delete gJetHists[iYear][iPer][2][1][2];
      gJetHists[iYear][iPer][2][1][2] = gJetHists[iYear][iPer][2][1][0];
     }
    }
   }
  }
  /**** End subtract hadronic background ****/


  /**** Save output ****/
  TFile* outFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "recreate");
  outFile->cd ();

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   for (short iPer = 0; iPer < 3; iPer++) {

    for (short iData = 0; iData < 3; iData++) {
     if (iData == 2 && skipSignalMC)
      continue; // ignore signal MC if not desired

     for (short iPCut = 0; iPCut < 3; iPCut++) {

      for (short iErr = 0; iErr < 3; iErr++) {
       if (iErr != 1 && iData != 0)
        continue; // ignore systematics for MC

       gJetHists[iYear][iPer][iData][iErr][iPCut]->Write ();
      }

      for (short iWgt = 0; iWgt < 2; iWgt++) {
       gJetCounts[iYear][iPer][iData][iPCut][iWgt]->Write ();
      }
     }
    }
   }
  }

  Delete5DArray (gJetHists, 2, 3, 3, 3, 3);
  Delete5DArray (gJetCounts, 2, 3, 3, 3, 2);

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }
  /**** End save output ****/

  return;
}

} // end namespace
