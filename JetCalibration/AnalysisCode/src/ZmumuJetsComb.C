#include "ZmumuJetsComb.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TH1D.h>
#include <TH3D.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <iostream>

namespace JetCalibration {

void ZmumuJetsComb () {

  SetupDirectories ("ZmumuJets/", "JetCalibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  TH3D***** zmumuJetHists = Get4DArray <TH3D*> (2, 3, 2, 3); // iYear, iPer, iData, iErr
  TH2D**** zmumuJetCounts = Get3DArray <TH2D*> (2, 3, 2);

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iPer = 0; iPer < 3; iPer++) {
    const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
     const TString dataType = (iData == 0 ? "data" : "mc_overlay");

     for (short iErr = 0; iErr < 3; iErr++) {
      if (iErr != 1 && iData != 0)
       continue;

      const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));
 
      const TString name = Form ("zmumuJetPtRatio_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data ());
      zmumuJetHists[iYear][iPer][iData][iErr] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
      zmumuJetHists[iYear][iPer][iData][iErr]->Sumw2 ();
     }

     const TString name = Form ("zmumuJetCounts_%s_%s_%s", year.Data (), dataType.Data (), period.Data ());
     zmumuJetCounts[iYear][iPer][iData] = new TH2D (name, "", numpbins, pbins, numetabins, etabins);
     zmumuJetCounts[iYear][iPer][iData]->Sumw2 ();
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
     if (debugStatements) cout << "Status: In ZmumuJetsComb.C: Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);

       for (short iErr = 0; iErr < 3; iErr++) {
        const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
        zmumuJetHists[iYear][iPer][0][iErr]->Add (temp3);
        zmumuJetHists[iYear][2][0][iErr]->Add (temp3);
       }

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("zmumuJetCounts_dataSet%i_data", runNumber));
       zmumuJetCounts[iYear][iPer][0]->Add (temp2);
       zmumuJetCounts[iYear][2][0]->Add (temp2);

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

       TH3D* temp = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_mc_overlay_stat", zmumuJetSampleId.Data ()));
       zmumuJetHists[iYear][iPer][1][1]->Add (temp);
       zmumuJetHists[iYear][2][1][1]->Add (temp);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("zmumuJetCounts_dataSet%s_mc_overlay", zmumuJetSampleId.Data ()));
       zmumuJetCounts[iYear][iPer][1]->Add (temp2);
       zmumuJetCounts[iYear][2][1]->Add (temp2);

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

  TFile* outFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "recreate");
  outFile->cd ();

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   for (short iPer = 0; iPer < 3; iPer++) {

    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC

     for (short iErr = 0; iErr < 3; iErr++) {
      if (iErr != 1 && iData != 0)
       continue;

      zmumuJetHists[iYear][iPer][iData][iErr]->Write ();
     }

     zmumuJetCounts[iYear][iPer][iData]->Write ();
    }
   }
  }

  Delete4DArray (zmumuJetHists, 2, 3, 2, 3);
  Delete3DArray (zmumuJetCounts, 2, 3, 2);

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }

  return;
}

} // end namespace
