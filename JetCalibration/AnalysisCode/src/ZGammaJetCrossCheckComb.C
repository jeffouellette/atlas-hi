#include "ZGammaJetCrossCheckComb.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <GlobalParams.h>
#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TGraphAsymmErrors.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace JetCalibration {

/**
 * Equivalent to adding h2 to h1, except with the weight defined by a TF1 at each x value.
 * Could be modified to take in a different weight at each x, y value or x, y, z value.
 * To ensure different histograms are reweighted correctly, they are normalized by default.
 */
void AddWx (TF1* fx, TH3D* h1, TH3D* h2, const bool normalize = true) {
  if (h2->Integral() == 0) {
   return; // true if there is nothing to add
  }
  for (int ix = 1; ix <= h1->GetNbinsX(); ix++) {
   const double wx = fx->Eval (h1->GetXaxis ()->GetBinCenter (ix)) / h2->Integral(); // this cannot = 0
   for (int iy = 1; iy <= h1->GetNbinsY(); iy++) {
    for (int iz = 1; iz <= h1->GetNbinsZ(); iz++) {
     h1->SetBinContent (ix, iy, iz, h1->GetBinContent (ix, iy, iz) + wx * h2->GetBinContent (ix, iy, iz));
     h1->SetBinError (ix, iy, iz, sqrt (pow (h1->GetBinError (ix, iy, iz), 2) + pow (wx * h2->GetBinError (ix, iy, iz), 2)));
    }
   }
  }
  return;
}


void ZGammaJetCrossCheckComb () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("", "JetCalibration/");

  vector<TString> subdirectories (0);
  subdirectories.push_back ("GammaJets");
  subdirectories.push_back ("ZeeJets");
  subdirectories.push_back ("ZmumuJets");
  
  TH3D****** vJetHists = Get5DArray <TH3D*> (2, 3, 3, 3, 4); // iYear, iPer, iData, iErr, iSpc
  TH2D***** vJetCounts = Get4DArray <TH2D*> (2, 3, 3, 4); // iYear, iPer, iData, iSpc

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue;
   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
    if (iData == 2 && skipSignalMC)
     continue; // ignore signal MC if not desired
    const TString dataType = (iData == 0 ? "data" : (iData == 1 ? "mc_overlay" : "mc_signalonly"));

    for (short iPer = 0; iPer < 3; iPer++) {
     const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

     for (short iSpc = 0; iSpc < 4; iSpc++) {
      const TString spc = (iSpc == 0 ? "g" : (iSpc == 1 ? "zee" : (iSpc == 2 ? "zmumu" : "v")));

      for (short iErr = 0; iErr < 3; iErr++) {
       if (iErr != 1 && iData != 0)
        continue; // ignore systematics for MC
       const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

       const TString name = Form ("%sJetPtRatio_%s_%s_%s_%s", spc.Data (), year.Data (), dataType.Data (), error.Data (), period.Data ());
       vJetHists[iYear][iPer][iData][iErr][iSpc] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
       vJetHists[iYear][iPer][iData][iErr][iSpc]->Sumw2();
      }

      const TString name = Form ("%sJetCounts_%s_%s_%s", spc.Data (), year.Data (), dataType.Data (), period.Data ());
      vJetCounts[iYear][iPer][iData][iSpc] = new TH2D (name, "", numpbins, pbins, numetabins, etabins);
      vJetCounts[iYear][iPer][iData][iSpc]->Sumw2();
     }
    }
   }
  }

  for (TString subdirectory : subdirectories) {
   for (short iYear = 0; iYear < 2; iYear++) {
    if (skipOldInsitu && iYear == 0)
     continue;
    const TString year = (iYear == 0 ? "2015" : "2016");

    const TString path = (iYear == 0 ? rootPath + "/2015factors/" : rootPath) + subdirectory + "/";

    TSystemDirectory dir (path.Data (), path.Data ());
    TList* sysfiles = dir.GetListOfFiles ();

    if (!sysfiles) {
     cout << "Cannot get list of files! Exiting." << endl;
     return;
    }
    TSystemFile *sysfile;
    TString fname;
    TIter next (sysfiles);
    while ( (sysfile= (TSystemFile*)next ())) {
     fname = sysfile->GetName ();
     if (fname != "outFile.root")
      continue;

     const short iSpc = (subdirectory == "GammaJets" ? 0 : (subdirectory == "ZeeJets" ? 1 : 2));
     const TString prefix = (subdirectory == "GammaJets" ? "gJet" : (subdirectory == "ZeeJets" ? "zeeJet" : "zmumuJet"));

     cout << "Reading in " << path+fname << endl;
     TFile* thisFile = new TFile (path + fname, "READ");

     for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
      if (iData == 2 && skipSignalMC)
       continue; // ignore signal MC if not desired
      const TString dataType = (iData == 0 ? "data" : (iData == 1 ? "mc_overlay" : "mc_signalonly"));

      for (short iPer = 0; iPer < 3; iPer++) {
       const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

       for (short iErr = 0; iErr < 3; iErr++) {
        if (iErr != 1 && iData != 0)
         continue; // ignore systematics for MC
        const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

        TString name = Form ("%sPtRatio_%s_%s_%s_%s", prefix.Data (), year.Data (), dataType.Data (), error.Data (), period.Data ());
        if (subdirectory == "GammaJets") {
         name = name + "_signal";
        }
        TH3D* temp3 = (TH3D*)thisFile->Get (name);
        vJetHists[iYear][iPer][iData][iErr][iSpc]->Add (temp3);
        vJetHists[iYear][iPer][iData][iErr][3]->Add (temp3);
       } // end loop over errors

       TString name = Form ("%sCounts_%s_%s_%s", prefix.Data (), year.Data (), dataType.Data (), period.Data ());
       if (subdirectory == "GammaJets") {
        name = name + "_tight_unweighted";
       }
       TH2D* temp2 = (TH2D*)thisFile->Get (name);
       vJetCounts[iYear][iPer][iData][iSpc]->Add (temp2);
       vJetCounts[iYear][iPer][iData][3]->Add (temp2);
      } // end loop over periods
     } // end loop over data types
    } // end loop over files
   } // end loop over insitu factors
  } // end loop over subdirectories

  //TF1* gWeight = new TF1 ("gWeight", "1 / (1 + exp (-[0] * (x - [1])))", 0, pbins[numpbins]);
  //gWeight->SetParNames ("g");
  //gWeight->SetParameter (0, 0.2); // guess- changes significantly over ~5 GeV?
  //gWeight->SetParameter (1, 70); // guess- centered at 70 GeV?

  //TF1* zWeight = new TF1 ("zWeight", "1-0.5*g", 0, pbins[numpbins]);

  //for (short iYear = 0; iYear < 2; iYear++) {
  // if (skipOldInsitu && iYear == 0) continue;

  // const TString path = (iYear == 0 ? rootPath + "/2015factors/" : rootPath);

  // TSystemDirectory dir (path.Data (), path.Data ());
  // TList* sysfiles = dir.GetListOfFiles ();
  // if (!sysfiles) {
  //  cout << "Cannot get list of files! Exiting." << endl;
  //  return;
  // }
  // TSystemFile *sysfile;
  // TString fname;
  // TString histName;
  // TIter next (sysfiles);
  // int numFiles = 0;
  // while ( (sysfile= (TSystemFile*)next ())) {
  //  fname = sysfile->GetName ();
  //  if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
  //   if (debugStatements) cout << "Status: In ZGammaJetCrossCheckHist.C: Found " << fname.Data () << endl;

  //   // do this if file is data
  //   for (int runNumber : runNumbers) { // check for data
  //    if (fname.Contains (to_string (runNumber))) { // if data, do this
  //     numFiles++;
  //     cout << "Reading in " << path+fname << endl;
  //     TFile* thisFile = new TFile (path + fname, "READ");
  //     const short iPer = (runNumber < 313500 ? 0 : 1);

  //     for (short iErr = 0; iErr < 3; iErr++) {
  //      TString error = "sys_lo";
  //      if (iErr == 1) error = "stat";
  //      else if (iErr == 2) error = "sys_hi";

  //      TH3D* temp3 = (TH3D*)thisFile->Get (Form ("vJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
  //      vJetHists[iYear][iPer][0][iErr]->Add (temp3);
  //      vJetHists[iYear][2][0][iErr]->Add (temp3);

  //      //TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
  //      //AddWx (gWeight, vJetHists[iYear][iPer][0][iErr], temp3);
  //      //AddWx (gWeight, vJetHists[iYear][2][0][iErr], temp3);

  //      //temp3 = (TH3D*)thisFile->Get (Form ("zeeJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
  //      //AddWx (zWeight, vJetHists[iYear][iPer][0][iErr], temp3);
  //      //AddWx (zWeight, vJetHists[iYear][2][0][iErr], temp3);

  //      //temp3 = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
  //      //AddWx (zWeight, vJetHists[iYear][iPer][0][iErr], temp3);
  //      //AddWx (zWeight, vJetHists[iYear][2][0][iErr], temp3);
  //     }

  //     TH2D* temp2 = (TH2D*)thisFile->Get (Form ("vJetCounts_dataSet%i_data", runNumber));
  //     vJetCounts[iYear][iPer][0]->Add (temp2);
  //     vJetCounts[iYear][2][0]->Add (temp2);

  //     //TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%i_data", runNumber));
  //     //vJetCounts[iYear][iPer][0]->Add (temp2);
  //     //vJetCounts[iYear][2][0]->Add (temp2);

  //     //temp2 = (TH2D*)thisFile->Get (Form ("zeeJetCounts_dataSet%i_data", runNumber));
  //     //vJetCounts[iYear][iPer][0]->Add (temp2);
  //     //vJetCounts[iYear][2][0]->Add (temp2);

  //     //temp2 = (TH2D*)thisFile->Get (Form ("zmumuJetCounts_dataSet%i_data", runNumber));
  //     //vJetCounts[iYear][iPer][0]->Add (temp2);
  //     //vJetCounts[iYear][2][0]->Add (temp2);

  //     thisFile->Close ();
  //     delete thisFile;
  //     break;
  //    }
  //   }
  //   // do this if gamma jet OVERLAY MC sample
  //   for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
  //    if (fname.Contains (gammaJetOverlaySampleId)) { // if gamma jet MC sample
  //     numFiles++;
  //     cout << "Reading in " << path+fname << endl;
  //     TFile* thisFile = new TFile (path + fname, "READ");
  //     const short iPer = (gammaJetOverlaySampleId.Contains ("pPb") ? 0 : 1);

  //     TH3D* temp3 = (TH3D*)thisFile->Get (Form ("vJetPtRatio_dataSet%s_mc_overlay_stat", gammaJetOverlaySampleId.Data ()));
  //     vJetHists[iYear][iPer][1][1]->Add (temp3);
  //     vJetHists[iYear][2][1][1]->Add (temp3);

  //     //TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_overlay_stat", gammaJetOverlaySampleId.Data ()));
  //     //AddWx (gWeight, vJetHists[iYear][iPer][1][1], temp3);
  //     //AddWx (gWeight, vJetHists[iYear][2][1][1], temp3);

  //     TH2D* temp2 = (TH2D*)thisFile->Get (Form ("vJetCounts_dataSet%s_mc_overlay", gammaJetOverlaySampleId.Data ()));
  //     vJetCounts[iYear][iPer][1]->Add (temp2);
  //     vJetCounts[iYear][2][1]->Add (temp2);

  //     //TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_overlay", gammaJetOverlaySampleId.Data ()));
  //     //vJetCounts[iYear][iPer][1]->Add (temp2);
  //     //vJetCounts[iYear][2][1]->Add (temp2);

  //     thisFile->Close ();
  //     delete thisFile;
  //     break;
  //    }
  //   }
  //   // do this if gamma jet SIGNAL MC sample
  //   for (TString gammaJetSignalSampleId : gammaJetSignalSampleIds) { // check for gamma jet MC
  //    if (fname.Contains (gammaJetSignalSampleId)) { // if gamma jet MC sample
  //     numFiles++;
  //     cout << "Reading in " << path+fname << endl;
  //     TFile* thisFile = new TFile (path + fname, "READ");
  //     const short iPer = (gammaJetSignalSampleId.Contains ("pPb") ? 0 : 1);

  //     TH3D* temp3 = (TH3D*)thisFile->Get (Form ("vJetPtRatio_dataSet%s_mc_signal_stat", gammaJetSignalSampleId.Data ()));
  //     vJetHists[iYear][iPer][2][1]->Add (temp3);
  //     vJetHists[iYear][2][2][1]->Add (temp3);

  //     //TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_signal_stat", gammaJetSignalSampleId.Data ()));
  //     //AddWx (gWeight, vJetHists[iYear][iPer][2][1], temp3);
  //     //AddWx (gWeight, vJetHists[iYear][2][2][1], temp3);

  //     TH2D* temp2 = (TH2D*)thisFile->Get (Form ("vJetCounts_dataSet%s_mc_signal", gammaJetSignalSampleId.Data ()));
  //     vJetCounts[iYear][iPer][2]->Add (temp2);
  //     vJetCounts[iYear][2][2]->Add (temp2);

  //     //TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_signal", gammaJetSignalSampleId.Data ()));
  //     //vJetCounts[iYear][iPer][2]->Add (temp2);
  //     //vJetCounts[iYear][2][2]->Add (temp2);

  //     thisFile->Close ();
  //     delete thisFile;
  //     break;
  //    }
  //   }
  //   // do this if Z->ee MC sample
  //   for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
  //    if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
  //     numFiles++;
  //     cout << "Reading in " << path+fname << endl;
  //     TFile* thisFile = new TFile (path + fname, "READ");
  //     const short iPer = (zeeJetSampleId.Contains ("pPb") ? 0 : 1);

  //     TH3D* temp3 = (TH3D*)thisFile->Get (Form ("vJetPtRatio_dataSet%s_mc_overlay_stat", zeeJetSampleId.Data ()));
  //     vJetHists[iYear][iPer][1][1]->Add (temp3);
  //     vJetHists[iYear][2][1][1]->Add (temp3);

  //     //TH3D* temp3 = (TH3D*)thisFile->Get (Form ("zeeJetPtRatio_dataSet%s_mc_overlay_stat", zeeJetSampleId.Data ()));
  //     //AddWx (zWeight, vJetHists[iYear][iPer][1][1], temp3);
  //     //AddWx (zWeight, vJetHists[iYear][2][1][1], temp3);

  //     TH2D* temp2 = (TH2D*)thisFile->Get (Form ("vJetCounts_dataSet%s_mc_overlay", zeeJetSampleId.Data ()));
  //     vJetCounts[iYear][iPer][1]->Add (temp2);
  //     vJetCounts[iYear][2][1]->Add (temp2);

  //     //TH2D* temp2 = (TH2D*)thisFile->Get (Form ("zeeJetCounts_dataSet%s_mc_overlay", zeeJetSampleId.Data ()));
  //     //vJetCounts[iYear][iPer][1]->Add (temp2);
  //     //vJetCounts[iYear][2][1]->Add (temp2);

  //     thisFile->Close ();
  //     delete thisFile;
  //     break;
  //    }
  //   }
  //   // do this if Z->mumu sample
  //   for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
  //    if (fname.Contains (zmumuJetSampleId)) { // if Z->mumu sample do this
  //     numFiles++;
  //     cout << "Reading in " << rootPath+fname << endl;
  //     TFile* thisFile = new TFile (rootPath + fname, "READ");
  //     const short iPer = (zmumuJetSampleId.Contains ("pPb") ? 0 : 1);

  //     TH3D* temp3 = (TH3D*)thisFile->Get (Form ("vJetPtRatio_dataSet%s_mc_overlay_stat", zmumuJetSampleId.Data ()));
  //     vJetHists[iYear][iPer][1][1]->Add (temp3);
  //     vJetHists[iYear][2][1][1]->Add (temp3);

  //     //TH3D* temp3 = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_mc_overlay_stat", zmumuJetSampleId.Data ()));
  //     //AddWx (zWeight, vJetHists[iYear][iPer][1][1], temp3);
  //     //AddWx (zWeight, vJetHists[iYear][2][1][1], temp3);

  //     TH2D* temp2 = (TH2D*)thisFile->Get (Form ("vJetCounts_dataSet%s_mc_overlay", zmumuJetSampleId.Data ()));
  //     vJetCounts[iYear][iPer][1]->Add (temp2);
  //     vJetCounts[iYear][2][1]->Add (temp2);

  //     //TH2D* temp2 = (TH2D*)thisFile->Get (Form ("zmumuJetCounts_dataSet%s_mc_overlay", zmumuJetSampleId.Data ()));
  //     //vJetCounts[iYear][iPer][1]->Add (temp2);
  //     //vJetCounts[iYear][2][1]->Add (temp2);

  //     thisFile->Close ();
  //     delete thisFile;
  //     break;
  //    }
  //   }
  //  }
  // }
  // cout << numFiles << " files read in." << endl;
  //}

  
  /**** Save output ****/
  TFile* outFile = new TFile (Form ("%s/ZGammaJetCrossCheck/outFile.root", rootPath.Data ()), "recreate");
  outFile->cd ();

  /**** End loop over input files ****/

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   for (short iPer = 0; iPer < 3; iPer++) {

    for (short iData = 0; iData < 3; iData++) {
     if (iData == 2 && skipSignalMC)
      continue; // ignore signal MC if not desired

     for (short iSpc = 0; iSpc < 4; iSpc++) {

      for (short iErr = 0; iErr < 3; iErr++) {
       if (iErr != 1 && iData != 0)
        continue; // ignore systematics for MC

       vJetHists[iYear][iPer][iData][iErr][iSpc]->Write ();
      }

      vJetCounts[iYear][iPer][iData][iSpc]->Write ();
     }
    }
   }
  }

  Delete5DArray (vJetHists, 2, 3, 3, 3, 4);
  Delete4DArray (vJetCounts, 2, 3, 3, 4);

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }

  return;
}

} // end namespace
