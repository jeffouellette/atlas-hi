#include "FCalDistributionHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>

#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasUtils.h>
#include <AtlasStyle.h>

namespace JetCalibration {

void FCalDistributionHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("FCalDistribution/", "JetCalibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back (TString ("Pbp") + (runValidation ? "_Signal":"_Overlay") + "_GammaJet_Slice" + to_string (i+1));
   gammaJetSampleIds.push_back (TString ("pPb") + (runValidation ? "_Signal":"_Overlay") + "_GammaJet_Slice" + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds (0);
  dijetSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  TH1D* fCal_p_et[2];
  fCal_p_et[0] = new TH1D ("fCal_p_et_data", "", 125, -50, 200);
  fCal_p_et[0]->Sumw2 ();
  fCal_p_et[1] = new TH1D ("fCal_p_et_mc", "", 125, -50, 200);
  fCal_p_et[1]->Sumw2 ();

  TH1D* fCal_Pb_et[2];
  fCal_Pb_et[0] = new TH1D ("fCal_Pb_et_data", "", 125, -50, 200);
  fCal_Pb_et[0]->Sumw2 ();
  fCal_Pb_et[1] = new TH1D ("fCal_Pb_et_mc", "", 125, -50, 200);
  fCal_Pb_et[1]->Sumw2 ();

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

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       fCal_p_et[0]->Add ( (TH1D*)thisFile->Get (Form ("fCal_p_et_dataSet%i", runNumber)));
       fCal_Pb_et[0]->Add ( (TH1D*)thisFile->Get (Form ("fCal_Pb_et_dataSet%i", runNumber)));
    
       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet MC sample
     for (TString gammaJetSampleId : gammaJetSampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       fCal_p_et[1]->Add ( (TH1D*)thisFile->Get (Form ("fCal_p_et_dataSet%s", gammaJetSampleId.Data ())));
       fCal_Pb_et[1]->Add ( (TH1D*)thisFile->Get (Form ("fCal_Pb_et_dataSet%s", gammaJetSampleId.Data ())));

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       fCal_p_et[1]->Add ( (TH1D*)thisFile->Get (Form ("fCal_p_et_dataSet%s", zeeJetSampleId.Data ())));
       fCal_Pb_et[1]->Add ( (TH1D*)thisFile->Get (Form ("fCal_Pb_et_dataSet%s", zeeJetSampleId.Data ())));

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

       fCal_p_et[1]->Add ( (TH1D*)thisFile->Get (Form ("fCal_p_et_dataSet%s", zmumuJetSampleId.Data ())));
       fCal_Pb_et[1]->Add ( (TH1D*)thisFile->Get (Form ("fCal_Pb_et_dataSet%s", zmumuJetSampleId.Data ())));

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



  /**** Canvas definitions ****/
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  canvas->Draw ();

  Style_t markers[2] = {kFullCircle, kOpenCircle};
  for (int iData = 0; iData < 2; iData++) {
   fCal_p_et[iData]->SetLineColor (kBlue);
   fCal_p_et[iData]->SetMarkerColor (kBlue);
   fCal_Pb_et[iData]->SetLineColor (kBlack);
   fCal_Pb_et[iData]->SetMarkerColor (kBlack);

   fCal_p_et[iData]->SetMarkerStyle (markers[iData]);
   fCal_Pb_et[iData]->SetMarkerStyle (markers[iData]);

   fCal_p_et[iData]->GetXaxis ()->SetTitle ("Forward #SigmaE_{T} #left[GeV#right]");
   fCal_Pb_et[iData]->GetXaxis ()->SetTitle ("Forward #SigmaE_{T} #left[GeV#right]");
   fCal_p_et[iData]->GetYaxis ()->SetTitle ("Counts / Total");
   fCal_Pb_et[iData]->GetYaxis ()->SetTitle ("Counts / Total");

   fCal_p_et[iData]->Scale (1./fCal_p_et[iData]->Integral ());
   fCal_Pb_et[iData]->Scale (1./fCal_Pb_et[iData]->Integral ());

   if (iData == 0)
    fCal_p_et[iData]->Draw ("e1");
   else
    fCal_p_et[iData]->Draw ("same e1");
   fCal_Pb_et[iData]->Draw ("same e1");
  }

  myMarkerText (0.55, 0.85, kBlack, kFullCircle, "Pb-going direction", 1.25, 0.05);
  myMarkerText (0.55, 0.78, kBlue, kFullCircle, "p-going direction", 1.25, 0.05);

  canvas->SaveAs (Form ("%s/fcal_et.pdf", plotPath.Data ()));

  return;
}

} // end namespace
