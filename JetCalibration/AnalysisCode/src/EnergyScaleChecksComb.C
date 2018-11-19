#include "EnergyScaleChecksComb.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TVectorT.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <iostream>

namespace JetCalibration {

const short plotMC = 0; // 0 plots overlay


TString GetMCType (const short iMC) {
  if (iMC == 0) return "overlay";
  else if (iMC == 1) return "signal";
  else if (iMC == 2) return "valid_signal";
  else return "";
}


void EnergyScaleChecksComb () {

  // Setup trigger vectors
  SetupDirectories ("EnergyScaleChecks/", "JetCalibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetOverlaySampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetOverlaySampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetOverlaySampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> gammaJetSignalSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSignalSampleIds.push_back (TString ("pPb_Signal_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds (0);
  dijetSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  vector<TString> dijetValidSampleIds (0);
  dijetValidSampleIds.push_back ("pPb_Signal_Valid_Dijet_Slice2");
  dijetValidSampleIds.push_back ("Pbp_Signal_Valid_Dijet_Slice2");


  /**** Initialize histograms ****/
  TH1D* electronEnergyScale = new TH1D ("electronEnergyScale", ";Electron #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 200, 0, 2.0);
  electronEnergyScale->Sumw2 ();

  TH1D****** jetEnergyResponseCalib = Get5DArray <TH1D*> (2, 3, numpbins+1, numetabins+1, numphibins); // DP slice, jet algo, MC type, pt bin, eta bin, phi bin
  TH1D****** jetEnergyResponseReco = Get5DArray <TH1D*> (2, 3, numpbins+1, numetabins+1, numphibins);
  TH1D**** photonEnergyResponse = Get3DArray <TH1D*> (3, numpbins+1, numetabins+1);

  //for (short iDP = 0; iDP < numdpbins; iDP++) {
  for (short iP = 0; iP <= numpbins; iP++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iMC = 0; iMC < 3; iMC++) {
     const TString mcType = GetMCType (iMC);
     for (short iPhi = 0; iPhi < numphibins; iPhi++) {
      for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
       const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
       jetEnergyResponseCalib[iAlgo][iMC][iP][iEta][iPhi] = new TH1D (Form ("%s_%s_jetEnergyResponseCalib_iP%i_iEta%i_iPhi%i", algo.Data (), mcType.Data (), iP, iEta, iPhi), "", 100, 0, 2);
       jetEnergyResponseCalib[iAlgo][iMC][iP][iEta][iPhi]->Sumw2 ();
       jetEnergyResponseReco[iAlgo][iMC][iP][iEta][iPhi] = new TH1D (Form ("%s_%s_jetEnergyResponseReco_iP%i_iEta%i_iPhi%i", algo.Data (), mcType.Data (), iP, iEta, iPhi), "", 100, 0, 2);
       jetEnergyResponseReco[iAlgo][iMC][iP][iEta][iPhi]->Sumw2 ();
      }
     }
     photonEnergyResponse[iMC][iP][iEta] = new TH1D (Form ("%s_photonEnergyResponse_iP%i_iEta%i", mcType.Data (), iP, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 100, 0, 2);
     photonEnergyResponse[iMC][iP][iEta]->Sumw2 ();
    }
   }
  }

  int***** nJet = Get5DArray <int> (2, 3, numpbins+1, numetabins+1, numphibins);
  int*** nGamma = Get3DArray <int> (3, numpbins+1, numetabins+1);

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
   TVectorD *nJetVec, *nGammaVec;
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data () << endl;

     // do this if gamma jet MC sample (OVERLAY)
     if (plotMC == 0) {
      for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
       if (fname.Contains (gammaJetOverlaySampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetOverlaySampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", gammaJetOverlaySampleId.Data ()));

        //int iDP = -1;
        //while (iDP < numdpbins && !gammaJetOverlaySampleId.Contains (Form ("Slice%i", iDP+1))) iDP++;
        //if (iDP == -1 || 6 <= iDP) continue;

        //const int iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iPhi = 0; iPhi < numphibins; iPhi++) {
           for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
             const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
             nJet[iAlgo][0][iP][iEta][iPhi] += (*nJetVec)[iPhi + (iP + (iEta + iAlgo*numetabins+1)*(numpbins+1))*numphibins];
             jetEnergyResponseCalib[iAlgo][0][iP][iEta][iPhi] = (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi));
             jetEnergyResponseReco[iAlgo][0][iP][iEta][iPhi] = (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi));
           }
           nGamma[0][iP][iEta] += (*nGammaVec)[iP + (iEta*(numpbins+1))];
           photonEnergyResponse[0][iP][iEta] = (TH1D*)thisFile->Get (Form ("photonEnergyResponse_iP%i_iEta%i", iP, iEta));
          }
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }
      // do this if gamma jet MC sample (SIGNAL)
      for (TString gammaJetSignalSampleId : gammaJetSignalSampleIds) { // check for gamma jet MC
       if (fname.Contains (gammaJetSignalSampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetSignalSampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", gammaJetSignalSampleId.Data ()));

        //int iDP = -1;
        //while (iDP < numdpbins && !gammaJetSignalSampleId.Contains (Form ("Slice%i", iDP+1))) iDP++;
        //if (iDP == -1 || 6 <= iDP) continue;

        //const int iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iPhi = 0; iPhi < numphibins; iPhi++) {
           for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iAlgo][1][iP][iEta][iPhi] += (*nJetVec)[iPhi + (iP + (iEta + iAlgo*numetabins+1)*(numpbins+1))*numphibins];
            jetEnergyResponseCalib[iAlgo][1][iP][iEta][iPhi] = (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi));
            jetEnergyResponseReco[iAlgo][1][iP][iEta][iPhi] = (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi));
           }
          }
          nGamma[1][iP][iEta] += (*nGammaVec)[iP + iEta*(numpbins+1)];
          photonEnergyResponse[1][iP][iEta] = (TH1D*)thisFile->Get (Form ("photonEnergyResponse_iP%i_iEta%i", iP, iEta));
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zeeJetSampleId.Data ()));

       electronEnergyScale->Add ( (TH1D*)thisFile->Get (Form ("electronEnergyScale", zeeJetSampleId.Data ())));

       //for (short iEta = 0; iEta <= numetabins; iEta++) {
       // for (short iP = 0; iP <= numpbins; iP++) {
       //  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
       //   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
       //   nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
       //   jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i", algo.Data (), iP, iEta)));
       //   jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i", algo.Data (), iP, iEta)));
       //  }
       // }
       //}

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     //// do this if Z->mumu sample
     //for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
     // if (fname.Contains (zmumuJetSampleId)) { // if Z->mumu sample do this
     //  numFiles++;
     //  cout << "Reading in " << rootPath+fname << endl;
     //  TFile* thisFile = new TFile (rootPath + fname, "READ");
     //  nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zmumuJetSampleId.Data ()));

     //  //for (short iEta = 0; iEta <= numetabins; iEta++) {
     //  // for (short iP = 0; iP <= numpbins; iP++) {
     //  //  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     //  //   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     //  //   nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
     //  //   jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i", algo.Data (), iP, iEta)));
     //  //   jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i", algo.Data (), iP, iEta)));
     //  //  }
     //  // }
     //  //}

     //  thisFile->Close ();
     //  delete thisFile;
     //  break;
     // }
     //}

     // do this if dijet MC sample (SIGNAL)
     if (plotMC == 1) {
      for (TString dijetSampleId : dijetSampleIds) { // check for gamma jet MC
       if (fname.Contains (dijetSampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", dijetSampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", dijetSampleId.Data ()));

        //const short iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iPhi = 0; iPhi < numphibins; iPhi++) {
           for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iAlgo][1][iP][iEta][iPhi] += (*nJetVec)[iPhi + (iP + (iEta + iAlgo*numetabins+1)*(numpbins+1))*numphibins];
            jetEnergyResponseCalib[iAlgo][1][iP][iEta][iPhi]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi)));
            jetEnergyResponseReco[iAlgo][1][iP][iEta][iPhi]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi)));
           }
          }
          //nGamma[1][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
          //photonEnergyResponse[1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("photonEnergyResponse_iP%i_iEta%i", iP, iEta)));
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }

      // do this if dijet MC sample (SIGNAL)
      for (TString dijetValidSampleId : dijetValidSampleIds) { // check for gamma jet MC
       if (fname.Contains (dijetValidSampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", dijetValidSampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", dijetValidSampleId.Data ()));

        //const short iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iPhi = 0; iPhi < numphibins; iPhi++) {
           for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iAlgo][2][iP][iEta][iPhi] += (*nJetVec)[iPhi + (iP + (iEta + iAlgo*numetabins+1)*(numpbins+1))*numphibins];
            jetEnergyResponseCalib[iAlgo][2][iP][iEta][iPhi]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi)));
            jetEnergyResponseReco[iAlgo][2][iP][iEta][iPhi]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_iP%i_iEta%i_iPhi%i", algo.Data (), iP, iEta, iPhi)));
           }
          }
          //nGamma[2][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
          //photonEnergyResponse[2][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("photonEnergyResponse_iP%i_iEta%i", iP, iEta)));
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }
     }
    }
   }
   cout << numFiles << " files read in." << endl;
  }
  /**** End loop over input files ****/

  /**** Save to output file ****/
  TFile* outFile = new TFile (Form ("%s/outFile.root", "read"));

  TVectorD nJetVec (2*3*(numpbins+1)*(numetabins+1)*numphibins);
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   for (short iMC = 0; iMC < 3; iMC++) {
    for (short iP = 0; iP <= numpbins; iP++) {
     for (short iEta = 0; iEta <= numetabins; iEta++) {
      for (short iPhi = 0; iPhi < numphibins; iPhi++) {
       nJetVec [iAlgo + (iMC + (iP + (iEta + (iPhi*(numetabins+1)))*(numpbins+1))*3)] = nJet[iAlgo][iMC][iP][iEta][iPhi];
       jetEnergyResponseCalib[iAlgo][iMC][iP][iEta][iPhi]->Write ();
       jetEnergyResponseReco[iAlgo][iMC][iP][iEta][iPhi]->Write ();
      }
     }
    }
   }
  }
  TVectorD nGammaVec (3*(numpbins+1)*(numetabins+1));
  for (short iMC = 0; iMC < 3; iMC++) {
   for (short iP = 0; iP <= numpbins; iP++) {
    for (short iEta = 0; iEta <= numetabins; iEta++) {
     nGammaVec [iMC + (iP + iEta*(numpbins*1))*3] = nGamma[iMC][iP][iEta];
     photonEnergyResponse[iMC][iP][iEta]->Write ();
    }
   }
  }
  nJetVec.Write ();
  nGammaVec.Write ();

  Delete5DArray (jetEnergyResponseCalib, 2, 3, numpbins+1, numetabins+1, numphibins);
  Delete5DArray (jetEnergyResponseReco, 2, 3, numpbins+1, numetabins+1, numphibins);
  Delete3DArray (photonEnergyResponse, 3, numpbins+1, numetabins+1);

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }

  return;
}

} // end namespace
