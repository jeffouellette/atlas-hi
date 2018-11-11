#include "PhotonAnalysis.h"
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

bool IsNonTight (const bool tight, const unsigned int isem) {
  return (!tight && !(0x45fc01 & isem));
}

bool IsIsolated (const float etcone30) {
  return etcone30 < 10;
}

bool IsNonIsolated (const float etcone30) {
  return etcone30 > 12;
}

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

float ShowerShapeBinsLow[12] = {-0.02, -0.04, 0, 0.7, 0.5, 0.2, 0.005, 0.5, -0.2, 0, -10, 0.6};
float ShowerShapeBinsHigh[12] = {0.02, 0.04, 100000, 1, 1.1, 1, 0.02, 4, 0.6, 0.9, 350, 1.1};

void PhotonAnalysis (const int dataSet, const bool isMC, const char* subdir) {

  SetAtlasStyle ();
  SetupDirectories ("", "OfflineAnalysis/");

  //TFile* inFile = GetFile ("Express/", dataSet, isMC);
  TFile* inFile = GetFile ((string(subdir)+"/").c_str(), dataSet, isMC);

  TTree* tree = (TTree*) inFile->Get ("jeffsjets");

  const int numEntries = tree->GetEntries ();

  vector<Trigger*> triggers = {};
  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
   Trigger* t = new Trigger (photonTrigNames[iTrig], photonTrigMinPtCuts[iTrig], -2.37, 2.37);
   tree->SetBranchAddress (t->name.c_str(), &(t->trigBool));
   tree->SetBranchAddress ((t->name+"_prescale").c_str(), &(t->trigPrescale));
   triggers.push_back (t);
  }

  unsigned int event_number = 0;
  float fcalA_et = 0;
  float fcalC_et = 0;

  int photon_n = 0;
  vector<float>* photon_pt = 0;
  vector<float>* photon_eta = 0;
  vector<float>* photon_phi = 0;
  vector<bool>* photon_tight = 0;
  vector<unsigned int>* photon_tight_isemvalue = 0;
  vector<float>* photon_etcone20 = 0;
  vector<float>* photon_etcone30 = 0;
  vector<float>* photon_etcone40 = 0;
  vector<float>* photon_rhad1 = 0;
  vector<float>* photon_rhad = 0;
  vector<float>* photon_e277 = 0;
  vector<float>* photon_reta = 0;
  vector<float>* photon_rphi = 0;
  vector<float>* photon_weta1 = 0;
  vector<float>* photon_weta2 = 0;
  vector<float>* photon_wtots1 = 0;
  vector<float>* photon_f1 = 0;
  vector<float>* photon_fracs1 = 0;
  vector<float>* photon_deltae = 0;
  vector<float>* photon_eratio = 0;

  int jet_n = 0;
  vector<float>* jet_pt = 0;
  vector<float>* jet_eta = 0;
  vector<float>* jet_phi = 0;
  vector<float>* jet_e = 0;

  tree->SetBranchAddress ("event_number", &event_number);
  tree->SetBranchAddress ("fcalA_et", &fcalA_et);
  tree->SetBranchAddress ("fcalC_et", &fcalC_et);
  tree->SetBranchAddress ("photon_n", &photon_n);
  tree->SetBranchAddress ("photon_pt", &photon_pt);
  tree->SetBranchAddress ("photon_eta", &photon_eta);
  tree->SetBranchAddress ("photon_phi", &photon_phi);
  tree->SetBranchAddress ("photon_tight", &photon_tight);
  tree->SetBranchAddress ("photon_tight_isemvalue", &photon_tight_isemvalue);
  tree->SetBranchAddress ("photon_etcone20", &photon_etcone20);
  tree->SetBranchAddress ("photon_etcone30", &photon_etcone30);
  tree->SetBranchAddress ("photon_etcone40", &photon_etcone40);
  tree->SetBranchAddress ("photon_rhad1", &photon_rhad1);
  tree->SetBranchAddress ("photon_rhad", &photon_rhad);
  tree->SetBranchAddress ("photon_e277", &photon_e277);
  tree->SetBranchAddress ("photon_reta", &photon_reta);
  tree->SetBranchAddress ("photon_rphi", &photon_rphi);
  tree->SetBranchAddress ("photon_w1", &photon_weta1);
  tree->SetBranchAddress ("photon_weta2", &photon_weta2);
  tree->SetBranchAddress ("photon_wtot", &photon_wtots1);
  tree->SetBranchAddress ("photon_f1", &photon_f1);
  tree->SetBranchAddress ("photon_fside", &photon_fracs1);
  tree->SetBranchAddress ("photon_deltae", &photon_deltae);
  tree->SetBranchAddress ("photon_eratio", &photon_eratio);

  tree->SetBranchAddress ("jet_n", &jet_n);
  tree->SetBranchAddress ("jet_pt", &jet_pt);
  tree->SetBranchAddress ("jet_eta", &jet_eta);
  tree->SetBranchAddress ("jet_phi", &jet_phi);
  tree->SetBranchAddress ("jet_e", &jet_e);

  TH1D**** ptSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 2);
  TH1D**** xJgDist = Get3DArray <TH1D*> (photonTrigN, 2, 2);

  TH2D** etaPhiMap = Get1DArray <TH2D*> (photonTrigN);

  TH1D**** etcone20 = Get3DArray <TH1D*> (photonTrigN, 2, 2);
  TH1D**** etcone30 = Get3DArray <TH1D*> (photonTrigN, 2, 2);
  TH1D**** etcone40 = Get3DArray <TH1D*> (photonTrigN, 2, 2);

  TH1D*** fcal_et = Get2DArray <TH1D*> (photonTrigN, 2);

  TH1D**** sidebandSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 4);
  TH1D**** showerShapeDists = Get3DArray <TH1D*> (photonTrigN, 2, 12);

  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
   etaPhiMap[iTrig] = new TH2D (Form ("etaPhiMap_%s", triggers[iTrig]->name.c_str()), "", 48, -2.4, 2.4, 48, -pi, pi);
   etaPhiMap[iTrig]->Sumw2 ();

   for (short iEta = 0; iEta < 2; iEta++) {
    const char* eta = iEta==0?"barrel":"endcap";

    for (short iCent = 0; iCent < 2; iCent++) {
     const char* cent = iCent==0?"cent":"periph";
     ptSpectrum[iTrig][iEta][iCent] = new TH1D (Form ("ptSpectrum_%s_%s_%s", triggers[iTrig]->name.c_str(), eta, cent), "", 150, 0, 300);
     ptSpectrum[iTrig][iEta][iCent]->Sumw2 ();

     xJgDist[iTrig][iEta][iCent] = new TH1D (Form ("xJgDist_%s_%s_%s", triggers[iTrig]->name.c_str(), eta, cent), "", 200, 0, 2);
     xJgDist[iTrig][iEta][iCent]->Sumw2 ();

     etcone20[iTrig][iEta][iCent] = new TH1D (Form ("etcone20_%s_%s_%s", triggers[iTrig]->name.c_str(), eta, cent), "", 90, -30, 60);
     etcone30[iTrig][iEta][iCent] = new TH1D (Form ("etcone30_%s_%s_%s", triggers[iTrig]->name.c_str(), eta, cent), "", 90, -30, 60);
     etcone40[iTrig][iEta][iCent] = new TH1D (Form ("etcone40_%s_%s_%s", triggers[iTrig]->name.c_str(), eta, cent), "", 90, -30, 60);

     etcone20[iTrig][iEta][iCent]->Sumw2 ();
     etcone30[iTrig][iEta][iCent]->Sumw2 ();
     etcone40[iTrig][iEta][iCent]->Sumw2 ();
    }

    fcal_et[iTrig][iEta] = new TH1D (Form ("fcal_et_%s_%s", triggers[iTrig]->name.c_str(), eta), "", 250, 0, 5);
    fcal_et[iTrig][iEta]->Sumw2 ();

    for (short iSide = 0; iSide < 4; iSide++) {
     char side = iSide==0?'a':(iSide==1?'b':(iSide==2?'c':'d'));
     sidebandSpectrum[iTrig][iEta][iSide] = new TH1D (Form ("sidebandSpectrum_%s_%s_%c", triggers[iTrig]->name.c_str(), eta, side), "", 250, 0, 500);
     sidebandSpectrum[iTrig][iEta][iSide]->Sumw2 ();
    }

    for (short iShower = 0; iShower < 12; iShower++) {
     showerShapeDists[iTrig][iEta][iShower] = new TH1D (Form ("showerShapeDist_%s_%s_%s", triggers[iTrig]->name.c_str(), eta, GetShowerShape (iShower)), "", 2000, ShowerShapeBinsLow[iShower], ShowerShapeBinsHigh[iShower]);
     showerShapeDists[iTrig][iEta][iShower]->Sumw2 ();
    }
   }
  }

  for (int entry = 0; entry < numEntries; entry++) {
   tree->GetEntry (entry);

   int iCent = (fcalA_et+fcalC_et)*1e-3 > 2 ? 0 : 1;

   for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
    if (! (triggers[iTrig]->trigBool))
     continue;

    for (int iP = 0; iP < photon_n; iP++) {

     int iEta = -1;
     if (fabs (photon_eta->at (iP)) <= 1.37)
      iEta = 0;
     else if (1.52 <= fabs (photon_eta->at (iP)) && fabs (photon_eta->at (iP)) <= 2.37)
      iEta = 1;

     if (photon_tight->at (iP) && photon_pt->at (iP) > 20)
      etaPhiMap[iTrig]->Fill (photon_eta->at (iP), photon_phi->at (iP));

     if (iEta == -1)
      continue;

     if (photon_tight->at (iP)) {
      // fill pT spectrum
      ptSpectrum[iTrig][iEta][iCent]->Fill (photon_pt->at (iP));

      // fill shower shape info
      if (photon_pt->at (iP) > 50) {
       showerShapeDists[iTrig][iEta][0]->Fill (photon_rhad1->at (iP));
       showerShapeDists[iTrig][iEta][1]->Fill (photon_rhad->at (iP));
       showerShapeDists[iTrig][iEta][2]->Fill (photon_e277->at (iP));
       showerShapeDists[iTrig][iEta][3]->Fill (photon_reta->at (iP));
       showerShapeDists[iTrig][iEta][4]->Fill (photon_rphi->at (iP));
       showerShapeDists[iTrig][iEta][5]->Fill (photon_weta1->at (iP));
       showerShapeDists[iTrig][iEta][6]->Fill (photon_weta2->at (iP));
       showerShapeDists[iTrig][iEta][7]->Fill (photon_wtots1->at (iP));
       showerShapeDists[iTrig][iEta][8]->Fill (photon_f1->at (iP));
       showerShapeDists[iTrig][iEta][9]->Fill (photon_fracs1->at (iP));
       showerShapeDists[iTrig][iEta][10]->Fill (photon_deltae->at (iP));
       showerShapeDists[iTrig][iEta][11]->Fill (photon_eratio->at (iP));
      }       

      if (IsIsolated (photon_etcone30->at (iP)) && photon_pt->at (iP) > 60) {
       int lJ = -1;
       for (int iJ = 0; iJ < jet_n; iJ++) {
        if (DeltaPhi (jet_phi->at (iJ), photon_phi->at (iP)) < 7*pi/8)
         continue;

        if (lJ == -1 || jet_pt->at (lJ) < jet_pt->at (iJ))
         lJ = iJ;
       }
       if (lJ != -1)
        xJgDist[iTrig][iEta][iCent]->Fill (jet_pt->at (lJ) / photon_pt->at (iP));
      }

      if (photon_pt->at (iP) > 20) {
       etcone20[iTrig][iEta][iCent]->Fill (photon_etcone20->at (iP));
       etcone30[iTrig][iEta][iCent]->Fill (photon_etcone30->at (iP));
       etcone40[iTrig][iEta][iCent]->Fill (photon_etcone40->at (iP));

       fcal_et[iTrig][iEta]->Fill ((fcalA_et+fcalC_et)*1e-3);
      }
 
     }

     if (photon_tight->at (iP) && IsIsolated (photon_etcone30->at (iP)))
      sidebandSpectrum[iTrig][iEta][0]->Fill (photon_pt->at (iP));
     else if (photon_tight->at (iP) && IsNonIsolated (photon_etcone30->at (iP)))
      sidebandSpectrum[iTrig][iEta][1]->Fill (photon_pt->at (iP));
     else if (IsNonTight (photon_tight->at (iP), photon_tight_isemvalue->at (iP)) && IsIsolated (photon_etcone30->at (iP)))
      sidebandSpectrum[iTrig][iEta][2]->Fill (photon_pt->at (iP));
     else if (IsNonTight (photon_tight->at (iP), photon_tight_isemvalue->at (iP)) && IsNonIsolated (photon_etcone30->at (iP)))
      sidebandSpectrum[iTrig][iEta][3]->Fill (photon_pt->at (iP));

     if (photon_pt->at (iP) > 250)
      cout << "Woah! Check out this event: run " << dataSet << ", tree entry " << entry << "!!" << endl;
    }
   }
  }

  TFile* outFile = new TFile (Form ("%s/PhotonAnalysis/dataSet_%i.root", rootPath.Data (), dataSet), "recreate");
  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
   etaPhiMap[iTrig]->Write ();

   for (short iEta = 0; iEta < 2; iEta++) {
    for (short iCent = 0; iCent < 2; iCent++) {
     ptSpectrum[iTrig][iEta][iCent]->Write ();
     xJgDist[iTrig][iEta][iCent]->Write ();

     etcone20[iTrig][iEta][iCent]->Write ();
     etcone30[iTrig][iEta][iCent]->Write ();
     etcone40[iTrig][iEta][iCent]->Write ();
    }

    fcal_et[iTrig][iEta]->Write ();

    for (short iSide = 0; iSide < 4; iSide++) {
     sidebandSpectrum[iTrig][iEta][iSide]->Write ();
    }

    for (short iShower = 0; iShower < 12; iShower++) {
     showerShapeDists[iTrig][iEta][iShower]->Write ();
    }
   }
  }
  outFile->Close ();

}

} // end namespace
