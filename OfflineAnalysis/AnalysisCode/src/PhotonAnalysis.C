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

void PhotonAnalysis (const int dataSet, const bool isMC, const char* subdir) {

  SetupDirectories ("", "OfflineAnalysis/");

  //TFile* inFile = GetFile ("Express/", dataSet, isMC);
  TFile* inFile = GetFile ((string(subdir)+"/").c_str (), dataSet, isMC);

  TTree* tree = (TTree*) inFile->Get ("jeffsjets");

  const int numEntries = tree->GetEntries ();

  int iRunNumber = 0;
  if (!isMC) {
    while (iRunNumber < numRunNumbers && runNumbers[iRunNumber] != dataSet)
      iRunNumber++;
  }

  vector<Trigger*> triggers = {};
  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
    Trigger* t = new Trigger (photonTrigNames[iTrig], photonTrigMinPtCuts[iTrig], -2.37, 2.37);
    tree->SetBranchAddress (t->name.c_str (), &(t->trigBool));
    tree->SetBranchAddress ((t->name+"_prescale").c_str (), &(t->trigPrescale));
    triggers.push_back (t);
  }

  unsigned int lumi_block = 0;
  unsigned int event_number = 0;
  float fcalA_et = 0;
  float fcalC_et = 0;
  float fcalA_et_Cos = 0;
  float fcalC_et_Cos = 0;
  float fcalA_et_Sin = 0;
  float fcalC_et_Sin = 0;

  int photon_n = 0;
  vector<float>* photon_pt = 0;
  vector<float>* photon_eta = 0;
  vector<float>* photon_phi = 0;
  vector<bool>* photon_tight = 0;
  vector<unsigned int>* photon_tight_isemvalue = 0;
  vector<float>* photon_etcone20 = 0;
  vector<float>* photon_etcone30 = 0;
  vector<float>* photon_etcone40 = 0;
  vector<float>* photon_Rhad1 = 0;
  vector<float>* photon_Rhad = 0;
  vector<float>* photon_e277 = 0;
  vector<float>* photon_Reta = 0;
  vector<float>* photon_Rphi = 0;
  vector<float>* photon_weta1 = 0;
  vector<float>* photon_weta2 = 0;
  vector<float>* photon_wtots1 = 0;
  vector<float>* photon_f1 = 0;
  vector<float>* photon_f3 = 0;
  vector<float>* photon_fracs1 = 0;
  vector<float>* photon_DeltaE = 0;
  vector<float>* photon_Eratio = 0;

  int jet_n = 0;
  vector<float>* jet_pt = 0;
  vector<float>* jet_eta = 0;
  vector<float>* jet_phi = 0;
  vector<float>* jet_e = 0;

  tree->SetBranchAddress ("lumi_block", &lumi_block);
  tree->SetBranchAddress ("event_number", &event_number);

  tree->SetBranchAddress ("fcalA_et", &fcalA_et);
  tree->SetBranchAddress ("fcalC_et", &fcalC_et);
  tree->SetBranchAddress ("fcalA_et_Cos", &fcalA_et_Cos);
  tree->SetBranchAddress ("fcalC_et_Cos", &fcalC_et_Cos);
  tree->SetBranchAddress ("fcalA_et_Sin", &fcalA_et_Sin);
  tree->SetBranchAddress ("fcalC_et_Sin", &fcalC_et_Sin);

  tree->SetBranchAddress ("photon_n", &photon_n);
  tree->SetBranchAddress ("photon_pt", &photon_pt);
  tree->SetBranchAddress ("photon_eta", &photon_eta);
  tree->SetBranchAddress ("photon_phi", &photon_phi);
  tree->SetBranchAddress ("photon_tight", &photon_tight);
  tree->SetBranchAddress ("photon_tight_isemvalue", &photon_tight_isemvalue);
  tree->SetBranchAddress ("photon_etcone20", &photon_etcone20);
  tree->SetBranchAddress ("photon_etcone30", &photon_etcone30);
  tree->SetBranchAddress ("photon_etcone40", &photon_etcone40);

  tree->SetBranchAddress ("photon_Rhad1", &photon_Rhad1);
  tree->SetBranchAddress ("photon_Rhad", &photon_Rhad);
  tree->SetBranchAddress ("photon_e277", &photon_e277);
  tree->SetBranchAddress ("photon_Reta", &photon_Reta);
  tree->SetBranchAddress ("photon_Rphi", &photon_Rphi);
  tree->SetBranchAddress ("photon_weta1", &photon_weta1);
  tree->SetBranchAddress ("photon_weta2", &photon_weta2);
  tree->SetBranchAddress ("photon_wtots1", &photon_wtots1);
  tree->SetBranchAddress ("photon_f1", &photon_f1);
  tree->SetBranchAddress ("photon_f3", &photon_f3);
  tree->SetBranchAddress ("photon_fracs1", &photon_fracs1);
  tree->SetBranchAddress ("photon_DeltaE", &photon_DeltaE);
  tree->SetBranchAddress ("photon_Eratio", &photon_Eratio);

  tree->SetBranchAddress ("akt4hi_jet_n", &jet_n);
  tree->SetBranchAddress ("akt4hi_jet_pt", &jet_pt);
  tree->SetBranchAddress ("akt4hi_jet_eta", &jet_eta);
  tree->SetBranchAddress ("akt4hi_jet_phi", &jet_phi);
  tree->SetBranchAddress ("akt4hi_jet_e", &jet_e);

  TH1D**** ptSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 2);

  TH1D***** xJgNgamma = Get4DArray <TH1D*> (photonTrigN, 2, 2, 3);
  TH1D***** xJgDist = Get4DArray <TH1D*> (photonTrigN, 2, 2, 3);

  TH2D** etaPhiMap = Get1DArray <TH2D*> (photonTrigN);

  TH1D***** etconeDists = Get4DArray <TH1D*> (3, photonTrigN, 2, 2);
  TH2D**** etconeFCalDists = Get3DArray <TH2D*> (3, photonTrigN, 2);

  TH1D***** etconeEventPlaneDists = Get4DArray <TH1D*> (3, photonTrigN, 2, 3);

  TH1D*** fcal_et = Get2DArray <TH1D*> (photonTrigN, 2);

  TH1D**** sidebandSpectrum = Get3DArray <TH1D*> (photonTrigN, 2, 4);
  TH1D**** showerShapeDists = Get3DArray <TH1D*> (photonTrigN, 2, 13);

  TH1D** photonCounts = Get1DArray <TH1D*> (photonTrigN);
  TH1D*** photonEventPlaneDist = Get2DArray <TH1D*> (photonTrigN, 2);

  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
    etaPhiMap[iTrig] = new TH2D (Form ("etaPhiMap_%s", triggers[iTrig]->name.c_str ()), "", 48, -2.4, 2.4, 48, -pi, pi);
    etaPhiMap[iTrig]->Sumw2 ();
 
    photonCounts[iTrig] = new TH1D (Form ("photonCounts_%s", triggers[iTrig]->name.c_str ()), "", numRunNumbers, -0.5, -0.5+numRunNumbers);
    photonCounts[iTrig]->Sumw2 ();

    for (short iEta = 0; iEta < 2; iEta++) {
      const char* eta = iEta==0?"barrel":"endcap";

      photonEventPlaneDist[iTrig][iEta] = new TH1D (Form ("photonEventPlaneDist_%s_%s", triggers[iTrig]->name.c_str (), eta), "", 90, 0, pi/2);
      photonEventPlaneDist[iTrig][iEta]->Sumw2 ();

      for (short iAngle = 0; iAngle < 3; iAngle++) {
        const char* angle = iAngle==0?"0_pi6":(iAngle==1?"pi6_pi3":"pi3_pi2");
        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          const char etcone = iEtcone==0?'2':(iEtcone==1?'3':'4');
          etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle] = new TH1D (Form ("etcone%c0_%s_%s_%s", etcone, triggers[iTrig]->name.c_str (), eta, angle), "", 90, -30, 60);
          etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle]->Sumw2 ();
        }
      }

      for (short iCent = 0; iCent < 2; iCent++) {
        const char* cent = iCent==0?"cent":"periph";
        ptSpectrum[iTrig][iEta][iCent] = new TH1D (Form ("ptSpectrum_%s_%s_%s", triggers[iTrig]->name.c_str (), eta, cent), "", 600, 50, 650);
        ptSpectrum[iTrig][iEta][iCent]->Sumw2 ();

        for (short iAngle = 0; iAngle < 3; iAngle++) {
          const char* angle = iAngle==0?"0_pi6":(iAngle==1?"pi6_pi3":"pi3_pi2");
          xJgNgamma[iTrig][iEta][iCent][iAngle] = new TH1D (Form ("xJgNgamma_%s_%s_%s_%s", triggers[iTrig]->name.c_str (), eta, cent, angle), "", 250, 0, 500);
          xJgNgamma[iTrig][iEta][iCent][iAngle]->Sumw2 ();
          xJgDist[iTrig][iEta][iCent][iAngle] = new TH1D (Form ("xJgDist_%s_%s_%s_%s", triggers[iTrig]->name.c_str (), eta, cent, angle), "", 200, 0, 2);
          xJgDist[iTrig][iEta][iCent][iAngle]->Sumw2 ();
        }

        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          const char etcone = iEtcone==0?'2':(iEtcone==1?'3':'4');
          etconeDists[iEtcone][iTrig][iEta][iCent] = new TH1D (Form ("etcone%c0_%s_%s_%s", etcone, triggers[iTrig]->name.c_str (), eta, cent), "", 90, -30, 60);
          etconeDists[iEtcone][iTrig][iEta][iCent]->Sumw2 ();
        }
      }
      for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
        const char etcone = iEtcone==0?'2':(iEtcone==1?'3':'4');
        etconeFCalDists[iEtcone][iTrig][iEta] = new TH2D (Form ("etcone%c0FCal_%s_%s", etcone, triggers[iTrig]->name.c_str (), eta), "", 5, 0, 5, 90, -30, 60);
        etconeFCalDists[iEtcone][iTrig][iEta]->Sumw2 ();
      }

      fcal_et[iTrig][iEta] = new TH1D (Form ("fcal_et_%s_%s", triggers[iTrig]->name.c_str (), eta), "", 250, 0, 5);
      fcal_et[iTrig][iEta]->Sumw2 ();

      for (short iSide = 0; iSide < 4; iSide++) {
        char side = iSide==0?'a':(iSide==1?'b':(iSide==2?'c':'d'));
        sidebandSpectrum[iTrig][iEta][iSide] = new TH1D (Form ("sidebandSpectrum_%s_%s_%c", triggers[iTrig]->name.c_str (), eta, side), "", numpbins, pbins);
        sidebandSpectrum[iTrig][iEta][iSide]->Sumw2 ();
      }

      for (short iShower = 0; iShower < 13; iShower++) {
        showerShapeDists[iTrig][iEta][iShower] = new TH1D (Form ("showerShapeDist_%s_%s_%s", triggers[iTrig]->name.c_str (), eta, GetShowerShape (iShower)), "", 2000, showerShapeBinsLow[iShower], showerShapeBinsHigh[iShower]);
        showerShapeDists[iTrig][iEta][iShower]->Sumw2 ();
      }
    }
  }

  for (int entry = 0; entry < numEntries; entry++) {
    tree->GetEntry (entry);

    // process FCal and event plane info
    int iCent = (fcalA_et+fcalC_et)*1e-3 > 2 ? 0 : 1;
    const double qx = fcalA_et_Cos + fcalC_et_Cos;
    const double qy = fcalA_et_Sin + fcalC_et_Sin;
    const double psi2 = 0.5 * atan2 (qy, qx);

    // check each trigger
    for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
      if (! (triggers[iTrig]->trigBool))
        continue;

      // loop over photons
      for (int iP = 0; iP < photon_n; iP++) {

        // check whether in barrel, endcaps, or cracks
        int iEta = -1;
        if (fabs (photon_eta->at (iP)) <= 1.37)
          iEta = 0;
        else if (1.52 <= fabs (photon_eta->at (iP)) && fabs (photon_eta->at (iP)) <= 2.37)
          iEta = 1;

        double delta_phi = DeltaPhi (psi2, photon_phi->at (iP));
        if (delta_phi > pi/2.)
          delta_phi = pi - delta_phi;
        int iAngle = -1; 
        if (delta_phi < pi/6.)
          iAngle = 0;
        else if (delta_phi < pi/3.)
          iAngle = 1;
        else
          iAngle = 2;

        // fill eta-phi map
        if (photon_tight->at (iP) && photon_pt->at (iP) > 20)
          etaPhiMap[iTrig]->Fill (photon_eta->at (iP), photon_phi->at (iP));

        if (iEta == -1)
          continue;

        // fill tight, isolated counts for this run
        if (photon_tight->at (iP) && IsIsolated (photon_etcone30->at (iP)) && photon_pt->at (iP) > 60) {
          photonCounts[iTrig]->Fill (iRunNumber);

          photonEventPlaneDist[iTrig][iEta]->Fill (delta_phi);
        }

        if (photon_tight->at (iP)) {
          // fill pT spectrum
          ptSpectrum[iTrig][iEta][iCent]->Fill (photon_pt->at (iP));

          // fill shower shape info
          if (photon_pt->at (iP) > 50 && IsIsolated (photon_etcone30->at (iP))) {
            showerShapeDists[iTrig][iEta][0]->Fill (photon_Rhad1->at (iP));
            showerShapeDists[iTrig][iEta][1]->Fill (photon_Rhad->at (iP));
            showerShapeDists[iTrig][iEta][2]->Fill (photon_e277->at (iP));
            showerShapeDists[iTrig][iEta][3]->Fill (photon_Reta->at (iP));
            showerShapeDists[iTrig][iEta][4]->Fill (photon_Rphi->at (iP));
            showerShapeDists[iTrig][iEta][5]->Fill (photon_weta1->at (iP));
            showerShapeDists[iTrig][iEta][6]->Fill (photon_weta2->at (iP));
            showerShapeDists[iTrig][iEta][7]->Fill (photon_wtots1->at (iP));
            showerShapeDists[iTrig][iEta][8]->Fill (photon_f1->at (iP));
            showerShapeDists[iTrig][iEta][9]->Fill (photon_f3->at (iP));
            showerShapeDists[iTrig][iEta][10]->Fill (photon_fracs1->at (iP));
            showerShapeDists[iTrig][iEta][11]->Fill (photon_DeltaE->at (iP));
            showerShapeDists[iTrig][iEta][12]->Fill (photon_Eratio->at (iP));
          }       

          // fill xJgamma distributions
          if (IsIsolated (photon_etcone30->at (iP)) && photon_pt->at (iP) > 60) {
            xJgNgamma[iTrig][iEta][iCent][iAngle]->Fill (photon_pt->at (iP));

            int lJ = -1;
            for (int iJ = 0; iJ < jet_n; iJ++) {
              if (DeltaPhi (jet_phi->at (iJ), photon_phi->at (iP)) < 7*pi/8)
                continue;

              if (lJ == -1 || jet_pt->at (lJ) < jet_pt->at (iJ))
                lJ = iJ;
            }
            if (lJ != -1) {
              xJgDist[iTrig][iEta][iCent][iAngle]->Fill (jet_pt->at (lJ) / photon_pt->at (iP));
              if (photon_pt->at (iP) > 250 && photon_tight->at (iP) && IsIsolated (photon_etcone30->at (iP)) && jet_pt->at(lJ) / photon_pt->at(iP) < 0.7)
                cout << "Woah! Check out this event: run " << dataSet << ", event " << event_number << ", lumiblock " << lumi_block << ", tree entry " << entry << "!!" << endl;
            }

          }

          // fill etcone distributions
          if (photon_pt->at (iP) > 20) {
            etconeDists[0][iTrig][iEta][iCent]->Fill (photon_etcone20->at (iP));
            etconeDists[1][iTrig][iEta][iCent]->Fill (photon_etcone30->at (iP));
            etconeDists[2][iTrig][iEta][iCent]->Fill (photon_etcone40->at (iP));

            etconeFCalDists[0][iTrig][iEta]->Fill ((fcalA_et+fcalC_et)*1e-3, photon_etcone20->at (iP));
            etconeFCalDists[1][iTrig][iEta]->Fill ((fcalA_et+fcalC_et)*1e-3, photon_etcone30->at (iP));
            etconeFCalDists[2][iTrig][iEta]->Fill ((fcalA_et+fcalC_et)*1e-3, photon_etcone40->at (iP));

            fcal_et[iTrig][iEta]->Fill ((fcalA_et+fcalC_et)*1e-3);
          }

          // fill etcone vs. angle w.r.t. event plane
          if (photon_pt->at (iP) > 20) {
            etconeEventPlaneDists[0][iTrig][iEta][iAngle]->Fill (photon_etcone20->at (iP));
            etconeEventPlaneDists[1][iTrig][iEta][iAngle]->Fill (photon_etcone30->at (iP));
            etconeEventPlaneDists[2][iTrig][iEta][iAngle]->Fill (photon_etcone40->at (iP));
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

      }
    }
  }

  TFile* outFile = new TFile (Form ("%s/PhotonAnalysis/dataSet_%i.root", rootPath.Data (), dataSet), "recreate");
  for (int iTrig = 0; iTrig < photonTrigN; iTrig++) {
    etaPhiMap[iTrig]->Write ();
    photonCounts[iTrig]->Scale (1. / lumis[iRunNumber]);
    photonCounts[iTrig]->Write ();

    for (short iEta = 0; iEta < 2; iEta++) {

      for (short iAngle = 0; iAngle < 3; iAngle++) {
        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          etconeEventPlaneDists[iEtcone][iTrig][iEta][iAngle]->Write ();
        }
      }

      for (short iCent = 0; iCent < 2; iCent++) {
        ptSpectrum[iTrig][iEta][iCent]->Write ();

        for (short iAngle = 0; iAngle < 3; iAngle++) {
          xJgNgamma[iTrig][iEta][iCent][iAngle]->Write ();
          xJgDist[iTrig][iEta][iCent][iAngle]->Write ();
        }

        for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
          etconeDists[iEtcone][iTrig][iEta][iCent]->Write ();
        }
      }
      for (short iEtcone = 0; iEtcone < 3; iEtcone++) {
        etconeFCalDists[iEtcone][iTrig][iEta]->Write ();
      }

      photonEventPlaneDist[iTrig][iEta]->Write ();

      fcal_et[iTrig][iEta]->Write ();

      for (short iSide = 0; iSide < 4; iSide++) {
        sidebandSpectrum[iTrig][iEta][iSide]->Write ();
      }

      for (short iShower = 0; iShower < 13; iShower++) {
        showerShapeDists[iTrig][iEta][iShower]->Write ();
      }
    }
  }
  outFile->Close ();

}

} // end namespace
