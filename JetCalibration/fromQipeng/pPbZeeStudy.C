#define pPbZeeStudy_cxx
#include "pPbZeeStudy.h"
#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void pPbZeeStudy::Loop()
{
//   In a ROOT session, you can do:
//      root> .L pPbZeeStudy.C
//      root> pPbZeeStudy t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile* fout = new TFile("hist_data.root","RECREATE");
   //TFile* fout = new TFile("hist_mc_overlay8.root","RECREATE");
   //TFile* fout = new TFile("hist_mc_overlay4.root","RECREATE");
   
   TH1F* h_ZeePt = new TH1F("h_ZeePt", "", 100, 0, 100);

   TH1F* h_ZeeMass_eta0     = new TH1F("h_ZeeMass_eta0",    "",60,60,120);
   TH1F* h_ZeeMass_eta0_raw = new TH1F("h_ZeeMass_eta0_raw","",60,60,120);

   TH1F* h_ZeeMass_eta1     = new TH1F("h_ZeeMass_eta1",    "",60,60,120);
   TH1F* h_ZeeMass_eta1_raw = new TH1F("h_ZeeMass_eta1_raw","",60,60,120);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (!HLT_e15_lhloose_nod0) continue; 
      for (size_t i=0; i<electron_e->size(); i++) {
        if (fabs(electron_etaBE2->at(i)) > 1.37 && fabs(electron_etaBE2->at(i)) < 1.52) continue;
        if (fabs(electron_d0Sig->at(i)) > 5) continue;
        if (fabs(electron_dz0SinTheta->at(i)) > 0.5) continue;
        if (!electron_isMedium->at(i)) continue;

        TLorentzVector _mom1;
        TLorentzVector _mom1_raw;
        _mom1.SetPtEtaPhiM(electron_pt->at(i), electron_eta->at(i), electron_phi->at(i), 0.511e-3);
        _mom1_raw.SetPtEtaPhiM(electron_pt_raw->at(i), electron_eta->at(i), electron_phi->at(i), 0.511e-3);

        for (size_t j=i+1; j<electron_e->size(); j++) {
          if (fabs(electron_etaBE2->at(j)) > 1.37 && fabs(electron_etaBE2->at(j)) < 1.52) continue;
          if (fabs(electron_d0Sig->at(j)) > 5) continue;
          if (fabs(electron_dz0SinTheta->at(j)) > 0.5) continue;
          if (!electron_isMedium->at(j)) continue;

          if (!electron_match_e15_nod0->at(i) && !electron_match_e15_nod0->at(j)) continue;
          if (electron_charge->at(i) == electron_charge->at(j)) continue;

          TLorentzVector _mom2;
          TLorentzVector _mom2_raw;
          _mom2.SetPtEtaPhiM(electron_pt->at(j), electron_eta->at(j), electron_phi->at(j), 0.511e-3);
          _mom2_raw.SetPtEtaPhiM(electron_pt_raw->at(j), electron_eta->at(j), electron_phi->at(j), 0.511e-3);
          TLorentzVector _mom = _mom1 + _mom2;
          TLorentzVector _mom_raw = _mom1_raw + _mom2_raw;
          if (_mom.M() > 50 && _mom.M() < 200) {
            h_ZeePt->Fill(_mom.Pt());
          }
          if (_mom.M() > 60 && _mom.M() < 120) {
            if (fabs(_mom.Rapidity()) <= 1.) h_ZeeMass_eta0->Fill(_mom.M());
            else h_ZeeMass_eta1->Fill(_mom.M());
          }
          if (_mom_raw.M() > 60 && _mom_raw.M() < 120) {
            if (fabs(_mom_raw.Rapidity()) <= 1.) h_ZeeMass_eta0_raw->Fill(_mom_raw.M());
            else h_ZeeMass_eta1_raw->Fill(_mom_raw.M());
          }
        }
      }
    }
    //h_ZeeMass->Scale(1./h_ZeeMass->Integral());
    //h_ZeeMass->Draw("E");

    //h_ZeeMass_raw->Scale(1./h_ZeeMass_raw->Integral());
    //h_ZeeMass_raw->SetMarkerStyle(24);
    //h_ZeeMass_raw->SetMarkerColor(2);
    //h_ZeeMass_raw->SetLineColor(2);
    //h_ZeeMass_raw->Draw("ESAME");

    h_ZeePt->Write();
    h_ZeeMass_eta0->Write();
    h_ZeeMass_eta0_raw->Write();
    h_ZeeMass_eta1->Write();
    h_ZeeMass_eta1_raw->Write();

}
