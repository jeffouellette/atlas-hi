#ifndef __ppMixingStudy_C__
#define __ppMixingStudy_C__

void ppMixingStudy () {

  TFile* f_mb = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/FCalCalibration/Nominal/data18hi.root", "read");
  TFile* f_ztag = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/pp_fcal_et_ztag.root", "read");

  TH1D* h_mb = (TH1D*) f_mb->Get ("h_pp_fcal_et");
  TH1D* h_ztag = (TH1D*) f_ztag->Get ("h_pp_fcal_et_ztag");

  TH1D* h_ztag_sim = new TH1D ("h_ztag_sim", "", 150, -100, 200);
  TH1D* h_mb_sim = new TH1D ("h_mb_sim", "", 150, -100, 200);
  TH1D* h_diff = new TH1D ("h_diff", "", 150, -100, 200);

  for (int i = 0; i < 1000000; i++) {
    float ztag = h_ztag->GetRandom ();
    float mb = h_mb->GetRandom ();

    h_ztag_sim->Fill (ztag);
    h_mb_sim->Fill (mb);
    h_diff->Fill (ztag-mb);
  }

  h_diff->Draw ();

  h_ztag_sim->SetLineStyle (2);
  h_mb_sim->SetLineStyle (2);

  h_ztag_sim->Draw ("hist same");
  h_mb_sim->Draw ("hist same");

}

#endif
