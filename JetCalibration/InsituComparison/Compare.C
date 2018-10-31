
void Compare () {
  TFile* f2015 = new TFile ("InsituCalibration_20_7_4EM_4LC_30Sep2016.root", "read");
  TFile* f2016 = new TFile ("InsituCalibration_2016data_19Dec2016.root", "read");

  TH1D* insitu2015 = (TH1D*)f2015->Get("AntiKt4EMTopo_InsituCalib");
  TH1D* insitu2016 = (TH1D*)f2016->Get("AntiKt4EMTopo_InsituCalib");

  TCanvas* c = new TCanvas ("c", "", 800, 600);
  FormatTH2Canvas (c, false);
  gPad->SetLogx();

  insitu2015->SetLineColor(kBlack);
  insitu2016->SetLineColor(kBlue);

  insitu2015->GetXaxis()->SetTitle ("jet #it{p}_{T} #left[GeV#right]");
  insitu2015->GetYaxis()->SetTitle ("<Insitu Factor>");

  insitu2015->Draw("hist");
  insitu2016->Draw("same hist");

  myMarkerText (0.22, 0.88, kBlack, kFullCircle, "2015 Insitu Factors", 1.25, 0.06);
  myMarkerText (0.22, 0.80, kBlue, kFullCircle, "2016 Insitu Factors", 1.25, 0.06);

  c->SaveAs("CompareSummary.pdf");

  TStyle* myStyle = AtlasStyle();
  myStyle->SetPalette(55);
  SetAtlasStyle();

  TH2D* etaInsitu2015 = (TH2D*)f2015->Get("AntiKt4EMTopo_EtaInterCalibration");
  TH2D* etaInsitu2016 = (TH2D*)f2016->Get("AntiKt4EMTopo_EtaInterCalibration");

  etaInsitu2015->Draw("col");
  c->SaveAs ("EtaInsitu2015.pdf");

  etaInsitu2016->Draw("col");
  c->SaveAs ("EtaInsitu2016.pdf");

  FormatTH2Canvas (c);

  TH2D* ratio = (TH2D*)etaInsitu2016->Clone();
  ratio->SetName("AntiKt4EMTopo_EtaInterCalibration_2016_over_2015");
  ratio->Divide(etaInsitu2015);

  ratio->GetYaxis()->SetTitleOffset (0.75);
  ratio->GetZaxis()->SetTitleOffset (1.3);
  ratio->GetZaxis()->SetTitle ("Ratio of Insitu Factors, 2016 / 2015");

  ratio->Draw("colz");

  c->SaveAs("CompareEtaInterCcalibration.pdf");

  
}
