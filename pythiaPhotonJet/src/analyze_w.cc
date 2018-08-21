#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>

#include <iostream>
#include <sstream>

//#include "eps09_cxx/eps09.h"

#include <AtlasStyle.h>
#include <AtlasUtils.h>

int main() {

  SetAtlasStyle();

  /*
  TFile *f[7];
  f[0] = new TFile("output/gammajet_pt10_10k.root","READ");
  f[1] = new TFile("output/gammajet_pt30_10k.root","READ");
  f[2] = new TFile("output/gammajet_pt45_10k.root","READ");
  f[3] = new TFile("output/gammajet_pt65_10k.root","READ");
  f[4] = new TFile("output/gammajet_pt135_10k.root","READ");
  f[5] = new TFile("output/gammajet_pt275_10k.root","READ");
  f[6] = new TFile("output/gammajet_pt495_10k.root","READ");

  float sigma[7] = {
    7.278e-04,
    2.395e-05,
    5.932e-06,
    1.573e-06,
    9.134e-08,
    4.049e-09,
    1.959e-10
  };
  */

  TFile *f[7];
  f[0] = new TFile("output/gammajet_pt12_100k.root","READ");
  f[1] = new TFile("output/gammajet_pt26_100k.root","READ");
  f[2] = new TFile("output/gammajet_pt37_100k.root","READ");
  f[3] = new TFile("output/gammajet_pt52_100k.root","READ");
  f[4] = new TFile("output/gammajet_pt105_100k.root","READ");
  f[5] = new TFile("output/gammajet_pt210_100k.root","READ");
  f[6] = new TFile("output/gammajet_pt375_100k.root","READ");

  float sigma[7] = {
    4.262e-04,
    3.832e-05,
    1.166e-05,
    3.533e-06,
    2.519e-07,
    1.391e-08,
    8.773e-10
  };


  //f[0]->ls();

  TTree *t[7];
  for (int n = 0; n < 7; n++) 
    t[n] = (TTree*) f[n]->Get("tree");
  
  //t[0]->Show(0);

  double ruv, rdv, ru, rd, rs, rc, rb, rg;

  //std::cout << ruv << " " << rdv << " " << ru << " " << rd << " " << rs << " " << rc << " " << rb << " " << rg << std::endl;

  double PTRANGES[ 8 ] = { 15.85, 35, 50, 70, 140, 280, 500, 1000 };

  //const int NPTBINS = 20;
  //double PTBINS[21];
  const int NPTBINS = 17;
  double PTBINS[] = {25,35,45,55,65,75,85,105,125,150,175,200,250,300,350,400,470,550};
  double alpha = pow( 10.0, 1 / 10.0 );
  //for (int n = 0; n < 21; n++) {
  //  PTBINS[ n ] = 10 * pow( alpha, n );
  //  //std::cout << PTBINS[ n ] << std::endl;
  //}

  const int NPTBINS2 = 170;
  double PTBINS2[ NPTBINS2+1 ];
  double alpha2 = pow( 10.0, 1 / 100.0 );
  for (int n = 0; n < NPTBINS2+1; n++) {
    PTBINS2[ n ] = 10 * pow( alpha2, n );
    std::cout << n << " " << PTBINS2[ n ] << std::endl;
  }

  TH1D *h1_pt_ETA[3];
  h1_pt_ETA[0] = new TH1D("h1_pt_ETA0","", NPTBINS, PTBINS );
  h1_pt_ETA[1] = new TH1D("h1_pt_ETA1","", NPTBINS, PTBINS );
  h1_pt_ETA[2] = new TH1D("h1_pt_ETA2","", NPTBINS, PTBINS );

  TH2D *h2_pt_xA_ETA[3];
  h2_pt_xA_ETA[0] = new TH2D("h2_pt_xA_ETA0","", NPTBINS, PTBINS, 10000, 0, 1 );
  h2_pt_xA_ETA[1] = new TH2D("h2_pt_xA_ETA1","", NPTBINS, PTBINS, 10000, 0, 1 );
  h2_pt_xA_ETA[2] = new TH2D("h2_pt_xA_ETA2","", NPTBINS, PTBINS, 10000, 0, 1 );

  TH2D *h2_pt_xp_ETA[3];
  h2_pt_xp_ETA[0] = new TH2D("h2_pt_xp_ETA0","", NPTBINS, PTBINS, 10000, 0, 1 );
  h2_pt_xp_ETA[1] = new TH2D("h2_pt_xp_ETA1","", NPTBINS, PTBINS, 10000, 0, 1 );
  h2_pt_xp_ETA[2] = new TH2D("h2_pt_xp_ETA2","", NPTBINS, PTBINS, 10000, 0, 1 );

  TProfile *h1_pt_meanxA_ETA[3];
  TProfile *h1_pt_meanxp_ETA[3];

  TH1D *h1_pt = new TH1D("h1_pt","",NPTBINS2,PTBINS2);
  h1_pt->Sumw2();

  TH1D *h1_pt_alt_ETA[3];
  h1_pt_alt_ETA[0] = new TH1D("h1_pt_alt_ETA0","", NPTBINS, PTBINS );
  h1_pt_alt_ETA[1] = new TH1D("h1_pt_alt_ETA1","", NPTBINS, PTBINS );
  h1_pt_alt_ETA[2] = new TH1D("h1_pt_alt_ETA2","", NPTBINS, PTBINS );


  TH1D *hFrame;
  hFrame = new TH1D("hFrame","", NPTBINS-1, PTBINS );
  TH1D *hFrame2;
  hFrame2 = new TH1D("hFrame2","", NPTBINS-1, PTBINS );
  

  int id2;
  float x2pdf;
  float x1pdf;
  float Q;

  int photon_n;
  float photon_pt[50];
  float photon_eta[50];
  float photon_iso[50];

  for (int n = 0; n < 7; n++) {

    // beam-1 is the nucleus
    // so boost is in the direction of beam-2...

    t[n]->SetBranchAddress("id2",&id2);
    t[n]->SetBranchAddress("x2pdf",&x2pdf);
    t[n]->SetBranchAddress("x1pdf",&x1pdf);
    t[n]->SetBranchAddress("Q",&Q);

    t[n]->SetBranchAddress("photon_n",&photon_n);
    t[n]->SetBranchAddress("photon_pt",photon_pt);
    t[n]->SetBranchAddress("photon_eta",photon_eta);
    t[n]->SetBranchAddress("photon_iso",photon_iso);

    for (int e = 0; e < t[n]->GetEntries(); e++) {

      t[n]->GetEntry( e );

      for (int g = 0; g < photon_n; g++) {

			photon_eta[ g ] += 0.465;

			int eta_bin = -1;
			if ( photon_eta[ g ] > -2.37 && photon_eta[ g ] < -1.56 ) eta_bin = 0;
			if ( photon_eta[ g ] > -1.37 && photon_eta[ g ] < +1.37 ) eta_bin = 1;
			if ( photon_eta[ g ] > +1.56 && photon_eta[ g ] < +2.37 ) eta_bin = 2;
			if (eta_bin == -1) continue;

			if ( photon_pt[ g ] < PTRANGES[ n ] || photon_pt[ g ] > PTRANGES[ n + 1 ] ) continue;

			if ( photon_iso[ g ] > 5 ) continue;

			//eps09(2, 1, 208, x2pdf, Q,
			//      ruv, rdv, ru, rd, rs, rc, rb, rg);

			float nPDF = 0.0;
			if (id2 == 21) nPDF = rg; 
			if ( abs(id2) == 1 ) nPDF = ru;
			if ( abs(id2) == 2 ) nPDF = rd;
			if ( abs(id2) == 3 ) nPDF = rs;
			if ( abs(id2) == 4 ) nPDF = rc;
			if ( abs(id2) == 5 ) nPDF = rb;

			if (nPDF < 0.01) std::cout << " ERROR: " << id2 << std::endl;

			h1_pt_ETA[ eta_bin ]->Fill( photon_pt[ g ], sigma[ n ] );
			h2_pt_xA_ETA[ eta_bin ]->Fill( photon_pt[ g ], x2pdf, sigma[ n ] );
			h2_pt_xp_ETA[ eta_bin ]->Fill( photon_pt[ g ], x1pdf, sigma[ n ] );
			h1_pt_alt_ETA[ eta_bin ]->Fill( photon_pt[ g ], sigma[ n ] * nPDF );

			h1_pt->Fill( photon_pt[ g ], sigma[ n ] );

      }

    }

  }

  TCanvas *tc = new TCanvas();

  TF1 *tf1_fit = new TF1("tf1_fit","[0]*pow(x,[1]+[2]*x)",10,600);
  tf1_fit->SetParameter(0, 1);
  tf1_fit->SetParameter(1, -4);
  tf1_fit->SetParameter(2, 0);
  tf1_fit->SetLineColor(kRed);

  gPad->SetLogy(1);
  gPad->SetLogx(1);

  h1_pt->Fit( tf1_fit, "R" );
  h1_pt->Draw("HIST");
  tf1_fit->Draw("same");

  tc->Print("plot/h1_pt.pdf");

  h1_pt->Divide( tf1_fit );

  h1_pt->Draw("HIST");
  tc->Print("plot/h1_pt_ratio.pdf");


  for (int eta = 0; eta < 3; eta++) {

    {
      std::ostringstream tA; tA << "h1_pt_meanxA_ETA" << eta;
      h1_pt_meanxA_ETA[ eta ] = (TProfile*) h2_pt_xA_ETA[ eta ]->ProfileX( tA.str().c_str(), 1, -1, "" );
      std::ostringstream tp; tp << "h1_pt_meanxp_ETA" << eta;
      h1_pt_meanxp_ETA[ eta ] = (TProfile*) h2_pt_xp_ETA[ eta ]->ProfileX( tp.str().c_str(), 1, -1, "" );
    }

    gPad->SetLogy(1);
    gPad->SetLogx(1);
    
    h1_pt_alt_ETA[ eta ]->SetLineColor( kRed );
    h1_pt_alt_ETA[ eta ]->SetMarkerColor( kRed );

    h1_pt_ETA[ eta ]->Draw();
    h1_pt_alt_ETA[ eta ]->Draw("same");
    {
      std::ostringstream t; t << "plot/h1_pt_ETA" << eta << ".pdf";
      tc->Print( t.str().c_str() );
    }
    
    h1_pt_alt_ETA[ eta ]->Divide( h1_pt_ETA[ eta ] );

    gPad->SetLogy(0);
    h1_pt_alt_ETA[ eta ]->GetYaxis()->SetRangeUser( 0.5, 1.5 );
    h1_pt_alt_ETA[ eta ]->Draw();

    {
      std::ostringstream t; t << "plot/h1_pt_ETA" << eta << "_ratio.pdf";
      tc->Print( t.str().c_str() );
    }

  }

  //

  gPad->SetLogy(1);


  TGraph *tg2A[3];
  tg2A[0] = new TGraph();
  tg2A[1] = new TGraph();
  tg2A[2] = new TGraph();

  for (int n = 0; n < h1_pt_meanxA_ETA[0]->GetNbinsX(); n++) {
    for (int a = 0; a < 3; a++) {
      h1_pt_meanxA_ETA[a]->SetBinError( n+1, 0.01 );
      
      if ( a == 0 && h1_pt_meanxA_ETA[a]->GetBinCenter( n + 1 ) > 300 ) continue;
      if ( a == 1 && h1_pt_meanxA_ETA[a]->GetBinCenter( n + 1 ) > 500 ) continue;
      if ( a == 2 && h1_pt_meanxA_ETA[a]->GetBinCenter( n + 1 ) > 300 ) continue;

      if (h1_pt_meanxA_ETA[a]->GetBinContent(n+1) > 0)
	      tg2A[a]->SetPoint( tg2A[a]->GetN(), h1_pt_meanxA_ETA[a]->GetBinCenter( n + 1 ),  h1_pt_meanxA_ETA[a]->GetBinContent( n + 1 ) );
      
    }
  }

  tg2A[0]->SetLineColor(kRed-6);
  tg2A[1]->SetLineColor(kBlue-6);
  tg2A[2]->SetLineColor(kGreen-6);

  tg2A[0]->SetLineWidth(5);
  tg2A[1]->SetLineWidth(5);
  tg2A[2]->SetLineWidth(5);

  hFrame2->GetYaxis()->SetRangeUser( 1e-3, 1);
  hFrame2->SetTitle(";#it{p}_{T}^{#gamma} [GeV];#LT #it{x}_{Pb} #GT");
  //hFrame2->SetTitle(";#it{p}_{T}^{#gamma} [GeV];#LT #it{x}_{p} #GT");
  hFrame2->Draw();

  tg2A[0]->Draw("C,same");
  tg2A[1]->Draw("C,same");
  tg2A[2]->Draw("C,same");


  if (0) {
    myText(0.25, 0.85, kRed-6, "-2.37 < #it{#eta}^{lab} < -1.56",0.06);
    myText(0.25, 0.78, kBlue-6, "-1.37 < #it{#eta}^{lab} < +1.37",0.06);
    myText(0.25, 0.71, kGreen-6, "+1.56 < #it{#eta}^{lab} < +2.37",0.06);
    
    myText( 0.65, 0.22, kBlack, "#it{p}+Pb, 8.16 TeV"  );
  } else {
    myText(0.60, 0.34, kRed-6, "-2.37 < #it{#eta}^{lab} < -1.56",0.06);
    myText(0.60, 0.27, kBlue-6, "-1.37 < #it{#eta}^{lab} < +1.37",0.06);
    myText(0.60, 0.20, kGreen-6, "+1.56 < #it{#eta}^{lab} < +2.37",0.06);
    
    myText( 0.25, 0.85, kBlack, "#it{p}+Pb, 8.16 TeV"  );
  }

  tc->Print("plot/mean_xA.pdf");
  //tc->Print("plot/mean_xp.pdf");

  TGraph *tg2p[3];
  tg2p[0] = new TGraph();
  tg2p[1] = new TGraph();
  tg2p[2] = new TGraph();

  for (int n = 0; n < h1_pt_meanxp_ETA[0]->GetNbinsX(); n++) {
    for (int a = 0; a < 3; a++) {
      h1_pt_meanxp_ETA[a]->SetBinError( n+1, 0.01 );
      
      if ( a == 0 && h1_pt_meanxp_ETA[a]->GetBinCenter( n + 1 ) > 300 ) continue;
      if ( a == 1 && h1_pt_meanxp_ETA[a]->GetBinCenter( n + 1 ) > 500 ) continue;
      if ( a == 2 && h1_pt_meanxp_ETA[a]->GetBinCenter( n + 1 ) > 300 ) continue;

      if (h1_pt_meanxp_ETA[a]->GetBinContent(n+1) > 0)
	      tg2p[a]->SetPoint( tg2p[a]->GetN(), h1_pt_meanxp_ETA[a]->GetBinCenter( n + 1 ),  h1_pt_meanxp_ETA[a]->GetBinContent( n + 1 ) );
      
    }
  }

  tg2p[0]->SetLineColor(kRed-6);
  tg2p[1]->SetLineColor(kBlue-6);
  tg2p[2]->SetLineColor(kGreen-6);

  tg2p[0]->SetLineWidth(5);
  tg2p[1]->SetLineWidth(5);
  tg2p[2]->SetLineWidth(5);

  hFrame2->GetYaxis()->SetRangeUser( 1e-3, 1);
  hFrame2->SetTitle(";#it{p}_{T}^{#gamma} [GeV];#LT #it{x}_{p} #GT");
  hFrame2->Draw();

  tg2p[0]->Draw("C,same");
  tg2p[1]->Draw("C,same");
  tg2p[2]->Draw("C,same");


  if (0) {
    myText(0.25, 0.85, kRed-6, "-2.37 < #it{#eta}^{lab} < -1.56",0.06);
    myText(0.25, 0.78, kBlue-6, "-1.37 < #it{#eta}^{lab} < +1.37",0.06);
    myText(0.25, 0.71, kGreen-6, "+1.56 < #it{#eta}^{lab} < +2.37",0.06);
    
    myText( 0.65, 0.22, kBlack, "#it{p}+Pb, 8.16 TeV"  );
  } else {
    myText(0.60, 0.34, kRed-6, "-2.37 < #it{#eta}^{lab} < -1.56",0.06);
    myText(0.60, 0.27, kBlue-6, "-1.37 < #it{#eta}^{lab} < +1.37",0.06);
    myText(0.60, 0.20, kGreen-6, "+1.56 < #it{#eta}^{lab} < +2.37",0.06);
    
    myText( 0.25, 0.85, kBlack, "#it{p}+Pb, 8.16 TeV"  );
  }

  //tc->Print("plot/mean_xA.pdf");
  tc->Print("plot/mean_xp.pdf");

  TFile* outf = new TFile("../meanXmaps.root","RECREATE");
  for(int i = 0; i < 3; i++)
  {
    tg2A[i]->SetName(Form("g_meanxA_%d",i));
    tg2p[i]->SetName(Form("g_meanxp_%d",i));
    tg2A[i]->Write();
    tg2p[i]->Write();
  }


  //

  gPad->SetLogy(0);

  TGraph *tg[3];
  tg[0] = new TGraph();
  tg[1] = new TGraph();
  tg[2] = new TGraph();

  for (int n = 0; n < h1_pt_alt_ETA[0]->GetNbinsX(); n++) {
    for (int a = 0; a < 3; a++) {
      h1_pt_alt_ETA[a]->SetBinError( n+1, 0.01 );
      
      if ( a == 0 && h1_pt_alt_ETA[a]->GetBinCenter( n + 1 ) > 300 ) continue;
      if ( a == 1 && h1_pt_alt_ETA[a]->GetBinCenter( n + 1 ) > 500 ) continue;
      if ( a == 2 && h1_pt_alt_ETA[a]->GetBinCenter( n + 1 ) > 300 ) continue;

      if (h1_pt_alt_ETA[a]->GetBinContent(n+1) > 0)
	tg[a]->SetPoint( tg[a]->GetN(), h1_pt_alt_ETA[a]->GetBinCenter( n + 1 ),  h1_pt_alt_ETA[a]->GetBinContent( n + 1 ) );

    }
  }

  tg[0]->SetLineColor(kRed-6);
  tg[1]->SetLineColor(kBlue-6);
  tg[2]->SetLineColor(kGreen-6);

  tg[0]->SetLineWidth(5);
  tg[1]->SetLineWidth(5);
  tg[2]->SetLineWidth(5);

  hFrame->GetYaxis()->SetRangeUser( 0.7, 1.3 );
  hFrame->SetTitle(";#it{p}_{T}^{#gamma} [GeV];#it{R}_{pPb}^{ EPS09}");

  hFrame->Draw();
  {TLine *tl = new TLine(10,1,500,1); tl->SetLineWidth(2); tl->SetLineStyle(2); tl->Draw("same"); }
  tg[0]->Draw("C,same");
  tg[1]->Draw("C,same");
  tg[2]->Draw("C,same");

  myText(0.25, 0.85, kRed-6, "-2.37 < #it{#eta}^{lab} < -1.52",0.06);
  myText(0.25, 0.78, kBlue-6, "-1.37 < #it{#eta}^{lab} < +1.37",0.06);
  myText(0.25, 0.71, kGreen-6, "+1.52 < #it{#eta}^{lab} < +2.37",0.06);

  myText(0.65, 0.78, kBlack, "EPS09 only");
  myText(0.65, 0.71, kBlack, "(no isospin)");

  myText( 0.65, 0.85, kBlack, "#it{p}+Pb, 8.16 TeV"  );


  tc->Print("plot/nPDF.pdf");

  return 0;
}
