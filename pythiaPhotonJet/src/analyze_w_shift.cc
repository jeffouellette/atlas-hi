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

  TFile *f[7];
  f[0] = new TFile("output/gammajet_pt12_50k.root","READ");
  f[1] = new TFile("output/gammajet_pt26_50k.root","READ");
  f[2] = new TFile("output/gammajet_pt37_50k.root","READ");
  f[3] = new TFile("output/gammajet_pt52_50k.root","READ");
  f[4] = new TFile("output/gammajet_pt105_50k.root","READ");
  f[5] = new TFile("output/gammajet_pt210_50k.root","READ");
  f[6] = new TFile("output/gammajet_pt375_50k.root","READ");

  float sigma[7] = {
    4.262e-04,
    3.832e-05,
    1.166e-05,
    3.533e-06,
    2.519e-07,
    1.391e-08,
    8.773e-10
  };

  TFile *fB[7];
  fB[0] = new TFile("output/gammajetB_pt12_100k.root","READ");
  fB[1] = new TFile("output/gammajetB_pt26_100k.root","READ");
  fB[2] = new TFile("output/gammajetB_pt37_100k.root","READ");
  fB[3] = new TFile("output/gammajetB_pt52_100k.root","READ");
  fB[4] = new TFile("output/gammajetB_pt105_100k.root","READ");
  fB[5] = new TFile("output/gammajetB_pt210_100k.root","READ");
  fB[6] = new TFile("output/gammajetB_pt375_100k.root","READ");

  float sigmaB[7] = {
    4.168e-04,
    3.722e-05,
    1.139e-05,
    3.424e-06,
    2.434e-07,
    1.331e-08,
    8.267e-10
  };


  TTree *t[7];
  for (int n = 0; n < 7; n++) 
    t[n] = (TTree*) fB[n]->Get("tree");
  TTree *tB[7];
  for (int n = 0; n < 7; n++) 
    tB[n] = (TTree*) fB[n]->Get("tree");
  

  double PTRANGES[ 8 ] = { 15.85, 35, 50, 70, 140, 280, 500, 1000 };

  const int NPTBINS = 20;
  double PTBINS[21];
  double alpha = pow( 10.0, 1 / 10.0 );
  for (int n = 0; n < 21; n++) {
    PTBINS[ n ] = 10 * pow( alpha, n );
    //std::cout << PTBINS[ n ] << std::endl;
  }

  const int NPTBINS2 = 170;
  double PTBINS2[ NPTBINS2+1 ];
  double alpha2 = pow( 10.0, 1 / 100.0 );
  for (int n = 0; n < NPTBINS2+1; n++) {
    PTBINS2[ n ] = 10 * pow( alpha2, n );
    std::cout << n << " " << PTBINS2[ n ] << std::endl;
  }

  TH1D *h1_pt_ETA[4];
  h1_pt_ETA[0] = new TH1D("h1_pt_ETA0","", NPTBINS, PTBINS );
  h1_pt_ETA[1] = new TH1D("h1_pt_ETA1","", NPTBINS, PTBINS );
  h1_pt_ETA[2] = new TH1D("h1_pt_ETA2","", NPTBINS, PTBINS );
  h1_pt_ETA[3] = new TH1D("h1_pt_ETA3","", NPTBINS, PTBINS );

  TH1D *h1B_pt_ETA[4];
  h1B_pt_ETA[0] = new TH1D("h1B_pt_ETA0","", NPTBINS, PTBINS );
  h1B_pt_ETA[1] = new TH1D("h1B_pt_ETA1","", NPTBINS, PTBINS );
  h1B_pt_ETA[2] = new TH1D("h1B_pt_ETA2","", NPTBINS, PTBINS );
  h1B_pt_ETA[3] = new TH1D("h1B_pt_ETA3","", NPTBINS, PTBINS );


  TH1D *h1_pt = new TH1D("h1_pt","",NPTBINS2,PTBINS2);
  h1_pt->Sumw2();

  TH1D *h1B_pt = new TH1D("h1B_pt","",NPTBINS2,PTBINS2);
  h1B_pt->Sumw2();

  int photon_n;
  float photon_pt[50];
  float photon_eta[50];
  float photon_iso[50];

  for (int n = 0; n < 7; n++) {

    // beam-1 is the nucleus
    // so boost is in the direction of beam-2...

    t[n]->SetBranchAddress("photon_n",&photon_n);
    t[n]->SetBranchAddress("photon_pt",photon_pt);
    t[n]->SetBranchAddress("photon_eta",photon_eta);
    t[n]->SetBranchAddress("photon_iso",photon_iso);

    for (int e = 0; e < t[n]->GetEntries(); e++) {

      t[n]->GetEntry( e );

      for (int g = 0; g < photon_n; g++) {

	photon_eta[ g ] -= 0.465;

	int eta_bin = -1;
	if ( photon_eta[ g ] > -2.37 && photon_eta[ g ] < -1.52 ) eta_bin = 0;
	if ( photon_eta[ g ] > -1.37 && photon_eta[ g ] < 0 ) eta_bin = 1;
	if ( photon_eta[ g ] > 0 && photon_eta[ g ] < +1.37 ) eta_bin = 2;
	if ( photon_eta[ g ] > +1.52 && photon_eta[ g ] < +2.37 ) eta_bin = 3;
	
	if (eta_bin == -1) continue;

	if ( photon_pt[ g ] < PTRANGES[ n ] || photon_pt[ g ] > PTRANGES[ n + 1 ] ) continue;

	if ( photon_iso[ g ] > 5 ) continue;

	h1_pt_ETA[ eta_bin ]->Fill( photon_pt[ g ], sigma[ n ] );

	h1_pt->Fill( photon_pt[ g ], sigma[ n ] );

      }

    }

  }


  for (int n = 0; n < 7; n++) {

    // beam-1 is the nucleus
    // so boost is in the direction of beam-2...

    tB[n]->SetBranchAddress("photon_n",&photon_n);
    tB[n]->SetBranchAddress("photon_pt",photon_pt);
    tB[n]->SetBranchAddress("photon_eta",photon_eta);
    tB[n]->SetBranchAddress("photon_iso",photon_iso);

    for (int e = 0; e < tB[n]->GetEntries(); e++) {

      tB[n]->GetEntry( e );

      for (int g = 0; g < photon_n; g++) {

	//photon_eta[ g ] += 0.465;

	int eta_bin = -1;
	if ( photon_eta[ g ] > -2.37 && photon_eta[ g ] < -1.52 ) eta_bin = 0;
	if ( photon_eta[ g ] > -1.37 && photon_eta[ g ] < 0 ) eta_bin = 1;
	if ( photon_eta[ g ] > 0 && photon_eta[ g ] < +1.37 ) eta_bin = 2;
	if ( photon_eta[ g ] > +1.52 && photon_eta[ g ] < +2.37 ) eta_bin = 3;

	
	if (eta_bin == -1) continue;

	if ( photon_pt[ g ] < PTRANGES[ n ] || photon_pt[ g ] > PTRANGES[ n + 1 ] ) continue;

	if ( photon_iso[ g ] > 5 ) continue;

	h1B_pt_ETA[ eta_bin ]->Fill( photon_pt[ g ], sigma[ n ] );

	h1B_pt->Fill( photon_pt[ g ], sigma[ n ] );

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

  //

  TF1 *tf1B_fit = new TF1("tf1B_fit","[0]*pow(x,[1]+[2]*x)",10,600);
  tf1B_fit->SetParameter(0, 1);
  tf1B_fit->SetParameter(1, -4);
  tf1B_fit->SetParameter(2, 0);
  tf1B_fit->SetLineColor(kRed);

  gPad->SetLogy(1);
  gPad->SetLogx(1);

  h1B_pt->Fit( tf1B_fit, "R" );
  h1B_pt->Draw("HIST");
  tf1_fit->Draw("same");

  tc->Print("plot/h1B_pt.pdf");

  h1B_pt->Divide( tf1B_fit );

  h1B_pt->Draw("HIST");
  tc->Print("plot/h1B_pt_ratio.pdf");

  //

  for (int eta = 0; eta < 4; eta++) {

    gPad->SetLogy(1);
    gPad->SetLogx(1);
    
    h1B_pt_ETA[ eta ]->SetLineColor( kRed );
    h1B_pt_ETA[ eta ]->SetMarkerColor( kRed );

    h1_pt_ETA[ eta ]->Draw();
    h1B_pt_ETA[ eta ]->Draw("same");
    {
      std::ostringstream t; t << "plot/h1_shift_pt_ETA" << eta << ".pdf";
      tc->Print( t.str().c_str() );
    }
    
    h1_pt_ETA[ eta ]->Divide( h1B_pt_ETA[ eta ] );

    gPad->SetLogy(0);
    h1_pt_ETA[ eta ]->GetYaxis()->SetRangeUser( 0.0, 2.0 );
    h1_pt_ETA[ eta ]->Draw();

    {
      std::ostringstream t; t << "plot/h1_shift_pt_ETA" << eta << "_ratio.pdf";
      tc->Print( t.str().c_str() );
    }

  }

  //

  gPad->SetLogy(0);

  for (int n = 0; n < h1B_pt_ETA[0]->GetNbinsX(); n++) {
    //for (int a = 0; a < 3; a++)
      //h1B_pt_ETA[a]->SetBinError( n+1, 0.01 );
  }

  h1_pt_ETA[0]->SetMarkerColor( kBlack );
  h1_pt_ETA[0]->SetLineColor( kBlack );

  h1_pt_ETA[1]->SetMarkerColor( kRed );
  h1_pt_ETA[1]->SetLineColor( kRed );

  h1_pt_ETA[2]->SetMarkerColor( kBlue );
  h1_pt_ETA[2]->SetLineColor( kBlue );

  h1_pt_ETA[3]->SetMarkerColor( kGreen+2 );
  h1_pt_ETA[3]->SetLineColor( kGreen+2 );

  h1_pt_ETA[0]->GetYaxis()->SetRangeUser( 0.0, 2.0 );
  h1_pt_ETA[0]->SetTitle(";#it{p}_{T}^{#gamma} [GeV];''#it{R}_{pPb}'' = #it{pp}(#Delta#it{y}=-0.465) / #it{pp}(#Delta#it{y}=0)");
  h1_pt_ETA[0]->Draw("");
  {TLine *tl = new TLine(10,1,1000,1); tl->SetLineWidth(2); tl->SetLineStyle(2); tl->Draw("same"); }
  h1_pt_ETA[1]->Draw("same");
  h1_pt_ETA[2]->Draw("same");
  h1_pt_ETA[3]->Draw("same");

  myText(0.2, 0.90, kBlack, "-2.37 < #it{#eta}^{lab} < -1.52");
  myText(0.2, 0.83, kRed, "-1.37 < #it{#eta}^{lab} < 0");
  myText(0.2, 0.76, kBlue, "0 < #it{#eta}^{lab} < +1.37");
  myText(0.2, 0.69, kGreen+2, "+1.52 < #it{#eta}^{lab} < +2.37");

  tc->Print("plot/h1_shift.pdf");

  return 0;
}
