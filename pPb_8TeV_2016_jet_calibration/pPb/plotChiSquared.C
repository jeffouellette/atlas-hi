#include <fstream>

void plotChiSquared() {

  ifstream infile;
  infile.open("chi-squared.dat");
  double chi=0;
  TH1F* chidist = new TH1F("chidist", "", 40, 0, 40);
  TH1F* chidistpt14 = new TH1F("chidistpt14", "", 40, 0, 40);
  for (int bin=0; bin<90; bin++) {
    infile >> chi;
    chidist->Fill(chi);
  }
  if (infile) infile.close();
  infile.open("chi-squared-pt14.dat");
  for (int bin=0; bin<90; bin++) { 
    infile >> chi;
    chidistpt14->Fill(chi);
  }
  chidist->SetMarkerColor(kBlack);
  chidist->SetLineColor(kBlack);
  chidist->Draw("e1");
  chidistpt14->SetMarkerColor(kBlue);
  chidistpt14->SetLineColor(kBlue);
  chidistpt14->Draw("same e1");
}
