#include <fstream>

void plotChiSquared() {

  ifstream infile;
  infile.open("chi-squared.dat");
  double chi=0;
  //TH2F* chidist = new TH2F("chidist", "", 40, 0, 40, -90, -4.5, 4.5);
  TGraphErrors* chidist = new TGraphErrors(90);
  for (int bin=0; bin<90; bin++) {
    infile >> chi;
    chidist->SetPoint(bin, -4.45+bin*0.1, chi);
    chidist->SetPointError(bin, 0.05, 0);
  }
  chidist->SetMarkerColor(kBlack);
  chidist->SetLineColor(kBlack);
  chidist->Draw("ab");
}
