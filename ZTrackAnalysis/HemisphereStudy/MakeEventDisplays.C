#include <Utilities.h>
#include <GlobalParams.h>

using namespace atlashi;

void MakeEventDisplays (TString _inFileName, TString _outputDir) {

  TFile* inFile = new TFile (_inFileName, "read");
  TTree* inTree = (TTree*) inFile->Get ("tree");

  TFile* outFile = new TFile (Form ("%s/missingEt.root", _outputDir.Data ()), "recreate");

  int rn, en;
  float zpt, zphi;
  int njet = 0;
  float jpt[10] = {};
  float jphi[10] = {};

  inTree->SetBranchAddress ("run_number", &rn);
  inTree->SetBranchAddress ("event_number", &en);
  inTree->SetBranchAddress ("z_pt", &zpt);
  inTree->SetBranchAddress ("z_phi", &zphi);
  inTree->SetBranchAddress ("njt", &njet);
  inTree->SetBranchAddress ("jtpt", jpt);
  inTree->SetBranchAddress ("jtphi", jphi);

  TCanvas* c1 = new TCanvas ("c1", "", 800, 800);
  c1->Draw ();
  c1->Divide (5, 5);

  const float rcirc = 0.46;
  TEllipse* circ = new TEllipse (0.5, 0.5, rcirc+0.008, rcirc+0.008);
  circ->SetLineWidth (2);
  //circ->Draw ();

  const int nEvt = fmin (25, inTree->GetEntries ());
  //const int nEvt = 48;

  TH1D* h_met = new TH1D ("h_met", ";#sqrt{#left(#Sigma#it{p}_{x}#right)^{2} + #left(#Sigma#it{p}_{y}#right)^{2}};Events", 25, 0, 200);
  for (int iEvt = 0; iEvt < nEvt; iEvt++) {
    inTree->GetEntry (iEvt);
    c1->cd (iEvt+1);

    gPad->Clear ();

    circ->Draw ();

    float maxpt = zpt;
    for (int iJ = 0; iJ < njet; iJ++)
      maxpt = fmax (maxpt, jpt[iJ]);

    float sum_px = 0, sum_py = zpt;
    for (int iJ = 0; iJ < njet; iJ++) {
      sum_px += jpt[iJ] * cos (DeltaPhi (zphi, jphi[iJ], true)+pi/2);
      sum_py += jpt[iJ] * sin (DeltaPhi (zphi, jphi[iJ], true)+pi/2);
    }

    float missing_pt = sqrt (sum_px*sum_px + sum_py*sum_py);
    h_met->Fill (missing_pt); 

    maxpt = fmax (maxpt, missing_pt);

    TArrow* arr;
    arr = new TArrow (0.5, 0.5, 0.5, rcirc * (zpt / maxpt) + 0.5, 0.04, "-|>");
    arr->SetLineColor (kBlack);
    arr->SetFillColor (kBlack);
    arr->SetLineWidth ( (int) (3 * (zpt/maxpt)));
    arr->SetLineWidth ( fmax ((int) (3 * zpt/maxpt), 1));
    arr->SetArrowSize (0.01 * zpt/maxpt);
    arr->Draw ();

    for (int iJ = 0; iJ < njet; iJ++) {
      arr = new TArrow (0.5, 0.5, rcirc * jpt[iJ] * cos (DeltaPhi (zphi, jphi[iJ], true)+pi/2) / maxpt + 0.5, rcirc * jpt[iJ] * sin (DeltaPhi (zphi, jphi[iJ], true)+pi/2) / maxpt + 0.5, 0.04, "-|>");
      arr->SetLineColor (kRed+1);
      arr->SetFillColor (kRed+1);
      arr->SetLineWidth ( fmax ((int) (3 * jpt[iJ]/maxpt), 1));
      arr->SetArrowSize (0.01 * jpt[iJ]/maxpt);
      arr->Draw ();
    }

    arr = new TArrow (0.5, 0.5, rcirc * (-sum_px) / maxpt + 0.5, rcirc * (-sum_py) / maxpt + 0.5, 0.04, "-|>");
    arr->SetLineColor (kGray+1);
    arr->SetFillColor (kGray+1);
    arr->SetLineStyle (9);
    arr->SetLineWidth ( fmax ((int) (3 * missing_pt/maxpt), 1));
    arr->SetArrowSize (0.01 * missing_pt/maxpt);
    arr->Draw ();

    myText (0.03, 0.965, kBlack, Form ("Run %i", rn), 0.03);
    myText (0.03, 0.925, kBlack, Form ("Event %i", en), 0.03);
    myText (0.80, 0.965, kBlack, "Event selection:", 0.03);
    myText (0.80, 0.925, kBlack, "#it{p}_{T}^{h} > 40 GeV", 0.03);
    myText (0.80, 0.885, kBlack, "&& #Delta#phi < #pi/2", 0.03);
    myText (0.03, 0.075, kBlack, Form ("#it{p}_{T}^{Z} = %.1f GeV", zpt), 0.03);
    myText (0.03, 0.035, kBlack, Form ("MET %.1f GeV", missing_pt), 0.03);

    //c1->SaveAs (Form ("%s/event%i.png", _outputDir.Data (), iEvt));
  }
  c1->SaveAs (Form ("%s/canvas.pdf", _outputDir.Data ()));

  inFile->Close ();

  h_met->Write ();
  outFile->Close ();
  

}
