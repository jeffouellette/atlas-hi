

void MakeZRunPlot () {

  double runs[45] = {365502, 365512, 365573, 365602, 365627, 365678, 365681, 365709, 365752, 365834, 365914, 365932, 366011, 366029, 366092, 366142, 366268, 366337, 366383, 366413, 366476, 366526, 366528, 366627, 366691, 366754, 366805, 366860, 366878, 366919, 366931, 366994, 367023, 367099, 367134, 367165, 367170, 367233, 367273, 367318, 367321, 367363, 367364, 367365, 367384};
  double numZs[45] = {6, 12, 120, 153, 302, 229, 394, 574, 557, 91, 101, 576, 663, 637, 632, 615, 801, 639, 1, 637, 635, 16, 472, 147, 494, 490, 503, 540, 394, 600, 482, 129, 557, 641, 734, 125, 684, 624, 520, 365, 406, 77, 420, 220, 41};
  double intZs[45] = {};
  double numMBs[45] = {142597, 195764, 157790, 188214, 201742, 100591, 279175, 276244, 878256, 236829, 380645, 411455, 408188, 375611, 411941, 414022, 655447, 571949, 7025, 391784, 410285, 297, 262080, 186607, 422185, 357036, 449550, 385276, 329717, 447512, 213376, 42223, 408918, 388807, 379964, 11356, 566336, 402547, 417474, 457606, 433393, 22005, 462025, 193131, 30043};
  double intMBs[45] = {};

  for (int i = 0; i < 45; i++) {
    for (int j = 0; j <= i; j++) {
      intZs[i] += numZs[j];
      intMBs[i] += numMBs[j];
    }
  }

  TGraph* g_numZs = new TGraph (45, runs, numZs);
  TGraph* g_numMBs = new TGraph (45, runs, numMBs);
  TGraph* g_intZs = new TGraph (45, runs, intZs);
  TGraph* g_intMBs = new TGraph (45, runs, intMBs);

  g_numZs->GetYaxis()->SetRangeUser (1e-2, 1e8);
  g_numZs->GetYaxis ()->SetTitle ("N_{evt}");
  g_numZs->GetXaxis ()->SetTitle ("Run number");
  g_numZs->Draw ("AP");
  g_intZs->SetMarkerStyle (kDot);
  g_intZs->Draw ("LP");
 
  g_numMBs->SetMarkerColor (kRed+1);
  g_numMBs->SetLineColor (kRed+1);
  g_numMBs->SetMarkerStyle (kOpenCircle); 
  g_numMBs->Draw ("P");

  //g_intMBs->SetMarkerColor (kRed+1);
  g_intMBs->SetLineColor (kRed+1);
  g_intMBs->SetMarkerStyle (kDot); 
  g_intMBs->Draw ("LP");

  myText (0.21, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.21, 0.84, kBlack, "2018 Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV");

  myText (0.42, 0.28, kBlack, "Z-tagged events", 0.036);
  myText (0.42, 0.22, kRed+1, "Mixing pool", 0.036);

  myText (0.19, 0.34, kBlack, "Per-run", 0.032);
  myText (0.28, 0.34, kBlack, "Cumulative sum", 0.032);
  myMarkerTextNoLine (0.24, 0.28, kBlack, kFullCircle, "", 1.5, 0.038);
  myMarkerTextNoLine (0.24, 0.22, kRed+1, kOpenCircle, "", 1.5, 0.038);
  myMarkerText (0.37, 0.28, kBlack, kDot, "", 1.5, 0.038);
  myMarkerText (0.37, 0.22, kRed+1, kDot, "", 1.5, 0.038);

  gPad->SetLogy ();

}
