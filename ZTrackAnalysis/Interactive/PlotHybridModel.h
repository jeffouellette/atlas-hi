#ifndef __PlotHybridModel_h__
#define __PlotHybridModel_h__

#include "Params.h"

void PlotHybridModel (const bool useTrkPt = true) {
  const int axisTextSize = 35;

  const short iCent = 3;

  //const string modelFileName = (useTrkPt ? "../HybridModel/IAAs/010_IAA_pt_wake_0_ignore_neg_0.dat" : "../HybridModel/IAAs/010_IAA_z_wake_0_ignore_neg_0.dat"); // means no medium response.data"
  const string modelFileName = (useTrkPt ? "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_0.dat" : "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_0.dat"); // means medium response including only the positive contribution from the wake
  //const string modelFileName = (useTrkPt ? "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_1.dat" : "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_1.dat"); // means full medium response, including also the negative contribution from the wake

  TGAE** g_hybridModel = Get1DArray <TGAE*> (nPtZBins);
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++)
    g_hybridModel[iPtZ] = new TGAE ();

  {
    ifstream f;
    f.open (modelFileName.c_str ());
    float dummy = 0, x = 0, y1 = 0, y2 = 0, y = 0, yerr = 0;
    while (f) {
      f >> x;
      for (int iPtZ = 2; iPtZ <= 4; iPtZ++) {
        f >> dummy >> dummy >> dummy >> dummy >> y1 >> y2;
        if (useTrkPt && x > zPtBins[iPtZ])
          continue;
        else if (!useTrkPt && x*zPtBins[iPtZ] < 1)
          continue;
        y = 0.5*(y1+y2);
        yerr = 0.5*(y1-y2); 
        g_hybridModel[iPtZ]->SetPoint (g_hybridModel[iPtZ]->GetN (), x, y);
        g_hybridModel[iPtZ]->SetPointEYhigh (g_hybridModel[iPtZ]->GetN ()-1, yerr);
        g_hybridModel[iPtZ]->SetPointEYlow (g_hybridModel[iPtZ]->GetN ()-1, yerr);
      }
    }
    f.close ();
  }

  //const char* canvasName = Form ("c_z_trk_zpt_%s_iaa_iCent%i", useTrkPt ? "pttrk" : "xhz", iCent);
  //const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  //TCanvas* c = nullptr;
  //if (canvasExists)
  //  c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  //else {
  //  c = new TCanvas (canvasName, "", 800, 800);
  //  gDirectory->Add (c);
  //}

  double xmin = 0, xmax = 0;
  //gPad->SetLogx ();
  //gPad->SetLogy ();

  for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    TGAE* g = g_hybridModel[iPtZ];

    g->GetYaxis ()->SetRangeUser (0.1, max_iaa);

    g->SetFillColorAlpha (modelFillColors[iPtZ-2], 0.3);
    //g->SetFillStyle (1001);

    //g->Draw (!canvasExists && iPtZ == 3 ? "A3" : "3");
    g->Draw ("3");
  }

  for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++)
    myOnlyBoxText (0.603, 0.868-0.05*(iPtZ-3), 1.2, modelFillColors[iPtZ-2], kBlack, 1, "", 0.040, 1001, 0.3);

  myText (0.54, 0.90, kBlack, "Hybrid", 0.030);
  myText (0.63, 0.90, kBlack, "Data", 0.030);


}

#endif
