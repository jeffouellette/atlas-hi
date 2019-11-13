#ifndef __PlotJewel_C__
#define __PlotJewel_C__

void PlotJewel () {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/Jewel/hists.root", rootPath.Data ()), "read");


}


#endif
