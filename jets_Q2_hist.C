void jets_Q2_hist(std::vector<int> runNumbers) {

	const float xbins_n200eta490[13] = {40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};              // Eta range: -4.9 < eta <= -2
        const float xbins_0eta200[13] = {40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};                 // Eta range: -2 < eta < 2
        const float xbins_p200eta320[15] = {40, 50, 55, 60, 70, 75, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};	// Eta range: 2 <= eta < 3.2
        const float xbins_p320eta490[15] = {25, 40, 50, 60, 65, 70, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};      // Eta range: 3.2 <= eta < 4.9

        int numHists = 8;
        TH1F* harr[numHists];
	harr[0] = new TH1F("eta0", "-4.9 < #eta < -3.2 (#times 0.005);|#it{Q}|| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_n200eta490)/sizeof(xbins_n200eta490[0])-1, xbins_n200eta490);
	harr[1] = new TH1F("eta1", "-3.2 < #eta < -2 (#times 0.03);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_n200eta490)/sizeof(xbins_n200eta490[0])-1, xbins_n200eta490);
	harr[2] = new TH1F("eta2", "-2 < #eta < -1 (#times 0.1);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[3] = new TH1F("eta3", "-1 < #eta < 0 (#times 0.5);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[4] = new TH1F("eta4", "0 < #eta < 1 (#times 1);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[5] = new TH1F("eta5", "1 < #eta < 2 (#times 0.3);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[6] = new TH1F("eta6", "2 < #eta < 3.2 (#times 0.05);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_p200eta320)/sizeof(xbins_p200eta320[0])-1, xbins_p200eta320);
	harr[7] = new TH1F("eta7", "3.2 < #eta < 4.9 (#times 0.01);|#it{Q}| [GeV];d#sigma/d|#it{Q}| dy [pb/GeV]", sizeof(xbins_p320eta490)/sizeof(xbins_p320eta490[0])-1, xbins_p320eta490);
	for (int i = 0; i< sizeof(harr)/sizeof(harr[0]); i++) {
		harr[i]->Sumw2(); // instruct each histogram to propagate errors
	}

        for (int runNumber : runNumbers) {
                TFile* thisfile = new TFile(Form("./Q2_data/run_%i.root", runNumber), "READ");
                for (int j = 0; j < numHists; j++) {
                        harr[j]->Add((TH1F*)thisfile->Get(Form("%ieta%i", runNumber, j)));
                }
        }

	TCanvas* c = new TCanvas("c", "", 1000, 800);	

	gPad->SetLogy();
        gPad->SetLogx();

        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(kFALSE);
	
	const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
        const Color_t mkcolors[8] = {kAzure-5, kOrange-5, kTeal-5, kPink-5, kSpring-5, kViolet-5, kGray+3, kRed-2};

        const int draw_order[8] = {4, 3, 5, 2, 6, 1, 7, 0};

	for (int i = 0; i < sizeof(harr)/sizeof(harr[0]); i++) {
		harr[draw_order[i]]->SetMarkerStyle(mkstyles[draw_order[i]]);
	        harr[draw_order[i]]->SetMarkerColor(mkcolors[draw_order[i]]);
                harr[draw_order[i]]->SetLineColor(mkcolors[draw_order[i]]);
        	harr[draw_order[i]]->Draw("same e1");
	}
        c->Draw();

        gPad->BuildLegend();

        c->SaveAs("./Plots/jets_Q2_8.16TeV.pdf");
}
