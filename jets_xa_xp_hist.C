void jets_xa_xp_hist() {

	TFile* output = new TFile("jets_xa_xp.root", "READ");
	TCanvas* c = new TCanvas("c", "", 1000, 800);	

	gPad->SetLogy();

        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(kFALSE);
	
        TH1F* harr[2] = {(TH1F*)output->Get("xp"), (TH1F*)output->Get("xa")};

	const Style_t mkstyles[2] = {kFullCircle, kFullTriangleUp};
        const Color_t mkcolors[2] = {kAzure-7, kOrange-7};

	for (int i = 0; i < sizeof(harr)/sizeof(harr[0]); i++) {
		harr[i]->SetMarkerStyle(mkstyles[i]);
	        harr[i]->SetMarkerColor(mkcolors[i]);
                harr[i]->SetLineColor(mkcolors[i]);
        	harr[i]->Draw("same e1");
	}
        c->Draw();

        gPad->BuildLegend();

        c->SaveAs("./Plots/jets_xa_xp_8.16TeV.pdf");
}
