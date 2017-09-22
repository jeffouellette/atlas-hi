void jets_xa_xp_hist(std::vector<int> runNumbers) {

    const float eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range

    const float ymin = 5e1;
    const float ymax = 5e8;

    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[8] = {kAzure-5, kOrange-5, kTeal-5, kPink-5, kSpring-5, kViolet-5, kGray+3, kRed-2};

    const float xbins[29] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.6, 0.68, 0.76, 0.84, 0.92, 1.00, 1.08, 1.16, 1.24, 1.32, 1.40, 1.48, 1.56, 1.64};
    const int nbins = sizeof(xbins)/sizeof(xbins[0]) - 1; 
    const int numhists = 16;
    TH1F* harr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        harr[i] = new TH1F(Form("eta%i", i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d#sigma^{2}/d#it{x}_{p} dy #left[pb#right]", eta_cuts[i], eta_cuts[i+2]), nbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        harr[i] = new TH1F(Form("eta%i", i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d#sigma^{2}/d#it{x}_{a} dy #left[pb#right]", eta_cuts[i%(numhists/2)], eta_cuts[(i%(numhists/2))+2]), nbins, xbins);
        harr[i]->Sumw2();
    }

    for (int runNumber : runNumbers) {
        TFile* thisfile = new TFile(Form("./xdata/run_%i.root", runNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1F*)thisfile->Get(Form("%ieta%i", runNumber, j)));
        }
    }

    TLegend* legend = new TLegend(0.56, 0.65, 0.9, 0.9);
    legend->SetHeader("Leading jet pseudorapidity domain selections", "C");
    for (int i = 0; i < numhists/2; i+=2) {
            legend->AddEntry(harr[i], "");
    }
    legend->SetTextSize(0.022);

    TCanvas* c1 = new TCanvas("c1", "", 1000, 800); 
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);    
    for (int i = 0; i < numhists/2; i+=2) {
        harr[i]->Add(harr[i+1]);
        harr[i]->Scale(0.5);
        harr[i]->SetAxisRange(ymin, ymax, "Y");
        harr[i]->SetMarkerStyle(mkstyles[i]);
        harr[i]->SetMarkerColor(mkcolors[i]);
        harr[i]->SetLineColor(mkcolors[i]);
        harr[i]->Draw("same e1");
    }
//        gPad->SetTitle("#it{x}_{p} distribution");
    c1->Draw();
    legend->Draw();
    c1->SaveAs("./Plots/jets_xp_8.16TeV.pdf");


    TCanvas* c2 = new TCanvas("c2", "", 1000, 800);
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    for (int i = numhists/2; i < numhists; i+=2) {
        harr[i]->Add(harr[i+1]);
        harr[i]->Scale(0.5);
        harr[i]->SetAxisRange(ymin, ymax, "Y");
        harr[i]->SetMarkerStyle(mkstyles[i%(numhists/2)]);
        harr[i]->SetMarkerColor(mkcolors[i%(numhists/2)]);
        harr[i]->SetLineColor(mkcolors[i%(numhists/2)]);
        harr[i]->Draw("same e1");
    }
    c2->Draw();
    legend->Draw();
    c2->SaveAs("./Plots/jets_xa_8.16TeV.pdf");
        
}
