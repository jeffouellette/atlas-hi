void jets_Q2_hist(std::vector<int> runNumbers) {

    const double xbins[17] = {25, 30, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};
    const float eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[8] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

    const float ymin = 1e-2;
    const float ymax = 6e7;
    const float xmin = 1e1;
    const float xmax = 6e3;

    const int numhists = 8;
    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #eta < %1.1f (#times %1.3f);#it{Q}^{dijet}_{12} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{Q}_{12}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[i+1], harr_scales[i]), sizeof(xbins)/sizeof(xbins[0])-1, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    for (int runNumber : runNumbers) {
        TFile* thisfile = new TFile(Form("./Q2_data/run_%i.root", runNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1F*)thisfile->Get(Form("%ieta%i", runNumber, j)));
        }
    }

    TCanvas* c = new TCanvas("c", "", 1000, 800);   
    TLegend* legend = new TLegend(0.6, 0.55, 0.9, 0.9);
    legend->SetHeader("Leading jet pseudorapidities", "C");
    for (int i = 0; i < numhists; i++) {
        legend->AddEntry(harr[i], "");
    }
    legend->SetTextSize(0.022);

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
        harr[draw_order[i]]->GetXaxis()->SetLimits(xmin, xmax);
        harr[draw_order[i]]->SetMinimum(ymin);
        harr[draw_order[i]]->SetMaximum(ymax);
        //harr[draw_order[i]]->SetAxisRange(ymin, ymax, "Y");
        //harr[draw_order[i]]->SetAxisRange(xmin, xmax, "X");
        harr[draw_order[i]]->Draw("same e1");
    }
    c->Draw();

    legend->Draw();

    c->SaveAs("./Plots/jets_Q2_8.16TeV.pdf");
}
