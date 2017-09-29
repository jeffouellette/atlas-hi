void jets_pt_hist(std::vector<int> runNumbers) {
    const int numhists = 8;

    const double xbins_n200eta490[16] = {25, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};           // Eta range: -4.9 < eta <= -2
    const double xbins_0eta200[16] = {25, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};              // Eta range: -2 < eta < 2
    const double xbins_p200eta320[18] = {25, 40, 50, 55, 60, 70, 75, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};   // Eta range: 2 <= eta < 3.2
    const double xbins_p320eta490[17] = {25, 40, 50, 60, 65, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};   // Eta range: 3.2 <= eta < 4.9
    const double* xbins[numhists] = {xbins_n200eta490, xbins_n200eta490, xbins_0eta200, xbins_0eta200, xbins_0eta200, xbins_0eta200, xbins_p200eta320, xbins_p320eta490};
    const int len_xbins[numhists] = {16, 16, 16, 16, 16, 16, 18, 17};

    const double ymin = 1e-2;
    const double ymax = 6e7;
    const double xmin = 1e1;
    const double xmax = 6e3;

    const double eta_cuts[numhists+1] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[numhists] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%g < #eta < %g (#times %g); #it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[i+1], harr_scales[i]), len_xbins[i]-1, xbins[i]);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    double integrated_luminosity = 0;
    for (int runNumber : runNumbers) {
        TFile* thisfile = new TFile(Form("./pt_data/run_%i.root", runNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1D*)thisfile->Get(Form("%ieta%i", runNumber, j)));
        }
        TVectorD* thisluminosityvec = (TVectorD*)(thisfile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
        integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
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
        harr[draw_order[i]]->GetXaxis()->SetLimits(xmin, xmax);
        harr[draw_order[i]]->SetMinimum(ymin);
        harr[draw_order[i]]->SetMaximum(ymax);
        harr[draw_order[i]]->Draw("same e1");
    }
    c->Draw();

    TLegend* legend = new TLegend(0.6, 0.55, 0.9, 0.9);
    legend->SetHeader("Leading jet pseudorapidities", "C");
    for (int i = 0; i < numhists; i++) {
        legend->AddEntry(harr[i], "");
    }
    legend->SetTextSize(0.024);
    legend->Draw();        

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} #it{p-Pb}");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.51, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.42, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", integrated_luminosity*1000)); 


    c->SaveAs("./Plots/jets_pt_8.16TeV.pdf");
}
