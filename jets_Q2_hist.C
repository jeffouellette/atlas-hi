#include "triggerUtil.C"

void jets_Q2_hist(std::vector<int> thisRunNumbers) {

    const double xbins[42] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 500., 600., 800., 1100., 1500., 2000., 2500., 6000.};
    const int numbins = sizeof(xbins)/sizeof(xbins[0]) - 1;
    const double eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[8] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

    const double ymin = 1e-6;
    const double ymax = 2e7;

    const int numhists = 8;
    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%g < #eta < %g (#times %g);#it{Q}^{avg}_{JJ} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{Q}_{JJ}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[i+1], harr_scales[i]), numbins, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    double integrated_luminosity = 0;
    for (int thisRunNumber : thisRunNumbers) {
        TFile* thisfile = new TFile(Form("./Q2_data/run_%i.root", thisRunNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1F*)thisfile->Get(Form("%ieta%i", thisRunNumber, j)));
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

    c->SaveAs("./Plots/jets_Q2_8.16TeV.pdf");
}
