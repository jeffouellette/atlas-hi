#include "triggerUtil.C"

void jets_pratio_hist(std::vector<int> runNumbers) {

    const int numbins = 40;
    const double* xbins = linspace(0, 20, numbins);
    const float eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[8] = {1, 1, 1, 1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap

    const float ymin = 1e-2;
    const float ymax = 6e7;

    const int numhists = 8;
    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%g < #eta < %g (#times %g);#it{p}_{1}/#it{p}_{2} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{1}/#it{p}_{2}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[(i+2)%numhists], harr_scales[i]), numbins, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    double integrated_luminosity = 0;
    for (int runNumber : runNumbers) {
        TFile* thisfile = new TFile(Form("./pratio_data/run_%i.root", runNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1F*)thisfile->Get(Form("%ieta%i", runNumber, j)));
        }
        TVectorD* thisluminosityvec = (TVectorD*)(thisfile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
        integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
    }

    TCanvas* c = new TCanvas("c", "", 1000, 800);   

    gPad->SetLogy();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    
    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[8] = {kAzure-5, kOrange-5, kTeal-5, kPink-5, kSpring-5, kViolet-5, kGray+3, kRed-2};

    for (int i = 0; i < numhists; i+=2) {
        harr[i]->Add(harr[i+1]);
        harr[i]->Scale(0.5);
        harr[i]->SetMarkerStyle(mkstyles[i]);
        harr[i]->SetMarkerColor(mkcolors[i]);
        harr[i]->SetLineColor(mkcolors[i]);
        harr[i]->Draw("same e1");
        harr[i]->SetMinimum(ymin);
        harr[i]->SetMaximum(ymax);
        harr[i]->Draw("same e1");
    }
    c->Draw();

    TLegend* legend = new TLegend(0.6, 0.55, 0.9, 0.9);
    legend->SetHeader("Leading jet pseudorapidities", "C");
    for (int i = 0; i < numhists; i+=2) {
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

    c->SaveAs("./Plots/jets_pratio_8.16TeV.pdf");
}
