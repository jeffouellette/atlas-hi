#include "../triggerUtil.C"

void jets_xa_xp_hist(std::vector<int> thisRunNumbers) {

    const int numbins = 40;
    const int numhists = 16;
    const double eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double d_eta[8] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};

    const double ymin = 1e0;
    const double ymax = 3e10;

    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[8] = {kAzure-5, kTeal-5, kOrange-5, kPink-5, kSpring-5, kViolet-5, kRed-2, kGray+3};

    const double* xbins = logspace(0, 1.6, numbins);
    TH1D* harr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d#sigma^{2}/d#it{x}_{p} dy #left[pb#right]", eta_cuts[i], eta_cuts[i+2]), numbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d#sigma^{2}/d#it{x}_{a} dy #left[pb#right]", eta_cuts[i%(numhists/2)], eta_cuts[(i%(numhists/2))+2]), numbins, xbins);
        harr[i]->Sumw2();
    }

    double integrated_luminosity = 0;
    for (int thisRunNumber : thisRunNumbers) {
        TFile* thisfile = new TFile(Form("./xdata/run_%i.root", thisRunNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1D*)thisfile->Get(Form("%ieta%i", thisRunNumber, j)));
        }
        TVectorD* thisluminosityvec = (TVectorD*)(thisfile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
        integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
    }

    TLegend* legend = new TLegend(0.60, 0.67, 0.9, 0.9);
    legend->SetHeader("Leading jet pseudorapidities", "C");
    for (int i = 0; i < numhists/2; i+=2) {
            legend->AddEntry(harr[i], "");
    }
    legend->SetTextSize(0.024);

    TCanvas* c1 = new TCanvas("c1", "", 1000, 800); 
    gPad->SetLogy();
    gPad->SetLogx();
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
    c1->Draw();
    legend->Draw();

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} #it{p-Pb}");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.51, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.42, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", integrated_luminosity*1000)); 
    
    c1->SaveAs("./Plots/jets_xp_8.16TeV.pdf");


    TCanvas* c2 = new TCanvas("c2", "", 1000, 800);
    gPad->SetLogy();
    gPad->SetLogx();
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

    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} #it{p-Pb}");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.51, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.42, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", integrated_luminosity*1000)); 
    
    c2->SaveAs("./Plots/jets_xa_8.16TeV.pdf");
        
}
