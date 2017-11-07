#include "TMathText.h"

void jets_mjj_hist(std::vector<int> thisRunNumbers) {

    const int numhists = 5;

    const float ymin = 1e-6;
    const float ymax = 2e7;
    const float xmin = 1e1;
    const float xmax = 6e3;

    const double xbins[42] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 500., 600., 800., 1100., 1500., 2000., 2500., 6000.};
    const int numbins = sizeof(xbins)/sizeof(xbins[0]) - 1;
    const int len_xbins = sizeof(xbins)/sizeof(xbins[0]);
    const float etastarcuts[numhists+1] = {0, 0.5, 1, 1.5, 2, 3};

    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%g #leq #left|#eta*#right| #leq %g; #it{M}_{JJ} #left[GeV/#it{c}^{2}#right];d^{2}#sigma/d#it{M}_{JJ}d#left|#eta*#right| #left[pb (GeV/#it{c}^{#it{2}})^{-1}#right]", etastarcuts[i], etastarcuts[i+1]), numbins, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    double integrated_luminosity = 0;
    for (int thisRunNumber : thisRunNumbers) {
        TFile* thisfile = new TFile(Form("./mjj_data/run_%i.root", thisRunNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1D*)thisfile->Get(Form("%ieta%i", thisRunNumber, j)));
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

    for (int i = 0; i < sizeof(harr)/sizeof(harr[0]); i++) {
        harr[i]->SetMarkerStyle(mkstyles[i]);
        harr[i]->SetMarkerColor(mkcolors[i]);
        harr[i]->SetLineColor(mkcolors[i]);
        harr[i]->Draw("same e1");
        harr[i]->GetXaxis()->SetLimits(xmin, xmax);
        harr[i]->SetMinimum(ymin);
        harr[i]->SetMaximum(ymax);
        //harr[i]->SetAxisRange(ymin, ymax, "Y");
        //harr[i]->SetAxisRange(xmin, xmax, "X");
        harr[i]->Draw("same e1");
    }
    c->Draw();

    TLegend* legend = new TLegend(0.57, 0.65, 0.9, 0.9);
    legend->SetHeader("Dijet center-of-mass pseudorapidities", "C");
    for (int i = 0; i < numhists; i++) {
        legend->AddEntry(harr[i], "");
    }
    legend->SetTextFont(42);
    legend->SetTextSize(0.024);
    legend->Draw();

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} #it{p-Pb}");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.61, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.52, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", integrated_luminosity*1000)); 
 
/*    TLatex* eta_star_text = new TLatex();
    eta_star_text->SetTextAlign(22);
    eta_star_text->SetTextFont(42);  
    eta_star_text->SetTextSize(0.036);
    eta_star_text->DrawLatexNDC(0.78, 0.55, "#eta* #equiv #frac{#eta_{1} #minus #eta_{2}}{2}");
*/


    c->SaveAs("./Plots/jets_mjj_8.16TeV.pdf");
}
