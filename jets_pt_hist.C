#include "triggerUtil.C"

void jets_pt_hist(std::vector<int> thisRunNumbers) {
    const int numhists = 8;

    const double harr_scales[numhists] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

    const double ymin = 1e-8;
    const double ymax = 5e5;

    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%g < #eta < %g (#times %g); #it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", etabins[i], etabins[i+1], harr_scales[i]), numpbins, pbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    double integrated_luminosity = 0;
    std::vector<double> luminosity(thisRunNumbers.size(), 0);
    TH1D* thishist;
    TFile* thisfile;
    for (int thisRunNumber : thisRunNumbers) {
        
        thisfile = new TFile(Form("./pt_data/run_%i.root", thisRunNumber), "READ");
        TVectorD* thisluminosityvec = (TVectorD*)(thisfile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
        integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
        
        for (int j = 0; j < numhists; j++) {
            thishist = (TH1D*)thisfile->Get(Form("%ieta%i", thisRunNumber, j));
            for (int pbin = 0; pbin < numpbins; pbin++) {
                if (thishist->GetBinContent(pbin+1) > 0) luminosity[pbin] += (*thisluminosityvec)[0];
            }
            harr[j]->Add(thishist);
        }
        thisfile->Close();
    }

    TCanvas* c = new TCanvas("c", "", 1000, 800);   
    
    gPad->SetLogy();
    gPad->SetLogx();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
 
    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[8] = {kAzure-5, kOrange-5, kTeal-5, kPink-5, kSpring-5, kViolet-5, kGray+3, kRed-2};

    const int draw_order[8] = {4, 3, 5, 2, 6, 1, 7, 0};

    for (int i = 0; i < numhists; i++) {
        thishist = harr[draw_order[i]];
        for (int pbin = 0; pbin < numpbins; pbin++) {
            if (luminosity[pbin] != 0) thishist->SetBinContent(pbin+1, thishist->GetBinContent(pbin+1) / luminosity[pbin]);
        }
        thishist->SetMarkerStyle(mkstyles[draw_order[i]]);
        thishist->SetMarkerColor(mkcolors[draw_order[i]]);
        thishist->SetLineColor(mkcolors[draw_order[i]]);
        thishist->Draw("same e1");
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->Draw("same e1");
    }
    c->Draw();

    TLegend* legend = new TLegend(0.6, 0.55, 0.9, 0.9);
    legend->SetHeader("Jet pseudorapidities", "C");
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
