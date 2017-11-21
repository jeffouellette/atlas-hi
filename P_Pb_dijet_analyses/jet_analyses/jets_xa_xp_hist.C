#include "../triggerUtil.C"

void jets_xa_xp_hist() {
    if (!runPeriodA && !runPeriodB) return;
    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numbins = 100;
    const int numhists = 2*numetabins;

    const double ymin = 5e-1;
    const double ymax = 3e12;

    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[8] = {kAzure-5, kTeal-5, kOrange-5, kPink-5, kSpring-5, kViolet-5, kRed-2, kGray+3};

    const double* xbins = logspace(8e-5, 1.6, numbins);
    TH1D* harr[numhists];
    for (int i = 0; i < numetabins; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #eta_{lab} < %1.1f;#it{x}_{p};d^{2}N/#it{L}_{int}d#it{x}_{p} dy #left[pb#right]", etabins[i], etabins[i+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numetabins; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #eta_{lab} < %1.1f;#it{x}_{a};d^{2}N/#it{L}_{int}d#it{x}_{a} dy #left[pb#right]", etabins[i%(numetabins)], etabins[(i%(numetabins))+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    TH2D* xaxpcorr = new TH2D("xaxpcorr", ";#it{x}_{a};#it{x}_{p};d^{2}N/#it{L}_{int}d#it{x}_{p}d#it{x}_{a}", numbins, xbins, numbins, xbins);
    TH2D* fcalhist = new TH2D("fcalhist", ";#it{x}_{p};FCAL energy deposited;", numbins, xbins, numpbins, pbins);

    double integrated_luminosity = 0;
    for (int thisRunNumber : (*thisRunNumbers)) {
        if (skipRun(thisRunNumber)) continue;
        TFile* thisfile = new TFile(Form("%sxdata/run_%i.root", rootPath.c_str(), thisRunNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            harr[j]->Add((TH1D*)thisfile->Get(Form("%ieta%i", thisRunNumber, j)));
        }
        xaxpcorr->Add((TH2D*)thisfile->Get(Form("xaxpcorr_run%i", thisRunNumber)));
        fcalhist->Add((TH2D*)thisfile->Get(Form("fcalhist_run%i", thisRunNumber)));
        TVectorD* thisluminosityvec = (TVectorD*)(thisfile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
        integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
    }

    TLegend* legend = new TLegend(0.77, 0.65, 0.95, 0.95);
    legend->SetBorderSize(0);
    legend->SetHeader("Leading jet #eta (lab frame)", "C");
    for (int i = 0; i < numetabins; i++) {
            legend->AddEntry(harr[i], "");
    }
//    legend->SetTextSize(0.019);

    TCanvas* c1 = new TCanvas("c1", "", 1000, 800);
//    c1->SetMargin(0.1,0.05,0.07,0.05);
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);    
    for (int i = 0; i < numetabins; i++) {
        harr[i]->SetAxisRange(ymin, ymax, "Y");
        harr[i]->SetMarkerStyle(mkstyles[i]);
        harr[i]->SetMarkerColor(mkcolors[i]);
        harr[i]->SetLineColor(mkcolors[i]);
        harr[i]->Draw("same e1");
        harr[i]->GetXaxis()->SetTickLength(0.02);
        harr[i]->GetYaxis()->SetTickLength(0.02);
        
//        harr[i]->GetXaxis()->SetTitleOffset(0.7);
//        harr[i]->GetXaxis()->SetLabelOffset(0.0025);
    }
    legend->Draw();

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
//    description->SetTextSize(0.036);
//    description->DrawLatexNDC(0.22, 0.9, "#bf{#it{ATLAS}} p-Pb");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.67, 0.9, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.67, 0.82, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", integrated_luminosity*1000)); 


    string hname = "jets_xp_8.16TeV";
    if (runPeriodA && runPeriodB) { 
        description->DrawLatexNDC(0.42, 0.91, "#it{x}_{p} Period A(-#eta) & B(#eta)");
    } else if (runPeriodA && !runPeriodB) {
        description->DrawLatexNDC(0.44, 0.91, "#it{x}_{p} Period A");
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        description->DrawLatexNDC(0.44, 0.91, "#it{x}_{p} Period B");
        hname = hname + "_periodB";
    }
    c1->Draw();
    c1->SaveAs((plotPath + hname + ".pdf").c_str());


    TCanvas* c2 = new TCanvas("c2", "", 1000, 800);
//    c2->SetMargin(0.1,0.05,0.07,0.05);
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    for (int i = numetabins; i < numhists; i++) {
        harr[i]->SetAxisRange(ymin, ymax, "Y");
        harr[i]->SetMarkerStyle(mkstyles[i%(numetabins)]);
        harr[i]->SetMarkerColor(mkcolors[i%(numetabins)]);
        harr[i]->SetLineColor(mkcolors[i%(numetabins)]);
        harr[i]->Draw("same e1");
        harr[i]->GetXaxis()->SetTickLength(0.02);
        harr[i]->GetYaxis()->SetTickLength(0.02);

//        harr[i]->GetXaxis()->SetTitleOffset(0.7);
//        harr[i]->GetXaxis()->SetLabelOffset(0.0025);
    }
    legend->Draw();

//    description->SetTextSize(0.036);
//    description->DrawLatexNDC(0.22, 0.9, "#bf{#it{ATLAS}} p-Pb");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.67, 0.9, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.67, 0.82, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", integrated_luminosity*1000)); 

    hname = "jets_xa_8.16TeV";
    if (runPeriodA && runPeriodB) { 
        description->DrawLatexNDC(0.42, 0.91, "#it{x}_{a} Period A(-#eta) & B(#eta)");
    } else if (runPeriodA && !runPeriodB) {
        description->DrawLatexNDC(0.44, 0.91, "#it{x}_{a} Period A");
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        description->DrawLatexNDC(0.44, 0.91, "#it{x}_{a} Period B");
        hname = hname + "_periodB";
    }
    c2->Draw();
    c2->SaveAs((plotPath + hname + ".pdf").c_str());


    TCanvas* c3 = new TCanvas("c3", "", 1000, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    xaxpcorr->Draw("colz");
    description->SetTextAlign(22);
//    description->SetTextSize(0.036);
//    description->DrawLatexNDC(0.77, 0.9, "#bf{#it{ATLAS}} p-Pb");

//    c3->SetMargin(0.06, 0.14, 0.07, 0.05);
//    xaxpcorr->GetXaxis()->SetTitleOffset(0.7);
    xaxpcorr->GetXaxis()->SetTickLength(0.02);
//    xaxpcorr->GetXaxis()->SetLabelOffset(0.0025);

//    xaxpcorr->GetYaxis()->SetTitleOffset(0.7);
    xaxpcorr->GetYaxis()->SetTickLength(0.02);
//    xaxpcorr->GetYaxis()->SetLabelOffset(0.0025);

//    xaxpcorr->GetZaxis()->SetLabelOffset(0.0025);
    xaxpcorr->GetZaxis()->SetTickLength(0.01);
//    xaxpcorr->GetZaxis()->SetTitleOffset(1.2);

    description->SetTextSize(0.032);
    description->SetTextAlign(22);
    description->DrawLatexNDC(0.19, 0.55, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.19, 0.47, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", 1000.*integrated_luminosity)); 
    
    hname = "jets_xa_xp_correlation_8.16TeV";
    if (runPeriodA && runPeriodB) { 
        description->DrawLatexNDC(0.19, 0.395, "Period A(-#eta) & B(#eta)");
    } else if (runPeriodA && !runPeriodB) {
        description->DrawLatexNDC(0.19, 0.395, "Period A");
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        description->DrawLatexNDC(0.19, 0.395, "Period B");
        hname = hname + "_periodB";
    }

    c3->Draw();
    c3->SaveAs((plotPath + hname + ".pdf").c_str());
    
    TCanvas* c4 = new TCanvas("c4", "", 1000, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    fcalhist->Draw("colz");
    description->SetTextAlign(22);
//    c4->SetMargin(0.06, 0.14, 0.07, 0.05);
//    fcalhist->GetXaxis()->SetTitleOffset(0.7);
    fcalhist->GetXaxis()->SetTickLength(0.02);
//    fcalhist->GetXaxis()->SetLabelOffset(0.0025);

//    fcalhist->GetYaxis()->SetTitleOffset(0.7);
    fcalhist->GetYaxis()->SetTickLength(0.02);
//    fcalhist->GetYaxis()->SetLabelOffset(0.0025);

//    fcalhist->GetZaxis()->SetLabelOffset(0.0025);
    fcalhist->GetZaxis()->SetTickLength(0.01);
//    fcalhist->GetZaxis()->SetTitleOffset(1.2);

    description->SetTextSize(0.032);
    description->SetTextAlign(22);
    description->DrawLatexNDC(0.19, 0.55, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.19, 0.47, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", 1000.*integrated_luminosity)); 
    
    hname = "jets_xp_fcal_8.16TeV";
    if (runPeriodA && runPeriodB) { 
        description->DrawLatexNDC(0.19, 0.395, "Period A(-#eta) & B(#eta)");
    } else if (runPeriodA && !runPeriodB) {
        description->DrawLatexNDC(0.19, 0.395, "Period A");
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        description->DrawLatexNDC(0.19, 0.395, "Period B");
        hname = hname + "_periodB";
    }

    c4->Draw();
    c4->SaveAs((plotPath + hname + ".pdf").c_str());
    
}
