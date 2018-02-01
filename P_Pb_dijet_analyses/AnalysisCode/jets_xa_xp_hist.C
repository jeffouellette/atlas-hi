#include "../triggerUtil.C"
//#include <AtlasUtils.C>

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
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #it{#eta}_{B} < %1.1f;#it{x}_{p};d^{2}N/#it{L}_{int}d#it{x}_{p} dy #left[pb#right]", etabins[i], etabins[i+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numetabins; i < numhists; i++) {
        harr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #it{#eta}_{B} < %1.1f;#it{x}_{a};d^{2}N/#it{L}_{int}d#it{x}_{a} dy #left[pb#right]", etabins[i%(numetabins)], etabins[(i%(numetabins))+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    TH2D* xaxpcorr = new TH2D("xaxpcorr", ";#it{x}_{a};#it{x}_{p};d^{2}N/#it{L}_{int}d#it{x}_{p}d#it{x}_{a}", numbins, xbins, numbins, xbins);

    const int numfcalbins = 60;
    const double* fcalbins = logspace(10, 500, numfcalbins);
    TH2D* fcalhist = new TH2D("fcalhist", ";#it{x}_{p};FCAL energy [GeV];", numbins, xbins, numfcalbins, fcalbins);

    double integrated_luminosity = 0;

    string period_dir;
    if (runPeriodA && !runPeriodB) period_dir = "periodA/";
    else if (!runPeriodA && runPeriodB) period_dir = "periodB/";
    else period_dir = "periodAB/";

    for (int thisRunNumber : (*thisRunNumbers)) {
        if (skipRun(thisRunNumber)) continue;
        TFile* thisfile = new TFile(Form("%sxdata/%srun_%i.root", rootPath.c_str(), period_dir.c_str(), thisRunNumber), "READ");
        for (int j = 0; j < numhists; j++) {
            int act_j = j;
            if (thisRunNumber < 313500) act_j = numhists - (j<numetabins)*numetabins - (j%numetabins) - 1;
            harr[j]->Add((TH1D*)thisfile->Get(Form("%ieta%i", thisRunNumber, act_j)));
        }
        xaxpcorr->Add((TH2D*)thisfile->Get(Form("xaxpcorr_run%i", thisRunNumber)));
        fcalhist->Add((TH2D*)thisfile->Get(Form("fcalhist_run%i", thisRunNumber)));
        TVectorD* thisluminosityvec = (TVectorD*)(thisfile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
        integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
    }

    TCanvas* c1 = new TCanvas("c1", "", 800, 600);
    gPad->SetLogy();
    gPad->SetLogx();
    for (int i = 0; i < numetabins; i++) {
        harr[i]->SetAxisRange(ymin, ymax, "Y");
        harr[i]->SetMarkerStyle(mkstyles[i]);
        harr[i]->SetMarkerColor(mkcolors[i]);
        harr[i]->SetLineColor(mkcolors[i]);
        harr[i]->Draw("same e1");
        harr[i]->GetXaxis()->SetTickLength(0.02);
        harr[i]->GetYaxis()->SetTickLength(0.02);
        myMarkerText (0.22 + (i>=(numetabins/2))*0.24, 0.39 - (i%(numetabins/2))*0.05*(i>=(numetabins/2)) - (-i-1+numetabins/2)*0.05*(i<(numetabins/2)), mkcolors[i], mkstyles[i], Form("%g < #it{#eta}_{B} < %g", etabins[i], etabins[i+1]));
        //myMarkerText (0.22+(i>=(numetabins/2))*0.28, 0.39-(i%(numetabins/2))*0.06, mkcolors[i], mkstyles[i], Form("%1.1f < #eta_{lab} < %1.1f", etabins[i], etabins[i+1]));
    }
    myText (0.2, 0.47, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s_{NN}^{avg}} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.22, 0.45, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));

    string hname = "jets_xp_8.16TeV";
    if (runPeriodA && !runPeriodB) {
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        hname = hname + "_periodB";
    }
    c1->Draw();
    c1->SaveAs((plotPath + hname + ".pdf").c_str());


    TCanvas* c2 = new TCanvas("c2", "", 800, 600);
    gPad->SetLogy();
    gPad->SetLogx();
    for (int i = numetabins; i < numhists; i++) {
        int j = i%numetabins;
        harr[i]->SetAxisRange(ymin, ymax, "Y");
        harr[i]->SetMarkerStyle(mkstyles[j]);
        harr[i]->SetMarkerColor(mkcolors[j]);
        harr[i]->SetLineColor(mkcolors[j]);
        harr[i]->Draw("same e1");
        harr[i]->GetXaxis()->SetTickLength(0.02);
        harr[i]->GetYaxis()->SetTickLength(0.02);
        myMarkerText (0.22 + (j>=(numetabins/2))*0.24, 0.39 - (j%(numetabins/2))*0.05*(j>=(numetabins/2)) - (-j-1+numetabins/2)*0.05*(j<(numetabins/2)), mkcolors[j], mkstyles[j], Form("%g < #it{#eta}_{B} < %g", etabins[j], etabins[j+1]));
        //myMarkerText (0.22+(j>=(numetabins/2))*0.28, 0.39-(j%(numetabins/2))*0.06, mkcolors[j], mkstyles[j], Form("%1.1f < #eta_{lab} < %1.1f", etabins[j], etabins[j+1]));
    }
    myText (0.2, 0.47, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s_{NN}^{avg}} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.22, 0.45, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));

    hname = "jets_xa_8.16TeV";
    if (runPeriodA && !runPeriodB) {
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        hname = hname + "_periodB";
    }
    c2->Draw();
    c2->SaveAs((plotPath + hname + ".pdf").c_str());


    TCanvas* c3 = new TCanvas("c3", "", 800, 600);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    xaxpcorr->Draw("col");

    xaxpcorr->GetXaxis()->SetTickLength(0.02);
    xaxpcorr->GetYaxis()->SetTickLength(0.02);
    xaxpcorr->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.31, kBlack, Form("z = d^{2}N/#it{L}_{int}d#it{x}_{a}d#it{x}_{p}"));
    myText (0.19, 0.25, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s_{NN}^{avg}} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.19, 0.47, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));
    
    hname = "jets_xa_xp_correlation_8.16TeV";
    if (runPeriodA && !runPeriodB) {
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        hname = hname + "_periodB";
    }

    c3->Draw();
    c3->SaveAs((plotPath + hname + ".pdf").c_str());
    
    TCanvas* c4 = new TCanvas("c4", "", 800, 600);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    fcalhist->Draw("col");
    fcalhist->GetXaxis()->SetTickLength(0.02);
    fcalhist->GetYaxis()->SetTickLength(0.02);
    fcalhist->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.85, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s_{NN}^{avg}} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.19, 0.47, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));
    
    hname = "jets_xp_fcal_8.16TeV";
    if (runPeriodA && !runPeriodB) {
        hname = hname + "_periodA";
    } else if (!runPeriodA && runPeriodB) {
        hname = hname + "_periodB";
    }

    c4->Draw();
    c4->SaveAs((plotPath + hname + ".pdf").c_str());
    
}
