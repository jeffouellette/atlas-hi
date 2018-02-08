#include "../triggerUtil.C"
//#include <AtlasUtils.C>

void DijetAnalysisHist() {
    if (!runPeriodA && !runPeriodB) return;
    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numxbins = 100;
    const int numhists = 2*numetabins;

    const double ymin = 5e-3;
    const double ymax = 3e9;

    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[8] = {kAzure-5, kTeal-5, kOrange-5, kPink-5, kSpring-5, kViolet-5, kRed-2, kGray+3};

    const double* xbins = logspace(8e-5, 1.6, numxbins);
    TH1D* histArr[numhists];
    for (int i = 0; i < numetabins; i++) {
        histArr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #it{#eta}_{B} < %1.1f;#it{x}_{p};d^{2}N/#it{L}_{int}d#it{x}_{p} dy #left[nb#right]", etabins[i], etabins[i+1]), numxbins, xbins);
        histArr[i]->Sumw2();
    }
    for (int i = numetabins; i < numhists; i++) {
        histArr[i] = new TH1D(Form("eta%i", i), Form("%1.1f < #it{#eta}_{B} < %1.1f;#it{x}_{a};d^{2}N/#it{L}_{int}d#it{x}_{a} dy #left[nb#right]", etabins[i%(numetabins)], etabins[(i%(numetabins))+1]), numxbins, xbins);
        histArr[i]->Sumw2();
    }
    TH2D* xaxpcorr = new TH2D("xaxpcorr", ";#it{x}_{a};#it{x}_{p};d^{2}N/#it{L}_{int}d#it{x}_{p}d#it{x}_{a}", numxbins, xbins, numxbins, xbins);

    const int numfcalbins = 60;
    const double* fcalbins = logspace(10, 500, numfcalbins);
    TH2D* fcalHist = new TH2D("fcalHist", ";#it{x}_{p};FCAL energy [GeV];", numxbins, xbins, numfcalbins, fcalbins);

    double integrated_luminosity = 0;

    {
        string period_dir;
        if (runPeriodA && !runPeriodB) period_dir = "periodA/";
        else if (!runPeriodA && runPeriodB) period_dir = "periodB/";
        else period_dir = "periodAB/";

        for (int thisRunNumber : (*thisRunNumbers)) {
            if (skipRun(thisRunNumber)) continue;
            TFile* thisFile = new TFile(Form("%s%srun_%i.root", xPath.c_str(), period_dir.c_str(), thisRunNumber), "READ");
            for (int j = 0; j < numhists; j++) {
                int act_j = j;
                if (thisRunNumber < 313500) act_j = numhists - (j<numetabins)*numetabins - (j%numetabins) - 1; // flips period A pseudorapidities
                histArr[j]->Add((TH1D*)thisFile->Get(Form("%ieta%i", thisRunNumber, act_j)));
            }
            xaxpcorr->Add((TH2D*)thisFile->Get(Form("xaxpcorr_run%i", thisRunNumber)));
            fcalHist->Add((TH2D*)thisFile->Get(Form("fcalhist_run%i", thisRunNumber)));
            TVectorD* thisluminosityvec = (TVectorD*)(thisFile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
            integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
            thisFile->Close();
            delete thisFile;
        }
    }

    TCanvas* c1 = new TCanvas("c1", "", 800, 600);
    c1->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    for (int i = 0; i < numetabins; i++) {
        histArr[i]->SetAxisRange(ymin, ymax, "Y");
        histArr[i]->SetMarkerStyle(mkstyles[i]);
        histArr[i]->SetMarkerColor(mkcolors[i]);
        histArr[i]->SetLineColor(mkcolors[i]);
        histArr[i]->Draw("same e1");
        histArr[i]->GetXaxis()->SetTickLength(0.02);
        histArr[i]->GetYaxis()->SetTickLength(0.02);
        myMarkerText (0.22 + (i>=(numetabins/2))*0.24, 0.39 - (i%(numetabins/2))*0.05*(i>=(numetabins/2)) - (-i-1+numetabins/2)*0.05*(i<(numetabins/2)), mkcolors[i], mkstyles[i], Form("%g < #it{#eta}_{B} < %g", etabins[i], etabins[i+1]));
        //myMarkerText (0.22+(i>=(numetabins/2))*0.28, 0.39-(i%(numetabins/2))*0.06, mkcolors[i], mkstyles[i], Form("%1.1f < #eta_{lab} < %1.1f", etabins[i], etabins[i+1]));
    }
    myText (0.2, 0.47, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s_{NN}^{avg}} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.22, 0.45, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));

    string histName = "jets_xp_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    c1->Draw();
    c1->SaveAs((plotPath + histName + ".pdf").c_str());
    for (int i = 0; i < numetabins; i++) delete histArr[i];
    delete c1;


    TCanvas* c2 = new TCanvas("c2", "", 800, 600);
    c2->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    for (int i = numetabins; i < numhists; i++) {
        int j = i%numetabins;
        histArr[i]->SetAxisRange(ymin, ymax, "Y");
        histArr[i]->SetMarkerStyle(mkstyles[j]);
        histArr[i]->SetMarkerColor(mkcolors[j]);
        histArr[i]->SetLineColor(mkcolors[j]);
        histArr[i]->Draw("same e1");
        histArr[i]->GetXaxis()->SetTickLength(0.02);
        histArr[i]->GetYaxis()->SetTickLength(0.02);
        myMarkerText (0.22 + (j>=(numetabins/2))*0.24, 0.39 - (j%(numetabins/2))*0.05*(j>=(numetabins/2)) - (-j-1+numetabins/2)*0.05*(j<(numetabins/2)), mkcolors[j], mkstyles[j], Form("%g < #it{#eta}_{B} < %g", etabins[j], etabins[j+1]));
        //myMarkerText (0.22+(j>=(numetabins/2))*0.28, 0.39-(j%(numetabins/2))*0.06, mkcolors[j], mkstyles[j], Form("%1.1f < #eta_{lab} < %1.1f", etabins[j], etabins[j+1]));
    }
    myText (0.2, 0.47, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.22, 0.45, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));

    histName = "jets_xa_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    c2->Draw();
    c2->SaveAs((plotPath + histName + ".pdf").c_str());
    for (int i = numetabins; i < numhists; i++) delete histArr[i];
    delete c2;


    TCanvas* c3 = new TCanvas("c3", "", 800, 600);
    c3->cd();
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    xaxpcorr->Draw("col");

    xaxpcorr->GetXaxis()->SetTickLength(0.02);
    xaxpcorr->GetYaxis()->SetTickLength(0.02);
    xaxpcorr->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.31, kBlack, Form("z = d^{2}N/#it{L}_{int}d#it{x}_{a}d#it{x}_{p}"));
    myText (0.19, 0.25, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.19, 0.47, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));
    
    histName = "jets_xa_xp_correlation_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    c3->Draw();
    c3->SaveAs((plotPath + histName + ".pdf").c_str());
    delete xaxpcorr;
    delete c3;
    
    TCanvas* c4 = new TCanvas("c4", "", 800, 600);
    c4->cd();
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    fcalHist->Draw("col");
    fcalHist->GetXaxis()->SetTickLength(0.02);
    fcalHist->GetYaxis()->SetTickLength(0.02);
    fcalHist->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.91, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity*1000));
    //myText (0.19, 0.47, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity*1000));
    
    histName = "jets_xp_fcal_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    c4->Draw();
    c4->SaveAs((plotPath + histName + ".pdf").c_str());
    delete c4;

    TCanvas* c5 = new TCanvas("c5", "", 800, 600);
    c5->cd();
    gPad->SetLogx();
    gPad->SetLogy();
    TProfile* fcalProfileX = fcalHist->ProfileX("fcalProfileX");
    fcalProfileX->SetMarkerStyle(kDot);
    fcalProfileX->Draw("e1");
    myText (0.19, 0.91, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity*1000));
    histName = "jets_xp_fcal_profilex_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    c5->Draw();
    c5->SaveAs((plotPath + histName + ".pdf").c_str());
    delete c5;
    delete fcalHist;
    delete fcalProfileX;
}
