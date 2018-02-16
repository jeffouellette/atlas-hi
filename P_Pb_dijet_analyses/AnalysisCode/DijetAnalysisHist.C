#include "DijetAnalysis.C"

const Int_t Number = 3;
Double_t Red[Number]    = { 0.00, 1.00, 1.00};
Double_t Green[Number]  = { 0.00, 1.00, 0.00};
Double_t Blue[Number]   = { 1.00, 0.00, 0.00};
Double_t Length[Number] = { 0.00, 0.50, 1.00 };
Int_t nb=100;

void DijetAnalysisHist() {
    if (!runPeriodA && !runPeriodB) return;
    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    double eta_min, eta_max, q_min, q_max, mjj_min, mjj_max;
    double* histArrScalesX_eta;
    double* histArrScalesM_eta;
    double* histArrScalesX_Q;
    if (scaleDijetAnalyses) {
        eta_min = 1e-4;
        eta_max = 1e12;
        q_min = 1e-7;
        q_max = 1e13;
        mjj_min = 1e-3;
        mjj_max = 1e9;
        histArrScalesX_eta = linspace(-3.5, 3.5, numetabins-1); // returns -3, -1, 1, 3
        histArrScalesM_eta = linspace(-1.5, 1.5, numetabins/2 - 1); // returns -1.5, -0.5, 0.5, 1.5
        histArrScalesX_Q = linspace(-4.2, 4.2, numqbins - 1);
    } else {
        eta_min = 1e-3;
        eta_max = 1e9;
        q_min = 1e-4;
        q_max = 1e9;
        mjj_min = 5e-2;
        mjj_max = 5e7;
        histArrScalesX_eta = linspace(0, 0, numetabins/2 - 1);
        histArrScalesM_eta = linspace(0, 0, numetabins/2 - 1);
        histArrScalesX_Q = linspace(0, 0, numqbins - 1);
    }

//  Old plotting style
//    const Style_t mkstyles[8] = {kFullCircle, kFullDiamond, kFullSquare, kFullFourTrianglesX, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
//    const Color_t mkcolors[8] = {kAzure-5, kTeal-5, kOrange-5, kPink-5, kSpring-5, kViolet-5, kRed-2, kGray+3};

//  New plotting style
//    const Style_t mkstyles[2] = {kFullCircle, kOpenCircle};
//    const Color_t mkcolors[8] = {kGreen, kOrange+4, kCyan, kBlue, kViolet, kGreen+3, kRed, kOrange-3};
    const Color_t mkcolors[8] = {kOrange+4, kOrange-3, kRed, kViolet, kBlue, kCyan-2, kGreen+3, kTeal+9};

    TH1D* xHistArr[2*numetabins];
    TH1D* xqHistArr[2*numqbins];
    TH1D* mHistArr[numetabins];

    for (int etabin = 0; etabin < numetabins; etabin++) {
        xHistArr[etabin] = new TH1D(Form("eta%i", etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d^{2}#sigma/d#it{x}_{p}d#eta #left[nb#right]", etabins[etabin], etabins[etabin+1]), numxbins, xbins);
        xHistArr[etabin]->Sumw2();
    }
    for (int etabin = numetabins; etabin < 2*numetabins; etabin++) {
        xHistArr[etabin] = new TH1D(Form("eta%i", etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d^{2}#sigma/d#it{x}_{a}d#eta #left[nb#right]", etabins[etabin%numetabins], etabins[(etabin%numetabins)+1]), numxbins, xbins);
        xHistArr[etabin]->Sumw2();
    }

    for (int qbin = 0; qbin < numqbins; qbin++) {
        xqHistArr[qbin] = new TH1D(Form("q%i", qbin), Form("%g < #it{Q} < %g;#it{x}_{p};d^{2}#sigma/d#it{x}_{p}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin], qbins[qbin+1]), numxbins, xbins);
        xqHistArr[qbin]->Sumw2();
    }
    for (int qbin = numqbins; qbin < 2*numqbins; qbin++) {
        xqHistArr[qbin] = new TH1D(Form("q%i", qbin), Form("%g < #it{Q} < %g;#it{x}_{a};d^{2}#sigma/d#it{x}_{a}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin%numqbins], qbins[(qbin%numqbins)+1]), numxbins, xbins);
        xqHistArr[qbin]->Sumw2();
    }

    for (int etabin = 0; etabin < numetabins; etabin++) {
        mHistArr[etabin] = new TH1D(Form("mjj_eta%i", etabin), Form("%1.1f < #eta < %1.1f;#it{m}_{JJ} #left[GeV#right];d^{2}#sigma/d#it{m}_{JJ}d#eta #left[nb GeV^{-1}#right]", etabins[etabin], etabins[etabin+1]), nummbins, mbins);
        mHistArr[etabin]->Sumw2();
    }

    TH2D* qxcorr = new TH2D("xqcorr", ";#it{x}_{a};#it{Q}_{rms}^{2} #left[GeV^{2}#right];d^{2}#sigma/d#it{x}_{a}d#it{Q}_{rms}^{2} #left[nb GeV^{-2}#right]", numq2xbins, q2xbins, numq2bins, q2bins);
    TH2D* xaxpcorr = new TH2D("xaxpcorr", ";#it{x}_{a};#it{x}_{p};d^{2}#sigma/d#it{x}_{p}d#it{x}_{a}", numxbins, xbins, numxbins, xbins);
    TH2D* fcalHist = new TH2D("fcalHist", ";#it{x}_{p};FCAL energy [GeV];d^{2}#sigma/d#it{x}_{p}d#it{E}_{FCAL}", numxbins, xbins, numfcalbins, fcalbins);

    double integrated_luminosity = 0; // units are nb^{-1}
    TH1D* thisHist;

    {
        string period_dir;
        if (runPeriodA && !runPeriodB) period_dir = "periodA/";
        else if (!runPeriodA && runPeriodB) period_dir = "periodB/";
        else period_dir = "periodAB/";

        for (int thisRunNumber : (*thisRunNumbers)) {
            if (skipRun(thisRunNumber)) continue;
            TFile* thisFile = new TFile(Form("%s%srun_%i.root", xPath.c_str(), period_dir.c_str(), thisRunNumber), "READ");
            for (int etabin = 0; etabin < 2*numetabins; etabin++) {
                xHistArr[etabin]->Add((TH1D*)thisFile->Get(Form("%ieta%i", thisRunNumber, etabin)));
            }
            for (int qbin = 0; qbin < 2*numqbins; qbin++) {
                xqHistArr[qbin]->Add((TH1D*)thisFile->Get(Form("%iq%i", thisRunNumber, qbin)));
            }
            for (int etabin = 0; etabin < numetabins; etabin++) {
                mHistArr[etabin]->Add((TH1D*)thisFile->Get(Form("mjj_%ieta%i", thisRunNumber, etabin)));
            }
            qxcorr->Add((TH2D*)thisFile->Get(Form("xqcorr_run%i", thisRunNumber)));
            xaxpcorr->Add((TH2D*)thisFile->Get(Form("xaxpcorr_run%i", thisRunNumber)));
            fcalHist->Add((TH2D*)thisFile->Get(Form("fcalhist_run%i", thisRunNumber)));
            TVectorD* thisluminosityvec = (TVectorD*)(thisFile->Get("lum_vec")); // Accesses luminosity for this run and creates a pointer to it
            integrated_luminosity += (*thisluminosityvec)[0];   // Dereferences the luminosity vector pointer to add the run luminosity
            thisFile->Close();
            delete thisFile;
        }
    }


    double histScale;
    /**** Plots xp spectrum, binned in eta ****/
    TCanvas* canvas = new TCanvas("xp_etabinned", "", 800, 600);
    canvas->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetErrorX(0);
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thisHist = xHistArr[etabin];

        histScale = histArrScalesX_eta[etabin];
    //    histScale = histArrScalesX_eta[(int)(numetabins/2 - 0.5 -TMath::Abs((etabin%numetabins) - numetabins/2 + 0.5))];
        thisHist->Scale(TMath::Power(10, histScale));
        
        thisHist->SetAxisRange(eta_min, eta_max, "Y");

        //Style_t mkstyle = mkstyles[etabin%2];
        Style_t mkstyle = kFullDotMedium;
        Color_t mkcolor = mkcolors[etabin];

        thisHist->SetMarkerStyle(mkstyle);
        thisHist->SetMarkerColor(mkcolor);
        thisHist->SetLineColor(mkcolor);
        thisHist->Draw("same e1");
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);

        const double textx = 0.22 + 0.27*(etabin>=(numetabins/2));
        const double texty = 0.35 - 0.05*((etabin%(numetabins/2))*(etabin >= (numetabins/2)) + (numetabins/2 - etabin - 1)*(etabin < (numetabins/2)));
        char* text;
        if (scaleDijetAnalyses) text = Form("%g < #it{#eta}_{JJ} < %g (#times10^{%1.1f})", etabins[etabin], etabins[etabin+1], histScale);
        else text = Form("%g < #it{#eta}_{JJ} < %g", etabins[etabin], etabins[etabin+1]);
        myMarkerText (textx, texty, mkcolor, kFullCircle, text);
    }
    myText (0.2, 0.9, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));

    string histName = "dijets_xp_pseudorapidity_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";
    if (scaleDijetAnalyses) histName = histName + "_scaled";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    for (int etabin = 0; etabin < numetabins; etabin++) delete xHistArr[etabin];
    delete canvas;
    /**** End plot xp spectrum, eta binned ****/


    /**** Plots xa spectrum, binned in eta ****/
    canvas = new TCanvas("xa_etabinned", "", 800, 600);
    canvas->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetErrorX(0);
    for (int etabin = numetabins; etabin < 2*numetabins; etabin++) {
        thisHist = xHistArr[etabin];
        int temp = etabin%numetabins;

//        histScale = histArrScalesX_eta[(int)(numetabins/2 - 0.5 -TMath::Abs((temp) - numetabins/2 + 0.5))];
        histScale = histArrScalesX_eta[numetabins - temp - 1];
        thisHist->Scale(TMath::Power(10, histScale));

        thisHist->SetAxisRange(eta_min, eta_max, "Y");

        //Style_t mkstyle = mkstyles[temp%2];
        Style_t mkstyle = kFullDotMedium;
        Color_t mkcolor = mkcolors[temp];

        thisHist->SetMarkerStyle(mkstyle);
        thisHist->SetMarkerColor(mkcolor);
        thisHist->SetLineColor(mkcolor);
        thisHist->Draw("same e1");
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);

        const double textx = 0.22 + 0.27*(temp>=(numetabins/2));
        const double texty = 0.35 - 0.05*((temp%(numetabins/2))*(temp>=(numetabins/2)) + (numetabins/2 - temp - 1)*(temp < (numetabins/2)));
        char* text;
        if (scaleDijetAnalyses) text = Form("%g < #it{#eta}_{JJ} < %g (#times10^{%1.1f})", etabins[temp], etabins[temp+1], histScale);
        else text = Form("%g < #it{#eta}_{JJ} < %g", etabins[temp], etabins[temp+1]);
        myMarkerText (textx, texty, mkcolor, kFullCircle, text);
    }
    myText (0.2, 0.9, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));

    histName = "dijets_xa_pseudorapidity_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";
    if (scaleDijetAnalyses) histName = histName + "_scaled";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    for (int etabin = numetabins; etabin < 2*numetabins; etabin++) delete xHistArr[etabin];
    delete canvas;
    /**** End plot xa spectrum, eta binned ****/


    /**** Plots xp spectrum, binned in hardness ****/
    canvas = new TCanvas("xp_qbinned", "", 800, 600);
    canvas->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetErrorX(0);
    for (int qbin = 0; qbin < numqbins; qbin++) {
        thisHist = xqHistArr[qbin];

        histScale = histArrScalesX_Q[numqbins - qbin - 1];
        thisHist->Scale(TMath::Power(10, histScale));

        thisHist->SetAxisRange(q_min, q_max, "Y");

        //Style_t mkstyle = mkstyles[qbin%2];
        Style_t mkstyle = kFullDotMedium;
        Color_t mkcolor = mkcolors[qbin];

        thisHist->SetMarkerStyle(mkstyle);
        thisHist->SetMarkerColor(mkcolor);
        thisHist->SetLineColor(mkcolor);
        thisHist->Draw("same e1");
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);

        const double textx = 0.22;
        const double texty = 0.55 - 0.05*qbin;
        char* text;
        if (scaleDijetAnalyses) text = Form("%.0f < #left|Q#right| < %.0f (#times10^{%g})", qbins[qbin], qbins[qbin+1], histScale);
        else text = Form("%.0f < #left|Q#right| < %.0f", qbins[qbin], qbins[qbin+1]);
        myMarkerText (textx, texty, mkcolor, kFullCircle, text);
    }
    myText (0.2, 0.9, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));

    histName = "dijets_xp_hardness_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";
    if (scaleDijetAnalyses) histName = histName + "_scaled";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    for (int qbin = 0; qbin < numqbins; qbin++) delete xqHistArr[qbin];
    delete canvas;
    /**** End plot xp spectrum, q binned ****/


    /**** Plots xa spectrum, binned in hardness ****/
    canvas = new TCanvas("xa_qbinned", "", 800, 600);
    canvas->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetErrorX(0);
    for (int qbin = numqbins; qbin < 2*numqbins; qbin++) {
        thisHist = xqHistArr[qbin];
        int temp = qbin%numqbins;

        histScale = histArrScalesX_Q[numqbins - temp - 1];
        thisHist->Scale(TMath::Power(10, histScale));

        thisHist->SetAxisRange(q_min, q_max, "Y");

        //Style_t mkstyle = mkstyles[temp%2];
        Style_t mkstyle = kFullDotMedium;
        Color_t mkcolor = mkcolors[temp];

        thisHist->SetMarkerStyle(mkstyle);
        thisHist->SetMarkerColor(mkcolor);
        thisHist->SetLineColor(mkcolor);
        thisHist->Draw("same e1");
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);

        const double textx = 0.22;
        const double texty = 0.55 - 0.05*temp;
        char* text;
        if (scaleDijetAnalyses) text = Form("%.0f < #left|Q#right| < %.0f (#times10^{%g})", qbins[temp], qbins[temp+1], histScale);
        else text = Form("%.0f < #left|Q#right| < %.0f", qbins[temp], qbins[temp+1]);
        myMarkerText (textx, texty, mkcolor, kFullCircle, text);
    }
    myText (0.2, 0.9, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));

    histName = "dijets_xa_hardness_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";
    if (scaleDijetAnalyses) histName = histName + "_scaled";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    for (int qbin = numqbins; qbin < 2*numqbins; qbin++) delete xqHistArr[qbin];
    delete canvas;
    /**** End plot xa spectrum, q binned ****/


    /**** Plots the dijet invariant mass (M_jj) spectrum, eta binned ****/
    canvas = new TCanvas("mjj_etabinned", "", 800, 600);
    canvas->cd();
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetErrorX(0);
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thisHist = mHistArr[etabin];

        histScale = histArrScalesM_eta[(int)(numetabins/2 - 0.5 -TMath::Abs(etabin - numetabins/2 + 0.5))];
        thisHist->Scale(TMath::Power(10, histScale));

        thisHist->SetAxisRange(mjj_min, mjj_max, "Y");

        //Style_t mkstyle = mkstyles[etabin%2];
        Style_t mkstyle = kFullDotMedium;
        Color_t mkcolor = mkcolors[etabin];

        thisHist->SetMarkerStyle(mkstyle);
        thisHist->SetMarkerColor(mkcolor);
        thisHist->SetLineColor(mkcolor);
        thisHist->Draw("same e1");
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);

        const float textx = 0.45 + (etabin>=(numetabins/2))*0.27;
        const float texty = 0.85 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
        char* text;
        if (scaleDijetAnalyses) text = Form("%g < #it{#eta}_{JJ} < %g (#times10^{%1.1f})", etabins[etabin], etabins[etabin+1], histScale);
        else text = Form("%g < #it{#eta}_{JJ} < %g", etabins[etabin], etabins[etabin+1]);
        myMarkerText (textx, texty, mkcolor, kFullCircle, text);
    }
    myText (0.43, 0.9, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));

    histName = "dijets_mjj_etabinned_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";
    if (scaleDijetAnalyses) histName = histName + "_scaled";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    for (int etabin = 0; etabin < numetabins; etabin++) delete mHistArr[etabin];
    delete canvas;
    /**** End plot M_jj spectrum, eta binned ****/


    /**** Plot q2-xa orrelation plot ****/  
    canvas = new TCanvas("qx_correlation", "", 800, 600);
    canvas->cd();
    canvas->SetRightMargin(0.18);
    //canvas->SetTopMargin(-0.02);
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    qxcorr->Draw("colz");

    //xaxpcorr->GetXaxis()->SetTickLength(0.02);
    //xaxpcorr->GetYaxis()->SetTickLength(0.02);
//    xaxpcorr->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.25, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));
    
    histName = "dijets_Q2_xa_correlation_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    delete qxcorr;
    delete canvas;
    /**** End plot q2-xa correlation ****/


    /**** Plot xa-xp cross correlation plot ****/  
    canvas = new TCanvas("xaxp_correlation", "", 800, 600);
    canvas->cd();
    canvas->SetRightMargin(0.18);
    //canvas->SetTopMargin(-0.02);
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    xaxpcorr->Draw("colz");

    //xaxpcorr->GetXaxis()->SetTickLength(0.02);
    //xaxpcorr->GetYaxis()->SetTickLength(0.02);
//    xaxpcorr->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.25, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));
    
    histName = "dijets_xa_xp_correlation_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    delete xaxpcorr;
    delete canvas;
    /**** End plot xa-xp correlation ****/


    /**** Plots FCAL energy in Pb-going direction as a function of xp ****/    
    canvas = new TCanvas("fcal_xp", "", 800, 600);
    canvas->cd();
    canvas->SetRightMargin(0.18);
    //canvas->SetTopMargin(-0.02);
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    fcalHist->Draw("colz");
    //fcalHist->GetXaxis()->SetTickLength(0.02);
    //fcalHist->GetYaxis()->SetTickLength(0.02);
    //fcalHist->GetZaxis()->SetTickLength(0.01);

    myText (0.19, 0.91, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));
    //myText (0.19, 0.47, kBlack, Form("#it{L}_{int} = %.1f nb^{-1}", integrated_luminosity));
    
    histName = "dijets_xp_fcal_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    delete canvas;
    /**** End plot FCAL-xp ****/


    /**** Plots the TProfileX of the FCAL energy deposition plot above ****/
    canvas = new TCanvas("fcal_tprofilex", "", 800, 600);
    canvas->cd();
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetErrorX(0);
    
    TProfile* fcalProfileX = fcalHist->ProfileX("fcalProfileX");
    fcalProfileX->SetMarkerStyle(kFullDotMedium);
    fcalProfileX->Draw("e1");
    myText (0.19, 0.91, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{s} = 8.16 TeV", integrated_luminosity));
    histName = "dijets_xp_fcal_profilex_8.16TeV";
    if (runPeriodA && !runPeriodB) histName = histName + "_periodA";
    else if (!runPeriodA && runPeriodB) histName = histName + "_periodB";

    canvas->Draw();
    canvas->SaveAs((plotPath + "dijets/" + histName + ".pdf").c_str());
    delete canvas;
    delete fcalHist;
    delete fcalProfileX;
    /**** End plot FCAL TProfileX ****/
}
