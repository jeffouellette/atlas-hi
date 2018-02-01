#include "TriggerEfficiencyAnalysis.C"

void TriggerEfficiencyAnalysisHist() {

    initialize(0, false, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();
    const int numhists = 2*numtrigs;
    const int numruns = (*thisRunNumbers).size();

    TH1F* histArr[numhists];
    TH1F* plotHistArr[numhists/2];
    int pbin, ebin, index, referenceIndex;
    double reduced_pbins[40] = {};
    const int num_reduced_pbins = sizeof(reduced_pbins)/sizeof(reduced_pbins[0]) - 1;
    for (int reduced_pbin = 0; reduced_pbin < num_reduced_pbins+1; reduced_pbin++) {
        reduced_pbins[reduced_pbin] = pbins[reduced_pbin];
    }

    for (Trigger* trig : triggerVec) {
        index = trig->index;
        // standard trigger firings
        TString histName = Form("%s_efficiency", trig->name.c_str());
        histArr[index] = new TH1F(histName, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", numpbins, pbins);
        histArr[index]->Sumw2(); // instruct each histogram to propagate errors
        // reference trigger firings
        histName = Form("%s_reference", trig->name.c_str());
        histArr[index+numtrigs] = new TH1F(histName, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", numpbins, pbins);
        histArr[index+numtrigs]->Sumw2();
        // plotted histograms - difference is just fewer pbins.
        histName = Form("%s_efficiency_plot", trig->name.c_str());
        plotHistArr[index] = new TH1F(histName, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", num_reduced_pbins, reduced_pbins);
        plotHistArr[index]->Sumw2(); // instruct each histogram to propagate errors
    }

    int thisRunNumber;
    double totalLuminosity = 0;
    TH1F* thisHist;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            if (debugStatements) cout << "Status: In trigger_efficiencies_hist.C (41): Skipping run " << thisRunNumber << endl;
            continue;
        }

        TFile* thisFile = new TFile(Form("%srun_%i.root", effPath.c_str(), thisRunNumber), "READ");
        TVectorD* thisLuminosityVec = (TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber));
        TVectorD* thisRunVec = (TVectorD*)thisFile->Get(Form("run_vec_%i", thisRunNumber));

        assert (thisRunVec[0] == thisRunNumber); // Check that the data matches what we want to analyze.
        assert (thisRunVec[1] == numetabins);
        assert (thisRunVec[2] == numtrigs);

        totalLuminosity += (*thisLuminosityVec)[0];

        for (Trigger* trig : triggerVec) {
            if (!(trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber)) continue;
            index = trig->index;
            thisHist = (TH1F*)thisFile->Get(Form("%s_efficiency_run%i", trig->name.c_str(), thisRunNumber));
            histArr[index]->Add(thisHist);
            thisHist = (TH1F*)thisFile->Get(Form("%s_reference_run%i", trig->name.c_str(), thisRunNumber));
            histArr[index+numtrigs]->Add(thisHist);
        }
        thisFile->Close();
    }

    /** Go from counts in each histogram to efficiencies **/
    for (Trigger* trig : triggerVec) {
        histArr[trig->index]->Divide(histArr[(trig->index)+numtrigs]);
    }
    for (Trigger* trig : triggerVec) {
        bootstrap(trig, histArr);
    }
    for (Trigger* trig : triggerVec) {
        thisHist = plotHistArr[trig->index];
        for (int pbin = 0; pbin < num_reduced_pbins; pbin++) {
            thisHist->SetBinContent(pbin+1, histArr[trig->index]->GetBinContent(pbin+1));
            thisHist->SetBinError(pbin+1, histArr[trig->index]->GetBinError(pbin+1));
        }
    }

    /** Plotting routines **/

    TCanvas* trigCanvas = new TCanvas("trigCanvas", "", 800, 600);
    TFile* output = new TFile((effPath+"allEfficiencyHistograms.root").c_str(), "RECREATE");
    TGraph* thisGraph;
    for (Trigger* trig : triggerVec) {
        gPad->SetLogx();
        gPad->SetTicks();
        int index = trig->index;
        thisHist = plotHistArr[index];

        Color_t kColor = 46;

        thisGraph = new TGraph(num_reduced_pbins);
        for (int pbin = 0; pbin < num_reduced_pbins; pbin++) thisGraph->SetPoint(pbin, 0.5*(reduced_pbins[pbin+1]+reduced_pbins[pbin]), thisHist->GetBinContent(pbin+1));
        thisGraph->SetMarkerColor(kColor);
        thisGraph->SetLineColor(kBlack);
        thisGraph->SetLineWidth(2*thisGraph->GetLineWidth());
        thisGraph->SetMarkerStyle(kDot);
        thisGraph->SetMaximum(1.3);
        thisGraph->Draw();

        thisHist->SetMarkerStyle(kDot);
        thisHist->SetMarkerColor(kColor);
        thisHist->SetLineColor(kColor);
        thisHist->SetMinimum(0);
        thisHist->SetMaximum(1.3);
        thisHist->GetYaxis()->SetTitleOffset(1.35);
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);
        
        thisHist->Draw("SAME E1");
        myText (0.18, 0.91, kBlack, ("Trig: " + trig->name).c_str());
        myText (0.18, 0.85, kBlack, ("Ref: " + trig->referenceTrigger->name).c_str());
        myText (0.65, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity*1000.));
        myText (0.65, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));
        trigCanvas->Draw();
        
        string histName = Form("%s_efficiency", trig->name.c_str());

        trigCanvas->SaveAs((plotPath + "triggerEfficiencies/" + histName + ".pdf").c_str());
        histArr[index]->Write();
        delete plotHistArr[index];
        delete histArr[index];
        trigCanvas->Clear();
    }
    output->Close();
    delete trigCanvas;
    delete output;

    if (debugStatements) cout << "Status: In trigger_efficiencies_hist.C (119): Finished plotting efficiency curves" << endl;
    return;
}
