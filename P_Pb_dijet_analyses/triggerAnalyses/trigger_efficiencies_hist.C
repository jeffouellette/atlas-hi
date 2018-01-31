#include "trigger_efficiencies.C"

void trigger_efficiencies_hist() {

    initialize(0, false, false, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();
    const int numhists = 2*numtrigs;
    const int numruns = (*thisRunNumbers).size();

    TH1F* histArr[numhists];
    int pbin, ebin, index, referenceIndex;
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
    }

    int thisRunNumber;
    double totalLuminosity = 0;
    TH1F* thisHist;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            if (debugStatements) cout << "Status: In trigger_efficiencies_hist.C (30): Skipping run " << thisRunNumber << endl;
            continue;
        }

        TFile* thisFile = new TFile(Form("%srun_%i.root", effPath.c_str(), thisRunNumber), "READ");
        TVectorD* thisluminosityvec = (TVectorD*)thisFile->Get("lum_vec");
        TVectorD* thisrunvec = (TVectorD*)thisFile->Get("run_vec");

        assert (thisrunvec[0] == thisRunNumber); // Check that the data matches what we want to analyze.
        assert (thisrunvec[1] == numetabins);
        assert (thisrunvec[2] == numtrigs);

        totalLuminosity += (*thisluminosityvec)[0];

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

    /** Plotting routines **/

    TCanvas* trigCanvas = new TCanvas("trigCanvas", "", 800, 600);
    TFile* output = new TFile((effPath+"allEfficiencyHistograms.root").c_str(), "RECREATE");
    for (Trigger* trig : triggerVec) {
        gPad->SetLogx();
        gPad->SetTicks();
        int index = trig->index;
        thisHist = histArr[index];
        thisHist->SetMarkerStyle(kDot);
        Color_t kColor = 46;
        thisHist->SetMarkerColor(kColor);
        thisHist->SetLineColor(kColor);
        thisHist->SetMinimum(0);
        thisHist->SetMaximum(1.2);
        thisHist->GetYaxis()->SetTitleOffset(1.35);
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);
        thisHist->Draw("e1");
        myText (0.7, 0.33, kBlack, trig->name.c_str());
        myText (0.7, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity*1000));
        myText (0.7, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));
        trigCanvas->Draw();
        
        string histName = Form("%s_efficiency", trig->name.c_str());

        trigCanvas->SaveAs((plotPath + "triggerEfficiencies/" + histName + ".pdf").c_str());
        thisHist->Write();
        delete histArr[index];
        trigCanvas->Clear();
    }
    output->Close();
    delete trigCanvas;
    delete output;

    return;
}
