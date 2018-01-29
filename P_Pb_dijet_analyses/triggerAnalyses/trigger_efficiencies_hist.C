#include "trigger_efficiencies.C"

void trigger_efficiencies_hist() {

    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();
    const int numhists = 2*numtrigs;
    const int numruns = (*thisRunNumbers).size();

    TH1F* harr[numhists];
    int pbin, ebin, index, referenceIndex;
    for (Trigger* trig : trigger_vec) {
        index = trig->index;
        // standard trigger firings
        TString histname = Form("%s_efficiency", trig->name.c_str());
        harr[index] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", numpbins, pbins);
        harr[index]->Sumw2(); // instruct each histogram to propagate errors
        // reference trigger firings
        histname = Form("%s_reference", trig->name.c_str());
        harr[index+numtrigs] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", numpbins, pbins);
        harr[index+numtrigs]->Sumw2();
    }

    int thisRunNumber;
    double total_luminosity = 0;
    TH1F* thishist;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            cout << "Skipping run " << thisRunNumber << endl;
            continue;
        }

        TFile* thisfile = new TFile(Form("%srun_%i.root", effPath.c_str(), thisRunNumber), "READ");
        TVectorD* thisluminosityvec = (TVectorD*)thisfile->Get("lum_vec");
        TVectorD* thisrunvec = (TVectorD*)thisfile->Get("run_vec");

        assert (thisrunvec[0] == thisRunNumber); // Check that the data matches what we want to analyze.
        assert (thisrunvec[1] == numetabins);
        assert (thisrunvec[2] == numtrigs);

        total_luminosity += (*thisluminosityvec)[0];

        for (Trigger* trig : trigger_vec) {
            if (!(trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber)) continue;
            index = trig->index;
            thishist = (TH1F*)thisfile->Get(Form("%s_efficiency_run%i", trig->name.c_str(), thisRunNumber));
            harr[index]->Add(thishist);
            thishist = (TH1F*)thisfile->Get(Form("%s_reference_run%i", trig->name.c_str(), thisRunNumber));
            harr[index+numtrigs]->Add(thishist);
        }
        thisfile->Close();
    }

    /** Go from counts in each histogram to efficiencies **/
    for (Trigger* trig : trigger_vec) {
        harr[trig->index]->Divide(harr[(trig->index)+numtrigs]);
    }
    for (Trigger* trig : trigger_vec) {
        bootstrap(trig, harr);
    }

    /** Plotting routines **/

    TCanvas* trig_canvas = new TCanvas("trig_canvas", "", 800, 600);
    for (Trigger* trig : trigger_vec) {
        gPad->SetLogx();
        gPad->SetTicks();
        int index = trig->index;
        thishist = harr[index];
        thishist->SetMarkerStyle(kDot);
        Color_t kColor = 46;
        thishist->SetMarkerColor(kColor);
        thishist->SetLineColor(kColor);
        thishist->SetMinimum(0);
        thishist->SetMaximum(1.2);
        thishist->GetYaxis()->SetTitleOffset(1.35);
        thishist->GetXaxis()->SetTickLength(0.02);
        thishist->GetYaxis()->SetTickLength(0.02);
        thishist->Draw("e1");
        myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", total_luminosity*1000));
        myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));
        trig_canvas->Draw();
        
        string hname = Form("%s_efficiency", trig->name.c_str());

        trig_canvas->SaveAs((plotPath + "triggerEfficiencies/" + hname + ".pdf").c_str());
        delete harr[index];
        trig_canvas->Clear();
    }
    delete trig_canvas;

    return;
}
