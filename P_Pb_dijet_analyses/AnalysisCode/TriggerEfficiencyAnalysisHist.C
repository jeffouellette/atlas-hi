#include "../triggerUtil.C"

void bootstrap(Trigger* trig, TGraphAsymmErrors** effGraphArr) {
    if (trig->referenceTrigger->index == trig->index) { // iff trigger is the minbias reference trigger
        for (int pbin = 0; pbin < numpbins; pbin++) {
            effGraphArr[trig->referenceTrigger->index]->SetPoint(pbin, 0.5*(pbins[pbin+1]+pbins[pbin]), 1);
            effGraphArr[trig->referenceTrigger->index]->SetPointError(pbin, 0.5*(pbins[pbin+1]-pbins[pbin]), 0.5*(pbins[pbin+1]-pbins[pbin]), 0, 0);
        }
        trig->isBootstrapped = true;
    }

    if (trig->isBootstrapped) return;

    int index = trig->index;
    int referenceIndex = trig->referenceTrigger->index;
    if (!(trig->referenceTrigger->isBootstrapped)) {
        if (debugStatements) cout << "Status: In TriggerEfficiencyAnalysis.C (16): Bootstrapping trigger " << trig->referenceTrigger->name << " for " << trig->name << endl;
        bootstrap(trig->referenceTrigger, effGraphArr);
    }
    TGraphAsymmErrors* temp = effGraphArr[index];
    effGraphArr[index] = multiplyTGraphAsymmErrors(effGraphArr[index], effGraphArr[referenceIndex]);
    delete temp;
    trig->isBootstrapped = true;
    return;
}


void TriggerEfficiencyAnalysisHist() {

    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();
    const int numhists = numtrigs;
    const int numruns = (*thisRunNumbers).size();

    TEfficiency* effArr[numhists];
    TGraphAsymmErrors* effGraphArr[numhists];
    int pbin, ebin, index, referenceIndex;
    double reduced_pbins[40] = {};
    const int num_reduced_pbins = sizeof(reduced_pbins)/sizeof(reduced_pbins[0]) - 1;
    for (int reduced_pbin = 0; reduced_pbin < num_reduced_pbins+1; reduced_pbin++) {
        reduced_pbins[reduced_pbin] = pbins[reduced_pbin];
    }

    for (Trigger* trig : triggerVec) {
        index = trig->index;
        effArr[index] = new TEfficiency(Form("%s_t_efficiency", trig->name.c_str()), ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #Epsilon", numpbins, pbins);
        effGraphArr[index] = new TGraphAsymmErrors(numpbins);        
    }

    int thisRunNumber;
    double totalLuminosity = 0;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            if (debugStatements) cout << "Status: In TriggerEfficiencyAnalysisHist.C (41): Skipping run " << thisRunNumber << endl;
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
            if (!((trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber) || trig->referenceTrigger == trig)) continue;
            index = trig->index;
            effArr[index]->Add(*((TEfficiency*)thisFile->Get(Form("%s_t_efficiency_run%i", trig->name.c_str(), thisRunNumber))));
        }
        thisFile->Close();
    }

    /** Calculate efficiencies relative to reference trigger **/
    TGraphAsymmErrors* effGraph;
    for (Trigger* trig : triggerVec) {
        index = trig->index;
        effGraph = effGraphArr[index];
        for (pbin = 0; pbin < numpbins; pbin++) {
            //if (pbins[pbin] > trig->threshold_pt && effArr[index]->GetEfficiency(pbin) == 0) break;
            effGraph->SetPoint(pbin, 0.5*(pbins[pbin+1]+pbins[pbin]), effArr[index]->GetEfficiency(pbin));
            effGraph->SetPointError(pbin, 0.5*(pbins[pbin+1]-pbins[pbin]), 0.5*(pbins[pbin+1]-pbins[pbin]), effArr[index]->GetEfficiencyErrorLow(pbin), effArr[index]->GetEfficiencyErrorUp(pbin));
        }
    }

    for (Trigger* trig : triggerVec) {
        bootstrap(trig, effGraphArr);
    }


    /** Plotting routines **/

    TCanvas* trigCanvas;
    TFile* output = new TFile((effPath+"allEfficiencyHistograms.root").c_str(), "RECREATE");
    TLine* lineDrawer;
    TF1* fittedFunc;

    for (Trigger* trig : triggerVec) {
        if (trig->referenceTrigger == trig) continue; // skip the min bias trigger

        index = trig->index;

        trigCanvas = new TCanvas(Form("%s_canvas", trig->name.c_str()), "", 800, 600);
        trigCanvas->cd();
        trigCanvas->Draw();

        lineDrawer = new TLine();
        lineDrawer->SetLineColor(kBlack);
        lineDrawer->SetLineStyle(7);
        lineDrawer->SetLineWidth(lineDrawer->GetLineWidth()*2);

        gPad->SetLogx();
        gPad->SetTicks();
        effGraph = effGraphArr[index];

        Color_t kColor = 46;

        //fittedFunc = new TF1(Form("%s_erf", trig->name.c_str()), "TMath::Erf((x-[0])/[1])", pbins[0], pbins[numpbins]);
        if (fittedFunctionType == "fermi_dirac") {
            fittedFunc = new TF1(Form("%s_fermi_dirac", trig->name.c_str()), "1/(1+TMath::Exp(([0]-x)/[1]))", pbins[0], pbins[numpbins]);
            fittedFunc->SetParameters(trig->threshold_pt, 5); // guess parameters
        }
        else if (fittedFunctionType == "erf") {
            fittedFunc = new TF1(Form("%s_erf", trig->name.c_str()), "TMath::Erf([0]*(x-[1]))", pbins[0], pbins[numpbins]);
            fittedFunc->SetParameters(0.2, trig->threshold_pt);
        }
        effGraph->Fit(fittedFunc, "", "", (double)((trig->threshold_pt)-15), (double)((trig->threshold_pt)+75));

        effGraph->SetMarkerStyle(kDot);
        effGraph->SetMarkerColor(kColor);
        effGraph->SetLineColor(kColor);
        effGraph->SetMinimum(0);
        effGraph->SetMaximum(1.3);
        effGraph->GetYaxis()->SetTitleOffset(1.35);
        effGraph->GetXaxis()->SetTickLength(0.02);
        effGraph->GetYaxis()->SetTickLength(0.02);
        
        effGraph->Draw();
        myText (0.18, 0.91, kBlack, ("Trig: " + trig->name).c_str());
        myText (0.18, 0.85, kBlack, ("Ref: " + trig->referenceTrigger->name).c_str());

        lineDrawer->DrawLine(trig->threshold_pt, 0, trig->threshold_pt, 1.08);
        lineDrawer->DrawLine(trig->min_pt, 0, trig->min_pt, 1.08); // draws a line at the imposed pt cut
        lineDrawer->DrawLine(reduced_pbins[0], 1, reduced_pbins[num_reduced_pbins], 1); // draws a horizontal line at 1

        string histName = Form("%s_efficiency", trig->name.c_str());

        trigCanvas->SaveAs((plotPath + "triggerEfficiencies/" + histName + ".pdf").c_str());
        effGraph->Write(Form("%s_t_graph", trig->name.c_str()));
        delete effGraphArr[index];
        delete effArr[index];
        delete lineDrawer;
        delete trigCanvas;
    }
    output->Close();
    delete output;

    if (debugStatements) cout << "Status: In TriggerEfficiencyAnalysisHist.C (119): Finished plotting efficiency curves" << endl;
    return;
}
