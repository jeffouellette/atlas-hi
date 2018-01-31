#include "../triggerUtil.C"

void triggers_pt_counts_hist() {

    initialize(0, true, false, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numruns = (*thisRunNumbers).size();
    const int numhists = numtrigs * numruns * numetabins;
    if (debugStatements) {
        cout << "Status: In triggers_pt_counts_hist.C (11): Building trigger pt histograms with " << numruns << " runs being used" << endl;
        cout << "Status: In triggers_pt_counts_hist.C (12): Numtrigs = " << numtrigs << endl;
        cout << "Status: In triggers_pt_counts_hist.C (13): Numetabins = " << numetabins << endl;
        cout << "Status: In triggers_pt_counts_hist.C (14): Numpbins = " << numpbins << endl;
        cout << "Status: In triggers_pt_counts_hist.C (15): ptPath = " << ptPath << endl;
    }

    double ymin = 5e-6;
    double ymax = 5e7;
    const Style_t mkstyles[2] = {kFullTriangleUp, kFullTriangleDown};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    double hscale, deta;
    TH1F* thisHist;

    if (debugStatements) cout << "Status: In triggers_pt_counts_hist.C (37): Initialized histArr histograms..." << endl;
    TH1F* histArr[numetabins];
    for (int etabin = 0; etabin < numetabins; etabin++) {
        histArr[etabin] = new TH1F(Form("best_statistics_etabin%i", etabin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb GeV^{-1}#right]", numpbins, pbins);
    }
    if (debugStatements) {
        cout << "Status: In triggers_pt_counts_hist.C (43): histArr histograms initialized." << endl;
        cout << "Status: In triggers_pt_counts_hist.C (44): Starting loop over triggers..." << endl;
    }

    vector<int>* runNumbers = getRunNumbers();

    /**** Fill summed histograms with results from event loops ****/
    {
        TSystemDirectory dir(ptPath.c_str(), ptPath.c_str());
        TList* sysfiles = dir.GetListOfFiles();
        if (sysfiles) {
            TSystemFile *sysfile;
            TString fname;
            TString histName;
            TIter next(sysfiles);

            int* rn;
            while ((sysfile=(TSystemFile*)next())) {
                fname = sysfile->GetName();
                if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
                    if (debugStatements) cout << "Status: In triggers_pt_counts.C (41): Found " << fname.Data() << endl; 
                    for (int thisRunNumber : *runNumbers) {
                        if (fname.Contains(to_string(thisRunNumber))) {
                            TFile* thisFile = new TFile(ptPath + fname, "READ");
                            int actetabin;
                            for (int etabin = 0; etabin < numetabins; etabin++) {
                                if (thisRunNumber < 313500) actetabin = numetabins - etabin - 1;
                                else actetabin = etabin;
                                histName = Form("trig_pt_counts_run%i_etabin%i", thisRunNumber, actetabin);
                                histArr[etabin]->Add((TH1F*)thisFile->Get(histName));
                            }
                            thisFile->Close();
                            delete thisFile;
                            break;
                        }
                    }
                }
            }
        }
    }
    /**** End fill summed histograms ****/

    // First combine trigger data from all runs into one histogram for each trigger. If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
    /*int thisRunNumber;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            cout << "Skipping run " << thisRunNumber << endl;
            continue;
        }
        
        TFile* thisFile = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "READ");
        TVectorD* thisLuminosityVec = (TVectorD*)thisFile->Get("lum_vec");
        TVectorD* thisRunVec = (TVectorD*)thisFile->Get("run_vec");

        assert (thisRunVec[0] == thisRunNumber); // Check that the data matches what we want to analyze.
        assert (thisRunVec[1] == numetabins);
        assert (thisRunVec[2] == numtrigs);

        totalLuminosity += (*thisLuminosityVec)[0];
        for (Trigger* trig : triggerVec) {
            double integral_deta = 0;
            int index = trig->index;
            for (int etabin = 0; etabin < numetabins; etabin++) {
                thisHist = (TH1F*)thisFile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, etabin));
                integral_deta += thisHist->Integral(); // Add the histogram's integral to find the total counts in this trigger over this run.
            }
            for (int etabin = 0; etabin < numetabins; etabin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) { // Set the trigger luminosity at this etabin in this run. If the trigger didn't fire, add 0, otherwise add the run luminosity.
                    if (integral_deta > 0 && trig->min_pt[etabin] <= pbins[pbin] && trig->lower_eta <= etabins[etabin] && etabins[etabin+1] <= trig->upper_eta) {
                        totalLumiVec[rnIndex + (index + (pbin + etabin*numpbins)*numtrigs)*numruns] = (*thisLuminosityVec)[0];
                    }
                }
            }
        }
        TVectorD* thisnumtrigfirings = (TVectorD*)thisFile->Get("trig_fire_vec");
        // Find the bin with the maximum number of trigger firings.
        for (int pbin = 0; pbin < numpbins; pbin++) {
            for (int etabin = 0; etabin < numetabins; etabin++) {
                double maxtrigfirings = 0;
//                int lastbestbin = 0;
                for (Trigger* trig : triggerVec) {
                    if (pbins[pbin] < trig->min_pt[etabin] || etabins[etabin] < trig->lower_eta || trig->upper_eta < etabins[etabin+1] || trig->disabled) continue;
                    int index = trig->index;
                    if ((*thisnumtrigfirings)[index + (pbin + etabin*numpbins)*numtrigs] >= maxtrigfirings) {
                        maxtrigfirings = (*thisnumtrigfirings)[index + (pbin + etabin*numpbins)*numtrigs];
//                        lastbestbin = bestBins[rnIndex + (pbin + etabin*numpbins)*numruns];
                        bestBins[rnIndex + (pbin + etabin*numpbins)*numruns] = index;
                    }
                }
//                bestBins[rnIndex + (pbin + etabin*numpbins)*numruns] = lastbestbin;
            }
        }
        thisFile->Close();
        delete thisFile;
    }
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        TFile* thisFile = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "READ");
        for (int etabin = 0; etabin < numetabins; etabin++) {
            thisHist = histArr[etabin];
            int actEtabin = etabin;
            if (thisRunNumber < 313603) actEtabin = numetabins-etabin-1; // Correct for period A kinematics (so bin 7 --> bin 0 if there are 8 bins, e.g.)
            for (int pbin = 0; pbin < numpbins; pbin++) {

                //int best_hist_index = rnIndex + (bestBins[rnIndex + (pbin + actEtabin*numpbins)*numruns] + actEtabin*numtrigs)*numruns;
                TH1F* hist_index = (TH1F*)thisFile->Get(Form("trig_pt_counts_run%i_trig%i_etabin%i", thisRunNumber, bestBins[rnIndex + (pbin + actEtabin*numpbins)*numruns], actEtabin));
                thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) + hist_index->GetBinContent(pbin+1));
                thisHist->SetBinError(pbin+1, TMath::Sqrt(TMath::Power(thisHist->GetBinError(pbin+1), 2) + TMath::Power(hist_index->GetBinError(pbin+1), 2)));
                kinematicLumiVec[pbin + etabin*numpbins] += totalLumiVec[rnIndex + (bestBins[rnIndex + (pbin + actEtabin*numpbins)*numruns] + (pbin + actEtabin*numpbins)*numtrigs)*numruns];
                 
            }
        }
        thisFile->Close();
        delete thisFile;
    }

    if (debugStatements) {
        for (int pbin = 0; pbin < numpbins; pbin++) {
            cout << Form("For etabin 4, pbin %i, effective luminosity is %.3f [1/nb]", pbin, kinematicLumiVec[pbin + 4*numpbins]*1000) << endl;
        }
    }

    // Scale best-selected histograms by first deta, then by the best integrated luminosity for that kinematic bin.
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thisHist = histArr[etabin];
        thisHist->Scale(1/(etabins[etabin+1] - etabins[etabin]), "width");
        for (int pbin = 0; pbin < numpbins; pbin++) {
            double lumi = kinematicLumiVec[pbin + etabin*numpbins];
            if (lumi > 0) {
                thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1)/lumi);
                thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1)/lumi);
            }
        }
    }
    for (Trigger* trig : triggerVec) {
        int index = trig->index;
        for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
            if (skipRun((*thisRunNumbers)[rnIndex])) continue;
            for (int etabin = 0; etabin < numetabins; etabin++) {
                int actEtabin = etabin;
                if (rnIndex < 313500) actEtabin = numetabins - etabin - 1;
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    int actEtabin = 0;
                    lumiVecIntegratedRuns[index + (pbin + etabin*numpbins)*numtrigs] += totalLumiVec[rnIndex + (index + (pbin + actEtabin*numpbins)*numtrigs)*numruns];
                }
            }
        }
    }
    delete [] totalLumiVec;
    */


    /** Plotting routines **/

    // Plot best-selected pt spectra
    double* histArrScales = linspace(-1.5, 1.5, numetabins/2 - 1);
//    double* histArrScales = linspace(0, 0, numetabins/2 - 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts
    TCanvas* trigCanvas = new TCanvas("trigCanvas", "", 800, 600);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicks();
    trigCanvas->Draw();
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thisHist = histArr[etabin];
        for (int pbin = 0; pbin < numpbins; pbin++) {
            if (kinematicLumiVec[pbin + etabin*numpbins] != 0) {
                thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / (kinematicLumiVec[pbin + etabin*numpbins]));
                thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / (kinematicLumiVec[pbin + etabin*numpbins]));
            }
        }
        thisHist->Scale(TMath::Power(10, histArrScales[(int)(numetabins/2 - 0.5 -TMath::Abs(etabin - numetabins/2 + 0.5))]), "width"); // separate different etabins
        thisHist->SetMarkerStyle(mkstyles[etabin < (numetabins/2)]);
//        thisHist->SetMarkerStyle(kDot);
        Color_t kColor = mkcolors[(147*etabin)%20];
        thisHist->SetMarkerColor(kColor);
        thisHist->SetLineColor(kColor);
        thisHist->SetMinimum(ymin);
        thisHist->SetMaximum(ymax);
        thisHist->GetYaxis()->SetTitleOffset(1.35);
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);
        if (etabin == 0) thisHist->Draw("e1");
        else thisHist->Draw("same e1");
        
        const float textx = 0.46 + (etabin>=(numetabins/2))*0.26;
        const float texty = 0.91 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
        const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[etabin], etabins[etabin+1], histArrScales[(int)((0.5*(numetabins-1))-TMath::Abs(etabin-(0.5*(numetabins-1))))]);
        myMarkerText (textx, texty, kColor, mkstyles[etabin < (numetabins/2)], text);
    }

    myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity));
    myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

    string histName;
    if (numetabins > 1) {
        histName = "ptSpectra_combinedTriggers_etabinned";
    }
    else histName = "ptSpectra_combinedTriggers_0eta490";

    if (runPeriodA && !runPeriodB) {
        myText (0.19, 0.33, kBlack, "Period A");
        histName = histName + "_periodA";
    }
    else if (!runPeriodA && runPeriodB) {
        myText (0.19, 0.33, kBlack, "Period B");
        histName = histName + "_periodB";
    }
    else {
        myText (0.19, 0.33, kBlack, "Period A & B");
    }

    if (runPeriodA || runPeriodB) trigCanvas->SaveAs((plotPath + "ptSpectra/" + histName + ".pdf").c_str());
    for (int etabin = 0; etabin < numetabins; etabin++) {
        delete histArr[etabin];
    }
    delete trigCanvas;


    // For each etabin, plot the trigger pt spectra.
/*    for (int etabin = 0; etabin < numetabins; etabin++) {
        TCanvas* allTrigsCanvas = new TCanvas(Form("allTrigsCanvas_%i", etabin), "", 800, 600);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetTicks();
        TH1F* allTrigHistArr[numtrigs];
        for (Trigger* trig : triggerVec) { 
            if (trig->disabled) continue;
            int index = trig->index;
            allTrigHistArr[index] = new TH1F(Form("all_trig%i_etabin%i", index, etabin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb GeV^{-1}#right]", numpbins, pbins);
            allTrigHistArr[index]->Sumw2();
            for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
                thisRunNumber = (*thisRunNumbers)[rnIndex];
                if (skipRun(thisRunNumber)) continue;
                TFile* thisFile = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "READ");
                int actEtabin = etabin;
                if (thisRunNumber < 313500) actEtabin = numetabins - etabin - 1;
                allTrigHistArr[index]->Add((TH1F*)thisFile->Get(Form("trig_pt_counts_run%i_trig%i_etabin%i", thisRunNumber, index, actEtabin)));
                thisFile->Close();
                delete thisFile;
            }
        }

        float legend_height = 0.91;
        float legend_ion_height = 0.91;
        for (Trigger* trig : triggerVec) {
            if (trig->disabled) continue;
            int index = trig->index;
            thisHist = allTrigHistArr[index];
            if (thisHist->Integral() == 0) continue;
            thisHist->SetMarkerStyle(kDot);
            //Color_t kColor = (trig->iontrigger)*(kOrange+10) + (!trig->iontrigger)*(kAzure+10);
            Color_t kColor = mkcolors[index%20];
            thisHist->SetMarkerColor(kColor);
            thisHist->SetLineColor(kColor);
            for (int pbin = 0; pbin < numpbins; pbin++) { // scale each bin by luminosity
                double thislumi = lumiVecIntegratedRuns[index + (pbin + etabin*numpbins)*numtrigs];
                if (thislumi <= 0) continue;
                thisHist->SetBinContent(pbin+1, (thisHist->GetBinContent(pbin+1))/thislumi);
                thisHist->SetBinError(pbin+1, (thisHist->GetBinError(pbin+1))/thislumi);
            }
            thisHist->Scale(1/(etabins[etabin+1]-etabins[etabin]), "width");
            thisHist->SetMaximum(ymax);
            thisHist->SetMinimum(ymin);
            if (legend_height == 0.91) thisHist->Draw("e1");
            else thisHist->Draw("same e1");
            float textx = 0.56 + (!trig->iontrigger)*0.2;
            float texty = legend_ion_height*(trig->iontrigger) + legend_height*(!trig->iontrigger);
            myMarkerText (textx, texty, kColor, kFullCircle, trig->name.c_str(), 0.75, 0.016);
            legend_height -= 0.02*(!trig->iontrigger);
            legend_ion_height -= 0.02*(trig->iontrigger);
        }
        myText (0.19, 0.33, kBlack, Form("%.1f < #it{#eta} < %.1f", etabins[etabin], etabins[etabin+1]));
        myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity*1000));
        myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

        if (etabins[etabin] < 0)  histName = Form("ptSpectra_combinedTriggers_n%ieta%i", (int)(-etabins[etabin]*100), (int)(TMath::Abs(etabins[etabin+1])*100));
        else histName = Form("ptSpectra_combinedTriggers_p%ieta%i", (int)(etabins[etabin]*100), (int)(etabins[etabin+1]*100));

        if (runPeriodA && !runPeriodB) {
            myText (0.19, 0.39, kBlack, "Period A");
            histName = histName + "_periodA";
        }
        else if (!runPeriodA && runPeriodB) {
            myText (0.19, 0.39, kBlack, "Period B");
            histName = histName + "_periodB";
        }
        else {
            myText (0.19, 0.39, kBlack, "Period A & B");
        }
        allTrigsCanvas->Draw();
        allTrigsCanvas->SaveAs((plotPath + "trigPtSpectra/" + histName + ".pdf").c_str());
        delete allTrigsCanvas;
        for (Trigger* trig : triggerVec) {
            if (trig->disabled) continue;
            delete allTrigHistArr[trig->index];
        }
    }*/
    return;
}
