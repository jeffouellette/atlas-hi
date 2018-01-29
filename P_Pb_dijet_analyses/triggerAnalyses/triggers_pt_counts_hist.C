#include "../triggerUtil.C"

void triggers_pt_counts_hist() {

    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numruns = (*thisRunNumbers).size();
    const int numhists = numtrigs * numruns * numetabins;
    if (debugStatements) {
        cout << "Building trigger pt histograms with " << numruns << " runs being used" << endl;
        cout << "Numtrigs = " << numtrigs << endl;
        cout << "Numetabins = " << numetabins << endl;
        cout << "Numpbins = " << numpbins << endl;
        cout << "ptPath = " << ptPath << endl;
    }

    double ymin = 5e-6;
    double ymax = 5e7;
    const Style_t mkstyles[2] = {kFullTriangleUp, kFullTriangleDown};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    double hscale, deta;
    TH1F* thishist;
    float* total_lumi_vec = new float[numruns * numpbins * numetabins * numtrigs];
    float* lumi_vec_integrated_runs = new float[numpbins * numetabins * numtrigs];
    float* kinematic_lumi_vec = new float[numpbins * numetabins];
    float total_luminosity = 0;
    int best_bins[numruns * numpbins * numetabins];
    for (int n = 0; n < numtrigs*numruns*numpbins*numetabins; n++) {
        total_lumi_vec[n] = 0;
        if (n < numpbins*numetabins*numtrigs) lumi_vec_integrated_runs[n] = 0;
        if (n < numpbins*numetabins) kinematic_lumi_vec[n] = 0;
        if (n < numruns*numpbins*numetabins) best_bins[n] = 0;
    }

    if (debugStatements) cout << "Initialized bestharr histograms..." << endl;
    TH1F* bestharr[numetabins];
    for (int ebin = 0; ebin < numetabins; ebin++) {
        bestharr[ebin] = new TH1F(Form("best_statistics_ebin%i", ebin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb GeV^{-1}#right]", numpbins, pbins);
    }
    if (debugStatements) cout << "bestharr histograms initialized." << endl;

    cout << "Starting loop over triggers..." << endl;

    // First combine trigger data from all runs into one histogram for each trigger. If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
    int thisRunNumber;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            cout << "Skipping run " << thisRunNumber << endl;
            continue;
        }
        
        TFile* thisfile = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "READ");
        TVectorD* thisluminosityvec = (TVectorD*)thisfile->Get("lum_vec");
        TVectorD* thisrunvec = (TVectorD*)thisfile->Get("run_vec");

        assert (thisrunvec[0] == thisRunNumber); // Check that the data matches what we want to analyze.
        assert (thisrunvec[1] == numetabins);
        assert (thisrunvec[2] == numtrigs);

        total_luminosity += (*thisluminosityvec)[0];
        for (Trigger* trig : trigger_vec) {
            double integral_deta = 0;
            int index = trig->index;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                thishist = (TH1F*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin));
                integral_deta += thishist->Integral(); // Add the histogram's integral to find the total counts in this trigger over this run.
            }
            for (int ebin = 0; ebin < numetabins; ebin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) { // Set the trigger luminosity at this ebin in this run. If the trigger didn't fire, add 0, otherwise add the run luminosity.
                    if (integral_deta > 0 && trig->min_pt[ebin] <= pbins[pbin] && trig->lower_eta <= etabins[ebin] && etabins[ebin+1] <= trig->upper_eta) {
                        total_lumi_vec[rnIndex + (index + (pbin + ebin*numpbins)*numtrigs)*numruns] = (*thisluminosityvec)[0];
                    }
                }
            }
        }
        TVectorD* thisnumtrigfirings = (TVectorD*)thisfile->Get("trig_fire_vec");
        // Find the bin with the maximum number of trigger firings.
        for (int pbin = 0; pbin < numpbins; pbin++) {
            for (int ebin = 0; ebin < numetabins; ebin++) {
                double maxtrigfirings = 0;
//                int lastbestbin = 0;
                for (Trigger* trig : trigger_vec) {
                    if (pbins[pbin] < trig->min_pt[ebin] || etabins[ebin] < trig->lower_eta || trig->upper_eta < etabins[ebin+1] || trig->disabled) continue;
                    int index = trig->index;
                    if ((*thisnumtrigfirings)[index + (pbin + ebin*numpbins)*numtrigs] >= maxtrigfirings) {
                        maxtrigfirings = (*thisnumtrigfirings)[index + (pbin + ebin*numpbins)*numtrigs];
//                        lastbestbin = best_bins[rnIndex + (pbin + ebin*numpbins)*numruns];
                        best_bins[rnIndex + (pbin + ebin*numpbins)*numruns] = index;
                    }
                }
//                best_bins[rnIndex + (pbin + ebin*numpbins)*numruns] = lastbestbin;
            }
        }
        thisfile->Close();
        delete thisfile;
    }
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        TFile* thisfile = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "READ");
        for (int ebin = 0; ebin < numetabins; ebin++) {
            thishist = bestharr[ebin];
            int act_ebin = ebin;
            if (thisRunNumber < 313603) act_ebin = numetabins-ebin-1; // Correct for period A kinematics (so bin 7 --> bin 0 if there are 8 bins, e.g.)
            for (int pbin = 0; pbin < numpbins; pbin++) {

                //int best_hist_index = rnIndex + (best_bins[rnIndex + (pbin + act_ebin*numpbins)*numruns] + act_ebin*numtrigs)*numruns;
                TH1F* hist_index = (TH1F*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, best_bins[rnIndex + (pbin + act_ebin*numpbins)*numruns], act_ebin));
                thishist->SetBinContent(pbin+1, thishist->GetBinContent(pbin+1) + hist_index->GetBinContent(pbin+1));
                thishist->SetBinError(pbin+1, TMath::Sqrt(TMath::Power(thishist->GetBinError(pbin+1), 2) + TMath::Power(hist_index->GetBinError(pbin+1), 2)));
                kinematic_lumi_vec[pbin + ebin*numpbins] += total_lumi_vec[rnIndex + (best_bins[rnIndex + (pbin + act_ebin*numpbins)*numruns] + (pbin + act_ebin*numpbins)*numtrigs)*numruns];
                 
            }
        }
        thisfile->Close();
        delete thisfile;
    }

    if (debugStatements) {
        for (int pbin = 0; pbin < numpbins; pbin++) {
            cout << Form("For ebin 4, pbin %i, effective luminosity is %.3f [1/nb]", pbin, kinematic_lumi_vec[pbin + 4*numpbins]*1000) << endl;
        }
    }

    // Scale best-selected histograms by first deta, then by the best integrated luminosity for that kinematic bin.
    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        thishist->Scale(1/(etabins[ebin+1] - etabins[ebin]), "width");
        for (int pbin = 0; pbin < numpbins; pbin++) {
            double lumi = kinematic_lumi_vec[pbin + ebin*numpbins];
            if (lumi > 0) {
                thishist->SetBinContent(pbin+1, thishist->GetBinContent(pbin+1)/lumi);
                thishist->SetBinError(pbin+1, thishist->GetBinError(pbin+1)/lumi);
            }
        }
    }
    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
            if (skipRun((*thisRunNumbers)[rnIndex])) continue;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                int act_ebin = ebin;
                if (rnIndex < 313500) act_ebin = numetabins - ebin - 1;
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    int act_ebin = 0;
                    lumi_vec_integrated_runs[index + (pbin + ebin*numpbins)*numtrigs] += total_lumi_vec[rnIndex + (index + (pbin + act_ebin*numpbins)*numtrigs)*numruns];
                }
            }
        }
    }
    delete [] total_lumi_vec;



    /** Plotting routines **/

    // Plot best-selected pt spectra
    double* bestharrscales = linspace(-1.5, 1.5, numetabins/2 - 1);
//    double* bestharrscales = linspace(0, 0, numetabins/2 - 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts
    TCanvas* trig_canvas = new TCanvas("trig_canvas", "", 800, 600);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicks();
    trig_canvas->Draw();
    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        thishist->Scale(TMath::Power(10, bestharrscales[(int)(numetabins/2 - 0.5 -TMath::Abs(ebin - numetabins/2 + 0.5))]));
        thishist->SetMarkerStyle(mkstyles[ebin < (numetabins/2)]);
//        thishist->SetMarkerStyle(kDot);
        Color_t kColor = mkcolors[(147*ebin)%20];
        thishist->SetMarkerColor(kColor);
        thishist->SetLineColor(kColor);
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->GetYaxis()->SetTitleOffset(1.35);
        thishist->GetXaxis()->SetTickLength(0.02);
        thishist->GetYaxis()->SetTickLength(0.02);
        if (ebin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        
        const float textx = 0.46 + (ebin>=(numetabins/2))*0.26;
        const float texty = 0.91 - (ebin%(numetabins/2))*0.05*(ebin>=(numetabins/2)) - (numetabins/2 - ebin - 1)*0.05*(ebin<(numetabins/2));
        const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[ebin], etabins[ebin+1], bestharrscales[(int)((0.5*(numetabins-1))-TMath::Abs(ebin-(0.5*(numetabins-1))))]);
        myMarkerText (textx, texty, kColor, mkstyles[ebin < (numetabins/2)], text);
    }

    myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", total_luminosity*1000));
    myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

    string hname;
    if (numetabins > 1) {
        hname = "ptSpectra_combinedTriggers_etabinned";
    }
    else hname = "ptSpectra_combinedTriggers_0eta490";

    if (runPeriodA && !runPeriodB) {
        myText (0.19, 0.33, kBlack, "Period A");
        hname = hname + "_periodA";
    }
    else if (!runPeriodA && runPeriodB) {
        myText (0.19, 0.33, kBlack, "Period B");
        hname = hname + "_periodB";
    }
    else {
        myText (0.19, 0.33, kBlack, "Period A & B");
    }

    if (runPeriodA || runPeriodB) trig_canvas->SaveAs((plotPath + "ptSpectra/" + hname + ".pdf").c_str());
    for (int ebin = 0; ebin < numetabins; ebin++) {
        delete bestharr[ebin];
    }
    delete trig_canvas;


    // For each etabin, plot the trigger pt spectra.
    for (int ebin = 0; ebin < numetabins; ebin++) {
        TCanvas* all_trigs_canvas = new TCanvas(Form("all_trigs_canvas_%i", ebin), "", 800, 600);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetTicks();
        TH1F* all_trig_harr[numtrigs];
        for (Trigger* trig : trigger_vec) { 
            if (trig->disabled) continue;
            int index = trig->index;
            all_trig_harr[index] = new TH1F(Form("all_trig%i_ebin%i", index, ebin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb GeV^{-1}#right]", numpbins, pbins);
            all_trig_harr[index]->Sumw2();
            for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
                thisRunNumber = (*thisRunNumbers)[rnIndex];
                if (skipRun(thisRunNumber)) continue;
                TFile* thisfile = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "READ");
                int act_ebin = ebin;
                if (thisRunNumber < 313500) act_ebin = numetabins - ebin - 1;
                all_trig_harr[index]->Add((TH1F*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, act_ebin)));
                thisfile->Close();
                delete thisfile;
            }
        }

        float legend_height = 0.91;
        float legend_ion_height = 0.91;
        for (Trigger* trig : trigger_vec) {
            if (trig->disabled) continue;
            int index = trig->index;
            thishist = all_trig_harr[index];
            if (thishist->Integral() == 0) continue;
            thishist->SetMarkerStyle(kDot);
            //Color_t kColor = (trig->iontrigger)*(kOrange+10) + (!trig->iontrigger)*(kAzure+10);
            Color_t kColor = mkcolors[index%20];
            thishist->SetMarkerColor(kColor);
            thishist->SetLineColor(kColor);
            for (int pbin = 0; pbin < numpbins; pbin++) { // scale each bin by luminosity
                double thislumi = lumi_vec_integrated_runs[index + (pbin + ebin*numpbins)*numtrigs];
                if (thislumi <= 0) continue;
                thishist->SetBinContent(pbin+1, (thishist->GetBinContent(pbin+1))/thislumi);
                thishist->SetBinError(pbin+1, (thishist->GetBinError(pbin+1))/thislumi);
            }
            thishist->Scale(1/(etabins[ebin+1]-etabins[ebin]), "width");
            thishist->SetMaximum(ymax);
            thishist->SetMinimum(ymin);
            if (legend_height == 0.91) thishist->Draw("e1");
            else thishist->Draw("same e1");
            float textx = 0.56 + (!trig->iontrigger)*0.2;
            float texty = legend_ion_height*(trig->iontrigger) + legend_height*(!trig->iontrigger);
            myMarkerText (textx, texty, kColor, kFullCircle, trig->name.c_str(), 0.75, 0.016);
            legend_height -= 0.02*(!trig->iontrigger);
            legend_ion_height -= 0.02*(trig->iontrigger);
        }
        myText (0.19, 0.33, kBlack, Form("%.1f < #it{#eta} < %.1f", etabins[ebin], etabins[ebin+1]));
        myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", total_luminosity*1000));
        myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

        if (etabins[ebin] < 0)  hname = Form("ptSpectra_combinedTriggers_n%ieta%i", (int)(-etabins[ebin]*100), (int)(TMath::Abs(etabins[ebin+1])*100));
        else hname = Form("ptSpectra_combinedTriggers_p%ieta%i", (int)(etabins[ebin]*100), (int)(etabins[ebin+1]*100));

        if (runPeriodA && !runPeriodB) {
            myText (0.19, 0.39, kBlack, "Period A");
            hname = hname + "_periodA";
        }
        else if (!runPeriodA && runPeriodB) {
            myText (0.19, 0.39, kBlack, "Period B");
            hname = hname + "_periodB";
        }
        else {
            myText (0.19, 0.39, kBlack, "Period A & B");
        }
        all_trigs_canvas->Draw();
        all_trigs_canvas->SaveAs((plotPath + "trigPtSpectra/" + hname + ".pdf").c_str());
        delete all_trigs_canvas;
        for (Trigger* trig : trigger_vec) {
            if (trig->disabled) continue;
            delete all_trig_harr[trig->index];
        }
    }
    return;
}
