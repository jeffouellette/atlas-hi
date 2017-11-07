#include "../triggerUtil.C"

void triggers_pt_counts_hist() {

    initialize(0, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numruns = (*thisRunNumbers).size();
    cout << "Building trigger pt histograms with " << numruns << " runs being used" << endl;
    const int numhists = numtrigs * numruns * numetabins;

    double ymin = 5e-7;
    double ymax = 5e7;
    const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    TCanvas* trig_canvas = new TCanvas("trig_canvas", "", 1000, 800);
    TLegend* legend;
    if (numetabins == 1) {
        legend = new TLegend(0.55, 0.55, 0.9, 0.9);
        legend->SetTextSize(0.016);
    } else {
        legend = new TLegend(0.63, 0.6, 0.9, 0.9);
        legend->SetTextSize(0.019);
    }
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);

    double hscale, deta;
    TH1D* thishist;
    double total_lumi_vec[numruns * numpbins * numetabins * numtrigs];
    double kinematic_lumi_vec[numpbins*numetabins];
    double total_luminosity = 0;
    int best_bins[numruns * numpbins * numetabins];
    for (int n = 0; n < numtrigs*numruns*numpbins*numetabins; n++) {
        total_lumi_vec[n] = 0;
        if (n < numpbins*numetabins) kinematic_lumi_vec[n] = 0;
        if (n < numruns*numpbins*numetabins) best_bins[n] = 0;
    }

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
    description->SetTextSize(0.026);

    TH1D* bestharr[numetabins];
    for (int ebin = 0; ebin < numetabins; ebin++) {
        bestharr[ebin] = new TH1D(Form("best_statistics_ebin%i", ebin), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
    }

    TH1D* harr[numhists];
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        for (Trigger* trig : trigger_vec) {
            int index = trig->index; 
            for (int ebin = 0; ebin < numetabins; ebin++) {
                harr[rnIndex + (index + ebin*numtrigs)*numruns] = new TH1D(Form("run%i_trig%i_ebin%i", rnIndex, index, ebin), ";#it{p}_{T}^{jet} #left[Gev/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
                harr[rnIndex + (index + ebin*numtrigs)*numruns]->Sumw2();
            }
        }
    }

    cout << "Starting loop over triggers..." << endl;

    // First combine trigger data from all runs into one histogram for each trigger. If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
    int thisRunNumber;
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) {
            cout << "Skipping run " << thisRunNumber << endl;
            continue;
        }
        
        TFile* thisfile = new TFile(Form("%srun_%i.root", trigPath.c_str(), thisRunNumber), "READ");
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
                thishist = (TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin));
                integral_deta += thishist->Integral(); // Add the histogram's integral to find the total counts in this trigger over this run.
                harr[rnIndex + (index + ebin*numtrigs)*numruns]->Add(thishist); // Add the histogram over pbins at this ebin for this trigger in this run.
            }
            for (int ebin = 0; ebin < numetabins; ebin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) { // Set the trigger luminosity at this ebin in this run. If the trigger didn't fire, add 0, otherwise add the run luminosity.
                    if (integral_deta > 0) total_lumi_vec[rnIndex + (index + (pbin + ebin*numpbins)*numtrigs)*numruns] = (*thisluminosityvec)[0];
                    //else total_lumi_vec[rnIndex + (index + (pbin + ebin*numpbins)*numtrigs)*numruns] = 0;
                }
            }
        }
        TVectorD* thisnumtrigfirings = (TVectorD*)thisfile->Get("trig_fire_vec");
        // Find the bin with the maximum number of trigger firings.
        for (int pbin = 0; pbin < numpbins; pbin++) {
            for (int ebin = 0; ebin < numetabins; ebin++) {
                double maxtrigfirings = 0;

                for (Trigger* trig : trigger_vec) {
                    if (pbins[pbin] < trig->min_pt || etabins[ebin] < trig->lower_eta || trig->upper_eta < etabins[ebin+1]) continue;
                    int index = trig->index;
                    if ((*thisnumtrigfirings)[index + (pbin + ebin*numpbins)*numtrigs] > maxtrigfirings) {
                        maxtrigfirings = (*thisnumtrigfirings)[index + (pbin + ebin*numpbins)*numtrigs];
                        best_bins[rnIndex + (pbin + ebin*numpbins)*numruns] = index;
                    }
                }

                thishist = bestharr[ebin];
                int best_hist_index = rnIndex + (best_bins[rnIndex + (pbin + ebin*numpbins)*numruns] + ebin*numtrigs)*numruns;
                thishist->SetBinContent(pbin+1, thishist->GetBinContent(pbin+1) + harr[best_hist_index]->GetBinContent(pbin+1));
                thishist->SetBinError(pbin+1, TMath::Sqrt(TMath::Power(thishist->GetBinError(pbin+1), 2) + TMath::Power(harr[best_hist_index]->GetBinError(pbin+1), 2)));
                kinematic_lumi_vec[pbin + ebin*numpbins] += total_lumi_vec[rnIndex + (best_bins[rnIndex + (pbin + ebin*numpbins)*numruns] + (pbin + ebin*numpbins)*numtrigs)*numruns];
            }
        }
        thisfile->Close();
    }
    for (int rnIndex = 0; rnIndex < numruns; rnIndex++) {
        int thisRunNumber = (*thisRunNumbers)[rnIndex];
        if (skipRun(thisRunNumber)) continue;

        for (Trigger* trig : trigger_vec) {
            if (trig->iontrigger == thisRunNumber > 313500) continue;
            int index = trig->index;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                thishist = harr[rnIndex + (index + ebin*numtrigs)*numruns];
                deta = std::min((trig->upper_eta - trig->lower_eta), (etabins[ebin+1]-etabins[ebin]));
                thishist->Scale(1./deta, "width");
            }
        }
    }
    for (int pbin = 0; pbin < numpbins; pbin++) {
        cout << Form("For ebin 7, pbin %i, effective luminosity is %.3f [1/pb]", pbin, kinematic_lumi_vec[pbin + 7*numpbins]) << endl;
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



    /** Plotting routines **/

    // Plot best-selected pt spectra
    double* bestharrscales = linspace(-1.5, 1.5, numetabins/2 - 1);
    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        if (numetabins > 1) {
            thishist->Scale(TMath::Power(10, bestharrscales[(int)(numetabins/2 - 0.5 -TMath::Abs(ebin - numetabins/2 + 0.5))]));
            thishist->SetMarkerStyle(mkstyles[ebin%7]);
            thishist->SetMarkerColor(mkcolors[(147*ebin)%20]);
            thishist->SetLineColor(mkcolors[(147*ebin)%20]);
            legend->AddEntry(thishist, Form("%g < #eta_{lab} < %g (#times 10^{%g})", etabins[ebin], etabins[ebin+1], bestharrscales[(int)((0.5*(numetabins-1))-TMath::Abs(ebin-(0.5*(numetabins-1))))]));
        } else {
            thishist->Scale(10);
            thishist->SetMarkerStyle(4);
            thishist->SetMarkerColor(12);
            thishist->SetLineColor(12);
            legend->AddEntry(thishist, Form("Best statistics (#times 10)"));
        }
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->GetYaxis()->SetTitleOffset(1.35);
        thishist->GetXaxis()->SetTickLength(0.02);
        thishist->GetYaxis()->SetTickLength(0.02);
        thishist->Draw("same e1");
    }

/*    // Plot pt spectra for individual triggers - ILLUSTRATIVE
    if (numetabins == 1) {
        for (Trigger* trig : trigger_vec) {
            int index = trig->index;
            thishist = harr[index];
            cout << Form("Plotting trigger %s...", trig->name.c_str()) << endl;

            deta = trig->upper_eta - trig->lower_eta;
            if (1.6 < deta && deta < 1.8) hscale = 0.1;
            else if (1.1 < deta && deta < 1.3) hscale = 1;
            else hscale = 10;
//            hscale = 1;

            thishist->Scale(hscale);

            thishist->SetMarkerStyle(mkstyles[index%7]);
            thishist->SetMarkerColor(mkcolors[(index)%20]);
            thishist->SetLineColor(mkcolors[(index)%20]);
            thishist->SetMinimum(ymin);
            thishist->SetMaximum(ymax);
            thishist->GetYaxis()->SetTitleOffset(1.35);
            thishist->GetXaxis()->SetTickLength(0.02);
            thishist->GetYaxis()->SetTickLength(0.02);
            thishist->Draw("same e1");
            legend->AddEntry(thishist, Form("%s (#times %g)", trig->name.c_str(), hscale));
        }
    }*/

//    description->SetTextSize(0.036);
//    if (numetabins == 1) description->DrawLatexNDC(0.42, 0.85, "#bf{#it{ATLAS}} p-Pb");
//    else description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} p-Pb");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.56, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.47, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", total_luminosity*1000)); 

    legend->Draw();
    trig_canvas->Draw();
    string hname;
    if (numetabins > 1) {
        hname = "ptSpectra_combinedTriggers_etabinned";
    }
    else hname = "ptSpectra_combinedTriggers_0eta490";

    if (runPeriodA && !runPeriodB) {
        hname = hname + "_periodA";
        description->DrawTextNDC(0.48, 0.75, "Period A");
    }
    else if (!runPeriodA && runPeriodB) {
        hname = hname + "_periodB";
        description->DrawTextNDC(0.48, 0.75, "Period B");
    }
    else if (runPeriodA && runPeriodB) {
        description->DrawLatexNDC(0.48, 0.75, "Period A(-#eta) & B(#eta)");
    }

    if (runPeriodA || runPeriodB) trig_canvas->SaveAs((plotPath + hname + ".pdf").c_str());


    // DEPRECATED
    /** Create number of trigger firings plot to illustrate procedure - ILLUSTRATIVE **/
/*
    TCanvas* num_canvas = new TCanvas("num_canvas", "", 1000, 800);
    TH1D* numtrigfiringsharr[numtrigs];
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    TLegend* numlegend;
    numlegend = new TLegend(0.6, 0.55, 0.9, 0.9);
    numlegend->SetTextSize(0.016);
    ymin = 5e-6;
    ymax = 2e5;

    Trigger* trig;
    for (int index = 0; index < numtrigs; index++) {
        for (Trigger* t : trigger_vec) {
            if (t->index == index) {
                trig = t;
                break;
            }
        }
        deta = trig->upper_eta - trig->lower_eta;
        if (1.6 < deta && deta < 1.8) hscale = 0.1;
        else if (1.1 < deta && deta < 1.3) hscale = 1;
        else hscale = 10;
        hscale = 1;
        numtrigfiringsharr[index] = new TH1D(Form("trig%i_numcounts", index), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];dN^{trig}_{counts}/d#it{p}_{T}dy #left[(GeV/#it{c})^{-1}#right]", numpbins, pbins);
        thishist = numtrigfiringsharr[index];
        for (int pbin = 0; pbin < numpbins; pbin++) {
            double val = 0;
            for (int ebin = 0; ebin < numetabins; ebin++) { // integrate the number of counts over eta
                val += numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs];
            }
            thishist->SetBinContent(pbin+1, val);
            thishist->SetBinError(pbin+1, TMath::Sqrt(val));
        }
        thishist->Scale(hscale/deta, "width");
        thishist->SetMarkerStyle(mkstyles[index%7]);
        thishist->SetMarkerColor(mkcolors[(index)%20]);
        thishist->SetLineColor(mkcolors[(index)%20]);
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->GetYaxis()->SetTitleOffset(1.35);
        thishist->GetYaxis()->SetTickLength(0.02);
        thishist->Draw("same e1");
        numlegend->AddEntry(thishist, Form("%s", trig->name.c_str()));
    }
    TH1D* bestnumhist = new TH1D("numcounts_best", ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];dN^{trig}_{counts}/d#it{p}_{T}dy #left[(GeV/#it{c})^{-1}#right]", numpbins, pbins);
    for (int pbin = 0; pbin < numpbins; pbin++) {
        double maxval = 0;
        double maxerr = 0;
        for (int index = 0; index < numtrigs; index++) {
            if (maxval < numtrigfiringsharr[index]->GetBinContent(pbin+1)) {
                maxval = numtrigfiringsharr[index]->GetBinContent(pbin+1);
                maxerr = numtrigfiringsharr[index]->GetBinError(pbin+1);
            }
        }
        bestnumhist->SetBinContent(pbin+1, maxval);
        bestnumhist->SetBinError(pbin+1, maxerr);
    }    
    bestnumhist->SetMarkerStyle(4);
    bestnumhist->SetMarkerColor(12);
    bestnumhist->SetLineColor(12);
    bestnumhist->SetMinimum(ymin);
    bestnumhist->SetMaximum(ymax);
    num_canvas->SetTitle("Number of trigger firings");
    bestnumhist->Draw("same e1");
    numlegend->AddEntry(bestnumhist, "Best statistics");
    
//    description->SetTextSize(0.036);
//    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} p-Pb");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.50, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.41, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", total_luminosity*1000)); 
    num_canvas->Draw();
    numlegend->Draw();
    hname = "../Plots/triggers/trigger_fire_count";
    if (runPeriodA && !runPeriodB) {
        hname = hname + "_periodA";
        description->DrawTextNDC(0.48, 0.75, "Period A");
    }
    else if (!runPeriodA && runPeriodB) {
        hname = hname + "_periodB";
        description->DrawTextNDC(0.48, 0.75, "Period B");
    }
    else if (runPeriodA && runPeriodB) {
        description->DrawLatexNDC(0.48, 0.75, "Period A(-#eta) & B(#eta)");
    }

    if (runPeriodA || runPeriodB) num_canvas->SaveAs((hname + ".pdf").c_str());*/
    
}
