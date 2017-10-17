#include "../triggerUtil.C"

void triggers_pt_counts_hist(std::vector<int> thisRunNumbers) {

    initialize(0, false);
    const int numhists = numtrigs * numetabins;

    const double ymin = 5e-7;
    const double ymax = 5e7;
    const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};
    //const Color_t mkcolors[8] = {kAzure-5, kOrange-5, kTeal-5, kPink-5, kSpring-5, kViolet-5, kGray+3, kRed-2};

    TH1D* harr[numhists];
    TCanvas* trig_canvas = new TCanvas("trig_canvas", "", 1000, 800);
    TLegend* legend;
    if (numetabins == 1) {
        legend = new TLegend(0.55, 0.55, 0.9, 0.9);
        legend->SetHeader("Trigger Selection", "C");
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
    double this_trig_integrated_luminosity_vec[numetabins] = {};
    double integrated_luminosity_vec[10 * numetabins] = {};
//    TVectorD numtrigfirings(numhists);
    TVectorD numtrigfirings(numtrigs*numpbins*numetabins);

    cout << "Starting loop over triggers..." << endl;

    // First combine trigger data from all runs into one histogram for each trigger. If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
    TH1D* thishist;
    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        for (int ebin = 0; ebin < numetabins; ebin++) {
            harr[index + ebin*numtrigs] = new TH1D(Form("trig%s_ebin%i", trig->name.c_str(), ebin), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
            harr[index + ebin*numtrigs]->Sumw2();
            this_trig_integrated_luminosity_vec[ebin] = 0;
        }
        for (int thisRunNumber : thisRunNumbers) {
            TFile* thisfile = new TFile(Form("../rootFiles/pt_data/trig_bin/run_%i.root", thisRunNumber), "READ");
            TVectorD* thisluminosityvec = (TVectorD*)thisfile->Get("lum_vec");

            double integral_deta = 0;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                thishist = (TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin));
                integral_deta += thishist->Integral();
            }

            for (int ebin = 0; ebin < numetabins; ebin++) {
                thishist = harr[index + ebin*numtrigs];
                thishist->Add((TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin)));
                if (integral_deta > 0) this_trig_integrated_luminosity_vec[ebin] += (*(TVectorD*)thisfile->Get("lum_vec")A)[0];
                TVectorD* thisnumtrigfirings = (TVectorD*)thisfile->Get("trig_fire_vec");
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs] += (*(thisnumtrigfirings))[index + (pbin + ebin*numpbins)*numtrigs]
                }
//                numtrigfirings[index + ebin*numtrigs] += (*(thisnumtrigfirings))[index];
            }
            thisfile->Close();
        }
        for (int ebin = 0; ebin < numetabins; ebin++) {
            integrated_luminosity_vec[index + ebin*numtrigs] = this_trig_integrated_luminosity_vec[ebin];
        }
    }


    // Calculate the best trigger to use for each bin, and be sure to scale by the correct deta, number of events, and luminosity.
    TH1D* bestharr[numetabins];
    for (int ebin = 0; ebin < numetabins; ebin++) {
        bestharr[ebin] = new TH1D(Form("best_statistics_ebin%i", ebin), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
    }
    int best_bins[numpbins * numetabins] = {};
    for (int pbin = 0; pbin < numpbins; pbin++) {
        for (int ebin = 0; ebin < numetabins; ebin++) {
            double maxtrigfirings = 0;
            for (Trigger* trig : trigger_vec) {
                if (trig->min_pt > pbins[pbin]) continue;
                int index = trig->index;
                if (numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs] > maxtrigfirings) {
                    maxtrigfirings = numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs];
                    best_bins[pbin + ebin*numpbins] = index + (pbin + ebin*numpbins)*numtrigs;
                }
            }
        }
    }

    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        for (int ebin = 0; ebin < numetabins; ebin++) {
            deta = std::min(trig->upper_eta - trig->lower_eta, etabins[ebin+1]-etabins[ebin]);
            harr[index + ebin*numtrigs]->Scale(1 / (integrated_luminosity_vec[index + ebin*numtrigs] * deta));
        }
    }

    for (int pbin = 0; pbin < numpbins; pbin++) {
        for (int ebin = 0; ebin < numetabins; ebin++) {
            bestharr[ebin]->SetBinContent(pbin+1, harr[best_bins[pbin + ebin*numpbins]]->GetBinContent(pbin+1));
            bestharr[ebin]->SetBinError(pbin+1, harr[best_bins[pbin + ebin*numpbins]]->GetBinError(pbin+1));
        }
    }

    double* bestharrscales;
    if (numetabins > 1) bestharrscales = linspace(-1.5, 1.5, numetabins/2 - 1);
    else bestharrscales = linspace(1, 1, 1);

    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        thishist->Scale(TMath::Power(10, bestharrscales[(int)(3.5-TMath::Abs(ebin-3.5))]));
        if (numetabins > 1) {
            thishist->SetMarkerStyle(mkstyles[ebin%7]);
            thishist->SetMarkerColor(mkcolors[(147*ebin)%20]);
            thishist->SetLineColor(mkcolors[(147*ebin)%20]);
            legend->AddEntry(thishist, Form("%g < #eta < %g (#times 10^{%g})", etabins[ebin], etabins[ebin+1], bestharrscales[(int)((0.5*(numetabins-1))-TMath::Abs(ebin-(0.5*(numetabins-1))))]));
        } else {
            thishist->SetMarkerStyle(4);
            thishist->SetMarkerColor(12);
            thishist->SetLineColor(12);
            legend->AddEntry(thishist, Form("Best statistics (#times 10^{%g})", bestharrscales[(int)((0.5*(numetabins-1))-TMath::Abs(ebin-(0.5*(numetabins-1))))]));
        }
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->GetYaxis()->SetTitleOffset(1.35);
    }

    if (numetabins == 1) {
        for (Trigger* trig : trigger_vec) {
            int index = trig->index;
            thishist = harr[index];
            cout << Form("Plotting trigger %s...", trig->name.c_str()) << endl;

            deta = trig->upper_eta - trig->lower_eta;
            if (1.6 < deta && deta < 1.8) hscale = 0.1;
            else if (1.1 < deta && deta < 1.3) hscale = 1;
            else hscale = 10;

            thishist->Scale(hscale);

            thishist->SetMarkerStyle(mkstyles[index%7]);
            thishist->SetMarkerColor(mkcolors[index%20]);
            thishist->SetLineColor(mkcolors[index%20]);
            thishist->SetMinimum(ymin);
            thishist->SetMaximum(ymax);
            thishist->GetYaxis()->SetTitleOffset(1.35);
            thishist->Draw("same e1");
            legend->AddEntry(thishist, Form("%s (#times %g)", trig->name.c_str(), hscale));
        }
    }
    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        thishist->Draw("same e1");
    }

    legend->Draw();        
    trig_canvas->Draw();
    if (numetabins > 1) {
        trig_canvas->SaveAs("../Plots/triggers/ptSpectra_combined_etabinned.pdf");
    }
    else trig_canvas->SaveAs("../Plots/triggers/ptSpectra_combinedTriggers.pdf");

}
