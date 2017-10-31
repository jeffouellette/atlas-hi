#include "../triggerUtil.C"

void triggers_pt_counts_hist(std::vector<int> thisRunNumbers) {

    initialize(0, false);
    const int numhists = numtrigs * numetabins;

    double ymin = 5e-7;
    double ymax = 5e7;
    const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    TH1D* harr[numhists];
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
    double this_trig_integrated_luminosity_vec[numetabins] = {};
    double integrated_luminosity_vec[32 * numetabins] = {};
    double total_luminosity = 0;
    TVectorD numtrigfirings(numtrigs*numpbins*numetabins);

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
    description->SetTextSize(0.026);

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
            if (!runPeriodA && thisRunNumber < 313500) continue;
            if (!runPeriodB && thisRunNumber > 313500) continue;
            TFile* thisfile = new TFile(Form("%srun_%i.root", trig_dir.c_str(), thisRunNumber), "READ");
            TVectorD* thisluminosityvec = (TVectorD*)thisfile->Get("lum_vec");
            if (index == 0) total_luminosity += (*thisluminosityvec)[0];

            double integral_deta = 0;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                thishist = (TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin));
                integral_deta += thishist->Integral();
                harr[index + ebin*numtrigs]->Add(thishist);
            }
            TVectorD* thisnumtrigfirings = (TVectorD*)thisfile->Get("trig_fire_vec");
            for (int ebin = 0; ebin < numetabins; ebin++) {
                if (integral_deta > 0) this_trig_integrated_luminosity_vec[ebin] += (*thisluminosityvec)[0];
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs] += (*(thisnumtrigfirings))[index + (pbin + ebin*numpbins)*numtrigs];
                }
            }
            thisfile->Close();
        }
        for (int ebin = 0; ebin < numetabins; ebin++) {
            integrated_luminosity_vec[index + ebin*numtrigs] = this_trig_integrated_luminosity_vec[ebin];
        }
    }


    // Calculate the best trigger to use for each bin, and be sure to scale by the correct deta, number of events, and luminosity.
    int best_bins[numpbins * numetabins] = {};
    for (int pbin = 0; pbin < numpbins; pbin++) {
        for (int ebin = 0; ebin < numetabins; ebin++) {
            double maxtrigfirings = 0;
            for (Trigger* trig : trigger_vec) {
                if (pbins[pbin] < trig->min_pt || etabins[ebin] < trig->lower_eta || trig->upper_eta < etabins[ebin+1]) continue;
                int index = trig->index;
                if (numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs] > maxtrigfirings) {
                    maxtrigfirings = numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs];
                    best_bins[pbin + ebin*numpbins] = index + ebin*numtrigs;
                }
            }
            if (ebin == 7) {
                int index = best_bins[pbin+ebin*numpbins] - ebin*numtrigs;
                cout << Form("For ebin 7, pbin %i, selecting trigger %i with effective luminosity %.3f [1/pb]", pbin, index, integrated_luminosity_vec[index + ebin*numtrigs]) << endl;
            }
        }
    }

    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        for (int ebin = 0; ebin < numetabins; ebin++) {
            deta = std::min(trig->upper_eta - trig->lower_eta, etabins[ebin+1]-etabins[ebin]);
            if (integrated_luminosity_vec[index+ebin*numtrigs] != 0) harr[index + ebin*numtrigs]->Scale(1 / (integrated_luminosity_vec[index + ebin*numtrigs] * deta), "width");
        }
    }

    TH1D* bestharr[numetabins];
    for (int ebin = 0; ebin < numetabins; ebin++) {
        bestharr[ebin] = new TH1D(Form("best_statistics_ebin%i", ebin), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
        thishist = bestharr[ebin];
        for (int pbin = 0; pbin < numpbins; pbin++) {
            thishist->SetBinContent(pbin+1, harr[best_bins[pbin + ebin*numpbins]]->GetBinContent(pbin+1));
            thishist->SetBinError(pbin+1, harr[best_bins[pbin + ebin*numpbins]]->GetBinError(pbin+1));
        }
    }

    double* bestharrscales = linspace(-1.5, 1.5, numetabins/2 - 1);

    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        if (numetabins > 1) {
            thishist->Scale(TMath::Power(10, bestharrscales[(int)(numetabins/2 - 0.5 -TMath::Abs(ebin - numetabins/2 + 0.5))]));
            thishist->SetMarkerStyle(mkstyles[ebin%7]);
            thishist->SetMarkerColor(mkcolors[(147*ebin)%20]);
            thishist->SetLineColor(mkcolors[(147*ebin)%20]);
            legend->AddEntry(thishist, Form("%g < #eta < %g (#times 10^{%g})", etabins[ebin], etabins[ebin+1], bestharrscales[(int)((0.5*(numetabins-1))-TMath::Abs(ebin-(0.5*(numetabins-1))))]));
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
    }
    for (int ebin = 0; ebin < numetabins; ebin++) {
        thishist = bestharr[ebin];
        thishist->GetXaxis()->SetTickLength(0.02);
        thishist->GetYaxis()->SetTickLength(0.02);
        thishist->Draw("same e1");
    }

    description->SetTextSize(0.036);
    if (numetabins == 1) description->DrawLatexNDC(0.42, 0.85, "#bf{#it{ATLAS}} p-Pb");
    else description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} p-Pb");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.50, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.78, 0.41, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", total_luminosity*1000)); 

    legend->Draw();
    trig_canvas->Draw();
    string hname;
    if (numetabins > 1) {
        hname = "../Plots/triggers/ptSpectra_combinedTriggers_etabinned";
    }
    else hname = "../Plots/triggers/ptSpectra_combinedTriggers_0eta490";

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

    if (runPeriodA || runPeriodB) trig_canvas->SaveAs((hname + ".pdf").c_str());


    /** Create number of trigger firings plot to illustrate procedure **/

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
    
    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} p-Pb");
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

    if (runPeriodA || runPeriodB) num_canvas->SaveAs((hname + ".pdf").c_str());
    
}
