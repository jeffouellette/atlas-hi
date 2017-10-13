#include "../triggerUtil.C"

void triggers_pt_counts_hist(std::vector<int> thisRunNumbers) {

    initialize(0, false);
    const int numhists = numtrigs;

    const double ymin = 1e-8;
    const double ymax = 1e8;
    const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};
    //const Color_t mkcolors[8] = {kAzure-5, kOrange-5, kTeal-5, kPink-5, kSpring-5, kViolet-5, kGray+3, kRed-2};

    TH1D* this_trig_hist[numtrigs];
    TCanvas* trig_canvas = new TCanvas("trig_canvas", "", 1000, 800);
    TLegend* legend = new TLegend(0.55, 0.55, 0.9, 0.9);
    legend->SetTextSize(0.024);
    legend->SetHeader("Trigger fired", "C");
    legend->SetTextSize(0.016);
    gPad->SetLogy();
    gPad->SetLogx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);

    double integrated_luminosity, hscale, deta;
    double integrated_luminosity_vec[10] = {};
    TVectorD numtrigfirings(numtrigs);

    cout << "Starting loop over triggers..." << endl;

    // First combine trigger data from all runs into one histogram for each trigger. If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
    TH1D* thishist;
    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        integrated_luminosity = 0;
        
        this_trig_hist[index] = new TH1D(Form("trig%s", trig->name.c_str()), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
        thishist = this_trig_hist[index];
        
        thishist->Sumw2();
        for (int thisRunNumber : thisRunNumbers) {
            TFile* thisfile = new TFile(Form("../rootFiles/pt_data/trig_bin/run_%i.root", thisRunNumber), "READ");
            TVectorD* thisluminosityvec = (TVectorD*)thisfile->Get("lum_vec");
            thishist->Add((TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i", thisRunNumber, index)));
            if (thishist->Integral() > 0) integrated_luminosity += (*(TVectorD*)thisfile->Get("lum_vec"))[0];
            numtrigfirings[index] += (*((TVectorD*)thisfile->Get("trig_fire_vec")))[index];
            thisfile->Close();
        }
        integrated_luminosity_vec[index] = integrated_luminosity;
    }


    // Calculate the best trigger to use for each bin, and be sure to scale by the correct deta, number of events, and luminosity.
    TH1D* besthist = new TH1D("best_statistics", ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];1/N d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
    //besthist->Sumw2();
    int best_bins[numpbins] = {};
    for (int pbin = 0; pbin < numpbins; pbin++) {
        double maxtrigfirings = 0;
        for (Trigger* trig : trigger_vec) {
            if (trig->min_pt >= pbins[pbin]) continue;
            int index = trig->index;
            if (numtrigfirings[index] > maxtrigfirings) {
                maxtrigfirings = numtrigfirings[index];
                best_bins[pbin] = index;
            }
        }
    }


    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        thishist = this_trig_hist[index];

        cout << Form("Plotting trigger %s...", trig->name.c_str()) << endl;

        deta = trig->upper_eta - trig->lower_eta;
        if (deta == 1.7) hscale = 100;
        else if (deta == 1.2) hscale = 1;
        else hscale = 0.01;

        hscale=1;

        thishist->Scale(hscale / (/*numtrigfirings[index] */ integrated_luminosity_vec[index] * deta));
        thishist->SetMarkerStyle(mkstyles[index%7]);
        thishist->SetMarkerColor(mkcolors[index%20]);
        thishist->SetLineColor(mkcolors[index%20]);
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->GetYaxis()->SetTitleOffset(1.35);
        thishist->Draw("same e1");
        legend->AddEntry(thishist, Form("%s (#times %g)", trig->name.c_str(), hscale));
    }

    for (int pbin = 0; pbin < numpbins; pbin++) {
        besthist->SetBinContent(pbin+1, this_trig_hist[best_bins[pbin]]->GetBinContent(pbin+1));
        besthist->SetBinError(pbin+1, this_trig_hist[best_bins[pbin]]->GetBinError(pbin+1));
    }

    besthist->SetMarkerStyle(4);
    besthist->SetMarkerColor(12);
    besthist->SetLineColor(12);
    besthist->SetMinimum(ymin);
    besthist->SetMaximum(ymax);
    besthist->GetYaxis()->SetTitleOffset(1.35);
    besthist->Draw("same e1");
    legend->AddEntry(besthist, "Combined triggers");

    legend->Draw();        
    trig_canvas->Draw();

    trig_canvas->SaveAs("../Plots/triggers/ptSpectra_combined.pdf");

}
