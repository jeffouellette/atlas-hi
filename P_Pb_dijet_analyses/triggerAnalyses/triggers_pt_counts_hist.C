#include "triggerUtil.C"

void triggers_pt_counts_hist(std::vector<int> thisRunNumbers) {

    initialize(0, false);
    const int numhists = numtrigs;

    const double harr_scales[3] = {1, 0.01, 0.0001};   // rescaling factors so the histograms don't all just overlap

    const double ymin = 1e-8;
    const double ymax = 5e6;
    const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
    const Color_t mkcolors[10] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48};
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
    
    cout << "Starting loop over triggers..." << endl;

    TH1D* thishist;
    for (Trigger* trig : trigger_vec) {
        int index = trig->index;
        integrated_luminosity = 0;
        
        this_trig_hist[index] = new TH1D(Form("trig%s", trig->name.c_str()), ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
        thishist = this_trig_hist[index];
        
        thishist->Sumw2();
        for (int thisRunNumber : thisRunNumbers) {
            TFile* thisfile = new TFile(Form("./pt_data/trig_bin/run_%i.root", thisRunNumber), "READ");
            TVectorD* thisluminosityvec = (TVectorD*)thisfile->Get("lum_vec");
            integrated_luminosity += (*(TVectorD*)thisfile->Get("lum_vec"))[0];
            thishist->Add((TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i", thisRunNumber, index)));
            thisfile->Close();
        }

        cout << Form("Plotting trigger %s...", trig->name.c_str()) << endl;

        deta = trig->upper_eta - trig->lower_eta;
        if (deta < 1.5) hscale = harr_scales[0];
        else if (deta < 4) hscale = harr_scales[1];
        else hscale = harr_scales[2];

        thishist->Scale(hscale / (integrated_luminosity));
        thishist->SetMarkerStyle(mkstyles[index%7]);
        thishist->SetMarkerColor(mkcolors[index%10]);
        thishist->SetLineColor(mkcolors[index%10]);
        thishist->SetMinimum(ymin);
        thishist->SetMaximum(ymax);
        thishist->Draw("same e1");
        legend->AddEntry(thishist, Form("%s (#times %g)", trig->name.c_str(), hscale));
    }
    legend->Draw();        
    trig_canvas->Draw();

    trig_canvas->SaveAs("./Plots/triggers/ptSpectra_combined.pdf");

}
