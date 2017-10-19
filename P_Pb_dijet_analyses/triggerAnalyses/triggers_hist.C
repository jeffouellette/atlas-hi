#include "../triggerUtil.C"

void triggers_hist(int thisRunNumber) {

    initialize(thisRunNumber, false);
    const double* trigbins = linspace(-0.5, numtrigs+0.5, numtrigs+1);
    const double* ybins = linspace(-0.5, numpbins+0.5, numpbins+1);

    double minval = 1e0;
    double maxval;

    TFile* thisfile = new TFile(Form("../rootFiles/trig_data/run_%i.root", thisRunNumber), "READ");
    TH1D* harr[numtrigs];
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        harr[i] = (TH1D*)thisfile->Get(Form("trig%i", i));
    }
    
    TCanvas* c1 = new TCanvas("c1", "", 1000, 800);   

    gPad->SetLogy();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    
    const Style_t mkstyles[8] = {kDot, kDot, kDot, kDot, kDot, kDot, kDot, kDot};
    const Color_t mkcolors[20] = {30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49};

    TLatex* description = new TLatex();
    description->SetTextAlign(22);
    description->SetTextFont(42);
    description->SetTextSize(0.026);

    double maxbincontent = 0;
    double thisbincontent;
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        thisbincontent = harr[i]->GetBinContent(i+1);
        if (maxbincontent < thisbincontent) maxbincontent = thisbincontent;
    }
    maxval = 6 * maxbincontent;

    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        harr[i]->SetMarkerStyle(mkstyles[i%8]);
        harr[i]->SetMarkerColor(mkcolors[i%20]);
        harr[i]->SetLineColor(mkcolors[i%20]);
        harr[i]->SetFillColor(mkcolors[i%20]);
        harr[i]->SetMinimum(minval);
        harr[i]->SetMaximum(maxval);
        harr[i]->GetXaxis()->SetLabelOffset(999);
        harr[i]->GetXaxis()->SetLabelSize(0);
        harr[i]->GetXaxis()->SetTickLength(0);

        harr[i]->GetYaxis()->SetTickLength(0.015);

        harr[i]->Draw("BAR, SAME, E2");

        TLatex* text = description->DrawLatex(i, TMath::Power(10, TMath::Log10(minval) + 0.05*(TMath::Log10(maxval)-TMath::Log10(minval))), trig->name.c_str());
        text->SetTextAngle(90);
        text->SetTextAlign(12);
    }
    c1->Draw();
    
    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.48, 0.85, "#bf{#it{ATLAS}} #it{p-Pb}");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.78, 0.84, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    TVectorD* lum_vec = (TVectorD*)thisfile->Get("lum_vec");
    description->DrawLatexNDC(0.78, 0.75, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", 1000.*(*lum_vec)[0])); 

    c1->SaveAs(Form("./Plots/triggers/run_trig_%i.pdf", thisRunNumber));
    cout << Form("Triggers for run number %i finished", thisRunNumber) << endl;

    TCanvas* c2 = new TCanvas("c2", "", 1000, 800);

    maxbincontent = 0;
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        thisbincontent = harr[i]->GetBinContent(i+1);
        if (maxbincontent < thisbincontent) maxbincontent = thisbincontent;
    }
    maxval = 6 * maxbincontent;

    minval = 1e0;

    gPad->SetLogz();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    description->SetTextAlign(32);
    description->SetTextFont(42);
    description->SetTextSize(0.012);

    TH2D* h2d = (TH2D*)thisfile->Get("pt_trig");
    h2d->GetXaxis()->SetLabelOffset(999);
    h2d->GetXaxis()->SetLabelSize(0);
    h2d->GetXaxis()->SetTickLength(0);
    h2d->GetXaxis()->SetTitleOffset(1.3);

    h2d->GetYaxis()->SetLabelOffset(999);
    h2d->GetYaxis()->SetLabelSize(0);    
    h2d->GetYaxis()->SetTickLength(0);
    h2d->GetYaxis()->SetTitleOffset(1.2);

    
    h2d->GetZaxis()->SetLabelSize(0.02);
    h2d->GetZaxis()->SetTickLength(0.01);
    h2d->GetZaxis()->SetTitleOffset(0.6);
    h2d->Draw("COLZ");
    
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        TLatex* text = description->DrawLatex(i+0.25, -0.65, trig->name.c_str());
        text->SetTextAngle(21);
    }
    description->SetTextSize(0.022);
    for (int pbin = 0; pbin < numpbins; pbin++) {
        description->DrawLatex(-0.55, ybins[pbin], Form("%i", (int)pbins[pbin]));
    }
    description->DrawLatex(-0.55, ybins[numpbins], Form("%i", (int)pbins[numpbins]));

    c2->Draw();
    
    description->SetTextSize(0.036);
    description->DrawLatexNDC(0.6, 0.85, "#bf{#it{ATLAS}} #it{p-Pb}");
    description->SetTextSize(0.032);
    description->DrawLatexNDC(0.84, 0.84, "#sqrt{s_{NN}^{avg}} = 8.16 TeV");
    description->DrawLatexNDC(0.84, 0.75, Form("#int#it{L}d#it{t} = %.3f nb^{-1}", 1000.*(*lum_vec)[0])); 

    c2->SaveAs(Form("./Plots/triggers/run_trigpt_%i.pdf", thisRunNumber));
    cout << Form("Triggers-pt joint plot for run %i finished", thisRunNumber) << endl;
}
