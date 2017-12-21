#include "electronOfflineAnalysis.C"
#include "/Users/jeffouellette/RootUtils/AtlasLabels.C"

const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
const Color_t mkcolors[10] = {kRed, kOrange+8, kBlue, kGreen+2, kTeal+4, kAzure+5, kMagenta, kSpring+5, kCyan+2, kPink+4};

bool display_counts = true;

/**
 * Plotting routine for the electron ptspectrum.
 */
void plot_electron_ptspectrum () {

    canvasName = "electron_ptspectrum_" + triggers[useTrigger];
    initialize_new_canvas(true);

    thishist = electron_ptspectrum;
    float new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
//    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);
    cout << "Electron pt spectrum integral = " << thishist->Integral() << endl;

    thishist->Draw("e1");

//    myText (0.64, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    ATLASLabel(0.6, 0.87, "Preliminary", kBlack);
    myText (0.6, 0.81, kBlack, Form("2017 #it{pp}, %.2f fb^{-1}", total_lumi/1000));
    myText (0.6, 0.75, kBlack, Form("#sqrt{#it{s}} = 5.02 TeV"));
//    myText (0.64, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    myText (0.6, 0.68, kBlack, Form("#||{#it{#eta}^{e}} < 2.47"));
    
    if(printStatementChecks) cout << "\nPlotting electron Pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the dielectron invariant mass.
 */
void plot_invariantMass () {

    canvasName = "electron_invariantMass_" + triggers[useTrigger];
    initialize_new_canvas(false);

    thishist = invariantMass;
    int counts = (int)thishist->Integral();
    cout << "Invariant mass integral = " << thishist->Integral() << endl;
//    thishist->Scale(1., "width");

    thishist = invariantMass_samesign;
    int counts_ss = (int)thishist->Integral();
//    thishist->Scale(1., "width");

    thishist = invariantMass_allsigns;
    int counts_as = (int)thishist->Integral();
//    thishist->Scale(1., "width");

    thishist = invariantMass;
    float new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    TH1F* background_hist = invariantMass_samesign;
    if (new_ymax < background_hist->GetBinContent(background_hist->GetMaximumBin())) {
        new_ymax = background_hist->GetBinContent(background_hist->GetMaximumBin());
    }
    background_hist = invariantMass_allsigns;
    if (new_ymax < background_hist->GetBinContent(background_hist->GetMaximumBin())) {
        new_ymax = background_hist->GetBinContent(background_hist->GetMaximumBin());
    }

    thishist = invariantMass;
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    thishist->SetMaximum(1.08 * new_ymax);    
    thishist->SetMinimum(2e-1);
    thishist->Draw("e1");
    if (display_counts) myMarkerText(0.22, 0.87, mkcolors[0], mkstyles[0], Form("OS (%i counts)", counts));

    thishist = invariantMass_samesign;
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[5]);
    thishist->SetLineColor(mkcolors[5]);
    thishist->SetMaximum(1.08 * new_ymax);    
    thishist->SetMinimum(2e-1);
    thishist->Draw("same e1");
    if (display_counts) myMarkerText(0.22, 0.78, mkcolors[5], mkstyles[0], Form("SS (%i counts)", counts_ss)); 
  
 
    myText (0.6, 0.89, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.595, 0.87, "Internal", kBlack);
    myText (0.6, 0.82, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.6, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    myText (0.6, 0.68, kBlack, Form("#||{#it{#eta}} < 2.47, #it{p}_{T} > %i GeV", ptcut));
    myText (0.6, 0.61, kBlack, Form("HLT_%s", triggers[useTrigger].substr(8).c_str()));
    myText (0.6, 0.54, kBlack, "Likelihood tight electrons");
   
    if(printStatementChecks) cout << "\nPlotting tight dielectrons invariant mass on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());

    thishist = invariantMass_allsigns;
    canvasName = "electron_invariantMass_allsigns_" + triggers[useTrigger];
    initialize_new_canvas(false);
     
    thishist = invariantMass_allsigns;
    thishist->SetMarkerStyle(kDot);
//    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(kBlack);
    thishist->SetLineColor(kBlack);
    thishist->SetMaximum(1.08 * new_ymax);    
    thishist->SetMinimum(2e-1);
    thishist->Draw("hist");
//    myMarkerText(0.22, 0.87, kBlack, mkstyles[0], Form("AS (%i counts)", counts));

//    myText (0.2, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    ATLASLabel(0.2, 0.87, "Preliminary", kBlack);
    myText (0.2, 0.81, kBlack, Form("2017 #it{pp}, %.2f fb^{-1}", total_lumi/1000));
    myText (0.2, 0.75, kBlack, "#sqrt{#it{s}} = 5.02 TeV");
//    myText (0.2, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    myText (0.2, 0.69, kBlack, Form("HLT_%s", triggers[useTrigger].substr(8).c_str()));
    myText (0.2, 0.63, kBlack, "Likelihood tight electrons");
    myText (0.2, 0.57, kBlack, Form("#||{#it{#eta}^{e}} < 2.47, #it{p}_{T}^{e} > %i GeV", ptcut));

    if(printStatementChecks) cout << "\nPlotting tight dielectrons (all signs) invariant mass on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the dielectron invariant (Z) pt spectrum.
 */
void plot_Z_ptspectrum () {

    canvasName = "Z_ptspectrum_" + triggers[useTrigger];
    initialize_new_canvas(true);

    float new_ymax = Z_ptspectrum->GetBinContent(Z_ptspectrum->GetMaximumBin());
    if (new_ymax < Z_ptspectrum_samesign->GetBinContent(Z_ptspectrum_samesign->GetMaximumBin())) {
        new_ymax = Z_ptspectrum_samesign->GetBinContent(Z_ptspectrum_samesign->GetMaximumBin());
    }

    
    thishist = Z_ptspectrum;
    int counts = (int)thishist->Integral();
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    cout << "Z pt spectrum integral = " << thishist->Integral() << endl;
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Draw("e1");
    if (display_counts) myMarkerText(0.25, 0.87, mkcolors[0], mkstyles[0], Form("OS (%i counts)", counts));

    thishist = Z_ptspectrum_samesign;
    counts = (int)thishist->Integral();
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[5]);
    thishist->SetLineColor(mkcolors[5]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Draw("same e1");
    if (display_counts) myMarkerText(0.25, 0.78, mkcolors[5], mkstyles[0], Form("SS (%i counts)", counts));
    
    myText (0.6, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.6, 0.87, "Internal", kBlack);
    myText (0.6, 0.81, kBlack, Form("2017 #it{pp}, %.2f fb^{-1}", total_lumi/1000));
    myText (0.6, 0.75, kBlack, "#sqrt{#it{s}} = 5.02 TeV");
//    myText (0.6, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    myText (0.6, 0.69, kBlack, Form("HLT_%s", triggers[useTrigger].substr(8).c_str()));
    myText (0.6, 0.63, kBlack, "Likelihood tight electrons");
    myText (0.6, 0.57, kBlack, Form("#||{#it{#eta}^{e}} < 2.47, #it{p}_{T}^{e} > %i GeV", ptcut));
    
    if(printStatementChecks) cout << "\nPlotting Z pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the ptspectrum, eta binned.
 */
void plot_electron_ptspectrum_etabinned () { 
 
    canvasName = "electron_ptspectrum_etabinned_" + triggers[useTrigger];
    initialize_new_canvas(true);

    float new_ymax = 0;
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = electron_ptspectrum_etabinned[etabin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    float ylabel = 0.48;
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = electron_ptspectrum_etabinned[etabin];
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*etabin)%10]);
        thishist->SetLineColor(mkcolors[(147*etabin)%10]);
        thishist->Scale(1./(etabins[etabin+1]-etabins[etabin]), "width");
        thishist->SetMaximum(2 * new_ymax);
        if (etabin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        myMarkerText(0.24, ylabel, mkcolors[(147*etabin)%10], mkstyles[0], Form("%g < #eta < %g", etabins[etabin], etabins[etabin+1]));
        ylabel -= 0.08;
    }

    myText (0.64, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.64, 0.87, "Internal", kBlack);
    myText (0.64, 0.81, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.64, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));

    if(printStatementChecks) cout << "\nPlotting etabinned electrons inclusive Pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the ptspectrum, phi binned.
 */
void plot_electron_ptspectrum_phibinned () {

    canvasName = "electron_ptspectrum_phibinned_" + triggers[useTrigger];
    initialize_new_canvas(true);

    float new_ymax = 0;
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    float ylabel = 0.48;
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*phibin)%10]);
        thishist->SetLineColor(mkcolors[(147*phibin)%10]);
        thishist->Scale(1./(phibins[phibin+1]-phibins[phibin]), "width");
        thishist->SetMaximum(2 * new_ymax);    
        if (phibin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        myMarkerText(0.24, ylabel, mkcolors[(147*phibin)%10], mkstyles[0], Form("%s < #phi < %s", DoubleRadiansToRationalRadians(phibins[phibin]).c_str(), DoubleRadiansToRationalRadians(phibins[phibin+1]).c_str()));
        ylabel -= 0.08;
    }

    myText (0.64, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.64, 0.87, "Internal", kBlack);
    myText (0.64, 0.81, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.64, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));

    if(printStatementChecks) cout << "\nPlotting phi binned electrons inclusive Pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the dielectron invariant mass, eta binned.
 */
void plot_invariantMass_etabinned () {

    canvasName = "electron_invariantMass_etabinned_" + triggers[useTrigger];
    initialize_new_canvas (true);

    float new_ymax = 0;
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = invariantMass_etabinned[etabin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    float ylabel = 0.88;
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = invariantMass_etabinned[etabin];
        thishist->Scale(1./(etabins[etabin+1]-etabins[etabin]), "width");
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*etabin)%10]);
        thishist->SetLineColor(mkcolors[(147*etabin)%10]);
        thishist->SetMaximum(2 * new_ymax);
        if (etabin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        myMarkerText(0.24, ylabel, mkcolors[(147*etabin)%10], mkstyles[0], Form("%g < #eta < %g", etabins[etabin], etabins[etabin+1]));
        ylabel -= 0.08;
    }

    myText (0.64, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.64, 0.87, "Internal", kBlack);
    myText (0.64, 0.81, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.64, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));

    if(printStatementChecks) cout << "\nPlotting tight dielectrons invariant mass on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for inclusive electron eta-phi distributions.
 */
void plot_eta_phi (bool cutonpt) {
    TH2F* this2hist;

    canvasName = "electron_eta_phi_" + triggers[useTrigger];
    if (cutonpt) {
        this2hist = eta_phi_hist;
        thishist = eta_phi_int_hist;
    } else {
        this2hist = eta_phi_hist_no_pt_cut;
        thishist = eta_phi_int_hist_no_pt_cut;
        canvasName += "_no_ptcut";
    }
    initialize_new_canvas (false);

    this2hist->Scale(1., "width"); 
    this2hist->Draw("colz");
    myText (0.64, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.64, 0.87, "Internal", kBlack);
    myText (0.64, 0.81, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.64, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    if (cutonpt) myText (0.64, 0.67, kBlack, Form("#it{p}_{T} > %i GeV", ptcut));
    else myText (0.64, 0.67, kBlack, Form("No #it{p}_{T} cut"));
    
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());

    canvasName = "electron_eta_phi_int_" + triggers[useTrigger];
    if (!cutonpt) canvasName += "_no_ptcut";
    initialize_new_canvas (false);
    thishist->Scale(1., "width");
    thishist->SetMinimum(0);
    cout << "Eta integral = " << thishist->Integral() << endl;

    thishist->Draw("hist e1");
    myText (0.62, 0.42, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.62, 0.42, "Internal", kBlack);
    myText (0.62, 0.36, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.62, 0.3, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    if (cutonpt) myText (0.62, 0.22, kBlack, Form("#it{p}_{T} > %i GeV", ptcut));
    else myText (0.62, 0.22, kBlack, Form("No #it{p}_{T} cut"));
    
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());

    return;
}


void plot_j_over_Z () {
    canvasName = "j_over_Z_" + triggers[useTrigger];
    initialize_new_canvas (false);

    thishist = j_over_Z_hist;
    int counts = (int)thishist->Integral();
    cout << "j/Z integral = " << counts << endl;
    thishist->Scale(1., "width");
    thishist->SetMinimum(0);
    thishist->Draw("hist e1");
    if (display_counts) myText (0.64, 0.49, kBlack, Form("(%i counts)", counts));
 
    myText (0.64, 0.87, kBlack, Form("#it{#bf{ATLAS}} Preliminary"));
    //ATLASLabel(0.64, 0.87, "Internal", kBlack);
    myText (0.64, 0.81, kBlack, "2017 #it{pp}, #sqrt{#it{s}} = 5.02 TeV");
    myText (0.64, 0.75, kBlack, Form("#it{L}_{int} = %.1f pb^{-1}", total_lumi));
    myText (0.64, 0.67, kBlack, Form("#it{p}_{T}^{Z} > %i GeV/#it{c}", Z_ptcut));
    myText (0.64, 0.61, kBlack, Form("#it{p}_{T}^{jet} > %i GeV/#it{c}", j_ptcut));
    myText (0.64, 0.55, kBlack, Form("#Delta#phi > %s", DoubleRadiansToRationalRadians(delta_phi_cut, false).c_str()));

    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
}

/**
 * Main electron offline analysis histogram macro.
 */
void electronOfflineAnalysisHist (int trig, char dataStream, bool dispcnts) {

    useTrigger = trig;
    display_counts = dispcnts;

    isExpress = false;
    isMain = false;
    isMinBias = false;
    switch (dataStream) {
        case 'e': {
            isExpress = true;
            dataPath = "./Data/express/";
            break;
        }
        case 'm': {
            isMain = true;
            dataPath = "./Data/main/";
            break;
        }
        case 'b': {
            isMinBias = true;
            dataPath = "./Data/minbias/";
            break;
        }
    }

    runNumber = 0;
    plotPath = "./Plots/";

    initialize_text();
    initialize_histograms ();

    for (int i = 0; i < fileRunNumbers.size(); i++) {
        int rn = fileRunNumbers[i];
        TFile* thisfile = new TFile(Form("./Data/RootOutput/%i_%s_out.root", rn, triggers[useTrigger].c_str()), "READ");
         
        electron_ptspectrum->Add((TH1F*)thisfile->Get(Form("run_%i_electron_ptspectrum_hist", rn)));
        for (int etabin = 0; etabin < numetabins; etabin++) {
            electron_ptspectrum_etabinned[etabin]->Add((TH1F*)thisfile->Get(Form("run_%i_electron_ptspectrum_hist_etabin%i", rn, etabin)));
            invariantMass_etabinned[etabin]->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_hist_etabin%i", rn, etabin)));
        }
        for (int phibin = 0; phibin < numphibins; phibin++) {
            electron_ptspectrum_phibinned[phibin]->Add((TH1F*)thisfile->Get(Form("run_%i_electron_ptspectrum_hist_phibin%i", rn, phibin)));
        }
        invariantMass->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_hist", rn)));
        invariantMass_samesign->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_samesign_hist", rn)));
        invariantMass_allsigns->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_allsigns_hist", rn)));
        Z_ptspectrum->Add((TH1F*)thisfile->Get(Form("run_%i_Z_ptspectrum_hist", rn)));
        Z_ptspectrum_samesign->Add((TH1F*)thisfile->Get(Form("run_%i_Z_ptspectrum_samesign_hist", rn)));
        j_over_Z_hist->Add((TH1F*)thisfile->Get(Form("run_%i_j_over_Z_hist", rn)));
        eta_phi_hist->Add((TH2F*)thisfile->Get(Form("run_%i_eta_phi_hist", rn)));
        eta_phi_hist_no_pt_cut->Add((TH2F*)thisfile->Get(Form("run_%i_eta_phi_hist_no_pt_cut", rn)));
        eta_phi_int_hist->Add((TH1F*)thisfile->Get(Form("run_%i_eta_phi_int_hist", rn)));
        eta_phi_int_hist_no_pt_cut->Add((TH1F*)thisfile->Get(Form("run_%i_eta_phi_int_hist_no_pt_cut", rn)));
    }

    plot_electron_ptspectrum ();
    plot_invariantMass ();
    plot_Z_ptspectrum (); 
    plot_j_over_Z ();
    plot_electron_ptspectrum_etabinned ();
    plot_electron_ptspectrum_phibinned ();
    plot_invariantMass_etabinned ();
    plot_eta_phi (true);
    plot_eta_phi (false);

    return;
}
