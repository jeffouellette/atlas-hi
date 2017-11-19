#include "electronOfflineAnalysis.C"
#include "/Users/jeffouellette/RootUtils/AtlasLabels.C"
#include "/Users/jeffouellette/RootUtils/AtlasUtils.C"

/**
 * Plotting routine for the electron ptspectrum.
 */
void plot_electron_ptspectrum () {

    canvasName = "electron_ptspectrum_" + triggers[useTrigger];
    initialize_new_canvas (true, 0.1, 0.03, 0.1, 0.05);

    thishist = electron_ptspectrum;
    float new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);
    thishist->Draw("e1");

    ATLASLabel(0.65, 0.87, "Internal", kBlack);
    myText (0.65, 0.78, kBlack, "pp #sqrt{s} = 5.02 TeV");
    myText (0.65, 0.69, kBlack, "#int#it{L}dt = 24 pb^{-1}");
    myText (0.65, 0.58, kBlack, "#left|#eta#right| < 1.37 or");
    myText (0.65, 0.49, kBlack, "1.56 < #left|#eta#right| < 2.37");
   // myMarkerText(0.65, 0.71, mkcolors[0], mkstyles[0], "e15_lhloose");
   // myMarkerText(0.65, 0.71, mkcolors[0], mkstyles[0], "e15_lhloose");
    
//    draw_title (0.5, 0.973, "Tight electron inclusive #it{p}_{T} spectrum");    
//    draw_information (0.38, 0.34);
//    draw_text (0.38, 0.34, "0 < |#eta| < 1.37 or 1.57 < #eta < 2.37");

    if(printStatementChecks) cout << "\nPlotting electron Pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the dielectron invariant mass.
 */
void plot_invariantMass () {

    canvasName = "electron_invariantMass_" + triggers[useTrigger];
    initialize_new_canvas (true, 0.1, 0.03, 0.1, 0.05);

//    TLegend* invariantMass_legend = new TLegend (0.75, 0.7, 0.97, 0.95);
//    invariantMass_legend->SetTextSize(0.028);
    // Add each histogram to the legend, then scale each histogram by bin width and draw
    thishist = invariantMass;
    float new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    if (new_ymax < invariantMass_samesign->GetBinContent(invariantMass_samesign->GetMaximumBin())) {
        new_ymax = invariantMass_samesign->GetBinContent(invariantMass_samesign->GetMaximumBin());
    }
    thishist = invariantMass;
//    invariantMass_legend->AddEntry (thishist, "#sum #it{q}_{i} = 0");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Draw("e1");

    thishist = invariantMass_samesign;
//    invariantMass_legend->AddEntry (thishist, "No charge cond.");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[1]);
    thishist->SetLineColor(mkcolors[1]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Draw("same e1");
   
    ATLASLabel(0.65, 0.87, "Internal", kBlack);
    myText (0.65, 0.78, kBlack, "pp #sqrt{s} = 5.02 TeV");
    myText (0.65, 0.69, kBlack, "#int#it{L}dt = 24 pb^{-1}");
    myText (0.65, 0.58, kBlack, "#left|#eta#right| < 1.37 or");
    myText (0.65, 0.49, kBlack, "1.56 < #left|#eta#right| < 2.37");
    
    myMarkerText(0.24, 0.86, mkcolors[0], mkstyles[0], "Opposite charges");
    myMarkerText(0.24, 0.77, mkcolors[1], mkstyles[0], "Same charges");

//    invariantMass_legend->Draw();

//    draw_title (0.5, 0.975, "Tight dielectrons invariant mass");
//    draw_information (0.43, 0.87);
//    draw_text (0.43, 0.795, Form("#it{p}_{T} #geq %i GeV/#it{c} and |#eta| < 1.37 or 1.57 < #eta < 2.37", ptcut));

    if(printStatementChecks) cout << "\nPlotting tight dielectrons invariant mass on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the dielectron invariant (Z) pt spectrum.
 */
void plot_Z_ptspectrum () {

    canvasName = "Z_ptspectrum_" + triggers[useTrigger];
    initialize_new_canvas (true, 0.1, 0.03, 0.1, 0.05); 

    // Add each histogram to the legend, then scale each histogram by bin width and draw    
    float new_ymax = Z_ptspectrum->GetBinContent(Z_ptspectrum->GetMaximumBin());
    if (new_ymax < Z_ptspectrum_samesign->GetBinContent(Z_ptspectrum_samesign->GetMaximumBin())) {
        new_ymax = Z_ptspectrum_samesign->GetBinContent(Z_ptspectrum_samesign->GetMaximumBin());
    }

//    TLegend* Z_ptspectrum_legend = new TLegend (0.75, 0.7, 0.97, 0.95);
//    Z_ptspectrum_legend->SetTextSize(0.028);
    
    thishist = Z_ptspectrum;
//    Z_ptspectrum_legend->AddEntry (thishist, "#sum q_{i} = 0");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Draw("e1");

    thishist = Z_ptspectrum_samesign;
//    Z_ptspectrum_legend->AddEntry (thishist, "No charge cond.");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[1]);
    thishist->SetLineColor(mkcolors[1]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Draw("same e1");

//    Z_ptspectrum_legend->Draw();

//    draw_title (0.5, 0.973, "#it{Z} #rightarrow #it{e}^{+}#it{e}^{-} inclusive #it{p}_{T} spectrum");    
//    draw_information(0.5, 0.855);
//    draw_text (0.5, 0.78, Form("#it{p}_{T} #geq %i GeV/#it{c} and |#eta| < 1.37 or 1.57 < #eta < 2.37", ptcut));

    if(printStatementChecks) cout << "\nPlotting Z pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the ptspectrum, eta binned.
 */
void plot_electron_ptspectrum_etabinned () { 
 
    canvasName = "electron_ptspectrum_etabinned_" + triggers[useTrigger];
    initialize_new_canvas (true, 0.1, 0.03, 0.1, 0.05);

    // Create the legend
//    TLegend* electron_ptspectrum_etabinned_legend = new TLegend(0.7, 0.7, 0.97, 0.95);
//    electron_ptspectrum_etabinned_legend->SetTextSize(0.028);
    // Add each histogram to the legend, then scale each histogram by bin width and draw    
    float new_ymax = 0;
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = electron_ptspectrum_etabinned[etabin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = electron_ptspectrum_etabinned[etabin];
//        electron_ptspectrum_etabinned_legend->AddEntry(thishist, Form("%g < #eta_{lab} < %g", etabins[etabin], etabins[etabin+1]));
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*etabin)%10]);
        thishist->SetLineColor(mkcolors[(147*etabin)%10]);
        thishist->Scale(1./(etabins[etabin+1]-etabins[etabin]), "width");
        thishist->SetMaximum(2 * new_ymax);
        if (etabin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
    }
//    draw_title (0.5, 0.973, "Tight electrons inclusive #it{p}_{T} spectrum");    
//    draw_information (0.38, 0.34);
//    electron_ptspectrum_etabinned_legend->Draw();
    if(printStatementChecks) cout << "\nPlotting etabinned electrons inclusive Pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the ptspectrum, phi binned.
 */
void plot_electron_ptspectrum_phibinned () {

    canvasName = "electron_ptspectrum_phibinned_" + triggers[useTrigger];
    initialize_new_canvas(true, 0.1, 0.03, 0.1, 0.05);

    // Create the legend
//    TLegend* electron_ptspectrum_phibinned_legend = new TLegend(0.7, 0.55, 0.97, 0.95);
//    electron_ptspectrum_phibinned_legend->SetTextSize(0.028);
    // Add each histogram to the legend, then scale each histogram by bin width and draw    
    float new_ymax = 0;
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
//        electron_ptspectrum_phibinned_legend->AddEntry(thishist, Form("%g#pi < #phi < %g#pi", phibins[phibin]/pi, phibins[phibin+1]/pi));
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*phibin)%10]);
        thishist->SetLineColor(mkcolors[(147*phibin)%10]);
        thishist->Scale(1./(phibins[phibin+1]-phibins[phibin]), "width");
        thishist->SetMaximum(2 * new_ymax);    
        if (phibin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
    }
//    electron_ptspectrum_phibinned_legend->Draw();
//    draw_title (0.5, 0.973, "Tight electrons inclusive #it{p}_{T} spectrum");
//    draw_information (0.38, 0.34);
    if(printStatementChecks) cout << "\nPlotting phi binned electrons inclusive Pt spectrum on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Plotting routine for the dielectron invariant mass, eta binned.
 */
void plot_invariantMass_etabinned () {

    canvasName = "electron_invariantMass_etabinned_" + triggers[useTrigger];
    initialize_new_canvas (true, 0.1, 0.03, 0.1, 0.05);

    // Create the legend
//    TLegend* invariantMass_legend = new TLegend(0.75, 0.7, 0.97, 0.95);
//    invariantMass_legend->SetTextSize(0.028);
    // Add each histogram to the legend, then scale each histogram by bin width and draw
    float new_ymax = 0;
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = invariantMass_etabinned[etabin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    for (int etabin = 0; etabin < numetabins; etabin++) {
        if (etabin == 1 || etabin == 4) continue;
        thishist = invariantMass_etabinned[etabin];
//        invariantMass_legend->AddEntry(thishist, Form("%g < #eta_{lab} < %g", etabins[etabin], etabins[etabin+1]));
        thishist->Scale(1./(etabins[etabin+1]-etabins[etabin]), "width");
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*etabin)%10]);
        thishist->SetLineColor(mkcolors[(147*etabin)%10]);
        thishist->SetMaximum(2 * new_ymax);
        if (etabin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
    }
//    invariantMass_legend->Draw();
//    draw_title (0.5, 0.975, "Tight dielectrons invariant mass");
//    draw_information (0.43, 0.87);
//    draw_text (0.43, 0.795, Form("#it{p}_{T} #geq %i GeV/#it{c} and |#eta| < 1.37 or 1.57 < #eta < 2.37", ptcut));
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

    initialize_new_canvas (false, 0.08, 0.12, 0.08, 0.05);
    this2hist->Scale(1., "width"); 
    this2hist->Draw("colz");
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());

    canvasName = "electron_eta_phi_int_" + triggers[useTrigger];
    if (!cutonpt) canvasName += "_no_ptcut";
    thishist->Scale(1., "width");
    thishist->SetMinimum(0);

    thishist->Draw("hist e1");
    
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());

    return;
}

/**
 * Main electron offline analysis histogram macro.
 */
void electronOfflineAnalysisHist (int trig, char dataStream) {

    useTrigger = trig;

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

    initialize_histograms ();

    const std::vector<char*> filesToProcess = {
        //Form("%suser.khill.pp5TeVAna_main.002.00340634.f895_m1902_myOutput_hadd.root", dataPath.c_str()),
        //Form("%suser.khill.pp5TeVAna_main.004.00340644.f895_m1902_myOutput_hadd.root", dataPath.c_str()),
        Form("%suser.khill.pp5TeVAna_main.004.00340683.f896_m1902_myOutput_hadd.root", dataPath.c_str()),
        Form("%suser.khill.pp5TeVAna_main.004.00340697.f896_m1902_myOutput_hadd.root", dataPath.c_str())
    };
    std::vector<int> fileRunNumbers = {
        //340634,
        //340644,
        340683,
        340697
    };
    
    for (int i = 0; i < filesToProcess.size(); i++) {
        int rn = fileRunNumbers[i];
        TFile* thisfile = new TFile(Form("./Data/%i_%s_out.root", rn, triggers[useTrigger].c_str()), "READ");
         
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
        Z_ptspectrum->Add((TH1F*)thisfile->Get(Form("run_%i_Z_ptspectrum_hist", rn)));
        Z_ptspectrum_samesign->Add((TH1F*)thisfile->Get(Form("run_%i_Z_ptspectrum_samesign_hist", rn)));
        eta_phi_hist->Add((TH2F*)thisfile->Get(Form("run_%i_eta_phi_hist", rn)));
        eta_phi_hist_no_pt_cut->Add((TH2F*)thisfile->Get(Form("run_%i_eta_phi_hist_no_pt_cut", rn)));
        eta_phi_int_hist->Add((TH1F*)thisfile->Get(Form("run_%i_eta_phi_int_hist", rn)));
        eta_phi_int_hist_no_pt_cut->Add((TH1F*)thisfile->Get(Form("run_%i_eta_phi_int_hist_no_pt_cut", rn)));
    }

    plot_electron_ptspectrum ();
    plot_invariantMass ();
    plot_Z_ptspectrum (); 
    plot_electron_ptspectrum_etabinned ();
    plot_electron_ptspectrum_phibinned ();
    plot_invariantMass_etabinned ();
    plot_eta_phi (true);
    plot_eta_phi (false);
    

    return;
}
