#include "electronOfflineAnalysis.C"

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
    thishist->Write();
    thishist->Draw("e1");

//    draw_title (0.5, 0.973, "Tight electron inclusive #it{p}_{T} spectrum");    
    draw_information (0.38, 0.34);
    draw_text (0.38, 0.34, "0 < |#eta| < 1.37 or 1.57 < #eta < 2.37");

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

    TLegend* invariantMass_legend = new TLegend (0.75, 0.7, 0.97, 0.95);
    invariantMass_legend->SetTextSize(0.028);
    // Add each histogram to the legend, then scale each histogram by bin width and draw
    thishist = invariantMass;
    float new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    if (new_ymax < invariantMass_nocharge->GetBinContent(invariantMass_nocharge->GetMaximumBin())) {
        new_ymax = invariantMass_nocharge->GetBinContent(invariantMass_nocharge->GetMaximumBin());
    }
    thishist = invariantMass;
    invariantMass_legend->AddEntry (thishist, "#sum #it{q}_{i} = 0");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Write();
    thishist->Draw("e1");

    thishist = invariantMass_nocharge;
    invariantMass_legend->AddEntry (thishist, "No charge cond.");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[1]);
    thishist->SetLineColor(mkcolors[1]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Write();
    thishist->Draw("same e1");

    invariantMass_legend->Draw();

//    draw_title (0.5, 0.975, "Tight dielectrons invariant mass");
    draw_information (0.43, 0.87);
    draw_text (0.43, 0.795, Form("#it{p}_{T} #geq %i GeV/#it{c} and |#eta| < 1.37 or 1.57 < #eta < 2.37", ptcut));

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
    if (new_ymax < Z_ptspectrum_nocharge->GetBinContent(Z_ptspectrum_nocharge->GetMaximumBin())) {
        new_ymax = Z_ptspectrum_nocharge->GetBinContent(Z_ptspectrum_nocharge->GetMaximumBin());
    }

    TLegend* Z_ptspectrum_legend = new TLegend (0.75, 0.7, 0.97, 0.95);
    Z_ptspectrum_legend->SetTextSize(0.028);
    
    thishist = Z_ptspectrum;
    Z_ptspectrum_legend->AddEntry (thishist, "#sum q_{i} = 0");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[0]);
    thishist->SetLineColor(mkcolors[0]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Write();
    thishist->Draw("e1");

    thishist = Z_ptspectrum_nocharge;
    Z_ptspectrum_legend->AddEntry (thishist, "No charge cond.");
    thishist->SetMarkerStyle(mkstyles[0]);
    thishist->SetMarkerColor(mkcolors[1]);
    thishist->SetLineColor(mkcolors[1]);
    thishist->Scale(1., "width");
    thishist->SetMaximum(2 * new_ymax);    
    thishist->Write();
    thishist->Draw("same e1");

    Z_ptspectrum_legend->Draw();

//    draw_title (0.5, 0.973, "#it{Z} #rightarrow #it{e}^{+}#it{e}^{-} inclusive #it{p}_{T} spectrum");    
    draw_information(0.5, 0.855);
    draw_text (0.5, 0.78, Form("#it{p}_{T} #geq %i GeV/#it{c} and |#eta| < 1.37 or 1.57 < #eta < 2.37", ptcut));

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
    TLegend* electron_ptspectrum_etabinned_legend = new TLegend(0.7, 0.7, 0.97, 0.95);
    electron_ptspectrum_etabinned_legend->SetTextSize(0.028);
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
        electron_ptspectrum_etabinned_legend->AddEntry(thishist, Form("%g < #eta_{lab} < %g", etabins[etabin], etabins[etabin+1]));
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*etabin)%10]);
        thishist->SetLineColor(mkcolors[(147*etabin)%10]);
        thishist->Scale(1./(etabins[etabin+1]-etabins[etabin]), "width");
        thishist->SetMaximum(2 * new_ymax);
        if (etabin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        thishist->Write();
    }
//    draw_title (0.5, 0.973, "Tight electrons inclusive #it{p}_{T} spectrum");    
    draw_information (0.38, 0.34);
    electron_ptspectrum_etabinned_legend->Draw();
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
    TLegend* electron_ptspectrum_phibinned_legend = new TLegend(0.7, 0.55, 0.97, 0.95);
    electron_ptspectrum_phibinned_legend->SetTextSize(0.028);
    // Add each histogram to the legend, then scale each histogram by bin width and draw    
    float new_ymax = 0;
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
        if (thishist->GetBinContent(thishist->GetMaximumBin()) > new_ymax) new_ymax = thishist->GetBinContent(thishist->GetMaximumBin());
    }
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
        electron_ptspectrum_phibinned_legend->AddEntry(thishist, Form("%g#pi < #phi < %g#pi", phibins[phibin]/pi, phibins[phibin+1]/pi));
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*phibin)%10]);
        thishist->SetLineColor(mkcolors[(147*phibin)%10]);
        thishist->Scale(1./(phibins[phibin+1]-phibins[phibin]), "width");
        thishist->SetMaximum(2 * new_ymax);    
        if (phibin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        thishist->Write();
    }
    electron_ptspectrum_phibinned_legend->Draw();
//    draw_title (0.5, 0.973, "Tight electrons inclusive #it{p}_{T} spectrum");
    draw_information (0.38, 0.34);
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
    TLegend* invariantMass_legend = new TLegend(0.75, 0.7, 0.97, 0.95);
    invariantMass_legend->SetTextSize(0.028);
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
        invariantMass_legend->AddEntry(thishist, Form("%g < #eta_{lab} < %g", etabins[etabin], etabins[etabin+1]));
        thishist->Scale(1./(etabins[etabin+1]-etabins[etabin]), "width");
        thishist->SetMarkerStyle(mkstyles[0]);
        thishist->SetMarkerColor(mkcolors[(147*etabin)%10]);
        thishist->SetLineColor(mkcolors[(147*etabin)%10]);
        thishist->SetMaximum(2 * new_ymax);
        if (etabin == 0) thishist->Draw("e1");
        else thishist->Draw("same e1");
        thishist->Write();
    }
    invariantMass_legend->Draw();
//    draw_title (0.5, 0.975, "Tight dielectrons invariant mass");
    draw_information (0.43, 0.87);
    draw_text (0.43, 0.795, Form("#it{p}_{T} #geq %i GeV/#it{c} and |#eta| < 1.37 or 1.57 < #eta < 2.37", ptcut));
    if(printStatementChecks) cout << "\nPlotting tight dielectrons invariant mass on canvas " << thiscanvas->GetName() << endl;
    thiscanvas->SaveAs((plotPath + canvasName + ".pdf").c_str());
    return;
}


/**
 * Main electron offline analysis histogram macro.
 */
void electronOfflineAnalysisHist.C (int trig, char dataStream) {

    runNumber = rn;

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

    plotPath = Form("./Plots/%i/", runNumber);

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
        TFile* thisfile = new TFile(Form("%s/%i_%s_out.root", dataPath.c_str(), fileRunNumbers[i], triggers[useTrigger].c_str()), "READ");
         
        electron_ptspectrum->Add((TH1F*)thisfile->Get(Form("run_%i_electron_ptspectrum_hist", runNumber)));
        for (int etabin = 0; etabin < numetabins; etabin++) {
            electron_ptspectrum_etabinned[numetabins]->Add((TH1F*)thisfile->Get(Form("run_%i_electron_ptspectrum_hist_etabin%i", etabin)));
            invariantMass_etabinned[numetabins]->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_hist_etabin%i", etabin)));
        }
        for (int phibin = 0; phibin < numphibins; phibin++) {
            electron_ptspectrum_phibinned[numphibins]->Add((TH1F*)thisfile->Get(Form("run_%i_electron_ptspectrum_hist_phibin%i", phibin)));
        }
        invariantMass->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_hist", runNumber)));
        invariantMass_nocharge->Add((TH1F*)thisfile->Get(Form("run_%i_invariantMass_nocharge_hist", runNumber)));
        Z_ptspectrum->Add((TH1F*)thisfile->Get(Form("run_%i_Z_ptspectrum_hist", runNumber)));
        Z_ptspectrum_nocharge->Add((TH1F*)thisfile->Get(Form("run_%i_Z_ptspectrum_nocharge_hist", runNumber)));
        eta_phi_hist->Add((TH2F*)thisfile->Get(Form("run_%i_eta_phi_hist", runNumber)));
        eta_phi_hist_no_pt_cut->Add((TH2F*)thisfile->Get(Form("run_%i_eta_phi_hist_no_pt_cut", runNumber)));
        eta_phi_int_hist->Add((TH1F*)thisfile->Get(Form("run_%i_eta_phi_int_hist", runNumber)));
        eta_phi_int_hist_no_pt_cut->Add((TH1F*)thisfile->Get(Form("run_%i_eta_phi_int_hist_no_pt_cut", runNumber)));
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
