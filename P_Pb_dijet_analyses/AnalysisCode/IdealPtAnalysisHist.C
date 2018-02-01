#include "../triggerUtil.C"

void IdealPtAnalysisHist() {

    initialize(0, true, false);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numruns = (*thisRunNumbers).size();
    const int numhists = numtrigs * numruns * numetabins;
    if (debugStatements) {
        cout << "Status: In triggers_pt_counts_hist.C (11): Building trigger pt histograms with " << numruns << " runs being used" << endl;
        cout << "Status: In triggers_pt_counts_hist.C (12): Numtrigs = " << numtrigs << endl;
        cout << "Status: In triggers_pt_counts_hist.C (13): Numetabins = " << numetabins << endl;
        cout << "Status: In triggers_pt_counts_hist.C (14): Numpbins = " << numpbins << endl;
        cout << "Status: In triggers_pt_counts_hist.C (15): ptPath = " << ptPath << endl;
    }

    double ymin = 5e-8;
    double ymax = 5e6;
    const Style_t mkstyles[2] = {kFullTriangleUp, kFullTriangleDown};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    double hscale, deta;
    TH1F* thisHist;

    if (debugStatements) cout << "Status: In triggers_pt_counts_hist.C (26): Initialized histArr histograms..." << endl;
    TH1F* histArr[numetabins];
    for (int etabin = 0; etabin < numetabins; etabin++) {
        histArr[etabin] = new TH1F(Form("best_statistics_etabin%i", etabin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/d#it{p}_{T}dy #left[nb GeV^{-1}#right]", numpbins, pbins);
    }
    if (debugStatements) {
        cout << "Status: In triggers_pt_counts_hist.C (32): histArr histograms initialized." << endl;
        cout << "Status: In triggers_pt_counts_hist.C (33): Starting loop over triggers..." << endl;
    }

    vector<int>* runNumbers = getRunNumbers();
    totalLuminosity = 0;

    /**** Fill summed histograms with results from event loops ****/
    {
        TSystemDirectory dir(ptPath.c_str(), ptPath.c_str());
        TList* sysfiles = dir.GetListOfFiles();
        if (sysfiles) {
            TSystemFile *sysfile;
            TString fname;
            TString histName;
            TIter next(sysfiles);
            TVectorD* run_vec;

            while ((sysfile=(TSystemFile*)next())) {
                fname = sysfile->GetName();
                if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
                    if (debugStatements) cout << "Status: In triggers_pt_counts.C (51): Found " << fname.Data() << endl; 
                    for (int thisRunNumber : *runNumbers) {
                        if (skipRun(thisRunNumber)) continue;
                        if (fname.Contains(to_string(thisRunNumber))) {
                            TFile* thisFile = new TFile(ptPath + fname, "READ");
                        //    totalLuminosity += (TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber))[0];

                            // quickly check the parameters stored in this root file
                            run_vec = (TVectorD*)thisFile->Get(Form("run_vec_%i", thisRunNumber));
                            assert (run_vec[0] == thisRunNumber);
                            assert (run_vec[1] == numetabins);
                            assert (run_vec[2] == numtrigs);
                            assert (run_vec[3] == numpbins);

                            int actetabin; // used to flip period A pseudorapidities
                            for (int etabin = 0; etabin < numetabins; etabin++) {
                                if (thisRunNumber < 313500) actetabin = numetabins - etabin - 1;
                                else actetabin = etabin;
                                histName = Form("trig_pt_counts_run%i_etabin%i", thisRunNumber, actetabin);
                                histArr[etabin]->Add((TH1F*)thisFile->Get(histName));
                            }
                            thisFile->Close();
                            delete thisFile;
                            break;
                        }
                    }
                }
            }
        }
    }
    /**** End fill summed histograms ****/


    /** Plotting routines **/

    // Plot best-selected pt spectra
    double* histArrScales = linspace(-1.5, 1.5, numetabins/2 - 1);
//    double* histArrScales = linspace(0, 0, numetabins/2 - 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts
    TCanvas* trigCanvas = new TCanvas("trigCanvas", "", 800, 600);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicks();
    trigCanvas->Draw();
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thisHist = histArr[etabin];
        for (int pbin = 0; pbin < numpbins; pbin++) {
            if (kinematicLumiVec[pbin + etabin*numpbins] != 0) {
                thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1) / (kinematicLumiVec[pbin + etabin*numpbins]));
                thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1) / (kinematicLumiVec[pbin + etabin*numpbins]));
            }
            else if (debugStatements) cout << "Warning: In triggers_pt_counts_hist.C (103): No exposed luminosity between pt= " << pbins[pbin] << ", " << pbins[pbin+1] << " and eta= " << etabins[etabin] << ", " << etabins[etabin+1] << endl;
        }
        deta = etabins[etabin+1] - etabins[etabin];
        thisHist->Scale(1e3*TMath::Power(10, histArrScales[(int)(numetabins/2 - 0.5 -TMath::Abs(etabin - numetabins/2 + 0.5))])/deta, "width"); // separate different etabins
        thisHist->SetMarkerStyle(mkstyles[etabin < (numetabins/2)]);
//        thisHist->SetMarkerStyle(kDot);
        Color_t kColor = mkcolors[(147*etabin)%20];
        thisHist->SetMarkerColor(kColor);
        thisHist->SetLineColor(kColor);
        thisHist->SetMinimum(ymin);
        thisHist->SetMaximum(ymax);
        thisHist->GetYaxis()->SetTitleOffset(1.35);
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);
        if (etabin == 0) thisHist->Draw("e1");
        else thisHist->Draw("same e1");
        
        const float textx = 0.46 + (etabin>=(numetabins/2))*0.26;
        const float texty = 0.91 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
        const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[etabin], etabins[etabin+1], histArrScales[(int)((0.5*(numetabins-1))-TMath::Abs(etabin-(0.5*(numetabins-1))))]);
        myMarkerText (textx, texty, kColor, mkstyles[etabin < (numetabins/2)], text);
    }

    myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity));
    myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

    string histName;
    if (numetabins > 1) {
        histName = "ptSpectra_combinedTriggers_etabinned";
    }
    else histName = "ptSpectra_combinedTriggers_0eta490";

    if (runPeriodA && !runPeriodB) {
        myText (0.19, 0.33, kBlack, "Period A");
        histName = histName + "_periodA";
    }
    else if (!runPeriodA && runPeriodB) {
        myText (0.19, 0.33, kBlack, "Period B");
        histName = histName + "_periodB";
    }
    else {
        myText (0.19, 0.33, kBlack, "Period A & B");
    }

    if (runPeriodA || runPeriodB) trigCanvas->SaveAs((plotPath + "ptSpectra/" + histName + ".pdf").c_str());

    /**** Free memory and quit ****/
    for (int etabin = 0; etabin < numetabins; etabin++) delete histArr[etabin];
    delete trigCanvas;
    delete[] histArrScales;
    delete runNumbers;

    if (debugStatements) cout << "Status: In triggers_pt_counts_hist.C (154): Finished plotting pt spectrum" << endl;
    return;
}
