#include "IdealRpPbAnalysis.C"

void IdealRpPbAnalysisHist() {

    TH1D** ppHistArr = setupPPConfiguration(); // creates the pp histograms and setups the relevant environment variables

    initialize(0, true);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numruns = (*thisRunNumbers).size();
    const int numhists = numtrigs * numruns * numetabins;
    if (debugStatements) {
        cout << "Status: In IdealRpPbAnalysisHist.C (11): Building trigger pt histograms with " << numruns << " runs being used" << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (12): Numtrigs = " << numtrigs << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (13): Numetabins = " << numppEtabins << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (14): Numpbins = " << numpbins << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (15): ptPath = " << ptPath << endl;
    }

    double ymin = 0;
    double ymax = 5;
    const Style_t mkstyles[2] = {kFullTriangleUp, kFullTriangleDown};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    double hscale, deta;
    TH1D* thisHist;

    if (debugStatements) cout << "Status: In IdealRpPbAnalysisHist.C (26): Initialized pPbHistArr histograms..." << endl;
    TH1D* pPbHistArr[numpPbCoMEtabins];
    for (int pPbCoMEtabin = 0; pPbCoMEtabin < numpPbCoMEtabins; pPbCoMEtabin++) {
        int ppEtabin_equiv = getppEtabin(TMath::Abs(0.5*(pPbCoMEtabins[pPbCoMEtabin]+pPbCoMEtabins[pPbCoMEtabin+1]) - etaCoM));
        pPbHistArr[pPbCoMEtabin] = new TH1D(Form("pPb_spectrum_pPbCoMEtabin%i", pPbCoMEtabin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}dy #left[nb GeV^{-1}#right]", numppPtEtabins[ppEtabin_equiv], ppPtbins);
        pPbHistArr[pPbCoMEtabin]->Sumw2();
    }
    TH1D* RpPbHistArr[numppEtabins];
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        RpPbHistArr[ppEtabin] = new TH1D(Form("R_pPb_ppEtabin%i", ppEtabin), ";#it{p}_{T}^{jet} #left[GeV#right];R_{#it{pPb}}", numppPtEtabins[ppEtabin], ppPtbins);
        RpPbHistArr[ppEtabin]->Sumw2();
    }

    if (debugStatements) {
        cout << "Status: In IdealRpPbAnalysisHist.C (32): pPbHistArr histograms initialized." << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (33): Starting loop over triggers..." << endl;
    }

    vector<int>* runNumbers = getRunNumbers();
    double totalpPbLuminosity = 0; // nb^-1
    double totalppLuminosity = 20.2; // fb^-1

    /**** Fill summed histograms with results from event loops ****/
    {
        TSystemDirectory dir(ptPath.c_str(), RpPbPath.c_str());
        TList* sysfiles = dir.GetListOfFiles();
        if (sysfiles) {
            TSystemFile *sysfile;
            TString fname;
            TString histName;
            TIter next(sysfiles);
            TVectorD* run_vec;
            TVectorD* lum_vec;

            while ((sysfile=(TSystemFile*)next())) {
                fname = sysfile->GetName();
                if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
                    if (debugStatements) cout << "Status: In triggers_pt_counts.C (51): Found " << fname.Data() << endl; 
                    for (int thisRunNumber : *runNumbers) {
                        if (skipRun(thisRunNumber)) continue;
                        if (fname.Contains(to_string(thisRunNumber))) {
                            TFile* thisFile = new TFile(RpPbPath + fname, "READ");
//                            lum_vec = (TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber));
  //                          cout << lum_vec[0] << endl;
 //                           totalpPbLuminosity += lum_vec[0]; 
 //                           totalpPbLuminosity += ((TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber)))[0];

                            // quickly check the parameters stored in this root file
                            run_vec = (TVectorD*)thisFile->Get(Form("run_vec_%i", thisRunNumber));
                            assert (run_vec[0] == thisRunNumber);
                            assert (run_vec[1] == numppEtabins);
                            assert (run_vec[2] == numpPbCoMEtabins);
                            assert (run_vec[3] == numtrigs);
                            assert (run_vec[4] == numppPtbins);

                            for (int pPbCoMEtabin = 0; pPbCoMEtabin < numpPbCoMEtabins; pPbCoMEtabin++) {
                                histName = Form("pPb_spectrum_run%i_pPbCoMEtabin%i", thisRunNumber, pPbCoMEtabin);
                                pPbHistArr[pPbCoMEtabin]->Add((TH1D*)thisFile->Get(histName));
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

    /**** Add histograms to match pp binning in eta, then divide the pPb spectrum by the pp spectrum ****/
    for (int pPbCoMEtabin = 0; pPbCoMEtabin < numpPbCoMEtabins; pPbCoMEtabin++) {
        int ppEtabin = getppEtabin(TMath::Abs(0.5*(pPbCoMEtabins[pPbCoMEtabin] +pPbCoMEtabins[pPbCoMEtabin+1]) - etaCoM)); // takes the bin center, then shifts "out" of the CoM frame so it looks like we're in the lab frame, then takes the absolute value and finds the pp eta bin (since the pp eta bins are in absolute value of eta)
        RpPbHistArr[ppEtabin]->Add(pPbHistArr[pPbCoMEtabin]);
        //cout << "eta_B: " << pPbCoMEtabins[pPbCoMEtabin] << ", " << pPbCoMEtabins[pPbCoMEtabin+1]; // verifies that we are mapping to the right pp eta bin
        //cout << "; eta_pp: " << ppEtabins[ppEtabin] << ", " << ppEtabins[ppEtabin+1] << endl;
    }
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) RpPbHistArr[ppEtabin]->Scale(1., "width");

    int n = 5;
    for (int ppPtbin = 0; ppPtbin < numppPtEtabins[n]; ppPtbin++) {
        cout << "ppPtbin="<< ppPtbin;
        if (ppHistArr[n]->GetBinContent(ppPtbin+1) != 0) cout << ", r_pPb= " << RpPbHistArr[n]->GetBinContent(ppPtbin+1)/ppHistArr[n]->GetBinContent(ppPtbin+1); 
        cout << ", y=" << RpPbHistArr[n]->GetBinContent(ppPtbin+1);
        cout << ", err=" << RpPbHistArr[n]->GetBinError(ppPtbin+1);
        if (RpPbHistArr[0]->GetBinContent(ppPtbin+1) != 0) cout << ", \% error=" << 100*RpPbHistArr[n]->GetBinError(ppPtbin+1)/RpPbHistArr[n]->GetBinContent(ppPtbin+1);
        cout << ", pp_y=" << ppHistArr[n]->GetBinContent(ppPtbin+1);
        cout << ", pp_y err=" << ppHistArr[n]->GetBinError(ppPtbin+1);
        cout << endl;
        //cout << "bin " << ppPtbin << ": " << ppPtbins[ppPtbin] << ", " << ppPtbins[ppPtbin+1] << endl;
    }

    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        //RpPbHistArr[ppEtabin]->Scale(1., "width");
        RpPbHistArr[ppEtabin]->Divide(ppHistArr[ppEtabin]); // now divide by the pp spectrum
    }

    /**** End add histograms ****/


    /** Plotting routines **/

    // Plot best-selected pt spectra
//    double* histArrScales = linspace(-1.5, 1.5, numppEtabins - 1);
//    double* histArrScales = linspace(1, 1, numppEtabins- 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts
    double histArrScales[numppEtabins] = {0, 0, 0, 0, 0, 0};
    TCanvas* canvas = new TCanvas("RpPbCanvas", "", 800, 600);
    gPad->SetLogx();
    gPad->SetTicks();
    canvas->Draw();
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        thisHist = RpPbHistArr[ppEtabin];

        deta = 2.*(ppEtabins[ppEtabin+1] - ppEtabins[ppEtabin]); // always = 1
        thisHist->Scale(TMath::Power(10, histArrScales[(int)(numppEtabins/2 - 0.5 -TMath::Abs(ppEtabin - numppEtabins/2 + 0.5))])/deta); // separate different ppEtabins

        thisHist->SetMarkerStyle(mkstyles[ppEtabin < (numppEtabins/2)]);
//        thisHist->SetMarkerStyle(kDot);
        Color_t kColor = mkcolors[(147*ppEtabin)%20];
        thisHist->SetMarkerColor(kColor);
        thisHist->SetLineColor(kColor);
        thisHist->SetMinimum(ymin);
        thisHist->SetMaximum(ymax);
        thisHist->GetYaxis()->SetTitleOffset(1.35);
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);
        if (ppEtabin == 0) thisHist->Draw("e1");
        else thisHist->Draw("same e1");
        
        const float textx = 0.46 + (ppEtabin>=(numppEtabins/2))*0.26;
        const float texty = 0.91 - (ppEtabin%(numppEtabins/2))*0.05*(ppEtabin>=(numppEtabins/2)) - (numppEtabins/2 - ppEtabin - 1)*0.05*(ppEtabin<(numppEtabins/2));
        const char* text = Form("%g < #left|#it{#eta}_{CoM}#right| < %g"/* (#times10^{%g})"*/, ppEtabins[ppEtabin], ppEtabins[ppEtabin+1]/*, histArrScales[(int)((0.5*(numppEtabins-1))-TMath::Abs(ppEtabin-(0.5*(numppEtabins-1))))]*/);
        myMarkerText (textx, texty, kColor, mkstyles[ppEtabin < (numppEtabins/2)], text);
    }

    myText (0.19, 0.85, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalpPbLuminosity));
    myText (0.19, 0.79, kBlack, Form("2012 #it{pp}, %.1f fb^{-1}", totalppLuminosity));
    myText (0.19, 0.73, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));

    string histName;
    histName = "R_pPb_combinedTriggers";

    if (runPeriodA && !runPeriodB) {
        myText (0.19, 0.91, kBlack, "Period A");
        histName = histName + "_periodA";
    }
    else if (!runPeriodA && runPeriodB) {
        myText (0.19, 0.91, kBlack, "Period B");
        histName = histName + "_periodB";
    }
    else {
        myText (0.19, 0.91, kBlack, "Period A & B");
    }

    if (runPeriodA || runPeriodB) canvas->SaveAs((plotPath + "R_pPb/" + histName + ".pdf").c_str());

    /**** Free memory and quit ****/
//    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) delete pPbHistArr[ppEtabin];
//    delete canvas;
 //   delete[] histArrScales;
    delete runNumbers;

    if (debugStatements) cout << "Status: In IdealRpPbAnalysisHist.C (154): Finished plotting pt spectrum" << endl;
    return;
}
