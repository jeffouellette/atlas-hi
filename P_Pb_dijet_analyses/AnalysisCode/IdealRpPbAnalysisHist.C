#include "IdealRpPbAnalysis.C"

void IdealRpPbAnalysisHist() {

    TH1D** ppHistArr = setupPPConfiguration(); // creates the pp histograms and setups the relevant environment variables

    initialize(0, true);
    std::vector<int>* thisRunNumbers = getRunNumbers();

    const int numruns = (*thisRunNumbers).size();
    const int numhists = numtrigs * numruns * numetabins;
    if (debugStatements) {
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint A): Building trigger pt histograms with " << numruns << " runs being used" << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint B): Numtrigs = " << numtrigs << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint C): Numetabins = " << numppEtabins << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint D): Numpbins = " << numpbins << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint E): ptPath = " << ptPath << endl;
    }

    const bool shiftBins = false;
    double ymin = 0;
    double ymax = 5;
    if (shiftBins) ymax = 15;
    const Style_t mkstyles[2] = {kFullCircle, kOpenCircle};
    const Color_t mkcolors[8] = {kOrange-3, kGreen, kCyan, kBlue, kViolet, kRed, kMagenta, kYellow};

    TH1D* thisHist;

    if (debugStatements) cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint F): Initialized pPbHistArr histograms..." << endl;
    TH1D* pPbHistArr[numpPbCoMEtabins];
    for (int pPbCoMEtabin = 0; pPbCoMEtabin < numpPbCoMEtabins; pPbCoMEtabin++) {
        int ppEtabin_equiv = getppEtabin(TMath::Abs(0.5*(pPbCoMEtabins[pPbCoMEtabin]+pPbCoMEtabins[pPbCoMEtabin+1]) - etaCoM));
        pPbHistArr[pPbCoMEtabin] = new TH1D(Form("pPb_spectrum_pPbCoMEtabin%i", pPbCoMEtabin), ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T} #left[nb GeV^{-1}#right]", numppPtEtabins[ppEtabin_equiv], ppPtbins);
        pPbHistArr[pPbCoMEtabin]->Sumw2();
    }
    TH1D* RpPbHistArr[numppEtabins];
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        RpPbHistArr[ppEtabin] = new TH1D(Form("R_pPb_ppEtabin%i", ppEtabin), ";#it{p}_{T}^{jet} [GeV];R_{#it{pPb}}", numppPtEtabins[ppEtabin], ppPtbins);
        RpPbHistArr[ppEtabin]->Sumw2();
    }

    if (debugStatements) {
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint G): pPbHistArr histograms initialized." << endl;
        cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint H): Starting loop over triggers..." << endl;
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
                    if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl; 
                    for (int thisRunNumber : *runNumbers) {
                        if (skipRun(thisRunNumber)) continue;
                        if (fname.Contains(to_string(thisRunNumber))) {
                            TFile* thisFile = new TFile(RpPbPath + fname, "READ");
//                            lum_vec = (TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber));
//                            cout << lum_vec[0] << endl;
//                            totalpPbLuminosity += lum_vec[0]; 
//                            totalpPbLuminosity += ((TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber)))[0];

                            // quickly check the parameters stored in this root file
                            run_vec = (TVectorD*)thisFile->Get(Form("run_vec_%i", thisRunNumber));
                            assert ((*run_vec)[0] == thisRunNumber);
                            assert ((*run_vec)[1] == numppEtabins);
                            assert ((*run_vec)[2] == numpPbCoMEtabins);
                            assert ((*run_vec)[3] == numtrigs);
                            assert ((*run_vec)[4] == numppPtbins);
                            totalpPbLuminosity += (*run_vec)[5];

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
    }
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) RpPbHistArr[ppEtabin]->Scale(1., "width");


    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        const double deta = 2.*(ppEtabins[ppEtabin+1] - ppEtabins[ppEtabin]); // always = 1
        RpPbHistArr[ppEtabin]->Scale(1/deta); 
        RpPbHistArr[ppEtabin]->Divide(ppHistArr[ppEtabin]); // now divide by the pp spectrum
    }

    /**** End add histograms ****/


    /** Plotting routines **/

    // Plot best-selected pt spectra
    const double histArrShifts[numppEtabins] = {0, 1.5, 3, 4.5, 6, 7.5};
    //const double histArrShifts[numppEtabins] = {7.5, 6, 4.5, 3, 1.5, 0};
    TCanvas* canvas = new TCanvas("RpPbCanvas", "", 800, 600);
    canvas->Divide(2, 3, 0.01, 0.01);
    TLine* lineDrawer = new TLine();  
    gPad->SetLogx();
    gPad->SetTicks();
    canvas->Draw();
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        canvas->cd(ppEtabin+1);
        thisHist = RpPbHistArr[ppEtabin];
        Style_t kStyle = mkstyles[ppEtabin%2];
        Color_t kColor = mkcolors[ppEtabin%8];

        if (shiftBins) {
            for (int ppPtbin = 0; ppPtbin < numppPtbins; ppPtbin++) {
                RpPbHistArr[ppEtabin]->SetBinContent(ppPtbin+1, RpPbHistArr[ppEtabin]->GetBinContent(ppPtbin+1) + histArrShifts[ppEtabin]); // shifts each histogram appropriately
            }
            lineDrawer->SetLineColor(kColor);
            lineDrawer->SetLineStyle(9);
            //lineDrawer->DrawLine(ppPtbins[0], histArrShifts[ppEtabin]+1, ppPtbins[numppPtbins], histArrShifts[ppEtabin]+1);
        }

        thisHist->SetMarkerStyle(kStyle);
        thisHist->SetMarkerColor(kColor);
        thisHist->SetLineColor(kColor);
        thisHist->SetMinimum(ymin);
        thisHist->SetMaximum(ymax);
        thisHist->GetYaxis()->SetLabelSize(0.07);
        thisHist->GetYaxis()->SetTitleSize(0.07);
        thisHist->GetXaxis()->SetLabelSize(0.07);
        thisHist->GetXaxis()->SetTitleSize(0.07);

        thisHist->GetYaxis()->SetTitleOffset(1.1);
        thisHist->GetXaxis()->SetTitleOffset(1.1);
        thisHist->GetXaxis()->SetTickLength(0.02);
        thisHist->GetYaxis()->SetTickLength(0.02);
        if (ppEtabin == 0) thisHist->Draw("e1");
        else thisHist->Draw("same e1");

        if (shiftBins) lineDrawer->DrawLine(ppPtbins[0], histArrShifts[ppEtabin]+1, ppPtbins[numppPtbins], histArrShifts[ppEtabin]+1);

        const float textx = 0.6; //- 0.04*shiftBins;
        const float texty = 0.85; // - (ppEtabin)*0.05;
        char* text;
        if (shiftBins) text = Form("%g < #left|#it{y}#right|_{CoM} < %g (+%g)", ppEtabins[ppEtabin], ppEtabins[ppEtabin+1], histArrShifts[ppEtabin]); /*histArrScales[(int)((0.5*(numppEtabins-1))-TMath::Abs(ppEtabin-(0.5*(numppEtabins-1))))]);*/
        else text = Form("%g < #left|#it{y}#right|_{CoM} < %g", ppEtabins[ppEtabin], ppEtabins[ppEtabin+1]);
        myText (textx, texty, kBlack, text,  0.1);
    }

    myText (0.19, 0.35, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}, #sqrt{#it{s}} = 8.16 TeV", totalpPbLuminosity), 0.1);
    myText (0.19, 0.225, kBlack, Form("2012 #it{pp}, %.1f fb^{-1}, #sqrt{#it{s}} = 8 TeV", totalppLuminosity), 0.1);

    string histName;
    histName = "R_pPb_combinedTriggers";

    if (shiftBins) histName = histName + "_shifted";

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

    if (debugStatements) cout << "Status: In IdealRpPbAnalysisHist.C (breakpoint J): Finished plotting pt spectrum" << endl;
    return;
}
