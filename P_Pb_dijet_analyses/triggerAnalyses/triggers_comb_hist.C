#include "../triggerUtil.C"
void triggers_comb_hist() {

    initialize(0, false, false);

    double ymin = 5e-8;
    double ymax = 5e6;
    const Style_t mkstyles[2] = {kFullTriangleUp, kFullTriangleDown};
    const Color_t mkcolors[20] = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49};

    double hscale, deta;

    vector<int>* runNumbers = getRunNumbers();
    const int numruns = runNumbers->size();
    double* triggerLuminosities = getTriggerLuminosities();
    TH1D* histArr[numtrigs*numetabins];
    for (Trigger* trig : triggerVec) {
        for (int etabin = 0; etabin < numetabins; etabin++) {
            TString histName = Form("counts_%s_etabin_%i", trig->name.c_str(), etabin);
            histArr[trig->index + etabin*numtrigs] = new TH1D(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}dy #left[nb GeV^{-1}#right]", numpbins, pbins);
            histArr[trig->index + etabin*numtrigs]->Sumw2();
        }
    }

    totalLuminosity = 0;

    /**** Loop over all runs gathering relevant trigger histograms ****/
    for (int thisRunNumber : *runNumbers) {
        TFile* thisFile = new TFile(Form("%srun_%i.root", trigPath.c_str(), thisRunNumber), "READ");
        TVectorD* lum_vec = (TVectorD*)thisFile->Get(Form("lum_vec_%i", thisRunNumber));
        totalLuminosity += (*lum_vec)[0];

        /**** Generate list of physics triggers ****/
        vector<Trigger*> triggerSubList(0);
        for (Trigger* trig : triggerVec) {
            if (trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber && trig->name != minbiasTriggerName) triggerSubList.push_back(trig);
        }
        if (debugStatements) {
            cout << "Status: In triggers_comb_hist.C (39): Processing run " << thisRunNumber << " with triggers:" << endl;
            for (Trigger* trig : triggerSubList) {
                cout << "\t" << trig->name << endl;
            }
        }
        /**** End generate list of physics triggers ****/

        /**** Add histograms from all physics triggers ****/
        for (Trigger* trig : triggerSubList) {
            int actetabin;
            for (int etabin = 0; etabin < numetabins; etabin++) {
                if (thisRunNumber < 313500) actetabin = numetabins - etabin - 1;
                else actetabin = etabin;
                histArr[trig->index + etabin*numtrigs]->Add((TH1D*)thisFile->Get(Form("counts_%s_etabin_%i_run%i", trig->name.c_str(), actetabin, thisRunNumber)));
            }
        }
        /**** End add histograms ****/
        thisFile->Close();
        delete thisFile;
    }
    /**** End loop over runs ****/


    /**** Calculate trigger luminosity rescaling factors ****/
    TH1D* thisHist;
    for (Trigger* trig : triggerVec) {
        for (int etabin = 0; etabin < numetabins; etabin++) {
            thisHist = histArr[trig->index + etabin*numtrigs];
            for (int pbin = 0; pbin < numpbins; pbin++) {
                double thisLumi = 0;
                for (int rn_itr = 0; rn_itr < numruns; rn_itr++) {
                    thisLumi += triggerLuminosities[rn_itr + (trig->index + (pbin + etabin*numpbins)*numtrigs)*numruns];
                }
                if (thisLumi == 0.) continue;
                thisHist->SetBinContent(pbin+1, thisHist->GetBinContent(pbin+1)/thisLumi);
                thisHist->SetBinError(pbin+1, thisHist->GetBinError(pbin+1)/thisLumi);
            }
        }
    }
    /**** End calculate trigger luminosities ****/


    /**** Plotting routines ****/

    const double* histArrScales = linspace(-1.5, 1.5, numetabins/2 - 1);
//    double* histArrScales = linspace(0, 0, numetabins/2 - 1); // for "un-unscaling" to see if one eta bin is particularly lacking in counts

    TCanvas* canvas;
    for (Trigger* trig : triggerVec) {
        canvas = new TCanvas((trig->name + "_canvas").c_str(), "", 800, 600);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetTicks();
        int index = trig->index;
        for (int etabin = 0; etabin < numetabins; etabin++) {
            thisHist = histArr[index + etabin*numtrigs];
            thisHist->Scale(TMath::Power(10, histArrScales[(int)(numetabins/2 - 0.5 -TMath::Abs(etabin - numetabins/2 + 0.5))])/((double)(A)), "width"); // separate different etabins

            thisHist->SetMarkerStyle(mkstyles[etabin < (numetabins/2)]);
//            thisHist->SetMarkerStyle(kDot);
            Color_t kColor = mkcolors[(147*etabin)%20];
            thisHist->SetMarkerColor(kColor);
            thisHist->SetLineColor(kColor);
            thisHist->SetMinimum(ymin);
            thisHist->SetMaximum(ymax);
            thisHist->GetYaxis()->SetTitleOffset(1.35);
            thisHist->GetXaxis()->SetTickLength(0.02);
            thisHist->GetYaxis()->SetTickLength(0.02);

            const float textx = 0.46 + (etabin>=(numetabins/2))*0.26;
            const float texty = 0.91 - (etabin%(numetabins/2))*0.05*(etabin>=(numetabins/2)) - (numetabins/2 - etabin - 1)*0.05*(etabin<(numetabins/2));
            const char* text = Form("%g < #it{#eta}_{B} < %g (#times10^{%g})", etabins[etabin], etabins[etabin+1], histArrScales[(int)((0.5*(numetabins-1))-TMath::Abs(etabin-(0.5*(numetabins-1))))]);
            myMarkerText (textx, texty, kColor, mkstyles[etabin < (numetabins/2)], text);

            if (etabin == 0) thisHist->Draw("e1");
            else thisHist->Draw("same e1");
        }

        canvas->Draw();

        myText (0.19, 0.33, kBlack, trig->name.c_str());
        myText (0.19, 0.27, kBlack, Form("2016 #it{p-Pb}, %.1f nb^{-1}", totalLuminosity*1000.));
        myText (0.19, 0.21, kBlack, Form("#sqrt{#it{s}} = 8.16 TeV"));
        canvas->SaveAs(Form("%striggerPtSpectra/%s.pdf", plotPath.c_str(), trig->name.c_str()));

        /**** Memory cleanup ****/
        delete canvas;
        for (int etabin = 0; etabin < numetabins; etabin++) delete histArr[index + etabin*numtrigs];
        /**** End memory cleanup ****/
    }
    /**** End plotting routines ****/

    if (debugStatements) cout << "Status: In triggers_comb_hist.C (128): Finished plotting trigger pt spectra" << endl;
    return;
}
