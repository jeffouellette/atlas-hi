/**
 * Macros for offline electron analyses during the Nov. 2017 5TeV p-p run.
 * @author Jeff Ouellette
 * @email jeff.ouellette@cern.ch
 */


/***** Variable declarations *****/

// General constants
const float pi = (float)(TMath::Pi());

// Run information
int runNumber; // current run being analyzed
bool isMain;
bool isExpress;
bool isMinBias;

bool considerCharge;

// Analysis information
const bool printStatementChecks = false; // for debugging
const int ptcut = 20;
int useTrigger = 5;
const int numTriggers = 10;
const string triggers[numTriggers] = {"trigger_e10_loose_L1EM7", "trigger_e10_lhloose_L1EM7", "trigger_e15_loose_L1EM7", "trigger_e15_lhloose_L1EM7", "trigger_e15_lhloose_nod0_L1EM7", "trigger_e15_lhloose_L1EM12", "trigger_e15_lhmedium_L1EM12", "trigger_e15_lhloose_nod0_L1EM12", "trigger_e20_loose_L1EM15", "trigger_2e15_lhloose_L12EM12"};

// Binning information
const double etabinshist[33] = {-2.53, -2.37, -2.21, -2.05, -1.89, -1.73, -1.56, -1.37, -1.22, -1.07, -0.92, -0.77, -0.62, -0.47, -0.32, -0.17, 0, 0.17, 0.32, 0.47, 0.62, 0.77, 0.92, 1.07, 1.22, 1.37, 1.56, 1.73, 1.89, 2.05, 2.21, 2.37, 2.53};
const float etabins[7] = {-2.37, -1.56, -1.37, 0, 1.37, 1.56, 2.37};
const float phibins[9] = {(float)(-pi), (float)(-0.75*pi), (float)(-0.5*pi), float(-0.25*pi), (float)(0), (float)(0.25*pi), (float)(0.5*pi), (float)(0.75*pi), (float)(pi)};
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
const int numphibins = sizeof(phibins)/sizeof(phibins[0]) - 1;

// Path information
TFile* outputFile;
string plotPath;
string dataPath;

// Initialize histograms for electron pt spectrum and dielectron invariant mass spectrum
TH1F* electron_ptspectrum;
TH1F* electron_ptspectrum_etabinned[numetabins];
TH1F* electron_ptspectrum_phibinned[numphibins];
TH1F* invariantMass;
TH1F* invariantMass_nocharge;
TH1F* invariantMass_etabinned[numetabins];
TH1F* Z_ptspectrum;
TH1F* Z_ptspectrum_nocharge;
TH2F* eta_phi_hist;
TH2F* eta_phi_hist_no_pt_cut;
TH1F* eta_phi_int_hist;
TH1F* eta_phi_int_hist_no_pt_cut;

// Plotting variables
TCanvas* thiscanvas;
TH1F* thishist; // placeholder
TLatex* latex = new TLatex();
string canvasName;
string run_string;
const Style_t mkstyles[7] = {kFullCircle, kFullDiamond, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCrossX, kFullFourTrianglesPlus};
const Color_t mkcolors[10] = {kRed, kAzure+5, kMagenta, kGreen+2, kBlue, kOrange+8, kTeal+4, kSpring+5, kCyan+2, kPink+4};

// Tree for the current data set being analyzed
TFile* currFile;
TTree* tree;

/***** End variable declarations *****/


/***** General purpose function declarations *****/

/**
 * Returns true if the eta provided is outside the detector limits.
 */
bool etaIsOutsideDetectableRange(float eta) {
    return (TMath::Abs(eta) >= 1.37 && TMath::Abs(eta) <= 1.56) || TMath::Abs(eta) > 2.37;
}


/**
 * Initializes text settings.
 */
void initialize_text () {
    run_string = Form("Using run %i", runNumber);
    if (isMain) run_string += " main stream";
    else if (isMinBias) run_string += " min bias";
    else if (isExpress) run_string += " express stream";
    return;
}


/**
 * Initializes a new canvas.
 */
void initialize_new_canvas (bool logy, float left_margin=0.1, float right_margin=0.1, float bottom_margin=0.1, float top_margin=0.1) {

    thiscanvas = new TCanvas (canvasName.c_str(), "", 1000, 800);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    if (logy) gPad->SetLogy();
    gPad->SetTicks();
    thiscanvas->SetMargin (left_margin, right_margin, bottom_margin, top_margin);
    thiscanvas->Draw();
    return;
}


/**
 * Draws run information and trigger selection on the plot.
 */
void draw_information (float ndc_x, float ndc_y) {
    latex->SetTextAlign(22);
    latex->SetTextFont(42);
    latex->SetTextSize(0.028);
    latex->DrawLatexNDC (ndc_x, ndc_y-0.025, Form("Tight electrons, filling on %s", triggers[useTrigger].c_str()));
    latex->DrawLatexNDC (ndc_x, ndc_y+0.025, run_string.c_str());
    return;
}


/**
 * Draws a plot title.
 */
void draw_title (float ndc_x, float ndc_y, string title) {
    latex->SetTextAlign(22);
    latex->SetTextFont(42);
    latex->SetTextSize(0.036);
    latex->DrawLatexNDC (ndc_x, ndc_y, title.c_str());
    return;
}


/**
 * Draws additional text.
 */
void draw_text (float ndc_x, float ndc_y, string text) {
    latex->SetTextAlign(22);
    latex->SetTextFont(42);
    latex->SetTextSize(0.028);
    latex->DrawLatexNDC (ndc_x, ndc_y, text.c_str());
    return;
}

/***** End general purpose functions *****/


/***** Histogram-specific plotting functions *****/

/**
 * Initializes the histograms to be generated while looping over the trees.
 */
void initialize_histograms() {

    eta_phi_hist = new TH2F (Form("run_%i_eta_phi_hist", runNumber), ";#eta_{lab};#phi;d^{2}N/d#phid#eta_{lab}", 32, etabinshist, 50, -pi, pi);
    eta_phi_hist->Sumw2();
    eta_phi_int_hist = new TH1F(Form("run_%i_eta_phi_int_hist", runNumber), ";#eta_{lab};dN/d#eta", 32, etabinshist);
    eta_phi_hist_no_pt_cut = new TH2F (Form("run_%i_eta_phi_hist_no_pt_cut", runNumber), ";#eta_{lab};#phi;d^{2}N/d#phid#eta_{lab}", 32, etabinshist, 50, -pi, pi);
    eta_phi_hist_no_pt_cut->Sumw2();
    eta_phi_int_hist_no_pt_cut = new TH1F(Form("run_%i_eta_phi_int_hist_no_pt_cut", runNumber), ";#eta_{lab};dN/d#eta", 32, etabinshist);

    electron_ptspectrum = new TH1F (Form("run_%i_electron_ptspectrum_hist", runNumber), ";#it{p}_{T}^{electron} #left[GeV/#it{c}#right];dN/d#it{p}_{T} #left[(GeV/#it{c})^{-1}#right]", 50, 20, 70);
    electron_ptspectrum->Sumw2();

    invariantMass = new TH1F (Form("run_%i_invariantMass_hist", runNumber), ";#it{M}_{ee} #left[GeV/#it{c}^{2}#right];dN/d#it{M}_{ee} #left[(GeV/#it{c}^{2})^{-1}#right]", 50, 50, 140);
    invariantMass->Sumw2();
    invariantMass_nocharge = new TH1F (Form("run_%i_invariantMass_nocharge_hist", runNumber), ";#it{M}_{ee} #left[GeV/#it{c}^{2}#right];dN/d#it{M}_{ee} #left[(GeV/#it{c}^{2})^{-1}#right]", 50, 50, 140);
    invariantMass_nocharge->Sumw2();

    Z_ptspectrum = new TH1F(Form("run_%i_Z_ptspectrum", runNumber), ";#it{p}_{T}^{Z} #left[GeV/#it{c}#right];dN/d#it{p}_{T} #left[(GeV/#it{c})^{-1}#right]", 50, 0, 250);
    Z_ptspectrum_nocharge = new TH1F(Form("run_%i_Z_ptspectrum_nocharge", runNumber), ";#it{p}_{T}^{Z} #left[GeV/#it{c}#right];dN/d#it{p}_{T} #left[(GeV/#it{c})^{-1}#right]", 50, 0, 250);

    for (int etabin = 0; etabin < numetabins; etabin++) {
        electron_ptspectrum_etabinned[etabin] = new TH1F (Form("run_%i_electron_ptspectrum_hist_etabin%i", runNumber, etabin), ";#it{p}_{T}^{electron} #left[GeV/#it{c}#right];d^{2}N/d#it{p}_{T}dy #left[(GeV/#it{c})^{-1}#right]", 50, 20, 70);
        electron_ptspectrum_etabinned[etabin]->Sumw2();
        invariantMass_etabinned[etabin] = new TH1F (Form("run_%i_invariantMass_hist_etabin%i", runNumber, etabin), ";#it{M}_{ee} #left[GeV/#it{c}^{2}#right];d^{2}N/d#it{M}_{ee}dy #left[(GeV/#it{c}^{2})^{-1}#right]", 50, 50, 140);
        invariantMass_etabinned[etabin]->Sumw2();
    }
    for (int phibin = 0; phibin < numphibins; phibin++) {
        electron_ptspectrum_phibinned[phibin] = new TH1F (Form("run_%i_electron_ptspectrum_hist_phibin%i", runNumber, phibin), ";#it{p}_{T}^{electron} #left[GeV/#it{c}#right];d^{2}N/d#it{p}_{T}d#phi #left[(GeV/#it{c})^{-1}#right]", 50, 20, 70);
        electron_ptspectrum_phibinned[phibin]->Sumw2();
    }
    return;
}


/**
 * Plotting routine for the electron ptspectrum.
 */
void save_electron_ptspectrum () {
    thishist = electron_ptspectrum;
    thishist->Write();
    return;
}


/**
 * Plotting routine for the dielectron invariant mass.
 */
void save_invariantMass () {
    thishist = invariantMass;
    thishist->Write();
    thishist = invariantMass_nocharge;
    thishist->Write();
    return;
}


/**
 * Plotting routine for the dielectron invariant (Z) pt spectrum.
 */
void save_Z_ptspectrum () {
    thishist = Z_ptspectrum;
    thishist->Write();
    thishist = Z_ptspectrum_nocharge;
    thishist->Write();
    return;
}


/**
 * Plotting routine for the ptspectrum, eta binned.
 */
void save_electron_ptspectrum_etabinned () { 
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thishist = electron_ptspectrum_etabinned[etabin];
        thishist->Write();
    }
    return;
}


/**
 * Plotting routine for the ptspectrum, phi binned.
 */
void save_electron_ptspectrum_phibinned () {
    for (int phibin = 0; phibin < numphibins; phibin++) {
        thishist = electron_ptspectrum_phibinned[phibin];
        thishist->Write();
    }
    return;
}


/**
 * Plotting routine for the dielectron invariant mass, eta binned.
 */
void save_invariantMass_etabinned () {
    for (int etabin = 0; etabin < numetabins; etabin++) {
        thishist = invariantMass_etabinned[etabin];
        thishist->Write();
    }
    return;
}


/**
 * Plotting routine for the eta, phi distribution of electrons.
 */
void save_eta_phi (bool boolptcut) {

    TH2F* this2hist;
    if (boolptcut) {
        this2hist = eta_phi_hist;
        thishist = eta_phi_int_hist;
    }
    else {
        this2hist = eta_phi_hist_no_pt_cut;
        thishist = eta_phi_int_hist_no_pt_cut;
    }
    
    // First calculate the integral over phi before we make any scalings.
    for (int i = 0; i < 50; i++) {
        float sum = 0;
        float variance = 0;
        for (int j = 0; j < 50; j++) {
            sum += this2hist->GetBinContent(i+1, j+1);
            variance += TMath::Power(this2hist->GetBinError(i+1, j+1), 2);
        }
        thishist->SetBinContent(i+1, sum);
        thishist->SetBinError(i+1, TMath::Sqrt(variance));
    }

    thishist->Write();
    this2hist->Write();

    return;
}

/***** End histogram-specific plotting functions *****/


/***** Main macro code *****/

/**
 * Loops over the events in the tree and fills the appropriate histograms under the specified conditions.
 */
void treeLoop () {

    // Tree dependent variables
    const int numentries = tree->GetEntries();

    // Branch dependent variables
    float electron_pt[100] = {};
    float electron_eta[100] = {};
    float electron_phi[100] = {};
    bool electron_tight[100] = {};
    bool electron_loose[100] = {};
    float electron_etcone40[100] = {};
    int electron_charge[100] = {};
    int electron_n = 0;
    bool trig_bool[numTriggers];
    float trig_prescale[numTriggers];
    int njet = 0;
    float j_pt[100] = {};
    float j_phi[100] = {};
    float j_eta[100] = {};
   
    // Branching address declarations
    tree->SetBranchAddress("electron_tight", electron_tight);
    tree->SetBranchAddress("electron_loose", electron_loose);
    tree->SetBranchAddress("electron_pt", electron_pt);
    tree->SetBranchAddress("electron_n", &electron_n);
    tree->SetBranchAddress("electron_phi", electron_phi);
    tree->SetBranchAddress("electron_eta", electron_eta);
    tree->SetBranchAddress("electron_etcone40", electron_etcone40);
    tree->SetBranchAddress("electron_charge", electron_charge);
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress((triggers[useTrigger]).c_str(), &(trig_bool[useTrigger])); 
//    tree->SetBranchAddress((triggers[useTrigger]+"_prescale").c_str(), &(trig_prescale[useTrigger]);
   
    // Loop objects
    TLorentzVector tvec1;
    TLorentzVector tvec2;
    TLorentzVector tvec3;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i);
        if (trig_bool[useTrigger]) { // Require some trigger to be satisfied. Generally, this is a dielectron or lhloose trigger.
            for (int e1 = 0; e1 < electron_n; e1++) {
                
                if (!electron_tight[e1]) continue; // Require TIGHT electron condition on the 1st electron 
                
                int etabin = 0; // Binning
                while (etabins[etabin] < electron_eta[e1]) etabin++;
                etabin--;
                int phibin = 0; // Binning
                while (phibins[phibin] < electron_phi[e1]) phibin++;
                phibin--;

                if (!etaIsOutsideDetectableRange(electron_eta[e1])) { // Only fill electron pt if eta falls within the detector
                    electron_ptspectrum->Fill(electron_pt[e1]);
                    electron_ptspectrum_etabinned[etabin]->Fill(electron_pt[e1]);
                    electron_ptspectrum_phibinned[phibin]->Fill(electron_pt[e1]);
                }
                
                if (electron_pt[e1] > ptcut) eta_phi_hist->Fill(electron_eta[e1], electron_phi[e1]); // Only fill electron eta,phi if electron meets kinematic (pt) cut
                eta_phi_hist_no_pt_cut->Fill(electron_eta[e1], electron_phi[e1]);

                tvec1.SetPtEtaPhiM(electron_pt[e1], electron_eta[e1], electron_phi[e1], 0.511e-3); // Store 1st electron information in a TLorentzVector
                for (int e2 = e1+1; e2 < electron_n; e2++) {
                    if (!electron_tight[e2]) continue; // Require TIGHT electron condition on 2nd electron

                    if (etaIsOutsideDetectableRange(electron_eta[e1]) || etaIsOutsideDetectableRange(electron_eta[e2])) continue; // Require etas to fall within the detector for dielectron (Z) studies

                    if (electron_pt[e1] < ptcut || electron_pt[e2] < ptcut) continue; // Require high pt condition

                    tvec2.SetPtEtaPhiM(electron_pt[e2], electron_eta[e2], electron_phi[e2], 0.511e-3); // Store 2nd electron information in a TLorentzVector

                    if (electron_pt[e1] > electron_pt[e2]) { // rebin in the trigger with higher pt
                        etabin = 0;
                        while (etabins[etabin] < electron_eta[e1]) etabin++;
                        etabin--;
                    } else {
                        etabin = 0;
                        while (etabins[etabin] < electron_eta[e2]) etabin++;
                        etabin--;
                    }

                    tvec3 = tvec1 + tvec2; // Form Lorentz invariant vector for dielectron system
                    if (electron_charge[e1] * electron_charge[e2] < 0) { // Z bosons are neutral so the electrons must have opposite charge
                        invariantMass->Fill(tvec3.M());
                        invariantMass_etabinned[etabin]->Fill(tvec3.M()); // Fill invariant mass (length of invariant vector)
                        Z_ptspectrum->Fill(tvec3.Pt()); // Fill Z Pt spectrum (Pt of invariant vector)
                    }
                    Z_ptspectrum_nocharge->Fill(tvec3.Pt());
                    invariantMass_nocharge->Fill(tvec3.M());
                    for (int j = 0; j < njet; j++) {
                        float smallAngleDiff = TMath::Abs(j_phi[j] - tvec3.Phi());
                        while (smallAngleDiff > pi) smallAngleDiff = TMath::Abs(smallAngleDiff - 2*pi);
                        if (electron_charge[e1] * electron_charge[e2] < 0 && j_pt[j] > 100 && tvec3.Pt() > 100 && (j_pt[j]-tvec3.Pt())/(j_pt[j]+tvec3.Pt()) < 0.1 && smallAngleDiff > 7.*pi/8. && tvec3.Eta() * j_eta[j] < 0) {
                            cout << Form("Found back-to-back Z + jet in run %i, event number %i! j_pt = %.1f GeV, Z_pt = %.1f GeV, j_phi-Z_phi = %.1f pi, j_eta = %.1f, Z_eta = %.1f", runNumber, i, j_pt[j], tvec3.Pt(), smallAngleDiff/pi, j_eta[j], tvec3.Eta()) << endl;
                        }
                    }
                }
            }
        }
    } 
    return; 
}


/**
 * Fills the tree with a new file.
 */
void getTree (std::string pathToTreeFile) {
    TFile* currFile = new TFile(pathToTreeFile.c_str(), "READ");
    tree = (TTree*)(currFile->Get("etree"));
    return;
}


/**
 * Main driving macro for analyzing offline electrons. dataStream can be 'e' for express, 'm' for main, or 'b' for minbias data streams.
 */
void electronOfflineAnalysis (int rn, int trig, char dataStream) {

    runNumber = rn;
    considerCharge = rn != 340634;

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

    initialize_text ();
    initialize_histograms ();    

    const std::vector<char*> filesToProcess = {
        Form("%suser.khill.pp5TeVAna_main.002.00340634.f895_m1902_myOutput_hadd.root", dataPath.c_str()),
        Form("%suser.khill.pp5TeVAna_main.004.00340644.f895_m1902_myOutput_hadd.root", dataPath.c_str()),
        Form("%suser.khill.pp5TeVAna_main.004.00340683.f896_m1902_myOutput_hadd.root", dataPath.c_str()),
        Form("%suser.khill.pp5TeVAna_main.004.00340697.f896_m1902_myOutput_hadd.root", dataPath.c_str())
    };
    std::vector<int> fileRunNumbers = {
        340634,
        340644,
        340683,
        340697
    };
    
    for (int i = 0; i < filesToProcess.size(); i++) {
        if (runNumber != fileRunNumbers[i]) continue;
        getTree (filesToProcess[i]);
        treeLoop ();
    }

    outputFile = new TFile (Form("./Data/%i_%s_out.root", runNumber, triggers[useTrigger].c_str()), "RECREATE");

    save_electron_ptspectrum ();
    save_invariantMass ();
    save_Z_ptspectrum (); 
    save_electron_ptspectrum_etabinned ();
    save_electron_ptspectrum_phibinned ();
    save_invariantMass_etabinned ();
    save_eta_phi (true);
    save_eta_phi (false);

    outputFile->Close();

    return;
}

/***** End main macro code *****/
