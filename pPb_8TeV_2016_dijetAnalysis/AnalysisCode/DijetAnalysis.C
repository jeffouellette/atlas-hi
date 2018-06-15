#include "../Util.C"

const int numxbins = 50;
const int numqbins = 8;
const int numq2bins = 120;
const int numq2xbins = 100;
const int nummbins = 50;
const int numfcalbins = 60;

const double* xbins = logspace(1.6e-4, 1.6, numxbins);
const double* qbins = logspace(20, 1200, numqbins);
const double* q2bins = logspace(1, 1000000, numq2bins);
const double* q2xbins = logspace(1.6e-4, 1.6, numq2xbins);
const double* mbins = logspace(20, 2500, nummbins);
const double* fcalbins = logspace(10, 500, numfcalbins);

void DijetAnalysis(const int dataSet, // Data set identifier. If not MC, this should be a run number. If MC, this should be whatever number follows "tid" in the MC file.
                   const double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
  initialize(dataSet, true);

  const bool isMC = useDataVersion == 0;

  // Check whether to skip this particular analysis
  if (!isMC && skipRun(dataSet)) return;
  if (isMC && skipMC(dataSet)) return;

  vector<TF1*>* triggerEfficiencyFunctions = NULL;
  if (!isMC) triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();

  const bool periodA = isPeriodA(dataSet);
  const bool useIonTriggers = !isMC && dataSet < 313629; // false for MC

  /**** Generate list of physics triggers ****/
  vector<Trigger*>* triggerSubvector = NULL;
  if (!isMC) {
    triggerSubvector = getTriggerSubvector(dataSet);
    if (debugStatements) {
      cout << "Status: In DijetAnalysis.C (breakpoint A): Processing run " << dataSet << " with triggers:" << endl;
      for (Trigger* trig : (*triggerSubvector)) {
        cout << "\t" << trig->name << endl;
      }
    }
  }
  /**** End generate list of physics triggers ****/


  /**** Find the relevant TTree for this run ****/
  //TTree* tree = (TTree*)(new TFile(Form("%srun_%i_raw.root", dataPath.c_str(), dataSet)))->Get("tree");
  TFile* inputFile = NULL;
  TTree* tree = NULL;
  const bool containsPeriodLetter = useIonTriggers;
  {
    TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
    TList* sysfiles = dir.GetListOfFiles();
    if (sysfiles) {
      TSystemFile* sysfile;
      TString fname;
      TIter next(sysfiles);

      while ((sysfile = (TSystemFile*)next())) {
        fname = sysfile->GetName();
        if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
          if (debugStatements) cout << "Status: In DijetAnalysis.C (breakpoint B): Found " << fname.Data() << endl; 
          if (fname.Contains(to_string(dataSet))) {
            //containsPeriodLetter = periodA || fname.Contains(".104");
            inputFile = new TFile(dataPath+fname, "READ");
            tree = (TTree*)(inputFile->Get("tree"));
            break;
          }
        }
      }
    }
  }
  if (tree == NULL) {
    cout << "Error: In DijetAnalysis.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
    return;
  }
  /**** End find TTree ****/


  /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
  {
    vector<string> interestingBranchNames = {"njet", "jet_pt", "jet_eta", "jet_phi", "jet_e", "eventNumber", "fcalA_et", "fcalC_et", "nvert", "vert_type", "ntrk", "trk_quality_4", "trk_d0", "trk_z0", "trk_theta", "trk_charge", "trk_pt", "trk_eta", "trk_phi", "vert_x", "vert_y", "vert_z"};
    TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
    bool interestingBranch;
    for (TObject* obj : *branches) {
      TString branchName = (TString)obj->GetName();
      if (debugStatements) cout << "Status: In DijetAnalysis.C (breakpoint D): Tree contains branch \"" << branchName.Data() << "\"" << endl;
      interestingBranch = false;
      for (string s : interestingBranchNames) {
        interestingBranch = interestingBranch || (branchName.Data() == s);
      }
      if (!interestingBranch && !isMC) {
        for (Trigger* trig : (*triggerSubvector)) {
          if (branchName == trig->name || (containsPeriodLetter && branchName == trig->name + "_prescale_A") || (!containsPeriodLetter && branchName == trig->name + "_prescale")) {
            interestingBranch = true;
            break;
          }
        }
      }
      if (!interestingBranch) {
        tree->SetBranchStatus(branchName, 0);
      }
    }
  }
  /**** End disable unimportant branches ****/


  /**** Calculate eta-phi rescaling information ****/
  TH2D* etaPhiScaleFactorsHist = NULL;
  if (!isMC) {
    /**** Load eta-phi jet correlation plot ****/
    TFile* etaPhiFile = new TFile((rootPath + "etaPhiHist.root").c_str(), "READ");
    TH2D* etaPhiHist = (TH2D*)etaPhiFile->Get("etaPhiHist");
    TH2D* subleadingEtaPhiHist = (TH2D*)etaPhiFile->Get("subleadingEtaPhiHist");
    /**** End load eta-phi correlation ****/


    /**** Evaluate the eta-phi rescaling factor at each point so we don't have to do this for each event ****/
    int nbins_x = etaPhiHist->GetNbinsX();
    int nbins_y = etaPhiHist->GetNbinsY();
    //double* etaPhiScaleFactors = new double[nbins_x*nbins_y];
    etaPhiScaleFactorsHist = new TH2D("etaPhiScaleFactorsHist", "", 98, -4.9, 4.9, 100, 0, 2*pi);

    // for each point in space where the leading jet COULD be, evaluate what the rescaling should be
    for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
      double x = etaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
      int etabin = getEtabin(x);
      for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
        double y = etaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);

        // if the leading jet is in the HEC region don't bother finding the scale factor since it won't be selected
        if (lowerEtaCut < x && x < upperEtaCut && lowerPhiCut < y && y < upperPhiCut) continue;

        // First calculate the number of subleading jets you missed given the position of the leading jet.
        // This requires knowing the number of subleading jets meeting their event selection cuts.
        double x_prime, y_prime, dy, content, numerator, denominator;

        numerator = 0;
        denominator = 0;
        for (int bin_y_prime = 0; bin_y_prime < nbins_y; bin_y_prime++) {
          y_prime = subleadingEtaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
          dy = TMath::Abs(y - y_prime);
          if (dy > pi) dy = 2*pi - dy;
          if (dy < 7.*pi/8.) continue; // this checks whether our y' coordinate is further away from y by at least 7pi/8 (up to pi). Otherwise it is not in our integration region, we skip it.
          for (int bin_x_prime = 0; bin_x_prime < nbins_x; bin_x_prime++) {
            x_prime = subleadingEtaPhiHist->GetXaxis()->GetBinCenter(bin_x_prime+1);
            // now check if x' meets the eta cut requirements
            content = subleadingEtaPhiHist->GetBinContent(bin_x_prime+1, bin_y_prime+1);
            // if the ' coordinate is outside the HEC region, then add the counts there to your integral in the denominator
            if (!(lowerEtaCut < x_prime && x_prime < upperEtaCut && lowerPhiCut < y_prime && y_prime < upperPhiCut)) denominator += content;
            // as long as the ' coordinate meets the dphi cut, add the counts there to your integral in the numerator
            numerator += content;
          }
        }
        if (denominator == 0. && debugStatements) cout << "Warning: In DijetAnalysis.C (breakpoint E): 0 jets meeting HEC cut!" << endl;
        else if (denominator != 0.) etaPhiScaleFactorsHist->SetBinContent(bin_x+1, bin_y+1, numerator/denominator);
        else etaPhiScaleFactorsHist->SetBinContent(bin_x+1, bin_y+1, 0);

        // Calculate the number of leading jets you missed in this etabin. Often this factor comes out to 1.
        numerator = 0;
        denominator = 0;
        for (int bin_x_prime = 0; bin_x_prime < nbins_x; bin_x_prime++) {
          x_prime = etaPhiHist->GetXaxis()->GetBinCenter(bin_x_prime+1);
          if(!(etabins[etabin] < x_prime && x_prime < etabins[etabin+1])) continue;
          for (int bin_y_prime = 0; bin_y_prime < nbins_y; bin_y_prime++) {
            y_prime = etaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
            content = etaPhiHist->GetBinContent(bin_x_prime+1, bin_y_prime+1);
            if (!(lowerEtaCut < x_prime && x_prime < upperEtaCut && lowerPhiCut < y_prime && y_prime < upperPhiCut)) denominator += content;
            numerator += content;
          }
        }
        if (denominator == 0. && debugStatements) cout << "Warning: In DijetAnalysis.C (breakpoint F): 0 jets meeting HEC cut!" << endl;
        else if (denominator != 0.) etaPhiScaleFactorsHist->SetBinContent(bin_x+1, bin_y+1, etaPhiScaleFactorsHist->GetBinContent(bin_x+1, bin_y+1) * numerator/denominator);
        else etaPhiScaleFactorsHist->SetBinContent(bin_x+1, bin_y+1, 0);
      }
    }
    etaPhiFile->Close();
    if (etaPhiFile) delete etaPhiFile;
  }
  if (etaPhiScaleFactorsHist == NULL) {
    cout << "Error: In DijetAnalysis.C (breakpoint C): etaPhiScaleFactorsHist not obtained for given run number. Quitting." << endl;
    return;
  }
  /**** End evaluate eta-phi rescaling ****/


  /**** Create output TTree ****/
  string output_name = xPath;
  if (runPeriodA && !runPeriodB) output_name = output_name + "periodA/";
  else if (!runPeriodA && runPeriodB) output_name = output_name + "periodB/";
  else output_name = output_name + "periodAB/";
  output_name = Form("%sdataset_%i.root", output_name.c_str(), dataSet);

  TFile* outputFile = new TFile(output_name.c_str(), "RECREATE");
  TTree* outTree = tree->CloneTree(0);//new TTree("trackTree", "trackTree");
  /**** End create output TTree ****/


  /**** Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a, as well as other histograms ****/
  TH1D* xHistArr[2*numetabins];
  TH1D* xqHistArr[2*numqbins];
  TH1D* mHistArr[numetabins];

  for (int etabin = 0; etabin < numetabins; etabin++) {
    xHistArr[etabin] = new TH1D(Form("%ieta%i", dataSet, etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#eta #left[nb#right]", etabins[etabin], etabins[etabin+1]), numxbins, xbins);
    xHistArr[etabin]->Sumw2();
    xHistArr[etabin]->SetDirectory(outputFile);
  }
  for (int etabin = numetabins; etabin < 2*numetabins; etabin++) {
    xHistArr[etabin] = new TH1D(Form("%ieta%i", dataSet, etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d^{2}N/L_{int}d#it{x}_{a}d#eta #left[nb#right]", etabins[etabin%numetabins], etabins[(etabin%numetabins)+1]), numxbins, xbins);
    xHistArr[etabin]->Sumw2();
    xHistArr[etabin]->SetDirectory(outputFile);
  }

  for (int qbin = 0; qbin < numqbins; qbin++) {
    xqHistArr[qbin] = new TH1D(Form("%iq%i", dataSet, qbin), Form("%g < #it{Q} < %g;#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin], qbins[qbin+1]), numxbins, xbins);
    xqHistArr[qbin]->Sumw2();
    xqHistArr[qbin]->SetDirectory(outputFile);
  }
  for (int qbin = numqbins; qbin < 2*numqbins; qbin++) {
    xqHistArr[qbin] = new TH1D(Form("%iq%i", dataSet, qbin), Form("%g < #it{Q} < %g;#it{x}_{a};d^{2}N/L_{int}d#it{x}_{a}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin%numqbins], qbins[(qbin%numqbins)+1]), numxbins, xbins);
    xqHistArr[qbin]->Sumw2();
    xqHistArr[qbin]->SetDirectory(outputFile);
  }

  for (int etabin = 0; etabin < numetabins; etabin++) {
    mHistArr[etabin] = new TH1D(Form("mjj_%ieta%i", dataSet, etabin), Form("%1.1f < #eta < %1.1f;#it{m}_{JJ} #left[GeV#right];d^{2}N/L_{int}d#it{m}_{JJ}d#eta #left[nb GeV^{-1}#right]", etabins[etabin], etabins[etabin+1]), nummbins, mbins);
    mHistArr[etabin]->Sumw2();
    mHistArr[etabin]->SetDirectory(outputFile);
  }

  TH2D* qxcorr = new TH2D(Form("xqcorr_dataset%i", dataSet), ";#it{x}_{a};#it{#bar{Q}}^{2} #left[GeV^{2}#right];d^{2}N/L_{int}d{x}_{a}d#it{#bar{Q}}^{2} #left[nb GeV^{-2}#right]", numq2xbins, q2xbins, numq2bins, q2bins);
  qxcorr->SetDirectory(outputFile);
  TH2D* xaxpcorr = new TH2D(Form("xaxpcorr_dataset%i", dataSet), ";#it{x}_{a};#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#it{x}_{a}", numxbins, xbins, numxbins, xbins);
  xaxpcorr->SetDirectory(outputFile);
  TH2D* fcalhist = new TH2D(Form("fcalhist_dataset%i", dataSet), ";#it{x}_{p};FCAL energy #left[GeV#right];", numxbins, xbins, numfcalbins, fcalbins);
  fcalhist->SetDirectory(outputFile);

  TH1D* leadingJetEtaHist = new TH1D(Form("leadingJetEtaHist_dataSet%i", dataSet), ";#eta;Counts", 98, -4.9, 4.9);
  leadingJetEtaHist->SetDirectory(outputFile);
  TH1D* subleadingJetEtaHist = new TH1D(Form("subleadingJetEtaHist_dataSet%i", dataSet), ";#eta;Counts", 98, -4.9, 4.9);
  subleadingJetEtaHist->SetDirectory(outputFile);
  /**** End histogram initialization ****/


  /**** Create variables, arrays to store data for each event, then set branches ****/
  int eventNumber = 0;
  float fcal_et = 0;
  int njet = 0;
  vector<float>* jet_pt = NULL;
  vector<float> out_jet_pt;
  vector<float>* jet_eta = NULL;
  vector<float> out_jet_eta;
  vector<float>* jet_phi = NULL;
  vector<float> out_jet_phi;
  vector<float>* jet_e = NULL;
  vector<float> out_jet_e;
  int nvert = 0;
  vector<int>* vert_type = NULL;
  vector<int> out_vert_type;
  vector<float>* vert_x = NULL;
  vector<float> out_vert_x;
  vector<float>* vert_y = NULL;
  vector<float> out_vert_y;
  vector<float>* vert_z = NULL;
  vector<float> out_vert_z;
  int ntrk = 0;
  vector<bool>* trk_quality_4 = NULL;
  vector<bool> out_trk_quality_4;
  vector<float>* trk_d0 = NULL;
  vector<float> out_trk_d0;
  vector<float>* trk_z0 = NULL;
  vector<float> out_trk_z0;
  vector<float>* trk_theta = NULL;
  vector<float> out_trk_theta;
  vector<float>* trk_charge = NULL;
  vector<float> out_trk_charge;
  vector<float>* trk_pt = NULL;
  vector<float> out_trk_pt;
  vector<float>* trk_eta = NULL;
  vector<float> out_trk_eta;
  vector<float>* trk_phi = NULL;
  vector<float> out_trk_phi;
  
  tree->SetBranchAddress("eventNumber", &eventNumber);
  if (!periodA) tree->SetBranchAddress("fcalC_et", &fcal_et);
  else tree->SetBranchAddress("fcalA_et", &fcal_et);

  tree->SetBranchAddress("njet", &njet);
  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_phi", &jet_phi);
  tree->SetBranchAddress("jet_e", &jet_e);
  tree->SetBranchAddress("nvert", &nvert);
  tree->SetBranchAddress("vert_type", &vert_type);
  tree->SetBranchAddress("vert_x", &vert_x);
  tree->SetBranchAddress("vert_y", &vert_y);
  tree->SetBranchAddress("vert_z", &vert_z);
  tree->SetBranchAddress("ntrk", &ntrk);
  tree->SetBranchAddress("trk_quality_4", &trk_quality_4);
  tree->SetBranchAddress("trk_d0", &trk_d0);
  tree->SetBranchAddress("trk_z0", &trk_z0);
  tree->SetBranchAddress("trk_theta", &trk_theta);
  tree->SetBranchAddress("trk_charge", &trk_charge);
  tree->SetBranchAddress("trk_pt", &trk_pt);
  tree->SetBranchAddress("trk_eta", &trk_eta);
  tree->SetBranchAddress("trk_phi", &trk_phi);

  /*outTree->Branch("eventNumber", &eventNumber);
  outTree->Branch("njet", &njet);
  outTree->Branch("jet_pt", &out_jet_pt);
  outTree->Branch("jet_eta", &out_jet_eta);
  outTree->Branch("jet_e", &out_jet_e);
  outTree->Branch("nvert", &nvert);
  outTree->Branch("vert_type", &out_vert_type);
  outTree->Branch("vert_x", &out_vert_x);
  outTree->Branch("vert_y", &out_vert_y);
  outTree->Branch("vert_z", &out_vert_z);
  outTree->Branch("ntrk", &ntrk);
  outTree->Branch("trk_quality_4", &out_trk_quality_4);
  outTree->Branch("trk_d0", &out_trk_d0);
  outTree->Branch("trk_z0", &out_trk_z0);
  outTree->Branch("trk_theta", &out_trk_theta);
  outTree->Branch("trk_charge", &out_trk_charge);
  outTree->Branch("trk_pt", &out_trk_pt);
  outTree->Branch("trk_eta", &out_trk_eta);
  outTree->Branch("trk_phi", &out_trk_phi);*/

  // Set branch addresses for triggers
  if (!isMC) {
    for (Trigger* trig : (*triggerSubvector)) {
      tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
      if (containsPeriodLetter) tree->SetBranchAddress(Form("%s_prescale_A", trig->name.c_str()), &(trig->m_trig_prescale));
      else tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->m_trig_prescale));
    }
  }
  /**** End create variables and set branches ****/

  // Iterate over each event
  const int numentries = tree->GetEntries();

  double leadingjpt, leadingjeta, subleadingjpt, subleadingjeta, leadingjphi, subleadingjphi, leadingje, subleadingje, subsubleadingjpt, deltaphi; // jet parameters
  double xp, xa, eff, lumi, scale, q_avg, mjj, etajj;
  int leadingj, subleadingj, subsubleadingj, etabin, actetabin, pbin, index, qbin;
  bool takeEvent;
  Trigger* bestTrigger;
  TLorentzVector leadingj_tlv;
  TLorentzVector subleadingj_tlv;
  TLorentzVector dijet_tlv;
  int numGoodEvents = 0;

  int numberOfEventsThatWillEarnMeFreeCoffee = 0;
  int numberOfEventsThatWontEarnMeFreeCoffee = 0;
  TH1I* eventSelectionHist = new TH1I(Form("eventSelectionHist_dataset%i", dataSet), ";Event selection combination;\"Dijet\" events", 6, -0.5, 5.5);
  eventSelectionHist->SetDirectory(outputFile);
  for (long long entry = 0; entry < numentries; entry++) {
    tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses
    // Basic event selection: require a primary vertex and there to be at least 2 reconstructed jets
    if ((nvert == 0) || (nvert > 0 && vert_type->at(0) != 1) || njet < 2) continue;

    leadingj = -1;
    subleadingj = -1;
    subsubleadingj = -1;
    for (int j = 0; j < njet; j++) {
      if (leadingj == -1 || jet_pt->at(leadingj) < jet_pt->at(j)) {
        subsubleadingj = subleadingj;
        subleadingj = leadingj;
        leadingj = j;
      } else if (subleadingj == -1 || jet_pt->at(subleadingj) < jet_pt->at(j)) {
        subsubleadingj = subleadingj;
        subleadingj = j;
      } else if (subsubleadingj == -1 || jet_pt->at(subsubleadingj) < jet_pt->at(j)) {
        subsubleadingj = j;
      }
    }

    /** Stores parameters of the leading dijet pair **/
    leadingjpt = (double)jet_pt->at(leadingj);
    subleadingjpt = (double)jet_pt->at(subleadingj);
    if (njet >= 3) subsubleadingjpt = (double)jet_pt->at(subsubleadingj);
    else subsubleadingjpt = 0;

    leadingjphi = (double)jet_phi->at(leadingj);
    subleadingjphi = (double)jet_phi->at(subleadingj);
    leadingjeta = (double)jet_eta->at(leadingj);
    subleadingjeta = (double)jet_eta->at(subleadingj);

    leadingje = (double)jet_e->at(leadingj);
    subleadingje = (double)jet_e->at(subleadingj);

    // make sure phi variables are in the range 0 to 2pi
    while (leadingjphi < 0) leadingjphi += 2*pi;
    while (subleadingjphi < 0) subleadingjphi += 2*pi;

    deltaphi = TMath::Abs(leadingjphi - subleadingjphi);
    if (deltaphi > pi) deltaphi = 2*pi - deltaphi;
    /** End find leading dijets **/

    /** Event selection **/
    eventSelectionHist->Fill(0);
    if (lowerPhiCut < leadingjphi && leadingjphi < upperPhiCut && lowerEtaCut < leadingjeta && leadingjeta < upperEtaCut) continue;
    if (lowerPhiCut < subleadingjphi && subleadingjphi < upperPhiCut && lowerEtaCut < subleadingjeta && subleadingjeta < upperEtaCut) continue; // select outside disabled HEC region
    eventSelectionHist->Fill(1);
    if (deltaphi < 7.*pi/8.) continue; //require deltaPhi gap of 7pi/8
    eventSelectionHist->Fill(2);
    if (leadingjpt < dijetMinimumPt) continue; // minimum pt cut on leading jet
    eventSelectionHist->Fill(3);
    if (subleadingjpt < dijetMinimumPt) continue; // minimum pt cut on subleading jet
    eventSelectionHist->Fill(4);
    if (subsubleadingjpt/leadingjpt > dijetPtRatioCut) continue; // maximum pt cut on subsubleading jet as a function of the leading jet pt
    eventSelectionHist->Fill(5);
    numGoodEvents++;
    /** End event selection **/


    if (leadingjpt > 1200) cout << Form("High pt (%.0f GeV) jet detected in run %i, event %i!", leadingjpt, dataSet, eventNumber) << endl;


    /** Find scaling information to get a cross section measurement **/
    etabin = getEtabin(leadingjeta);
    pbin = getPbin(leadingjpt);
    if (pbin == -1 || pbin == numpbins || etabin == -1 || etabin == numetabins) continue; // make sure we are in a valid kinematic bin

    if (!isMC) {
      bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins]; // kinematicTriggerVec is created per run, so we do not need to flip etas
      if (bestTrigger == NULL) continue; // make sure its not a null trigger
      takeEvent = bestTrigger->m_trig_bool && bestTrigger->m_trig_prescale > 0 && bestTrigger->min_pt <= leadingjpt;

      if (periodA) actetabin = numetabins - etabin - 1;
      else actetabin = etabin;

      index = bestTrigger->index;
      lumi = kinematicLumiVec[pbin + actetabin*numpbins]; // kinematicLumiVec is created for all runs, so we DO need to flip etas!!!!
      eff = (*triggerEfficiencyFunctions)[index]->Eval(leadingjpt); // same with trigger efficiency, but there is no eta dependence
      takeEvent = takeEvent && lumi != 0 && eff != 0.;
      if (!takeEvent) continue; // make sure the trigger fired, and that we're not going to be dividing by 0
    } else {
      lumi = 1; // assume 1 ub of luminosity for simplicity right now if MC
      eff = 1; // no triggers for MC, so efficiency is just 1
    }

    scale = 1e3/(eff*lumi); // set the scale to convert counts -> "efficiency corrected cross-section"-esque measurement
    //if (etaPhiScaleFactorsHist != NULL) scale *= etaPhiScaleFactorsHist->GetBinContent(etaPhiScaleFactorsHist->FindBin(leadingjeta, leadingjphi));

    if (periodA) {
      leadingjeta *= -1;
      subleadingjeta *= -1;
    }
    /** End find scaling info **/

    xp = get_xp(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, false);
    xa = get_xa(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, false);

    if (xp >= 0.1) {
      out_jet_pt.clear();         for (auto element : *jet_pt) out_jet_pt.push_back(element);
      out_jet_eta.clear();        for (auto element : *jet_eta) out_jet_eta.push_back(element);
      out_jet_phi.clear();        for (auto element : *jet_phi) out_jet_phi.push_back(element);
      out_jet_e.clear();          for (auto element : *jet_e) out_jet_e.push_back(element);
      out_vert_type.clear();      for (auto element : *vert_type) out_vert_type.push_back(element);
      out_vert_x.clear();         for (auto element : *vert_x) out_vert_x.push_back(element);
      out_vert_y.clear();         for (auto element : *vert_y) out_vert_y.push_back(element);
      out_vert_z.clear();         for (auto element : *vert_z) out_vert_z.push_back(element);
      out_trk_quality_4.clear();  for (auto element : *trk_quality_4) out_trk_quality_4.push_back(element);
      out_trk_theta.clear();      for (auto element : *trk_theta) out_trk_theta.push_back(element);
      out_trk_d0.clear();         for (auto element : *trk_d0) out_trk_d0.push_back(element);
      out_trk_z0.clear();         for (auto element : *trk_z0) out_trk_z0.push_back(element);
      out_trk_charge.clear();     for (auto element : *trk_charge) out_trk_charge.push_back(element);
      out_trk_pt.clear();         for (auto element : *trk_pt) out_trk_pt.push_back(element);
      out_trk_eta.clear();        for (auto element : *trk_eta) out_trk_eta.push_back(element);
      out_trk_phi.clear();        for (auto element : *trk_phi) out_trk_phi.push_back(element);

      outTree->Fill();
      numberOfEventsThatWillEarnMeFreeCoffee++;
      leadingJetEtaHist->Fill(leadingjeta);
      subleadingJetEtaHist->Fill(subleadingjeta);
    }

    if (xp >= 0.01) numberOfEventsThatWontEarnMeFreeCoffee++;

    // Fill hardness (Q^2) plots
    q_avg = TMath::Sqrt(0.5*(get_q2(xp, leadingje, leadingjpt) + get_q2(xp, subleadingje, subleadingjpt)));
    {
      qbin = 0;
      while (qbin < numqbins+1 && qbins[qbin] < q_avg) qbin++;
      qbin -= 1;
      if (!(qbin == -1 || qbin >= numqbins)) {
        xqHistArr[qbin]->Fill(xp, scale);
        xqHistArr[qbin+numqbins]->Fill(xa, scale);
      }
    }

    // Fill q2 xa correlation plot
    qxcorr->Fill(xa, q_avg*q_avg, scale);

    // Fill xa xp correlation plot
    xaxpcorr->Fill(xa, xp, scale);

    // Fill Pb-going FCAL distribution plot
    fcalhist->Fill(xp, fcal_et, scale);

    leadingj_tlv.SetPtEtaPhiE(leadingjpt, leadingjeta, leadingjphi, leadingje);
    subleadingj_tlv.SetPtEtaPhiE(subleadingjpt, subleadingjeta, subleadingjphi, subleadingje);
    dijet_tlv = leadingj_tlv + subleadingj_tlv;

    mjj = dijet_tlv.Mag();
    etajj = dijet_tlv.Eta();

    etabin = getEtabin(etajj);
    if (etabin == -1 || etabin == numetabins) continue;
    
    mHistArr[etabin]->Fill(mjj, scale);
    // Fill xa xp distribution plots
    xHistArr[etabin]->Fill(xp, scale); // by filling with the period B etabin, we avoid having to flip the addition later on
    xHistArr[etabin+numetabins]->Fill(xa, scale);

  }

  // Save to root file
  for (int etabin = 0; etabin < 2*numetabins; etabin++) {
    scale = 1. / (etabins[(etabin%numetabins)+1] - etabins[etabin%numetabins]);
    xHistArr[etabin]->Scale(scale, "width");
    xHistArr[etabin]->Write();
  }
  for (int qbin = 0; qbin < 2*numqbins; qbin++) {
    scale = 1. / (qbins[(qbin%numqbins)+1] - qbins[qbin%numqbins]);
    xqHistArr[qbin]->Scale(scale, "width");
    xqHistArr[qbin]->Write();
  }

  for (int etabin = 0; etabin < numetabins; etabin++) {
    scale = 1. / (etabins[etabin+1] - etabins[etabin]);
    mHistArr[etabin]->Scale(scale, "width");
    mHistArr[etabin]->Write();
  }
  qxcorr->Scale(1., "width");
  xaxpcorr->Scale(1., "width");
  fcalhist->Scale(1., "width");

  qxcorr->Write();
  xaxpcorr->Write();
  fcalhist->Write();

  eventSelectionHist->Write();

  leadingJetEtaHist->Write();
  subleadingJetEtaHist->Write();

  TVectorD infoVec(4);
  infoVec[0] = luminosity;
  infoVec[1] = numGoodEvents;
  infoVec[2] = numberOfEventsThatWillEarnMeFreeCoffee;
  infoVec[3] = numberOfEventsThatWontEarnMeFreeCoffee;
  infoVec.Write("infoVec");

  outTree->Write();
  outputFile->Close();
  if (outputFile) delete outputFile;
  //if (outTree) delete outTree;
}
