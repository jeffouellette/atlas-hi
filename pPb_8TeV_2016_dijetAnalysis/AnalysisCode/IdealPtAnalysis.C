#include "../Util.C"

void IdealPtAnalysis(const int dataSet, // Run number identifier.
                     const double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
  initialize(dataSet, true);

  const bool isMC = useDataVersion == 0;
  if (!isMC && skipRun(dataSet)) return;
  if (isMC && skipMC(dataSet)) return;
 
  vector<TF1*>* triggerEfficiencyFunctions = NULL;
  vector<Trigger*>* triggerSubvector = NULL;
  if (!isMC) {
    triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();

    /**** Generate list of physics triggers ****/
    triggerSubvector = getTriggerSubvector(dataSet);
    if (debugStatements) {
      cout << "Status: In IdealPtAnalysis.C (breakpoint A): Processing run " << dataSet << " with triggers:" << endl;
      for (Trigger* trig : (*triggerSubvector)) {
        cout << "\t" << trig->name << endl;
      }
    }
    /**** End generate list of physics triggers ****/
  }


  /**** Find the relevant TTree for this run ****/
  TTree* tree = NULL;
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
          if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint B): Found " << fname.Data() << endl; 
          if (fname.Contains(to_string(dataSet))) {
            tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
            break;
          }
        }
      }
    }
  }
  if (tree == NULL) {
    cout << "Error: In IdealPtAnalysis.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
    return;
  }
  /**** End find TTree ****/


  /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
  {
    vector<string> interestingBranchNames = {"njet", "jet_pt", "jet_eta", "jet_phi", "jet_e", "nvert", "vert_type"};
    TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
    bool interestingBranch;
    for (TObject* obj : *branches) {
      TString branchName = (TString)obj->GetName();
      if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint D): Tree contains branch \"" << branchName.Data() << "\"" << endl;
      interestingBranch = false;
      for (string s : interestingBranchNames) {
        interestingBranch = interestingBranch || (branchName.Data() == s);
      }
      if (!isMC && !interestingBranch) {
        for (Trigger* trig : (*triggerSubvector)) {
          if (branchName == trig->name) {
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


  /**** Set branching addresses ****/
  vector<float>* jet_pt = NULL;
  vector<float>* jet_eta = NULL;
  vector<float>* jet_phi = NULL;
  vector<float>* jet_e = NULL;
  /*float jet_pt[60] = {};
  float jet_eta[60] = {};
  float jet_phi[60] = {};
  float jet_e[60] = {};*/
  int njet = 0;
  vector<int>* vert_type = NULL;
  //int vert_type[60] = {};
  int nvert = 0;
  

  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_phi", &jet_phi);
  tree->SetBranchAddress("jet_e", &jet_e);
  tree->SetBranchAddress("njet", &njet);
  tree->SetBranchAddress("nvert", &nvert);
  tree->SetBranchAddress("vert_type", &vert_type);
  if (!isMC) {
    for (Trigger* trig : (*triggerSubvector)) {
      tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
    }
  }
  /**** End set branch addresses ****/


  /**** Histogram initialization ****/
  TH1D* histArr[numetabins];
  TH2D* etaPhiHist;
  TH2D* subleadingEtaPhiHist;
  TH2D* yPhiHist;
  int pbin, etabin;
  for (etabin = 0; etabin < numetabins; etabin++) {
    TString histName = Form("trig_pt_counts_dataset%i_etabin%i", dataSet, etabin);
    histArr[etabin] = new TH1D(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}d#eta #left[nb GeV^{-1}#right]", numpbins, pbins);
    histArr[etabin]->Sumw2(); // instruct each histogram to propagate errors
  }
  {
    TString histName = Form("etaPhiHist_dataset%i", dataSet);
    etaPhiHist = new TH2D(histName, ";#eta;#phi;", 98, -4.9, 4.9, 100, 0, 2*pi);
    histName = Form("subleadingEtaPhiHist_dataset%i", dataSet);
    subleadingEtaPhiHist = new TH2D(histName, ";#eta;#phi;", 98, -4.9, 4.9, 100, 0, 2*pi);
    histName = Form("yPhiHist_dataset%i", dataSet);
    yPhiHist = new TH2D(histName, ";#it{y};#phi;", 138, -3.5, 3.5, 100, 0, 2*pi); // bins are spaced by 0.05 in y to match CoM boost exactly
  }

  /**** Iterate over each event ****/
  const int numentries = tree->GetEntries();
  int leadingj;
  double jpt, jeta, jphi, je, eff, scale;
  TLorentzVector tlv;
  Trigger* bestTrigger = NULL;
  for (long long entry = 0; entry < numentries; entry++) {
    tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

    if ((nvert == 0) || (nvert > 0 && vert_type->at(0) != 1)) continue; // basic event selection: require there to be a primary vertex

    for (int j = 0; j < njet; j++) {
      jpt = (double)jet_pt->at(j);
      jeta = (double)jet_eta->at(j);
      jphi = (double)jet_phi->at(j);
      je = (double)jet_e->at(j);
      while (jphi < 0) jphi += 2.*pi; // converts phi to 0, 2*pi range

      etabin = getEtabin(jeta);
      pbin = getPbin(jpt);
      if (pbin == -1 || etabin == -1 || pbin == numpbins || etabin == numetabins) continue; // this checks that the jets fall within the pt, eta bins

      if (!isMC) {
        bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
        if (bestTrigger == NULL) continue; // make sure we're not trying to look at a null trigger

        eff = (*triggerEfficiencyFunctions)[bestTrigger->index]->Eval(jpt);
        if (eff == 0) continue; // avoid dividing by 0 
        scale = 1./eff;
      }
      else scale = 1.;

      if (isMC || (!isMC && bestTrigger->m_trig_bool)) {
        if (!(lowerPhiCut < jphi && jphi < upperPhiCut && lowerEtaCut < jeta && jeta < upperEtaCut)) histArr[etabin]->Fill(jpt, scale);
        if (!highPtJetsOnly || jpt > 70) etaPhiHist->Fill(jeta, jphi);
        tlv.SetPtEtaPhiE(jpt, jeta, jphi, je);
        yPhiHist->Fill(tlv.Rapidity(), jphi);
      }
    }
    
    leadingj = 0;
    if (njet < 2) continue; // specialize to events with 2+ jets
    for (int j = 1; j < njet; j++) {
      if (jet_pt->at(leadingj) < jet_pt->at(j)) leadingj = j;
    }
    // fill eta-phi correlation with only non-leading jets
    for (int j = 0; j < njet; j++) {
      if (j == leadingj) continue;
      jpt = (double)jet_pt->at(j);
      jeta = (double)jet_eta->at(j);
      jphi = (double)jet_phi->at(j);
      while (jphi < 0) jphi += 2.*pi;

      etabin = getEtabin(jeta);
      pbin = getPbin(jpt);
      if (pbin == -1 || etabin == -1 || pbin == numpbins || etabin == numetabins) continue; // this checks that the jets fall within the pt, eta bins

      if (!isMC) {
        bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
        if (bestTrigger == NULL || !bestTrigger->m_trig_bool) continue; // make sure we're not trying to look at a null trigger, and if not, that it fired
      }

      subleadingEtaPhiHist->Fill(jeta, jphi);
    }
  }
  /**** End event iteration ****/


  /**** Write output histograms to a root file ****/
  TFile* output = new TFile(Form("%sdataset_%i.root", ptPath.c_str(), dataSet), "RECREATE");
  for (etabin = 0; etabin < numetabins; etabin++) {
    histArr[etabin]->Scale(1/A); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
    histArr[etabin]->Write();
  }
  etaPhiHist->Write();
  subleadingEtaPhiHist->Write();
  yPhiHist->Write();

  TVectorD infoVec(5);
  infoVec[0] = dataSet;
  infoVec[1] = numetabins;
  infoVec[2] = numtrigs;
  infoVec[3] = numpbins;
  infoVec[4] = luminosity;
  infoVec.Write(Form("infoVec_%i", dataSet));

  output->Close();
  /**** End write output ****/

  if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint E): Finished calculating pt spectrum for run " << dataSet << endl;
  return;
}
