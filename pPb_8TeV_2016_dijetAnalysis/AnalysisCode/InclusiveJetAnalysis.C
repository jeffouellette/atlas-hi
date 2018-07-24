#include "../Util.C"
#include "../TreeVariables.C"

void InclusiveJetAnalysis(const int dataSet, // Run number identifier.
                          const double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                          const bool isMC,
                          const bool isPeriodA)
{
  Initialize(dataSet, isMC, true);
  const TString fileIdentifier = GetIdentifier (dataSet, isMC, isPeriodA);

  if (!isMC && SkipRun(isPeriodA)) return;
  if (isMC && SkipMC(isPeriodA)) return;
 
  std::vector<TF1*>* triggerEfficiencyFunctions = NULL;
  std::vector<Trigger*>* triggerSubvector = NULL;
  if (!isMC) {
    triggerEfficiencyFunctions = GetTriggerEfficiencyFunctions();

    /**** Generate list of physics triggers ****/
    triggerSubvector = GetTriggerSubvector(dataSet);
    if (debugStatements) {
      cout << "Status: In IdealPtAnalysis.C (breakpoint A): Processing run " << dataSet << " with triggers:" << endl;
      for (Trigger* trig : *triggerSubvector) {
        cout << "\t" << trig->name << endl;
      }
    }
    /**** End generate list of physics triggers ****/
  }


  /**** Find the relevant TTree for this run ****/
  TFile* inFile = GetTFile(fileIdentifier);
  TTree* tree = (TTree*)inFile->Get("bush");

  /**** Set branching addresses ****/
  TreeVariables* t = new TreeVariables (tree);
  std::set<TString> activeBranches = {"njet", "jet_pt", "jet_eta", "jet_phi", "jet_e", "truth_njet", "truth_jet_pt", "truth_jet_eta", "truth_jet_phi", "truth_jet_e", "nvert", "vert_type", "fcalA_et", "fcalC_et"};
  t->SetBranchAddresses(isMC, activeBranches);
  if (!isMC) {
    for (Trigger* trig : (*triggerSubvector)) {
      tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->trigBool));
      tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->trigPrescale));
    }
  }
  /**** End set branch addresses ****/


  /**** Histogram initialization ****/
  TString outFileName = Form("dataSet_%i_%s.root", dataSet, (isPeriodA?"Pbp":"pPb"));
  TFile* outFile = new TFile(Form("%s%s", ptPath.c_str(), outFileName.Data()), "RECREATE");
  outFile->cd();
  TH1F* jetPtHistArr[numetabins];
  TH1F* fcalSumEt[2];
  for (short etabin = 0; etabin < numetabins; etabin++) {
   const char* histName = Form("jetPtSpectrum_dSet%i_%s_etabin%i", dataSet, (isPeriodA?"Pbp":"pPb"), etabin);
   const char* labels = Form(";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}d#eta #left[nb GeV^{-1}#right]");
   jetPtHistArr[etabin] = new TH1F(histName, labels, numpbins, pbins);
   jetPtHistArr[etabin]->Sumw2(); // instruct each histogram to propagate errors
  }
  for (short dir = 0; dir < 2; dir++) {
   fcalSumEt[dir] = new TH1F(Form("fcal%sSumEt_dSet%i_%s", (dir==0?"A":"C"), dataSet, (isPeriodA?"Pbp":"pPb")), ";#Sigma #it{E}_{T} #left[GeV#right];Counts", 100, -50, 250);
   fcalSumEt[dir]->Sumw2();
  }
  TH2F* jetEtaPhiHist = new TH2F(Form("etaPhiHist_dSet%i_%s", dataSet, (isPeriodA?"Pbp":"pPb")), ";#eta;#phi;", 98, -4.9, 4.9, 80, 0, 2*pi);
  TH2F* subJetEtaPhiHist = new TH2F(Form("subJetEtaPhiHist_dSet%i_%s", dataSet, (isPeriodA?"Pbp":"pPb")), ";#eta;#phi;", 98, -4.9, 4.9, 80, 0, 2*pi);
  TH2F* jetYPhiHist = new TH2F(Form("jetYPhiHist_dSet%i_%s", dataSet, (isPeriodA?"Pbp":"pPb")), ";#it{y};#phi;", 138, -3.5, 3.5, 80, 0, 2*pi); // bins are spaced by 0.05 in y to match CoM boost exactly
  TH1F* jetEnergyResponse = new TH1F(Form("jetEnergyResponse_dSet%i_%s", dataSet, (isPeriodA?"Pbp":"pPb")), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 50, 0, 2.0);
  jetEnergyResponse->Sumw2();

  /**** Iterate over each event ****/
  const int numentries = tree->GetEntries();
  int leadingj;
  double jpt, jeta, jphi, je, eff, scale;
  TLorentzVector tlv;
  Trigger* bestTrigger = NULL;
  for (long long entry = 0; entry < numentries; entry++) {
   tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

   if ((t->nvert <= 0) || (t->nvert >= 1 && t->vert_type->at(0) != 1)) continue; // basic event selection: require there to be a primary vertex

   fcalSumEt[0]->Fill(t->fcalA_et);
   fcalSumEt[1]->Fill(t->fcalC_et); 

   for (int j = 0; j < t->njet; j++) {
    jpt = (double)t->jet_pt->at(j);
    jeta = (double)t->jet_eta->at(j);
    jphi = (double)InTwoPi(t->jet_phi->at(j));
    je = (double)t->jet_e->at(j);

    const short etabin = getEtabin(jeta);
    const short pbin = getPbin(jpt);
    if (pbin == -1 || etabin == -1 || pbin == numpbins || etabin == numetabins) continue; // this checks that the jets fall within the pt, eta bins

    if (!isMC) {
     bestTrigger = kinematicTriggerVec[pbin][etabin];
     if (bestTrigger == NULL) continue; // make sure we're not trying to look at a null trigger

     if (!bestTrigger->trigBool) continue; // make sure the trigger fired

     eff = (*triggerEfficiencyFunctions)[bestTrigger->index]->Eval(jpt);
     if (eff == 0) continue; // avoid dividing by 0 
     scale = 1./eff;
    }
    //else scale = t->eventWeight;
    else scale = 1;

    //if (!InDisabledHEC (jeta, jphi)) // reject jets in the HEC
    jetPtHistArr[etabin]->Fill(jpt, scale);
    if (!highPtJetsOnly || jpt > 70)
     jetEtaPhiHist->Fill(jeta, jphi);

    tlv.SetPtEtaPhiE(jpt, jeta, jphi, je);
    jetYPhiHist->Fill(tlv.Rapidity(), jphi);

    double minDeltaR = 100.;
    int truthMatch = -1;
    for (int tj = 0; tj < t->truth_njet; tj++) {
     const dR = deltaR (jeta, t->truth_jet_eta->at(tj), jphi, t->truth_jet_phi->at(tj));
     if (dR < minDeltaR) {
      truth_jet = tj;
      minDeltaR = dR;
     }
    }
    if (0 <= truth_jet && truth_jet < t->truth-njet && minDeltaR < 0.3) {
     jetEnergyResponse->Fill(jpt / t->truth_jet_pt->at(truth_jet), scale);
    }
    
   }
   
   leadingj = 0;
   if (t->njet < 2) continue; // specialize to events with 2+ jets
   for (int j = 1; j < t->njet; j++) {
    if (t->jet_pt->at(leadingj) < t->jet_pt->at(j)) leadingj = j;
   }
   // fill eta-phi correlation with only non-leading jets
   for (int j = 0; j < t->njet; j++) {
    if (j == leadingj) continue;
    jpt = (double)t->jet_pt->at(j);
    jeta = (double)t->jet_eta->at(j);
    jphi = (double)InTwoPi(t->jet_phi->at(j));

    const short etabin = getEtabin(jeta);
    const short pbin = getPbin(jpt);
    if (pbin == -1 || etabin == -1 || pbin == numpbins || etabin == numetabins) continue; // this checks that the jets fall within the pt, eta bins

    if (!isMC) {
     bestTrigger = kinematicTriggerVec[pbin][etabin];
     if (bestTrigger == NULL || !bestTrigger->trigBool) continue; // make sure we're not trying to look at a null trigger, and if not, that it fired
    }

    subJetEtaPhiHist->Fill(jeta, jphi);
   }
  }
  /**** End event iteration ****/
  inFile->Close();
  if (inFile) delete inFile;

  /**** Write output histograms to a root file ****/
  for (short etabin = 0; etabin < numetabins; etabin++) {
   jetPtHistArr[etabin]->Scale(1./A); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
   jetPtHistArr[etabin]->Write();
   if (jetPtHistArr[etabin]) delete jetPtHistArr[etabin];
  }
  for (short dir = 0; dir < 2; dir++) {
   fcalSumEt[dir]->Write();
   if (fcalSumEt[dir]) delete fcalSumEt[dir];
  }
  jetEtaPhiHist->Write();
  if (jetEtaPhiHist) delete jetEtaPhiHist;
  subJetEtaPhiHist->Write();
  if (subJetEtaPhiHist) delete subJetEtaPhiHist;
  jetYPhiHist->Write();
  if (jetYPhiHist) delete jetYPhiHist;
  jetEnergyResponse->Write();
  if (jetEnergyResponse) delete jetEnergyResponse;

  TVectorD infoVec(5);
  infoVec[0] = dataSet;
  infoVec[1] = numetabins;
  infoVec[2] = numtrigs;
  infoVec[3] = numpbins;
  infoVec[4] = luminosity;
  infoVec.Write(Form("infoVec_%i_per%s", dataSet, (isPeriodA?"A":"B")));

  outFile->Close();
  if (outFile) delete outFile;
  /**** End write output ****/

  if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint E): Finished calculating pt spectrum for run " << dataSet << endl;
  return;
}
