#include "DijetAnalysis.h"
#include "DijetUtils.h"
#include "Params.h"

#include <Utils.h>
#include <TreeVariables.h>
#include <GlobalParams.h>

#include <TH2D.h>
#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TVectorT.h>

using namespace atlashi;

namespace JetAnalysis {

void DijetAnalysis(const int dataSet, // Data set identifier. If not MC, this should be a run number. If MC, this should be whatever number follows "tid" in the MC file.
                   const double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                   const bool isMC, // MC flag
                   const bool isPeriodA) // period A flag
{
  Initialize (dataSet, isMC, true);

  const TString fileIdentifier = GetIdentifier (dataSet, "", isMC, false, isPeriodA);
  cout << "File identifier: " << fileIdentifier << endl;

  // Check whether to skip this particular analysis
  if (!isMC && SkipRun(isPeriodA)) return;
  if (isMC && SkipMC(isPeriodA)) return;

  //vector<TF1*>* triggerEfficiencyFunctions = NULL;
  //if (!isMC) triggerEfficiencyFunctions = GetTriggerEfficiencyFunctions();

  const bool useIonTriggers = !isMC && dataSet < 313629; // false for MC

  /**** Generate list of physics triggers ****/
  vector<Trigger*>* triggerSubvector = NULL;
  if (!isMC) {
   triggerSubvector = GetTriggerSubvector(dataSet);
   if (debugStatements) {
    cout << "Status: In DijetAnalysis.C (breakpoint A): Processing run " << dataSet << " with triggers:" << endl;
    for (Trigger* trig : (*triggerSubvector)) {
     cout << "\t" << trig->name << endl;
    }
   }
  }
  /**** End generate list of physics triggers ****/


  /**** Find the relevant TTree for this run ****/
  //TFile* inFile = GetTFile (fileIdentifier);
  TFile* inFile = GetFile ("", dataSet, isMC, "");
  if (!inFile || !inFile->IsOpen ()) {
   cout << "Error: In DijetAnalysis.C: TFile for given data set not identified. Quitting." << endl;
   return;
  }
  TTree* tree = (TTree*)inFile->Get("tree");
  if (tree == NULL) {
   cout << "Error: In DijetAnalysis.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/


  /**** Set branching addresses ****/
  TreeVariables* t = new TreeVariables (tree, isMC);
  //std::set<TString> activeBranches = {"jet_n", "jet_pt", "jet_eta", "jet_phi", "jet_e", "truth_jet_n", "truth_jet_pt", "truth_jet_eta", "truth_jet_phi", "truth_jet_e", "nvert", "vert_type", "fcalA_et", "fcalC_et", "ntrk", "trk_quality4", "trk_d0", "trk_z0", "trk_theta", "trk_charge", "trk_pt", "trk_eta", "trk_phi", "vert_x", "vert_y", "vert_z"};
  //t->SetGetVertices();
  t->SetGetFCals();
  //t->SetGetTracks();
  t->SetGetHIJets(false);
  t->SetGetSimpleJets(false);
  t->SetBranchAddresses ();
  tree->SetBranchAddress ("njet", &t->jet_n);
  tree->SetBranchAddress ("jet_pt", &t->jet_pt);
  tree->SetBranchAddress ("jet_eta", &t->jet_eta);
  tree->SetBranchAddress ("jet_phi", &t->jet_phi);
  tree->SetBranchAddress ("jet_e", &t->jet_e);

  //t->SetBranchAddresses(activeBranches);
  if (!isMC) {
   for (Trigger* trig : (*triggerSubvector)) {
    if (trig->name == minbiasTriggerName) {
     trig->trigBool = false;
     trig->trigPrescale = 1;
     continue;
    }

    tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->trigBool));
    if (useIonTriggers)
     tree->SetBranchAddress(Form("%s_prescale_A", trig->name.c_str()), &(trig->trigPrescale));
    else
     tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->trigPrescale));
   }
  }
  /**** End set branch addresses ****/


  /**** Calculate eta-phi rescaling information ****/
  TH2D* etaPhiScaleFactorsHist = NULL;
  if (!isMC) {
   /**** Load eta-phi jet correlation plot ****/
   TFile* etaPhiFile = new TFile((rootPath + "etaPhiHist.root").Data(), "READ");
   TH2D* etaPhiHist = (TH2D*)etaPhiFile->Get("etaPhiHist");
   TH2D* subJetEtaPhiHist = (TH2D*)etaPhiFile->Get("subleadingEtaPhiHist");
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
      y_prime = subJetEtaPhiHist->GetYaxis()->GetBinCenter(bin_y_prime+1);
      dy = TMath::Abs(y - y_prime);
      if (dy > pi) dy = 2*pi - dy;
      if (dy < 7.*pi/8.) continue; // this checks whether our y' coordinate is further away from y by at least 7pi/8 (up to pi). Otherwise it is not in our integration region, we skip it.
      for (int bin_x_prime = 0; bin_x_prime < nbins_x; bin_x_prime++) {
       x_prime = subJetEtaPhiHist->GetXaxis()->GetBinCenter(bin_x_prime+1);
       // now check if x' meets the eta cut requirements
       content = subJetEtaPhiHist->GetBinContent(bin_x_prime+1, bin_y_prime+1);
       // if the ' coordinate is outside the HEC region, then add the counts there to your integral in the denominator
       if (!InDisabledHEC(x_prime, y_prime)) denominator += content;
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
       if (!InDisabledHEC(x_prime, y_prime)) denominator += content;
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


  /**** Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a, as well as other histograms ****/
  TH1D* xHistArr[2*numetabins];
  TH1D* xqHistArr[2*numqbins];
  TH1D* mHistArr[numetabins];

  for (int etabin = 0; etabin < numetabins; etabin++) {
   xHistArr[etabin] = new TH1D(Form("%ieta%i", dataSet, etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#eta #left[nb#right]", etabins[etabin], etabins[etabin+1]), numxbins, xbins);
   xHistArr[etabin]->Sumw2();
  }
  for (int etabin = numetabins; etabin < 2*numetabins; etabin++) {
   xHistArr[etabin] = new TH1D(Form("%ieta%i", dataSet, etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d^{2}N/L_{int}d#it{x}_{a}d#eta #left[nb#right]", etabins[etabin%numetabins], etabins[(etabin%numetabins)+1]), numxbins, xbins);
   xHistArr[etabin]->Sumw2();
  }

  for (int qbin = 0; qbin < numqbins; qbin++) {
   xqHistArr[qbin] = new TH1D(Form("%iq%i", dataSet, qbin), Form("%g < #it{Q} < %g;#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin], qbins[qbin+1]), numxbins, xbins);
   xqHistArr[qbin]->Sumw2();
  }
  for (int qbin = numqbins; qbin < 2*numqbins; qbin++) {
   xqHistArr[qbin] = new TH1D(Form("%iq%i", dataSet, qbin), Form("%g < #it{Q} < %g;#it{x}_{a};d^{2}N/L_{int}d#it{x}_{a}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin%numqbins], qbins[(qbin%numqbins)+1]), numxbins, xbins);
   xqHistArr[qbin]->Sumw2();
  }

  for (int etabin = 0; etabin < numetabins; etabin++) {
   mHistArr[etabin] = new TH1D(Form("mjj_%ieta%i", dataSet, etabin), Form("%1.1f < #eta < %1.1f;#it{m}_{JJ} #left[GeV#right];d^{2}N/L_{int}d#it{m}_{JJ}d#eta #left[nb GeV^{-1}#right]", etabins[etabin], etabins[etabin+1]), nummbins, mbins);
   mHistArr[etabin]->Sumw2();
  }

  TH2D* qxcorr = new TH2D(Form("xqcorr_dataSet%i", dataSet), ";#it{x}_{a};#it{#bar{Q}}^{2} #left[GeV^{2}#right];d^{2}N/L_{int}d{x}_{a}d#it{#bar{Q}}^{2} #left[nb GeV^{-2}#right]", numq2xbins, q2xbins, numq2bins, q2bins);
  TH2D* xaxpcorr = new TH2D(Form("xaxpcorr_dataSet%i", dataSet), ";#it{x}_{a};#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#it{x}_{a}", numxbins, xbins, numxbins, xbins);
  TH2D* fcalhist = new TH2D(Form("fcalhist_dataSet%i", dataSet), ";#it{x}_{p};FCAL energy #left[GeV#right];", numxbins, xbins, numfcalbins, fcalbins);

  TH1D* leadingJetEtaHist = new TH1D(Form("leadingJetEtaHist_dataSet%i", dataSet), ";#eta;Counts", 98, -4.9, 4.9);
  TH1D* subleadingJetEtaHist = new TH1D(Form("subleadingJetEtaHist_dataSet%i", dataSet), ";#eta;Counts", 98, -4.9, 4.9);
  /**** End histogram initialization ****/


  /**** Create variables, arrays to store data for each event, then set branches ****/
  //int out_eventNumber = 0;
  //float out_fcalA_et = 0;
  //float out_fcalC_et = 0;
  //int out_njet = 0;
  //vector<float> out_jet_pt;
  //vector<float> out_jet_eta;
  //vector<float> out_jet_phi;
  //vector<float> out_jet_e;
  //int out_nvert = 0;
  //vector<int> out_vert_type;
  //vector<float> out_vert_x;
  //vector<float> out_vert_y;
  //vector<float> out_vert_z;
  //int out_ntrk = 0;
  //vector<bool> out_trk_quality_4;
  //vector<float> out_trk_d0;
  //vector<float> out_trk_z0;
  //vector<float> out_trk_theta;
  //vector<float> out_trk_charge;
  //vector<float> out_trk_pt;
  //vector<float> out_trk_eta;
  //vector<float> out_trk_phi;
  
  //if (fillTracksTree) {
  // outTree->Branch("eventNumber", &out_eventNumber);
  // outTree->Branch("njet", &out_njet);
  // outTree->Branch("jet_pt", &out_jet_pt);
  // outTree->Branch("jet_eta", &out_jet_eta);
  // outTree->Branch("jet_e", &out_jet_e);
  // outTree->Branch("nvert", &out_nvert);
  // outTree->Branch("vert_type", &out_vert_type);
  // outTree->Branch("vert_x", &out_vert_x);
  // outTree->Branch("vert_y", &out_vert_y);
  // outTree->Branch("vert_z", &out_vert_z);
  // outTree->Branch("ntrk", &out_ntrk);
  // outTree->Branch("trk_quality_4", &out_trk_quality_4);
  // outTree->Branch("trk_d0", &out_trk_d0);
  // outTree->Branch("trk_z0", &out_trk_z0);
  // outTree->Branch("trk_theta", &out_trk_theta);
  // outTree->Branch("trk_charge", &out_trk_charge);
  // outTree->Branch("trk_pt", &out_trk_pt);
  // outTree->Branch("trk_eta", &out_trk_eta);
  // outTree->Branch("trk_phi", &out_trk_phi);
  //}

  /**** End create variables and set branches ****/

  // Iterate over each event
  const int numentries = tree->GetEntries();

  double leadingjpt, leadingjeta, subleadingjpt, subleadingjeta, leadingjphi, subleadingjphi, leadingje, subleadingje, subsubleadingjpt, deltaphi; // jet parameters
  double xp, xa, eff, lumi, scale, q_avg, mjj, etajj;
  int leadingj, subleadingj, subsubleadingj, etabin, actetabin, pbin, index, qbin;
  bool takeEvent;
  Trigger* bestTrigger = NULL;
  TLorentzVector leadingj_tlv, subleadingj_tlv, dijet_tlv;
  int numGoodEvents = 0;

  TH1I* eventSelectionHist = new TH1I(Form("eventSelectionHist_dataSet%i", dataSet), ";Event selection combination;\"Dijet\" events", 6, -0.5, 5.5);
  for (long long entry = 0; entry < numentries; entry++) {
   tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

   //// Basic event selection: require a primary vertex and there to be at least 2 reconstructed jets
   //if ((t->nvert == 0) || (t->nvert > 0 && t->vert_type->at(0) != 1) || t->jet_n < 2) continue;

   leadingj = -1;
   subleadingj = -1;
   subsubleadingj = -1;
   for (int j = 0; j < t->jet_n; j++) {
    if (leadingj == -1 || t->jet_pt->at(leadingj) < t->jet_pt->at(j)) {
     subsubleadingj = subleadingj;
     subleadingj = leadingj;
     leadingj = j;
    } else if (subleadingj == -1 || t->jet_pt->at(subleadingj) < t->jet_pt->at(j)) {
     subsubleadingj = subleadingj;
     subleadingj = j;
    } else if (subsubleadingj == -1 || t->jet_pt->at(subsubleadingj) < t->jet_pt->at(j)) {
     subsubleadingj = j;
    }
   }
   if (leadingj == -1 || subleadingj == -1)
    continue;

   /** Stores parameters of the leading dijet pair **/
   leadingjpt = t->jet_pt->at(leadingj);
   subleadingjpt = t->jet_pt->at(subleadingj);
   if (t->jet_n >= 3) subsubleadingjpt = t->jet_pt->at(subsubleadingj);
   else subsubleadingjpt = 0;

   leadingjphi = InTwoPi (t->jet_phi->at(leadingj));
   subleadingjphi = InTwoPi (t->jet_phi->at(subleadingj));
   leadingjeta = t->jet_eta->at(leadingj);
   subleadingjeta = t->jet_eta->at(subleadingj);

   leadingje = t->jet_e->at(leadingj);
   subleadingje = t->jet_e->at(subleadingj);

   deltaphi = DeltaPhi (leadingjphi, subleadingjphi);
   /** End find leading dijets **/

   /** Event selection **/
   eventSelectionHist->Fill(0);
   if (InDisabledHEC (leadingjeta, leadingjphi))
    continue;
   if (InDisabledHEC (subleadingjeta, subleadingjphi))
    continue; // select outside disabled HEC region

   eventSelectionHist->Fill(1);
   if (deltaphi < 7.*pi/8.)
    continue; //require deltaPhi gap of 7pi/8

   eventSelectionHist->Fill(2);
   if (leadingjpt < dijetMinimumPt)
    continue; // minimum pt cut on leading jet

   eventSelectionHist->Fill(3);
   if (subleadingjpt < dijetMinimumPt)
    continue; // minimum pt cut on subleading jet

   eventSelectionHist->Fill(4);
   if (subsubleadingjpt/leadingjpt > dijetPtRatioCut)
    continue; // maximum pt cut on subsubleading jet as a function of the leading jet pt

   eventSelectionHist->Fill(5);
   numGoodEvents++;
   /** End event selection **/

   if (leadingjpt > 1200)
    cout << Form("High pt (%.0f GeV) jet detected in run %i, event %i!", leadingjpt, dataSet, t->eventNumber) << endl;


   /** Find scaling information to get a cross section measurement **/
   etabin = getEtabin(leadingjeta);
   pbin = getPbin(leadingjpt);
   if (pbin == -1 || pbin == numpbins || etabin == -1 || etabin == numetabins)
    continue; // make sure we are in a valid kinematic bin


   if (!isMC) {
    bestTrigger = kinematicTriggerVec[pbin][etabin]; // kinematicTriggerVec is created per run, so we do not need to flip etas
    if (bestTrigger == NULL)
     continue; // make sure its not a null trigger
    takeEvent = bestTrigger->trigBool && bestTrigger->trigPrescale > 0 && bestTrigger->minPt <= leadingjpt;
    if (!takeEvent)
     continue;

    if (isPeriodA) actetabin = numetabins - etabin - 1;
    else actetabin = etabin;

    index = bestTrigger->index;
    lumi = kinematicLumiVec[pbin][actetabin]; // kinematicLumiVec is created for all runs, so we DO need to flip etas!!!!
    //eff = (*triggerEfficiencyFunctions)[index]->Eval(leadingjpt); // same with trigger efficiency, but there is no eta dependence
    eff = 1.;
    if (lumi ==0  || eff == 0.)
     continue; // make sure the trigger fired, and that we're not going to be dividing by 0
   }
   else {
    lumi = 1; // assume 1 ub of luminosity for simplicity right now if MC
    eff = 1; // no triggers for MC, so efficiency is just 1
   }

   scale = 1e3/(eff*lumi); // set the scale to convert counts -> "efficiency corrected cross-section"-esque measurement
   //if (etaPhiScaleFactorsHist != NULL) scale *= etaPhiScaleFactorsHist->GetBinContent(etaPhiScaleFactorsHist->FindBin(leadingjeta, leadingjphi));

   if (isPeriodA) {
    leadingjeta *= -1;
    subleadingjeta *= -1;
   }
   /** End find scaling info **/

   xp = get_xp(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, false);
   xa = get_xa(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, false);

   //if (fillTracksTree && xp >= 0.1) {
   // out_njet = t->jet_n;
   // out_jet_pt.clear();         for (auto element : *t->jet_pt) out_jet_pt.push_back(element);
   // out_jet_eta.clear();        for (auto element : *t->jet_eta) out_jet_eta.push_back(element);
   // out_jet_phi.clear();        for (auto element : *t->jet_phi) out_jet_phi.push_back(element);
   // out_jet_e.clear();          for (auto element : *t->jet_e) out_jet_e.push_back(element);
   // out_nvert = t->nvert;
   // out_vert_type.clear();      for (auto element : *t->vert_type) out_vert_type.push_back(element);
   // out_vert_x.clear();         for (auto element : *t->vert_x) out_vert_x.push_back(element);
   // out_vert_y.clear();         for (auto element : *t->vert_y) out_vert_y.push_back(element);
   // out_vert_z.clear();         for (auto element : *t->vert_z) out_vert_z.push_back(element);
   // out_ntrk = t->ntrk;
   // out_trk_quality_4.clear();  for (auto element : *t->trk_quality_4) out_trk_quality_4.push_back(element);
   // out_trk_theta.clear();      for (auto element : *t->trk_theta) out_trk_theta.push_back(element);
   // out_trk_d0.clear();         for (auto element : *t->trk_d0) out_trk_d0.push_back(element);
   // out_trk_z0.clear();         for (auto element : *t->trk_z0) out_trk_z0.push_back(element);
   // out_trk_charge.clear();     for (auto element : *t->trk_charge) out_trk_charge.push_back(element);
   // out_trk_pt.clear();         for (auto element : *t->trk_pt) out_trk_pt.push_back(element);
   // out_trk_eta.clear();        for (auto element : *t->trk_eta) out_trk_eta.push_back(element);
   // out_trk_phi.clear();        for (auto element : *t->trk_phi) out_trk_phi.push_back(element);

   // outTree->Fill();
   // leadingJetEtaHist->Fill(leadingjeta);
   // subleadingJetEtaHist->Fill(subleadingjeta);
   //}


   // Fill hardness (Q^2) plots
   //q_avg = TMath::Sqrt(0.5*(get_q2(xp, leadingje, leadingjpt) + get_q2(xp, subleadingje, subleadingjpt)));
   q_avg = get_q (leadingjpt, subleadingjpt);
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
   fcalhist->Fill(xp, (isPeriodA ? t->fcalA_et : t->fcalC_et), scale);

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
  // end event loop.

  //if (inFile) delete inFile;


  /**** Create output TTree ****/
  TString outFileName = xPath;
  if (runPeriodA && !runPeriodB) outFileName = outFileName + "periodA/";
  else if (!runPeriodA && runPeriodB) outFileName = outFileName + "periodB/";
  else outFileName = outFileName + "periodAB/";
  outFileName = Form("%sdataSet_%i.root", outFileName.Data(), dataSet);

  TFile* outputFile = new TFile(outFileName, "RECREATE");
  //TTree* outTree = tree->CloneTree(0);//new TTree("trackTree", "trackTree");
  /**** End create output TTree ****/

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

  TVectorD infoVec(2);
  infoVec[0] = luminosity;
  infoVec[1] = numGoodEvents;
  infoVec.Write("infoVec");

  //outTree->Write();
  outputFile->Close();
  if (outputFile) delete outputFile;
}

} // end namespace
