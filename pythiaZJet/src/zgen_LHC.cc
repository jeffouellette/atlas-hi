#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <sstream>

#include "Pythia8/Pythia.h"
#include "Pythia8/FJcore.h"

#include <Utilities.h>

using namespace Pythia8;
using namespace atlashi;

int main (int argc, char *argv[]) {

  if (argc != 4) {

    std::cout << " usage: ./zgen MINPT NEVT FILENAMEOUT" << std::endl;

    return 0;
  }

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  pythia.readString ("Beams:eCM = 5020.");

  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 11 13");

  pythia.readString ("WeakZ0:gmZmode = 2"); // set to Z's
  pythia.readString ("WeakSingleBoson:ffbar2gmZ = on");       // code 221
  pythia.readString ("WeakDoubleBoson:ffbar2gmZgmZ = on");    // code 231
  pythia.readString ("WeakDoubleBoson:ffbar2ZW = on");        // code 232
  pythia.readString ("WeakBosonAndParton:qqbar2gmZg = on");   // code 241
  pythia.readString ("WeakBosonAndParton:qg2gmZq = on");      // code 242
  pythia.readString ("WeakBosonAndParton:ffbar2gmZgm = on");  // code 243
  pythia.readString ("WeakBosonAndParton:fgm2gmZf = on");     // code 244

  ostringstream ss; ss << "PhaseSpace:pTHatMin = " << argv[1] << ".";
  //pythia.readString("PhaseSpace:pTHatMin = 10.");
  pythia.readString (ss.str ().c_str ());
  pythia.readString ("PhaseSpace:mHatMin = 60");

  pythia.init ();

  SlowJet *antikT4 = new SlowJet (-1, 1.0, 1, 5, 2, 1);

  int NEVT = atoi (argv[2]);

  TFile *f = new TFile (argv[3] , "RECREATE");

  int b_code;
  int b_id1;
  int b_id2;
  float b_x1pdf;
  float b_x2pdf;
  float b_Q;
  bool b_isValence1;
  bool b_isValence2;

  int b_z_n, b_jet_n, b_l_n;
  vector<float> b_z_pt, b_z_eta, b_z_phi, b_z_m;
  vector<float> b_jet_pt, b_jet_eta, b_jet_phi, b_jet_e;
  vector<float> b_l_pt, b_l_eta, b_l_phi, b_l_m;

  TTree *t = new TTree("tree","a shambling vine tree");

  t->Branch ("code",  &b_code);
  t->Branch ("id1",   &b_id1);
  t->Branch ("id2",   &b_id2);
  t->Branch ("x1pdf", &b_x1pdf);
  t->Branch ("x2pdf", &b_x2pdf);
  t->Branch ("Q",     &b_Q);

  t->Branch ("z_n",   &b_z_n);
  t->Branch ("z_pt",  &b_z_pt);
  t->Branch ("z_eta", &b_z_eta);
  t->Branch ("z_phi", &b_z_phi);
  t->Branch ("z_m",   &b_z_m);

  t->Branch ("l_n",   &b_l_n);
  t->Branch ("l_pt",  &b_l_pt);
  t->Branch ("l_eta", &b_l_eta);
  t->Branch ("l_phi", &b_l_phi);
  t->Branch ("l_m",   &b_l_m);

  t->Branch ("jet_n",   &b_jet_n);
  t->Branch ("jet_pt",  &b_jet_pt);
  t->Branch ("jet_eta", &b_jet_eta);
  t->Branch ("jet_phi", &b_jet_phi);
  t->Branch ("jet_e",   &b_jet_e);

  TLorentzVector l1, l2;
  
  for (int iEvent = 0; iEvent < NEVT; iEvent++) {
    if (!pythia.next ())
      continue;

    b_z_n = 0;
    b_z_pt.clear ();
    b_z_eta.clear ();
    b_z_phi.clear ();
    b_z_m.clear ();

    for (int i = 0; i < pythia.event.size (); i++) {

      //if (!pythia.event[i].isFinal()) continue; // check if in final state

      //if (abs (pythia.event[i].id ()) != 11 && abs (pythia.event[i].id ()) != 13) continue; // check if electron or muon, resp.

      //l1.SetPtEtaPhiM (pythia.event[i].pT (), pythia.event[i].eta (), pythia.event[i].phi (), pythia.event[i].m ());

      //for (int j = 0; j < i; j++) {

      //  if (!pythia.event[j].isFinal()) continue; // check if in final state

      //  if (pythia.event[i].id () != -pythia.event[j].id ()) continue; // check if anti-particle of first particle

      //  l2.SetPtEtaPhiM (pythia.event[j].pT (), pythia.event[j].eta (), pythia.event[j].phi (), pythia.event[j].m ());

      //  if ((l1+l2).M () < 40) continue; // loose invariant mass cut to make sure these are from Z decays

      //  // reconstruct Z
      //  b_z_pt.push_back ((l1+l2).Pt ());
      //  b_z_eta.push_back ((l1+l2).Eta ());
      //  b_z_phi.push_back ((l1+l2).Phi ());
      //  b_z_m.push_back ((l1+l2).M ());
      //  b_z_n++;
      //}

      if (abs (pythia.event[i].id ()) != 23) continue; // check if Z

      b_z_pt.push_back (pythia.event[i].pT ());
      b_z_eta.push_back (pythia.event[i].eta ());
      b_z_phi.push_back (pythia.event[i].phi ());
      b_z_m.push_back (pythia.event[i].m ());
      b_z_n++;
    }

    if (b_z_n == 0) {
      iEvent--;
      continue;
    }

    antikT4->analyze (pythia.event);

    b_code = pythia.info.code ();
    b_id1 = pythia.info.id1pdf ();
    b_id2 = pythia.info.id2pdf ();
    b_x1pdf = pythia.info.x1pdf ();
    b_x2pdf = pythia.info.x2pdf ();
    b_Q =  pythia.info.QFac ();
    
    b_isValence1 = pythia.info.isValence1 ();
    b_isValence2 = pythia.info.isValence2 ();

    b_l_n = 0;
    b_l_pt.clear ();
    b_l_eta.clear ();
    b_l_phi.clear ();
    b_l_m.clear ();

    for (int i = 0; i < pythia.event.size (); i++) {

      if (!pythia.event[i].isFinal()) continue; // check if in final state

      if (abs (pythia.event[i].id ()) != 11 && abs (pythia.event[i].id ()) != 13) continue; // check if electron or muon, resp.

      b_l_pt.push_back (pythia.event[i].pT ());
      b_l_eta.push_back (pythia.event[i].eta ());
      b_l_phi.push_back (pythia.event[i].phi ());
      b_l_m.push_back (pythia.event[i].m ());
      b_l_n++;
    }

    b_jet_n = 0;
    b_jet_pt.clear ();
    b_jet_eta.clear ();
    b_jet_phi.clear ();
    b_jet_e.clear ();

    for (int i = 0; i < antikT4->sizeJet (); i++) {

      bool matchesLepton = false;
      for (int j = 0; !matchesLepton && j < b_l_n; j++) {
        matchesLepton = DeltaR (antikT4->p (i).eta (), b_l_eta.at (j), antikT4->phi (i), b_l_phi.at (j)) < 0.2;
      }
      if (matchesLepton) continue;

      b_jet_pt.push_back (antikT4->pT (i));
      b_jet_eta.push_back (antikT4->p (i).eta ());
      b_jet_phi.push_back (antikT4->phi (i));
      b_jet_e.push_back (antikT4->p (i).e ());
      b_jet_n++;
    }
    
    t->Fill();

    if (iEvent % (NEVT/100) == 0)
      std::cout << iEvent / (NEVT/100) << "\% done...\r" << std::flush;
  }

  pythia.stat();
  
  f->Write();
  f->Close();

  return 0;
}
