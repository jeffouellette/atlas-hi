//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 19 15:31:28 2018 by ROOT version 6.14/02
// from TTree tree/tree
// found on file: pPbZee.sub1.root
//////////////////////////////////////////////////////////

#ifndef pPbZeeStudy_h
#define pPbZeeStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class pPbZeeStudy {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNumber;
   Int_t           EventNumber;
   Int_t           LumiBlock;
   Bool_t          HLT_e15_lhloose_nod0;
   vector<float>   *electron_e;
   vector<float>   *electron_pt;
   vector<float>   *electron_e_raw;
   vector<float>   *electron_pt_raw;
   vector<float>   *electron_eta;
   vector<float>   *electron_etaBE2;
   vector<float>   *electron_d0Sig;
   vector<float>   *electron_dz0SinTheta;
   vector<float>   *electron_phi;
   vector<float>   *electron_charge;
   vector<float>   *electron_author;
   vector<bool>    *electron_isLoose;
   vector<bool>    *electron_isMedium;
   vector<bool>    *electron_isTight;
   vector<bool>    *electron_isMedium_HI;
   vector<bool>    *electron_isLoose_HI;
   vector<bool>    *electron_match_e15_nod0;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_HLT_e15_lhloose_nod0;   //!
   TBranch        *b_electron_e;   //!
   TBranch        *b_electron_pt;   //!
   TBranch        *b_electron_e_raw;   //!
   TBranch        *b_electron_pt_raw;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_etaBE2;   //!
   TBranch        *b_electron_d0Sig;   //!
   TBranch        *b_electron_dz0SinTheta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_charge;   //!
   TBranch        *b_electron_author;   //!
   TBranch        *b_electron_isLoose;   //!
   TBranch        *b_electron_isMedium;   //!
   TBranch        *b_electron_isTight;   //!
   TBranch        *b_electron_isMedium_HI;   //!
   TBranch        *b_electron_isLoose_HI;   //!
   TBranch        *b_electron_match_e15_nod0;   //!

   pPbZeeStudy(TTree *tree=0);
   virtual ~pPbZeeStudy();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pPbZeeStudy_cxx
pPbZeeStudy::pPbZeeStudy(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pPbZee.sub1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pPbZee.sub1.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

pPbZeeStudy::~pPbZeeStudy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pPbZeeStudy::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pPbZeeStudy::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pPbZeeStudy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   electron_e = 0;
   electron_pt = 0;
   electron_e_raw = 0;
   electron_pt_raw = 0;
   electron_eta = 0;
   electron_etaBE2 = 0;
   electron_d0Sig = 0;
   electron_dz0SinTheta = 0;
   electron_phi = 0;
   electron_charge = 0;
   electron_author = 0;
   electron_isLoose = 0;
   electron_isMedium = 0;
   electron_isTight = 0;
   electron_isMedium_HI = 0;
   electron_isLoose_HI = 0;
   electron_match_e15_nod0 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("HLT_e15_lhloose_nod0", &HLT_e15_lhloose_nod0, &b_HLT_e15_lhloose_nod0);
   fChain->SetBranchAddress("electron_e", &electron_e, &b_electron_e);
   fChain->SetBranchAddress("electron_pt", &electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron_e_raw", &electron_e_raw, &b_electron_e_raw);
   fChain->SetBranchAddress("electron_pt_raw", &electron_pt_raw, &b_electron_pt_raw);
   fChain->SetBranchAddress("electron_eta", &electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron_etaBE2", &electron_etaBE2, &b_electron_etaBE2);
   fChain->SetBranchAddress("electron_d0Sig", &electron_d0Sig, &b_electron_d0Sig);
   fChain->SetBranchAddress("electron_dz0SinTheta", &electron_dz0SinTheta, &b_electron_dz0SinTheta);
   fChain->SetBranchAddress("electron_phi", &electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron_charge", &electron_charge, &b_electron_charge);
   fChain->SetBranchAddress("electron_author", &electron_author, &b_electron_author);
   fChain->SetBranchAddress("electron_isLoose", &electron_isLoose, &b_electron_isLoose);
   fChain->SetBranchAddress("electron_isMedium", &electron_isMedium, &b_electron_isMedium);
   fChain->SetBranchAddress("electron_isTight", &electron_isTight, &b_electron_isTight);
   fChain->SetBranchAddress("electron_isMedium_HI", &electron_isMedium_HI, &b_electron_isMedium_HI);
   fChain->SetBranchAddress("electron_isLoose_HI", &electron_isLoose_HI, &b_electron_isLoose_HI);
   fChain->SetBranchAddress("electron_match_e15_nod0", &electron_match_e15_nod0, &b_electron_match_e15_nod0);
   Notify();
}

Bool_t pPbZeeStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pPbZeeStudy::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pPbZeeStudy::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pPbZeeStudy_cxx
