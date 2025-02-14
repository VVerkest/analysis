//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 18 12:58:26 2024 by ROOT version 6.24/04
// from TTree T/MyJetAnalysis Tree
// found on file: out/JetVal_run2pp_ana437_2024p007_00049270_0000.root
//////////////////////////////////////////////////////////

/*
 root
 .L T.C
 T t
 t.Loop();
 
 */

#ifndef T_h
#define T_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class T {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Int_t           m_event;
  Int_t           nJet;
  Int_t           cent;
  Float_t         zvtx;
  Float_t         b;
  vector<int>     *id;
  vector<int>     *nComponent;
  vector<bool>    *triggerVector;
  vector<float>   *eta;
  vector<float>   *phi;
  vector<float>   *e;
  vector<float>   *pt;
  vector<float>   *UE_energy;
  vector<float>   *UE_eta;
  vector<float>   *UE_phi;
  vector<int>   *UE_caloID;
  vector<float>   *cleta;
  vector<float>   *clphi;
  vector<float>   *cle;
  vector<float>   *clecore;
  vector<float>   *clpt;
  vector<float>   *clprob;
  vector<float>   *pt_unsub;
  vector<float>   *subtracted_et;
  
  // List of branches
  TBranch        *b_event;   //!
  TBranch        *b_nJet;   //!
  TBranch        *b_cent;   //!
  TBranch        *b_zvtx;   //!
  TBranch        *b_b;   //!
  TBranch        *b_id;   //!
  TBranch        *b_nComponent;   //!
  TBranch        *b_triggerVector;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_e;   //!
  TBranch        *b_pt;   //!
  TBranch        *b_UE_energy;   //!
  TBranch        *b_UE_eta;   //!
  TBranch        *b_UE_phi;   //!
  TBranch        *b_UE_caloID;   //!
  TBranch        *b_cleta;   //!
  TBranch        *b_clphi;   //!
  TBranch        *b_cle;   //!
  TBranch        *b_clecore;   //!
  TBranch        *b_clpt;   //!
  TBranch        *b_clprob;   //!
  TBranch        *b_pt_unsub;
  TBranch        *b_subtracted_et;
  
  T(TTree *tree=0);
  virtual ~T();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef T_cxx
T::T(TTree *tree) : fChain(0) 
{
  /*
   root
   .L T.C
   T t
   t.Loop();
   */
  
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.  const char *fileName = "JetVal_run2pp_ana437_2024p007_00049270.root";
  //  const char *fileName = "JetVal_run2pp_ana437_2024p007_00049219_calib.root";
  //  const char *fileName = "JetVal_run2pp_ana437_2024p007_00049219.root";
//  const char *fileName = "JetVal_run2pp_ana446_2024p007_00049219_allUE.root";
//  const char *fileName = "JetVal_run2pp_ana462_2024p007_00049219_jetBkgCut.root";
//  const char *fileName = "JetVal_run2pp_ana462_2024p007_00049219_noJetCut.root";
  const char *fileName = "JetVal_pythia8_Jet10_0.root";
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
    if (!f || !f->IsOpen()) {
      f = new TFile(fileName);
    }
    f->GetObject("T",tree);
    
  }
  Init(tree);
}

T::~T()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t T::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t T::LoadTree(Long64_t entry)
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

void T::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set object pointer
  id = 0;
  nComponent = 0;
  triggerVector = 0;
  eta = 0;
  phi = 0;
  e = 0;
  pt = 0;
  UE_energy = 0;
  cleta = 0;
  clphi = 0;
  cle = 0;
  clecore = 0;
  clpt = 0;
  clprob = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("m_event", &m_event, &b_event);
  fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
  fChain->SetBranchAddress("cent", &cent, &b_cent);
  fChain->SetBranchAddress("zvtx", &zvtx, &b_zvtx);
  fChain->SetBranchAddress("b", &b, &b_b);
  fChain->SetBranchAddress("id", &id, &b_id);
  fChain->SetBranchAddress("nComponent", &nComponent, &b_nComponent);
  fChain->SetBranchAddress("triggerVector", &triggerVector, &b_triggerVector);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
  fChain->SetBranchAddress("e", &e, &b_e);
  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("UE_energy", &UE_energy, &b_UE_energy);
  fChain->SetBranchAddress("UE_eta", &UE_eta, &b_UE_eta);
  fChain->SetBranchAddress("UE_phi", &UE_phi, &b_UE_phi);
  fChain->SetBranchAddress("UE_caloID", &UE_caloID, &b_UE_caloID);
  fChain->SetBranchAddress("cleta", &cleta, &b_cleta);
  fChain->SetBranchAddress("clphi", &clphi, &b_clphi);
  fChain->SetBranchAddress("cle", &cle, &b_cle);
  fChain->SetBranchAddress("clecore", &clecore, &b_clecore);
  fChain->SetBranchAddress("clpt", &clpt, &b_clpt);
  fChain->SetBranchAddress("clprob", &clprob, &b_clprob);
//  fChain->SetBranchAddress("pt_unsub", &pt_unsub, &b_pt_unsub);
//  fChain->SetBranchAddress("subtracted_et", &subtracted_et, &b_subtracted_et);

  Notify();
}

Bool_t T::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void T::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t T::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef T_cxx
