// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HCALJETPHISHIFT_H
#define HCALJETPHISHIFT_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TTree;
class TFile;

class HCalJetPhiShift : public SubsysReco
{
public:
  
  HCalJetPhiShift(const std::string &name = "HCalJetPhiShift", const std::string &outputFile = "HCalJetPhiShift_test.root");
  
  ~HCalJetPhiShift() override;
  
  /** Called during initialization.
   Typically this is where you can book histograms, and e.g.
   register them to Fun4AllServer (so they can be output to file
   using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;
  
  /** Called for first event when run number is known.
   Typically this is where you may want to fetch data from
   database, because you know the run number. A place
   to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;
  
  /** Called for each event.
   This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;
  
  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;
  
  /// Called at the end of each run.
  int EndRun(const int runnumber) override;
  
  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;
  
  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;
  
  void Print(const std::string &what = "ALL") const override;
  
  void SetEventNumber(int event_number)
  {
    m_event = event_number;
  };
  
  void SetSingleParticleMode()
  {
    m_single_particle = true;
  };
  
  void SetJetMode()
  {
    m_single_particle = false;
  };
  
  bool IsSingleParticleMode()
  {
    return m_single_particle;
  }
    
private:
  std::string m_outputFileName;
  
  bool m_single_particle;
  
  //! Output Tree variables
  TTree *m_T;
  
  //! eventwise quantities
  int m_event;
  int m_nTow_in;
  int m_nTow_out;
  int m_nTow_emc;
  float m_eta;
  float m_phi;
  float m_e;
  float m_pt;
  float m_vx;
  float m_vy;
  float m_vz;

  //! towers and G4 hits
  std::vector<float> m_eta_in;
  std::vector<float> m_phi_in;
  std::vector<float> m_e_in;
  std::vector<float> m_eta_in_g4hit;
  std::vector<float> m_phi_in_g4hit;
  std::vector<float> m_e_in_g4hit;
  
  std::vector<float> m_eta_out;
  std::vector<float> m_phi_out;
  std::vector<float> m_e_out;
  std::vector<float> m_eta_out_g4hit;
  std::vector<float> m_phi_out_g4hit;
  std::vector<float> m_e_out_g4hit;
  
  std::vector<float> m_eta_emc;
  std::vector<float> m_phi_emc;
  std::vector<float> m_e_emc;
  std::vector<float> m_eta_emc_g4hit;
  std::vector<float> m_phi_emc_g4hit;
  std::vector<float> m_e_emc_g4hit;
  
  //! jets
  std::vector<float> m_truth_pt;
  std::vector<float> m_truth_eta;
  std::vector<float> m_truth_phi;
  std::vector<float> m_truth_mass;

  std::vector<float> m_jet_pt;
  std::vector<float> m_jet_eta;
  std::vector<float> m_jet_phi;
  std::vector<float> m_jet_mass;


  int FillTTree(PHCompositeNode *topNode);
};

#endif // HCALJETPHISHIFT_H
