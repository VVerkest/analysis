//_________________________________________________________________________________..
//  Veronica Verkest, June 2023
//  module for performing an analysis on the phi shift between truth and calo jets
//_________________________________________________________________________________..

#include "HCalJetPhiShift.h"

#include "g4eval/CaloEvalStack.h"
#include "g4eval/CaloRawClusterEval.h"
#include "g4eval/CaloRawTowerEval.h"
#include "g4eval/CaloTruthEval.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <jetbase/JetMap.h>
#include <jetbase/Jetv1.h>
#include <jetbase/JetContainerv1.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

//____________________________________________________________________________..
HCalJetPhiShift::HCalJetPhiShift(const std::string &name, const std::string &outputFile):
SubsysReco(name),
m_outputFileName(outputFile),
m_T(nullptr),
m_event(-1),
//m_single_particle(false),
m_nTow_in(0),
m_nTow_out(0),
m_nTow_emc(0),
m_eta(),
m_phi(),
m_e(),
m_pt(),
m_vx(),
m_vy(),
m_vz(),
m_eta_in(),
m_phi_in(),
m_e_in(),
m_eta_in_g4hit(),
m_phi_in_g4hit(),
m_e_in_g4hit(),
m_eta_out(),
m_phi_out(),
m_e_out(),
m_eta_out_g4hit(),
m_phi_out_g4hit(),
m_e_out_g4hit(),
m_eta_emc(),
m_phi_emc(),
m_e_emc(),
m_eta_emc_g4hit(),
m_phi_emc_g4hit(),
m_e_emc_g4hit(),
m_truth_pt(),
m_truth_eta(),
m_truth_phi(),
m_truth_mass(),
m_jet_pt(),
m_jet_eta(),
m_jet_phi(),
m_jet_mass()

{
  std::cout << "HCalJetPhiShift::HCalJetPhiShift(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
HCalJetPhiShift::~HCalJetPhiShift()
{
  std::cout << "HCalJetPhiShift::~HCalJetPhiShift() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int HCalJetPhiShift::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "HCalJetPhiShift::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  std::cout << "HCalJetPhiShift::Init - Output to " << m_outputFileName << std::endl;
  
  // configure Tree
  m_T = new TTree("T", "HCalJetPhiShift Tree");
  m_T->Branch("event", &m_event);
  if (m_single_particle)
  {
    m_T->Branch("eta", &m_eta);
    m_T->Branch("phi", &m_phi);
    m_T->Branch("e", &m_e);
    m_T->Branch("pt", &m_pt);
  }
  m_T->Branch("vx", &m_vx);
  m_T->Branch("vy", &m_vy);
  m_T->Branch("vz", &m_vz);
  if (m_single_particle)
  {
    m_T->Branch("nTow_in", &m_nTow_in);
    m_T->Branch("eta_in", &m_eta_in);
    m_T->Branch("phi_in", &m_phi_in);
    m_T->Branch("e_in", &m_e_in);
    m_T->Branch("eta_in_g4hit", &m_eta_in_g4hit);
    m_T->Branch("phi_in_g4hit", &m_phi_in_g4hit);
    m_T->Branch("e_in_g4hit", &m_e_in_g4hit);
    m_T->Branch("nTow_out", &m_nTow_out);
    m_T->Branch("eta_out", &m_eta_out);
    m_T->Branch("phi_out", &m_phi_out);
    m_T->Branch("e_out", &m_e_out);
    m_T->Branch("eta_out_g4hit", &m_eta_out_g4hit);
    m_T->Branch("phi_out_g4hit", &m_phi_out_g4hit);
    m_T->Branch("e_out_g4hit", &m_e_out_g4hit);
    m_T->Branch("nTow_emc", &m_nTow_emc);
    m_T->Branch("eta_emc", &m_eta_emc);
    m_T->Branch("phi_emc", &m_phi_emc);
    m_T->Branch("e_emc", &m_e_emc);
    m_T->Branch("eta_emc_g4hit", &m_eta_emc_g4hit);
    m_T->Branch("phi_emc_g4hit", &m_phi_emc_g4hit);
    m_T->Branch("e_emc_g4hit", &m_e_emc_g4hit);
  }
  if (!m_single_particle)
  {
    m_T->Branch("truth_pt", &m_truth_pt);
    m_T->Branch("truth_eta", &m_truth_eta);
    m_T->Branch("truth_phi", &m_truth_phi);
    m_T->Branch("truth_mass", &m_truth_mass);
    m_T->Branch("jet_pt", &m_jet_pt);
    m_T->Branch("jet_eta", &m_jet_eta);
    m_T->Branch("jet_phi", &m_jet_phi);
    m_T->Branch("jet_mass", &m_jet_mass);
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCalJetPhiShift::InitRun(PHCompositeNode *topNode)
{
  std::cout << "HCalJetPhiShift::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCalJetPhiShift::process_event(PHCompositeNode *topNode)
{
  //  std::cout << "HCalJetPhiShift::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  ++m_event;
  
  //fill the tree
  FillTTree(topNode);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCalJetPhiShift::ResetEvent(PHCompositeNode *topNode)
{
  //  std::cout << "HCalJetPhiShift::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  m_nTow_in = 0;
  m_nTow_out = 0;
  m_nTow_emc = 0;

  m_eta_in.clear();
  m_phi_in.clear();
  m_e_in.clear();
  
  m_eta_out.clear();
  m_phi_out.clear();
  m_e_out.clear();
  
  m_eta_emc.clear();
  m_phi_emc.clear();
  m_e_emc.clear();
  
  m_eta_in_g4hit.clear();
  m_phi_in_g4hit.clear();
  m_e_in_g4hit.clear();
  
  m_eta_out_g4hit.clear();
  m_phi_out_g4hit.clear();
  m_e_out_g4hit.clear();
  
  m_eta_emc_g4hit.clear();
  m_phi_emc_g4hit.clear();
  m_e_emc_g4hit.clear();
  
  m_truth_pt.clear();
  m_truth_eta.clear();
  m_truth_phi.clear();
  m_truth_mass.clear();
  
  m_jet_pt.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_mass.clear();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCalJetPhiShift::EndRun(const int runnumber)
{
  std::cout << "HCalJetPhiShift::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCalJetPhiShift::End(PHCompositeNode *topNode)
{
  std::cout << "HCalJetPhiShift::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);
  
  m_T->Write();
  std::cout << "HCalJetPhiShift::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCalJetPhiShift::Reset(PHCompositeNode *topNode)
{
  std::cout << "HCalJetPhiShift::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void HCalJetPhiShift::Print(const std::string &what) const
{
  std::cout << "HCalJetPhiShift::Print(const std::string &what) const Printing info for " << what << std::endl;
}

//____________________________________________________________________________..
int HCalJetPhiShift::FillTTree(PHCompositeNode *topNode)
{  
  //  std::cout << "HCalJetPhiShift::FillTTree(PHCompositeNode *topNode)" << std::endl;
  
    // fill truth info from nodetree
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!truthinfo)
    {
      std::cout << PHWHERE << " ERROR: Can't find G4TruthInfo" << std::endl;
      exit(-1);
    }
    
    PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second;
         ++iter)
    {
      PHG4Particle* primary = iter->second;
      
      if (primary->get_e() < 0.)
      {
        continue;
      }
      
      float gpx = primary->get_px();
      float gpy = primary->get_py();
      float gpz = primary->get_pz();
      m_e = primary->get_e();
      
      m_pt = std::sqrt(gpx * gpx + gpy * gpy);
      m_eta = NAN;
      if (m_pt != 0.0)
      {
        m_eta = std::asinh(gpz / m_pt);
      }
      m_phi = std::atan2(gpy, gpx);
      
      PHG4VtxPoint* gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      m_vx = gvertex->get_x();
      m_vy = gvertex->get_y();
      m_vz = gvertex->get_z();
    }
  
  if (m_single_particle)
  {
    //calorimeter towers
    TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    RawTowerGeomContainer *tower_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(!towersIH3 || !towersOH3 || !towersEM3){
      std::cout
      <<"HCalJetPhiShift::process_event - Error cannot find raw tower node "
      
      << std::endl;
      exit(-1);
    }
    
    if(!tower_geomIH || !tower_geomOH || !tower_geomEM){
      std::cout
      <<"HCalJetPhiShift::process_event - Error cannot find raw tower geometry "
      
      << std::endl;
      exit(-1);
    }
    
    TowerInfo *tower_in, *tower_out;//, *tower_emc;
    
    const int n_channels_IH = (int) towersIH3->size();
    
    // Inner HCal
    for (int i_chan=0; i_chan<n_channels_IH; ++i_chan) {
      tower_in = towersIH3->get_tower_at_channel(i_chan);
      if (tower_in->get_energy()>0.) {
        ++m_nTow_in;
        
        unsigned int calokey = towersIH3->encode_key(i_chan);
        int ieta = towersIH3->getTowerEtaBin(calokey);
        int iphi = towersIH3->getTowerPhiBin(calokey);
        
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
        float tower_phi = tower_geomIH->get_tower_geometry(key)->get_phi();
        float tower_eta = tower_geomIH->get_tower_geometry(key)->get_eta();
        
        m_eta_in.push_back(tower_eta);
        m_phi_in.push_back(tower_phi);
        m_e_in.push_back(tower_in->get_energy());
      }
      
      // Outer HCal
      tower_out = towersOH3->get_tower_at_channel(i_chan);
      if (tower_out->get_energy()>0.) {
        ++m_nTow_out;
        
        unsigned int calokey = towersOH3->encode_key(i_chan);
        int ieta = towersOH3->getTowerEtaBin(calokey);
        int iphi = towersOH3->getTowerPhiBin(calokey);
        
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
        float tower_phi = tower_geomOH->get_tower_geometry(key)->get_phi();
        float tower_eta = tower_geomOH->get_tower_geometry(key)->get_eta();
        
        m_eta_out.push_back(tower_eta);
        m_phi_out.push_back(tower_phi);
        m_e_out.push_back(tower_out->get_energy());
      }
      
    }
    
    const int n_channels_EM = (int) towersEM3->size();
    for (int i_chan=0; i_chan<n_channels_EM; ++i_chan) {
      // EMCal
      //    tower_emc = towersEM3->get_tower_at_channel(i_chan);
      TowerInfo *tower_emc = towersEM3->get_tower_at_channel(i_chan);
      if (tower_emc->get_energy()>0.) {
        ++m_nTow_emc;
        
        unsigned int calokey = towersEM3->encode_key(i_chan);
        int ieta = towersEM3->getTowerEtaBin(calokey);
        int iphi = towersEM3->getTowerPhiBin(calokey);
        
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
        RawTowerGeom *tower_geom = tower_geomEM->get_tower_geometry(key);
        float tower_phi = tower_geom->get_phi();
        float tower_eta = tower_geom->get_eta();
        
        //      unsigned int calokey = towersEM3->encode_key(i_chan);
        //      int ieta = towersEM3->getTowerEtaBin(calokey);
        //      int iphi = towersEM3->getTowerPhiBin(calokey);
        //
        //      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
        //      float tower_phi = tower_geomEM->get_tower_geometry(key)->get_phi();
        //      float tower_eta = tower_geomEM->get_tower_geometry(key)->get_eta();
        
        m_eta_emc.push_back(tower_eta);
        m_phi_emc.push_back(tower_phi);
        m_e_emc.push_back(tower_emc->get_energy());
      }
      
    }
    
    PHG4HitContainer* _hcalout_hit_container = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_HCALOUT");
    PHG4HitContainer* _hcalin_hit_container = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_HCALIN");
    PHG4HitContainer* _cemc_hit_container = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_CEMC");
    // https://github.com/VVerkest/analysis/blob/master/EMCal-calibration/EMCalCalib/EMCalCalib.C
    
    // INNER HCAL
    if (_hcalin_hit_container){
      PHG4HitContainer::ConstRange hcalin_hit_range = _hcalin_hit_container->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hcalin_hit_range.first; hit_iter != hcalin_hit_range.second; hit_iter++){
        PHG4Hit *this_hit = hit_iter->second;
        assert(this_hit);
        const TVector3 hit(this_hit->get_avg_x(), this_hit->get_avg_y(), this_hit->get_avg_z());
        m_eta_in_g4hit.push_back(hit.Eta());
        m_phi_in_g4hit.push_back(hit.Phi());
        m_e_in_g4hit.push_back(this_hit->get_edep());
      }
    }
    
    // OUTER HCAL
    if (_hcalout_hit_container){
      PHG4HitContainer::ConstRange hcalout_hit_range = _hcalout_hit_container->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hcalout_hit_range.first; hit_iter != hcalout_hit_range.second; hit_iter++){
        PHG4Hit *this_hit = hit_iter->second;
        assert(this_hit);
        const TVector3 hit(this_hit->get_avg_x(), this_hit->get_avg_y(), this_hit->get_avg_z());
        m_eta_out_g4hit.push_back(hit.Eta());
        m_phi_out_g4hit.push_back(hit.Phi());
        m_e_out_g4hit.push_back(this_hit->get_edep());
      }
    }
    
    // EMCAL
    if (_cemc_hit_container){
      PHG4HitContainer::ConstRange cemc_hit_range = _cemc_hit_container->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = cemc_hit_range.first; hit_iter != cemc_hit_range.second; hit_iter++){
        PHG4Hit *this_hit = hit_iter->second;
        assert(this_hit);
        const TVector3 hit(this_hit->get_avg_x(), this_hit->get_avg_y(), this_hit->get_avg_z());
        m_eta_emc_g4hit.push_back(hit.Eta());
        m_phi_emc_g4hit.push_back(hit.Phi());
        m_e_emc_g4hit.push_back(this_hit->get_edep());
      }
    }
    
  }
  
  //  topNode->print();
  
  if (!m_single_particle)
  {
    
    // interface to reco jets
    std::string calojetNodeName = "AntiKt_Tower_r04";
    JetContainerv1* jets = (JetContainerv1*)findNode::getClass<JetContainerv1>(topNode, calojetNodeName);
    if (!jets)
    {
      std::cout
      << "HCalJetPhiShift::FillTTree - Error can not find DST Reco JetMap node "
      << calojetNodeName << std::endl;
      exit(-1);
    }
    
    if (jets->size()>0) {
      for (int iter = 0; iter < (int)jets->size(); ++iter)
      {
        Jet* jet = jets->get_jet(iter);
        if (jet->get_pt()<5.) { continue; }
        m_jet_pt.push_back(jet->get_pt());
        m_jet_eta.push_back(jet->get_eta());
        m_jet_phi.push_back(jet->get_phi());
        m_jet_mass.push_back(jet->get_mass());
      }
    }
//    topNode->print();
    // interface to truth jets
    std::string jetNodeName = "AntiKt_Truth_r04";
    JetMap* truth_jets = (JetMap*)findNode::getClass<JetMap>(topNode, jetNodeName);
    if (!truth_jets)
    {
      std::cout
      << "HCalJetPhiShift::FillTTree - Error can not find DST Truth JetMap node "
      << jetNodeName << std::endl;
      exit(-1);
    }
    
    for (JetMap::Iter iter = truth_jets->begin(); iter != truth_jets->end(); ++iter)
    {
      Jet* jet = iter->second;
      if (jet->get_pt()>=5.) {
        m_truth_pt.push_back(jet->get_pt());
        m_truth_eta.push_back(jet->get_eta());
        m_truth_phi.push_back(jet->get_phi());
        m_truth_mass.push_back(jet->get_mass());
      }
    }
  }
  
  m_T->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
