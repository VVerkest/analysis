#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <ffamodules/CDBInterface.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4centrality/PHG4CentralityReco.h>

#include <calotrigger/TriggerRunInfoReco.h>

#include <Calo_Calib.C>
#include <HIJetReco.C>
#include <JetValidation.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libJetValidation.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4dst.so)


#endif
void Fun4All_JetVal(const char *filelistjetcalo = "dst_jet_calo.list",
                    const char *filelistjet = "dst_jet.list",
		    const char *filelistcalo = "dst_calo.list",
		    const char *outname = "outputest.root")
{

  
  Fun4AllServer *se = Fun4AllServer::instance();
  int verbosity = 0;

  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  
//  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(filelistjet);
//    int runnumber = runseg.first;
  
  rc->set_StringFlag("CDB_GLOBALTAG", "2024p009");
  rc->set_uint64Flag("TIMESTAMP",49219);

//  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(filelistjetcalo);
//  int runnumber = runseg.first;
//  if (runnumber != 0)
    rc->set_IntFlag("RUNNUMBER",49219);

//  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
//  rc->set_uint64Flag("TIMESTAMP", 6);
  
  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(verbosity);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem( cent );


  Enable::VERBOSITY = verbosity;
  
  Process_Calo_Calib();  //  DST_JETCALO has raw towers, so you can call process_calib() to produce the calibrated tower on the fly
  
//  HIJetReco();

  JetValidation *myJetVal = new JetValidation("AntiKt_Tower_r04_Sub1", "AntiKt_Truth_r04", outname);

  myJetVal->setPtRange(5, 100);
  myJetVal->setEtaRange(-0.7, 0.7);
  myJetVal->doUnsub(0);
  myJetVal->doTruth(0);
  myJetVal->doSeeds(0);
  se->registerSubsystem(myJetVal);

  TriggerRunInfoReco *triggerruninforeco = new TriggerRunInfoReco();
  se->registerSubsystem(triggerruninforeco);
 
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTjetcalo");
//  in1->AddListFile(filelistjetcalo,1);
  in1->AddFile(filelistjetcalo);
  se->registerInputManager(in1);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTjet");
//  in2->AddListFile(filelistjet,1);
  in2->AddFile(filelistjet);
  se->registerInputManager(in2);

  Fun4AllInputManager *ingeom = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  ingeom->AddFile(geoLocation);
  se->registerInputManager(ingeom);

//  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTgeom");
////  in3->AddListFile(filelistgeom,1);
//  in3->AddFile(filelistgeom);
//  se->registerInputManager(in3);

//  Fun4AllInputManager *in4 = new Fun4AllDstInputManager("DSTcalo");
//  in4->AddFile(filelistcalo);
//  se->registerInputManager(in4);
  

  
  se->run(-1);
  se->End();

  gSystem->Exit(0);
  return 0;

}
