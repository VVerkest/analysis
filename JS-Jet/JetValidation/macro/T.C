#define T_cxx
#include "T.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double delta_phi(double phi1, double phi2){
  double dPhi = phi1 - phi2;
  while (dPhi>=M_PI) { dPhi -= 2.*M_PI; }
  while (dPhi<=-M_PI) { dPhi += 2.*M_PI; }
  return dPhi;
};

//DST_JET_run2pp_ana450_2024p009-00049219
// total sum E: no cut
// N towers with cut

void T::Loop()
{
//   In a ROOT session, you can do:
  /*
root
.L T.C
T t
t.Loop();

hLeadPtVsUEenergy->Draw("COLZ")
new TCanvas;
hLeadingJetSpectrum->Draw("E");

   */
  
//  double minE = 0.1;
  double minE = 0.2;
//  double minE = -10.;
//  TFile *outfile = new TFile("analysis_output_ana462_2024p007_00049219_noJetCut_noMinE.root","RECREATE");
//  TFile *outfile = new TFile("analysis_output_ana462_2024p007_00049219_jetBkgCut_100MeV.root","RECREATE");
//  TFile *outfile = new TFile("analysis_output_ana462_2024p007_00049219_jetBkgCut_noMinE.root","RECREATE");
  TFile *outfile = new TFile("analysis_output_JetVal_pythia8_Jet10_0_200MeV.root","RECREATE");
//  TFile *outfile = new TFile("analysis_output_JetVal_pythia8_Jet10_0_noMinE.root","RECREATE");

  //      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries

  
//  Int_t           m_event;
//  Int_t           nJet;
//  Int_t           cent;
//  Float_t         zvtx;
//  Float_t         b;
//  vector<int>     *id;
//  vector<int>     *nComponent;
//  vector<bool>    *triggerVector;
//  vector<float>   *eta;
//  vector<float>   *phi;
//  vector<float>   *e;
//  vector<float>   *pt;
//  vector<float>   *UE_energy;
//  vector<float>   *UE_eta;
//  vector<float>   *UE_phi;
//  vector<int>   *UE_caloID;
//  vector<float>   *cleta;
//  vector<float>   *clphi;
//  vector<float>   *cle;
//  vector<float>   *clecore;
//  vector<float>   *clpt;
//  vector<float>   *clprob;
//  vector<float>   *pt_unsub;
//  vector<float>   *subtracted_et;

// WIDE AXIS RANGES
  TH1D *hLeadingJetSpectrum = new TH1D("hLeadingJetSpectrum",";leading jet p_{T} [GeV]",125,0.,125.);
  TH1D *hTriggeredLeadJetSpectrum = new TH1D("hTriggeredLeadJetSpectrum",";leading jet p_{T} [GeV]",125,0.,125.);
  TH1D *hUE_spectrum = new TH1D("hUE_spectrum",";UE E [GeV]",60,-15.,15.);
  TH2D *hLeadingJetEtaPhi = new TH2D("hLeadingJetEtaPhi",";leading jet #eta;leading jet #phi",100,-1.,1.,100,-M_PI,M_PI);
  TH2D *hLeadingJetEtaPhi_pTweighted = new TH2D("hLeadingJetEtaPhi_pTweighted",";leading jet #eta;leading jet #phi",100,-1.,1.,100,-M_PI,M_PI);
  TH2D *hLeadPtVsUEenergy = new TH2D("hLeadPtVsUEenergy",";leading jet p_{T} [GeV];UE E [GeV]",125,0.,125.,80,-15.,25.);
  TH2D *hLeadPtVsUE_spectrum = new TH2D("hLeadPtVsUE_spectrum",";leading jet p_{T} [GeV];UE E spectrum [GeV]",125,0.,125.,50,-5.,20.);
  TH2D *hTotalClusterEnergyVsUEenergy = new TH2D("hTotalClusterEnergyVsUEenergy",";EMCal cluster E [GeV];UE E [GeV]",100,0.,50,50,-20.,30.);
  TH2D *h_nUEtowersVSjetPt = new TH2D("h_nUEtowersVSjetPt",";leading jet p_{T} [GeV];# UE towers",125,0.,125.,50,0.,50.);
  TH2D *h_nUEclustersVSjetPt = new TH2D("h_nUEclustersVSjetPt",";leading jet p_{T} [GeV];# UE clusters",125,0.,125.,120,0.,120.);

  TH1D *hUEcluster_spectrum = new TH1D("hUEcluster_spectrum",";UE cluster E [GeV]",30,0.,15.);
  TH2D *hLeadPtVsUEclusterEnergy = new TH2D("hLeadPtVsUEclusterEnergy",";leading jet p_{T} [GeV];UE cluster E [GeV]",125,0.,125.,80,-10.,30.);
  TH2D *hLeadPtVsUEcluster_spectrum = new TH2D("hLeadPtVsUEcluster_spectrum",";leading jet p_{T} [GeV];UE cluster E spectrum [GeV]",125,0.,125.,30,0.,15.);

  TH3D *hLeadJetPhi_caloE = new TH3D("hLeadJetPhi_caloE",";leading jet #phi;#Sigma calo E;caloID",100,-M_PI,M_PI,200,-100.,100.,3,0.,3.);
  TH3D *hLeadJetPtPhi_nTowers = new TH3D("hLeadJetPtPhi_nTowers",";leading jet p_{T} [GeV];leading jet #phi;number of towers",125,0.,125.,100,-M_PI,M_PI,1000,0.,1000.);
  
//  // NARROWER AXIS RANGES
//  TH1D *hLeadingJetSpectrum = new TH1D("hLeadingJetSpectrum",";leading jet p_{T} [GeV]",125,0.,125.);
//  TH1D *hTriggeredLeadJetSpectrum = new TH1D("hTriggeredLeadJetSpectrum",";leading jet p_{T} [GeV]",125,0.,125.);
//  TH1D *hUE_spectrum = new TH1D("hUE_spectrum",";UE E [GeV]",60,0.,15.);
//  TH2D *hLeadingJetEtaPhi = new TH2D("hLeadingJetEtaPhi",";leading jet #eta;leading jet #phi",100,-1.,1.,100,-M_PI,M_PI);
//  TH2D *hLeadingJetEtaPhi_pTweighted = new TH2D("hLeadingJetEtaPhi_pTweighted",";leading jet #eta;leading jet #phi",100,-1.,1.,100,-M_PI,M_PI);
//  TH2D *hLeadPtVsUEenergy = new TH2D("hLeadPtVsUEenergy",";leading jet p_{T} [GeV];UE E [GeV]",125,0.,125.,40,0.,20.);
//  TH2D *hLeadPtVsUE_spectrum = new TH2D("hLeadPtVsUE_spectrum",";leading jet p_{T} [GeV];UE E spectrum [GeV]",125,0.,125.,40,0.,20.);
//  TH2D *hTotalClusterEnergyVsUEenergy = new TH2D("hTotalClusterEnergyVsUEenergy",";EMCal cluster E [GeV];UE E [GeV]",100,0.,50.,100,0.,50.);
//  TH2D *h_nUEtowersVSjetPt = new TH2D("h_nUEtowersVSjetPt",";leading jet p_{T} [GeV];# UE towers",125,0.,125.,1000,0.,1000.);
//  TH2D *h_nUEclustersVSjetPt = new TH2D("h_nUEclustersVSjetPt",";leading jet p_{T} [GeV];# UE clusters",125,0.,125.,100,0.,100.);
//
//  TH1D *hUEcluster_spectrum = new TH1D("hUEcluster_spectrum",";UE cluster E [GeV]",60,0.,15.);
//  TH2D *hLeadPtVsUEclusterEnergy = new TH2D("hLeadPtVsUEclusterEnergy",";leading jet p_{T} [GeV];UE cluster E [GeV]",125,0.,125.,100,0.,50.);
//  TH2D *hLeadPtVsUEcluster_spectrum = new TH2D("hLeadPtVsUEcluster_spectrum",";leading jet p_{T} [GeV];UE cluster E spectrum [GeV]",125,0.,125.,30,0.,15.);
//
//  TH3D *hLeadJetPhi_caloE = new TH3D("hLeadJetPhi_caloE",";leading jet #phi;#Sigma calo E;caloID",100,-M_PI,M_PI,100,0.,50.,3,0.,3.);
//  TH3D *hLeadJetPtPhi_nTowers = new TH3D("hLeadJetPtPhi_nTowers",";leading jet p_{T} [GeV];leading jet #phi;number of towers",125,0.,125.,100,-M_PI,M_PI,500,0.,500.);

  
/*
  const int regions = 3;
  const char *regionName = {"transverse","towards","away"};
 
  TH1D *hUE_spectrum[regions];
  TH1D *hUEcluster_spectrum[regions];

  TH2D *hLeadPtVsUEenergy[regions];
  TH2D *hLeadPtVsUE_spectrum[regions];
  TH2D *hTotalClusterEnergyVsUEenergy[regions];
  TH2D *h_nUEtowersVSjetPt[regions];
  TH2D *hLeadPtVsUEclusterEnergy[regions];
  TH2D *hLeadPtVsUEcluster_spectrum[regions];
*/
  
  const int layers = 3;
  const char *caloLabel[layers] = {"oHCal","iHCal","EMCal"};

  TH2D *hLeadPtVsUEenergy_layer[layers];
  TH2D *hLeadPtVsUE_spectrum_layer[layers];
  TH2D *h_nUEtowersVSjetPt_layer[layers];
  TH1D *hUE_spectrum_layer[layers];
  
  for (int i=0; i<layers; ++i) {
    hLeadPtVsUEenergy_layer[i] = (TH2D*) hLeadPtVsUEenergy->Clone(Form("hLeadPtVsUEenergy_%s",caloLabel[i]));
    hLeadPtVsUE_spectrum_layer[i] = (TH2D*) hLeadPtVsUE_spectrum->Clone(Form("hLeadPtVsUE_spectrum_%s",caloLabel[i]));
    h_nUEtowersVSjetPt_layer[i] = (TH2D*) h_nUEtowersVSjetPt->Clone(Form("h_nUEtowersVSjetPt_%s",caloLabel[i]));
    hUE_spectrum_layer[i] = (TH1D*) hUE_spectrum->Clone(Form("hUE_spectrum_%s",caloLabel[i]));
  }
  
  int trig_number[6] = {17, 18, 19, 21, 22, 23};
  
  const char *trigger_string[6] = {
    "(Trig17) Jet 8 GeV + MBD NS >= 1",
    "(Trig18) Jet 10 GeV + MBD NS >= 1",
    "(Trig19) Jet 12 GeV + MBD NS >= 1",
    "(Trig21) Jet 8 GeV",
    "(Trig22) Jet 10 GeV",
    "(Trig23) Jet 12 GeV"
  };
  const char *trigger_label[6] = {
    "Jet_8_GeV_MBDNS",
    "Jet_10_GeV_MBDNS",
    "Jet_12_GeV_MBDNS",
    "Jet_8_GeV",
    "Jet_10_GeV",
    "Jet_12_GeV"
  };
  TH1D *hLeadingJetSpectrum_trigger[6];
  for (int i=0; i<6; ++i) {
    hLeadingJetSpectrum_trigger[i] = new TH1D(Form("hLeadingJetSpectrum_%s",trigger_label[i]),Form("%s;leading jet p_{T} [GeV]",trigger_string[i]),100,0.,50.);
  }
  TH2D *hEMCalClusterEtaPhi = new TH2D("hEMCalClusterEtaPhi",";tower #eta;tower #phi",220,-1.1,1.1,156,-M_PI,M_PI);
  TH2D *hEMCalClusterEtaPhi_Eweighted = new TH2D("hEMCalClusterEtaPhi_Eweighted",";tower #eta;tower #phi",220,-1.1,1.1,156,-M_PI,M_PI);
    
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
  cout<<nentries<<" entries in tree"<<endl;
//  nentries = 10000;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
     if (ientry%10000==0) { cout<<"entry "<<ientry<<endl;}
     
//     if ( triggerVector->at(17)==0 && triggerVector->at(18)==0 && triggerVector->at(19)==0 && triggerVector->at(21)==0 && triggerVector->at(22)==0 && triggerVector->at(23)==0 ) {continue;}

//Very roughly I think the noise thresholds for the EMCal, IHCal and OHCal should be around 90 MeV for the EMCal, 15 MeV for the IHCal and 90 MeV for the OHCal so I feel like 0.1 GeV is a good cut?
     if (fabs(zvtx)>50.) {continue;}

     float lead_pt=0;
     float lead_eta = -9999;
     float lead_phi = -9999;
     float lead_e = -9999;
     
     for (int i=0; i<pt->size(); ++i) {
       if (e->at(i)<=0.) { continue; }

       if ( pt->at(i)>lead_pt ) {
         lead_pt = pt->at(i);
         lead_eta = eta->at(i);
         lead_phi = phi->at(i);
         lead_e = e->at(i);
       }
     } // end jet loop

//     if (lead_pt<1.) {continue;}
     if (lead_pt<10.) {continue;}
//     if (lead_pt<14.) {continue;}
     if (fabs(lead_eta)>0.7) {continue;}

     hLeadingJetSpectrum->Fill(lead_pt);
     hLeadingJetEtaPhi->Fill(lead_eta,lead_phi);
     hLeadingJetEtaPhi_pTweighted->Fill(lead_eta,lead_phi,lead_pt);
     if (triggerVector->size()>0) {
       bool trig8 = triggerVector->at(21);
       bool trig10 = triggerVector->at(22);
       bool trig12 = triggerVector->at(23);
       if ( (trig8) || (trig10 && !trig8) || (trig12 && !trig8 && !trig10) ) {
         hTriggeredLeadJetSpectrum->Fill(lead_pt);
       }
     }
     int nUEtowers[layers] = {0,0,0};
     int nTowers = 0;

     float totalCaloE[layers] = {0,0,0};
     float totalUE_E = 0;
     float totalUE_E_layer[layers] = {0.,0.,0.};
     for (int i=0; i<UE_energy->size(); ++i) {
       
       if (UE_energy->at(i)<minE) {continue;}
       
       int caloID = UE_caloID->at(i);
       
       totalCaloE[caloID] += UE_energy->at(i);
       ++nTowers;
       
       double dPhi = lead_phi - UE_phi->at(i);
       while (dPhi>=M_PI) { dPhi -= 2.*M_PI; }
       while (dPhi<=-M_PI) { dPhi += 2.*M_PI; }
       if ( dPhi < M_PI/3. || dPhi > 2.*M_PI/3. ) {continue;}  // only transverse to jet in phi

//       if (fabs(UE_eta->at(i))>.9) {continue;}
       
       totalUE_E += UE_energy->at(i);
       totalUE_E_layer[caloID] += UE_energy->at(i);
       hLeadPtVsUE_spectrum->Fill(lead_pt,UE_energy->at(i));
       hLeadPtVsUE_spectrum_layer[caloID]->Fill(lead_pt,UE_energy->at(i));
       hUE_spectrum->Fill(UE_energy->at(i));
       hUE_spectrum_layer[caloID]->Fill(UE_energy->at(i));
       ++nUEtowers[caloID];
     }
     hLeadPtVsUEenergy->Fill(lead_pt,totalUE_E);
     for (int i=0; i<layers; ++i){
       hLeadPtVsUEenergy_layer[i]->Fill(lead_pt,totalUE_E_layer[i]);
       h_nUEtowersVSjetPt_layer[i]->Fill(lead_pt,nUEtowers[i]);
     }
     h_nUEtowersVSjetPt->Fill(lead_pt,nUEtowers[0]+nUEtowers[1]+nUEtowers[2]);
     
     int nClusters = 0;
     float totalCluster_E = 0;
     for (int i=0; i<cle->size(); ++i) {
       
       if (cle->at(i)<minE) {continue;}

       double dPhi = lead_phi - clphi->at(i);
       while (dPhi>=M_PI) { dPhi -= 2.*M_PI; }
       while (dPhi<=-M_PI) { dPhi += 2.*M_PI; }
       if ( dPhi < M_PI/3. || dPhi > 2.*M_PI/3. ) {continue;}  // only transverse to jet in phi

       
       totalCluster_E += cle->at(i);
       hEMCalClusterEtaPhi->Fill(cleta->at(i),clphi->at(i));
       hEMCalClusterEtaPhi_Eweighted->Fill(cleta->at(i),clphi->at(i),cle->at(i));
       hLeadPtVsUEcluster_spectrum->Fill(lead_pt,cle->at(i));
       hUEcluster_spectrum->Fill(cle->at(i));
       ++nClusters;
     }
     hTotalClusterEnergyVsUEenergy->Fill(totalCluster_E,totalUE_E);
     hLeadPtVsUEclusterEnergy->Fill(lead_pt,totalCluster_E);
     h_nUEclustersVSjetPt->Fill(lead_pt,nClusters);

     for (int i=0; i<layers; ++i) {
       hLeadJetPhi_caloE->Fill(lead_phi,totalCaloE[i],i+0.5);
     }
     hLeadJetPtPhi_nTowers->Fill(lead_pt,lead_phi,nTowers);

     if (triggerVector->size()>0) {
       for (int j=0; j<6; ++j) {
         int trig = trig_number[j];
         if (triggerVector->at(trig)==1) {
           hLeadingJetSpectrum_trigger[j]->Fill(lead_pt);
         }
       }
     }

   }
  
  hEMCalClusterEtaPhi->Write();
  hEMCalClusterEtaPhi_Eweighted->Write();
  hTotalClusterEnergyVsUEenergy->Write();
  hLeadingJetSpectrum->Write();
  hTriggeredLeadJetSpectrum->Write();
  hUE_spectrum->Write();
  hUEcluster_spectrum->Write();
  hLeadingJetEtaPhi->Write();
  hLeadingJetEtaPhi_pTweighted->Write();
  hLeadPtVsUEenergy->Write();
  hLeadPtVsUE_spectrum->Write();
  hLeadPtVsUEcluster_spectrum->Write();
  hLeadPtVsUEclusterEnergy->Write();
  h_nUEtowersVSjetPt->Write();
  h_nUEclustersVSjetPt->Write();
  hLeadJetPhi_caloE->Write();
  hLeadJetPtPhi_nTowers->Write();
  for (int i=0; i<6; ++i) {
    hLeadingJetSpectrum_trigger[i]->Write();
  }
  for (int i=0; i<layers; ++i) {
    hUE_spectrum_layer[i]->Write();
    hLeadPtVsUEenergy_layer[i]->Write();
    hLeadPtVsUE_spectrum_layer[i]->Write();
    h_nUEtowersVSjetPt_layer[i]->Write();
  }
  
//  hEMCalClusterEtaPhi->Write();
//  hEMCalClusterEtaPhi_Eweighted->Write();
  
//hLeadingJetSpectrum
//hLeadingJetEtaPhi
//hLeadingJetEtaPhi_pTweighted
//hEMCalClusterEtaPhi
//hEMCalClusterEtaPhi_Eweighted
  
  cout<<"done!"<<endl;
}
