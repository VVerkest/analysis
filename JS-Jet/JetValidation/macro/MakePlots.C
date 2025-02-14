
void MakePlots(){
  
//  const char *infilename = "analysis_output_fullEta_200MeV.root";
//  const char *plotdir = "plots/200MeV/";
//
//  const char *infilename = "analysis_output_fullEta_100MeV.root";
//  const char *plotdir = "plots/100MeV/";
//
//  const char *infilename = "analysis_output_fullEta_150MeV.root";
//  const char *plotdir = "plots/150MeV/";
//
//  const char *infilename = "analysis_output_fullEta_50MeV.root";
//  const char *plotdir = "plots/50MeV/";
  
//  const char *infilename = "analysis_output_ana462_2024p007_00049219_noJetCut_100MeV.root";
//  const char *plotdir = "plots/100MeV/";
//  const char *infilename = "analysis_output_ana462_2024p007_00049219_jetBkgCut_100MeV.root";
//  const char *plotdir = "plots/jetBkgCut/100MeV/";

  const char *infilename = "analysis_output_JetVal_pythia8_Jet10_0_100MeV.root";
  const char *plotdir = "plots/sim/100MeV/";

  TFile *infile = new TFile(infilename,"READ");
  
  
  TH1D *hLeadingJetSpectrum = (TH1D*)infile->Get("hLeadingJetSpectrum");
  TH1D *hUE_spectrum = (TH1D*)infile->Get("hUE_spectrum");

  TH2D *hEMCalClusterEtaPhi = (TH2D*)infile->Get("hEMCalClusterEtaPhi");
  TH2D *hEMCalClusterEtaPhi_Eweighted = (TH2D*)infile->Get("hEMCalClusterEtaPhi_Eweighted");
  TH2D *hTotalClusterEnergyVsUEenergy = (TH2D*)infile->Get("hTotalClusterEnergyVsUEenergy");
  TH2D *hUEcluster_spectrum = (TH2D*)infile->Get("hUEcluster_spectrum");
  TH2D *hLeadingJetEtaPhi = (TH2D*)infile->Get("hLeadingJetEtaPhi");
  TH2D *hLeadingJetEtaPhi_pTweighted = (TH2D*)infile->Get("hLeadingJetEtaPhi_pTweighted");
  TH2D *hLeadPtVsUEenergy = (TH2D*)infile->Get("hLeadPtVsUEenergy");
  TH2D *hLeadPtVsUE_spectrum = (TH2D*)infile->Get("hLeadPtVsUE_spectrum");
  TH2D *hLeadPtVsUEcluster_spectrum = (TH2D*)infile->Get("hLeadPtVsUEcluster_spectrum");
  TH2D *hLeadPtVsUEclusterEnergy = (TH2D*)infile->Get("hLeadPtVsUEclusterEnergy");
  TH2D *h_nUEtowersVSjetPt = (TH2D*)infile->Get("h_nUEtowersVSjetPt");
  TH2D *h_nUEclustersVSjetPt = (TH2D*)infile->Get("h_nUEclustersVSjetPt");

  
  const int layers = 3;
  const char *caloLabel[layers] = {"oHCal","iHCal","EMCal"};
  int trig_number[6] = {17, 18, 19, 21, 22, 23};
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
    hLeadingJetSpectrum_trigger[i] = (TH1D*)infile->Get(Form("hLeadingJetSpectrum_%s",trigger_label[i]));
  }
  
  TH2D *hLeadPtVsUEenergy_layer[layers];
  TH2D *hLeadPtVsUE_spectrum_layer[layers];
  TH2D *h_nUEtowersVSjetPt_layer[layers];
  TH1D *hUE_spectrum_layer[layers];

  for (int i=0; i<layers; ++i) {
    hLeadPtVsUEenergy_layer[i] = (TH2D*) infile->Get(Form("hLeadPtVsUEenergy_%s",caloLabel[i]));
    hLeadPtVsUE_spectrum_layer[i] = (TH2D*) infile->Get(Form("hLeadPtVsUE_spectrum_%s",caloLabel[i]));
    h_nUEtowersVSjetPt_layer[i] = (TH2D*) infile->Get(Form("h_nUEtowersVSjetPt_%s",caloLabel[i]));
    hUE_spectrum_layer[i] = (TH1D*) infile->Get(Form("hUE_spectrum_%s",caloLabel[i]));
  }
  
  
  TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  hLeadingJetSpectrum->Draw("E");
  c1->SaveAs(Form("%sLeadingJetSpectrum.pdf",plotdir),"PDF");

  TCanvas *c2 = new TCanvas("c2");
  c2->SetLogz();
  hLeadingJetEtaPhi->Draw("COLZ");
  c2->SaveAs(Form("%sLeadingJetEtaPhi.pdf",plotdir),"PDF");

  TCanvas *c3 = new TCanvas("c3");
  c3->SetLogz();
  hLeadingJetEtaPhi_pTweighted->Draw("COLZ");
  c3->SaveAs(Form("%sLeadingJetEtaPhi_pTweighted.pdf",plotdir),"PDF");
  
  TH1D *hLeadPtVsUEenergy_pfx = (TH1D*) hLeadPtVsUEenergy->ProfileX("hLeadPtVsUEenergy_pfx");
  TCanvas *c4 = new TCanvas("c4");
  c4->SetLogz();
  hLeadPtVsUEenergy->Draw("COLZ");
  hLeadPtVsUEenergy_pfx->Draw("ESAME");
  hLeadPtVsUEenergy_pfx->Fit("pol1","d","",5,35);
  c4->SaveAs(Form("%sLeadPtVsUEenergy.pdf",plotdir),"PDF");

  TH1D *hLeadPtVsUE_spectrum_pfx = (TH1D*) hLeadPtVsUE_spectrum->ProfileX("hLeadPtVsUE_spectrum_pfx");
  TCanvas *c5 = new TCanvas("c5");
  c5->SetLogz();
  hLeadPtVsUE_spectrum->Draw("COLZ");
  hLeadPtVsUE_spectrum_pfx->Draw("ESAME");
  hLeadPtVsUE_spectrum_pfx->Fit("pol1","d","",5,33);
  c5->SaveAs(Form("%sLeadPtVsUEspectrum.pdf",plotdir),"PDF");

  TCanvas *c6 = new TCanvas("c6");
  c6->SetLogz();
  hTotalClusterEnergyVsUEenergy->Draw("COLZ");
  c6->SaveAs(Form("%sTotalClusterEnergyVsUEenergy.pdf",plotdir),"PDF");

  TCanvas *c7 = new TCanvas("c7");
  c7->SetLogy();
  hUE_spectrum->Draw("E");
  c7->SaveAs(Form("%sUEspectrum.pdf",plotdir),"PDF");

  TCanvas *c8 = new TCanvas("c8");
  c8->SetLogy();
  hLeadingJetSpectrum->Draw("EPLCPMC");
  for (int i=0; i<6; ++i) {
    hLeadingJetSpectrum_trigger[i]->Draw("EPLCPMCSAME");
  }
  c8->BuildLegend();
  c8->SaveAs(Form("%sLeadingJetSpectrum_trigger.pdf",plotdir),"PDF");

  TH1D *h_nUEtowersVSjetPt_pfx = (TH1D*) h_nUEtowersVSjetPt->ProfileX("h_nUEtowersVSjetPt_pfx");
  TCanvas *c9 = new TCanvas("c9");
  c9->SetLogz();
  h_nUEtowersVSjetPt->Draw("COLZ");
  h_nUEtowersVSjetPt_pfx->Draw("ESAME");
  c9->SaveAs(Form("%snUEtowersVSjetPt.pdf",plotdir),"PDF");

  TCanvas *c10 = new TCanvas("c10");
  c10->SetLogy();
  hUE_spectrum_layer[2]->SetLineColor(kRed);
  hUE_spectrum_layer[0]->SetLineColor(kBlack);
  hUE_spectrum_layer[1]->SetLineColor(kGreen+2);
  hUE_spectrum_layer[2]->SetMarkerColor(kRed);
  hUE_spectrum_layer[0]->SetMarkerColor(kBlack);
  hUE_spectrum_layer[1]->SetMarkerColor(kGreen+2);
  hUE_spectrum_layer[2]->Draw("ESAME");
  hUE_spectrum_layer[0]->Draw("ESAME");
  hUE_spectrum_layer[1]->Draw("ESAME");
  c10->BuildLegend();
  c10->SaveAs(Form("%sUE_spectrum_layers.pdf",plotdir),"PDF");
  
  TH1D *hLeadPtVsUEcluster_spectrum_pfx = (TH1D*) hLeadPtVsUEcluster_spectrum->ProfileX("hLeadPtVsUEcluster_spectrum_pfx");
  TCanvas *c11 = new TCanvas("c11");
  c11->SetLogz();
  hLeadPtVsUEcluster_spectrum->Draw("COLZ");
  hLeadPtVsUEcluster_spectrum_pfx->Draw("ESAME");
  hLeadPtVsUEcluster_spectrum_pfx->Fit("pol1","d","",5,35);
  c11->SaveAs(Form("%sLeadPtVsUEcluster_spectrum.pdf",plotdir),"PDF");

  TH1D *hLeadPtVsUEclusterEnergy_pfx = (TH1D*) hLeadPtVsUEclusterEnergy->ProfileX("hLeadPtVsUEclusterEnergy_pfx");
  TCanvas *c12 = new TCanvas("c12");
  c12->SetLogz();
  hLeadPtVsUEclusterEnergy->Draw("COLZ");
  hLeadPtVsUEclusterEnergy->ProfileX("hLeadPtVsUEclusterEnergy_pfx");
  hLeadPtVsUEclusterEnergy_pfx->Draw("ESAME");
  hLeadPtVsUEclusterEnergy_pfx->Fit("pol1","d","",5,33);
  c12->SaveAs(Form("%sLeadPtVsUEclusterEnergy.pdf",plotdir),"PDF");
  
  TCanvas *c13 = new TCanvas("c13");
  c13->SetLogy();
  hUEcluster_spectrum->Draw("E");
  c13->SaveAs(Form("%sUEcluster_spectrum.pdf",plotdir),"PDF");
  
  TH1D *h_nUEclustersVSjetPt_pfx = (TH1D*) h_nUEclustersVSjetPt->ProfileX("h_nUEclustersVSjetPt_pfx");
  TCanvas *c14 = new TCanvas("c14");
  c14->SetLogz();
  h_nUEclustersVSjetPt->Draw("COLZ");
  h_nUEclustersVSjetPt_pfx->Draw("ESAME");
  c14->SaveAs(Form("%snUEclustersVSjetPt.pdf",plotdir),"PDF");

}
