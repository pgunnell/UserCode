#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TPad.h"
#include "TLegend.h"
#include "TVirtualPad.h"
#include <iostream>
#include "TCanvas.h"

void SetHistoDATA(TH1F* Histo){
  
  Histo->SetLineWidth(2);
  Histo->SetMarkerStyle(20);
  Histo->SetMarkerSize(1.3);
  Histo->SetMarkerColor(1);
  Histo->SetLineColor(1);
}


void SetHistoTTbar(TH1F* Histo){
  
  Histo->SetLineWidth(3);
  Histo->SetMarkerStyle(0);
  Histo->SetMarkerSize(0);
  Histo->SetMarkerColor(0);
  Histo->SetLineColor(6);
  Histo->SetLineStyle(1);

}

void SetHistoQCD(TH1F* Histo){
  
  Histo->SetLineWidth(3);
  Histo->SetMarkerStyle(2);
  Histo->SetMarkerSize(0.7);
  Histo->SetMarkerColor(2);
  Histo->SetFillColor(0);
  Histo->SetLineColor(2);
  Histo->SetLineStyle(3);

}

void Detector_Level_MC(char* HistoLQModel, char* HistoQCD, char* HistoTop, char* LabelName, char* SaveName) { 

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);

  TCanvas *c1=new TCanvas("c1","c1",100,100,600,600);

  TFile* LQModel = TFile::Open("Yourfile.root","READ");
  TFile* TTBar = TFile::Open("AnotherFile.root","READ");
  TFile* QCD = TFile::Open("AnotherAnotherFile.root","READ");

  QCD->cd("analysis");
  
  TH1F *HistoPlotQCD  = (TH1F*)gDirectory->Get(HistoQCD);

  TTBar->cd("analysis");
    
  TH1F *HistoPlotTTBar  = (TH1F*)gDirectory->Get(HistoTop);
  
  LQModel->cd("analysis");

  TH1F *HistoPlotLQModel  = (TH1F*)gDirectory->Get(HistoLQModel);

  ////////////////////////////////////////////////////Normalization
    
  HistoPlotLQModel->Scale(1./HistoPlotLQModel->Integral());
  HistoPlotTTBar->Scale(1./HistoPlotTTBar->Integral());
  HistoPlotQCD->Scale(1./HistoPlotQCD->Integral());
  
  SetHistoDATA(HistoPlotLQModel);
  SetHistoQCD(HistoPlotQCD);
  SetHistoTTbar(HistoPlotTTBar);

  //TLegend
  TLegend* leg20;
  leg20 = new TLegend(0.68,0.57,0.97,0.90);

  leg20->SetBorderSize(1);
  leg20->SetTextSize(0.058);
  leg20->SetLineColor(0);
  leg20->SetLineStyle(1);
  leg20->SetLineWidth(1);
  leg20->SetFillColor(0);
  leg20->SetFillStyle(1001);

  TH1F *HistoPlotQCDRatio=(TH1F*)HistoPlotQCD->Clone();
  TH1F *HistoPlotTTBarRatio=(TH1F*)HistoPlotTTBar->Clone();
  
  HistoPlotTTBarRatio->Divide(HistoPlotLQModel);
  HistoPlotQCDRatio->Divide(HistoPlotLQModel);

  HistoPlotQCDRatio->GetXaxis()->SetTitle(LabelName);
  HistoPlotQCDRatio->GetYaxis()->SetTitle("MC/LQModel  ");

  HistoPlotQCDRatio->SetMaximum(2.6);
  HistoPlotQCDRatio->SetMinimum(0.02);

  TPaveText *pt = new TPaveText(0.2,0.94,0.8,1, "NDC");
  pt->SetFillColor(0); // text is black on white
  pt->SetTextSize(0.049); 
  pt->SetTextAlign(12);
  pt->AddText("CMS Internal Simulation (13 TeV)");

  leg20->AddEntry(HistoPlotLQModel,"LQ model");
  leg20->AddEntry(HistoPlotQCD,"QCD");
  leg20->AddEntry(HistoPlotTTBar,"TTBar");

  c1->cd();

  TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.32);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.01);
  c1_1->SetBottomMargin(0.22);
  c1_1->SetRightMargin(0.01);
  c1_1->SetLeftMargin(0.145);  
 
  HistoPlotQCDRatio->GetYaxis()->SetLabelFont(42);
  HistoPlotQCDRatio->GetXaxis()->SetLabelFont(42);
  HistoPlotQCDRatio->GetYaxis()->SetTitleFont(42);
  HistoPlotQCDRatio->GetXaxis()->SetTitleFont(42);

  HistoPlotQCDRatio->GetYaxis()->SetTitleSize(0.12);
  HistoPlotQCDRatio->GetYaxis()->SetTitleOffset(0.6);
  HistoPlotQCDRatio->GetYaxis()->SetLabelSize(0.12);
  HistoPlotQCDRatio->GetXaxis()->SetLabelSize(0.11);
  HistoPlotQCDRatio->GetXaxis()->SetTitleSize(0.1);
  HistoPlotQCDRatio->GetXaxis()->SetTitleOffset(0.9);

  HistoPlotQCDRatio->Draw("hist");  
  HistoPlotTTBarRatio->Draw("hist same");  
  
  c1->cd();

  TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.99,0.99);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTopMargin(0.07);
  c1_2->SetBottomMargin(0.005);
  c1_2->SetRightMargin(0.01);
  c1_2->SetLeftMargin(0.145);
  c1_2->SetLogy(1);

  HistoPlotLQModel->GetYaxis()->SetTitle("A.U.");
  HistoPlotLQModel->GetYaxis()->SetTitleSize(0.08);
  HistoPlotLQModel->GetYaxis()->SetTitleOffset(0.7);
  HistoPlotLQModel->GetYaxis()->SetLabelSize(0.06);

  HistoPlotLQModel->SetMaximum(5);

  HistoPlotLQModel->Draw("e2");
  HistoPlotQCD->Draw("hist same");
  HistoPlotTTBar->Draw("hist same");
  leg20->Draw("same");
  pt->Draw("same");

  c1->SaveAs(SaveName);

  c1->Clear();

}

void PlottingFunction(){

  //here you specify the name of the histos you want to plot
  Detector_Level_MC("ptGENJet","ptGENJet","ptGENJet","Jet p_{T} [GeV]","jetpT.pdf");
  Detector_Level_MC("yGENJet","yGENJet","yGENJet","Jet y","jetY.pdf");

}
