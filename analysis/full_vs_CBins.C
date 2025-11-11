/*
Illustrate centrality bins on the 8TeV dataset.
This work is part of the CMS Collaboration and uses CMS Preliminary Data.
Contact: lucasdoriadecarvalho@gmail.com
*/

/*
--- Energy cutoff values used by the CMS ---
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.
*/

// --- Headers ---
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>

double delta = 1e-6; // GeV
double low_eta = -1.5; // Low edge of pseudorapidity window.
double high_eta = 1.5; // High edge of pseudorapidity window. The p_T distribution will be integrated over (low_eta, high_eta).
TString plot_extension = ".pdf"; // .png to test .pdf to final result.
TString output_name = "full_vs_CBins";
TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

//const std::string myFile = "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root";
const std::string myFile = "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root";

TH1D* cut_centrality(double EHFmin = 0.0, double EHFmax = 250.0);
void draw_CMS_Header(TString latex_text = "#bf{CMS} #it{Work in Progress}", double x = 0.10, double y = 0.93, double text_size = 0.04, int align = 11);

void full_vs_CBins(){

    // Get histograms for each centrality bin
    TH1D *hist_full = cut_centrality(0.0, 250.0);      // Full sample
    TH1D *hist_30_80 = cut_centrality(2.5, 14.5);      // 30–80%
    TH1D *hist_1_30  = cut_centrality(14.5, 44.0);     // 1–30%
    TH1D *hist_0_1   = cut_centrality(44.0, 250.0);    // 0–1%

    // Create canvas
    TCanvas *c = new TCanvas("c", "pT Distributions vs Centrality", 900, 700);
    gStyle->SetOptStat(0);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);

    // first draw settings.
    hist_full->SetTitle("");
    hist_full->GetXaxis()->CenterTitle(true);
    hist_full->GetYaxis()->CenterTitle(true);
    hist_full->GetXaxis()->SetTitleOffset(1.);
    hist_full->GetYaxis()->SetTitleOffset(1.);
    hist_full->GetXaxis()->SetTitleFont(42);
    hist_full->GetYaxis()->SetTitleFont(42);
    hist_full->GetXaxis()->SetLabelFont(42);
    hist_full->GetYaxis()->SetLabelFont(42);
    hist_full->GetXaxis()->SetTitleSize(0.05);
    hist_full->GetYaxis()->SetTitleSize(0.044);
    hist_full->GetXaxis()->SetLabelSize(0.036);
    hist_full->GetYaxis()->SetLabelSize(0.036);
    hist_full->GetXaxis()->SetTitle("p_{T} [GeV]");
    hist_full->GetYaxis()->SetTitle("N of tracks");
    //hist_full->GetYaxis()->SetTitle("N of tracks (#times 10^{9})");

    // Style and draw histograms
    hist_full->SetMarkerStyle(21);
    hist_full->SetMarkerColor(kBlack);
    hist_full->SetMarkerSize(0.8);
    hist_full->SetLineColor(kBlack);
    hist_full->SetLineWidth(2);
    hist_full->SetTitle("");

    hist_30_80->SetMarkerStyle(21);
    hist_30_80->SetMarkerColor(kBlue);
    hist_30_80->SetMarkerSize(0.8);
    hist_30_80->SetLineColor(kBlue);
    hist_30_80->SetLineWidth(2);

    hist_1_30->SetMarkerStyle(21);
    hist_1_30->SetMarkerColor(kRed);
    hist_1_30->SetMarkerSize(0.8);
    hist_1_30->SetLineColor(kRed);
    hist_1_30->SetLineWidth(2);

    hist_0_1->SetMarkerStyle(21);
    hist_0_1->SetMarkerColor(kGreen+2);
    hist_0_1->SetMarkerSize(0.8);
    hist_0_1->SetLineColor(kGreen+2);
    hist_0_1->SetLineWidth(2);

    //hist_full->Scale(1./1e9);
    //hist_30_80->Scale(1./1e9);
    //hist_1_30->Scale(1./1e9);
    //hist_0_1->Scale(1./1e9);
    hist_full->GetYaxis()->SetRangeUser(1e7, 1e11);

    hist_full->Draw("E1");
    hist_30_80->Draw("E1 SAME");
    hist_1_30->Draw("E1 SAME");
    hist_0_1->Draw("E1 SAME");

    // Add legend
    auto legend = new TLegend(0.17, 0.20, 0.38, 0.40); //x = 22, y = 23
    legend->AddEntry(hist_full, "All events", "lep");
    legend->AddEntry(hist_30_80, "30-80%", "lep");
    legend->AddEntry(hist_1_30, "1-30%", "lep");
    legend->AddEntry(hist_0_1, "0-1%", "lep");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetMargin(0.2);
    legend->SetEntrySeparation(0.04);
    legend->Draw();

    draw_CMS_Header();
    draw_CMS_Header("pPb (186.0 nb^{#minus1}) 8.16 TeV", 0.95, 0.93, 0.038, 31);

    c->SetLogy();
    c->Update();
    TString full_path = base_output_path + output_name + plot_extension;
    c->SaveAs(full_path);
}

TH1D* cut_centrality(double EHFmin, double EHFmax){
    TFile *f = TFile::Open(myFile.c_str(), "READ");

    TDirectory *dir = (TDirectory*)f->Get("Analysis_histograms");
    
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    f->Close();

    // Setting pseudorapidity window [low_eta, high_eta].
    int z_min = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(low_eta + delta);
    int z_max = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(high_eta - delta);
    hist_HFSumPb_vs_pt_eta->GetZaxis()->SetRange(z_min, z_max);

    // Projection to TH2 by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *hist_HFSumPb_vs_pt = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");

    // Setting centrality class defined by HF energy cutoffs EHFmin, EHFmax.
    int bin_min = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmax - delta);

    // Projects TH2 on TH1 for the defined centrality class.
    TH1D *hist_pT = (TH1D*)hist_HFSumPb_vs_pt->ProjectionY("", bin_min, bin_max);
    hist_pT->SetDirectory(0);
    hist_pT->GetXaxis()->SetRangeUser(0.,3.);

    return hist_pT;
}

void draw_CMS_Header(TString latex_text, double x, double y, double text_size, int align){

TLatex latex;
latex.SetNDC(); // use normalized coordinates
latex.SetTextSize(text_size);
latex.SetTextFont(42);
latex.SetTextAlign(align);

// Alignment code format: 10 × vertical + horizontal
// Horizontal Value Meaning
// 1    Left aligned    
// 2    Center aligned  
// 3    Right aligned
// Vertical Value
// 1   Bottom  
// 2   Middle  
// 3   Top
latex.DrawLatex(x, y, latex_text);
}