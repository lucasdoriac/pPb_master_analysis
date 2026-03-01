// Macro developed to illustrate the centrality class determination method used in the work.
// This work is part of the CMS Collaboration and uses CMS Preliminary Data.

/*
--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.

---------------Lucas Carvalho---------------*/

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

void centrality_visualization()
{
    TFile *f = TFile::Open("../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root");
    if (!f || f->IsZombie()) {
        Error("plotCentralityBands", "Cannot open input file");
        return;
    }

    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*) f->Get("Analysis_histograms/hist_HFSumPb_vs_pt_eta");

    // Define pseudorapidity window.
    double low_eta = -1.5;
    double high_eta = 1.5;
    double delta = 1e-6;

    int z_min = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(low_eta + delta);
    int z_max = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(high_eta - delta);
    hist_HFSumPb_vs_pt_eta->GetZaxis()->SetRange(z_min, z_max);

    // Projection to TH2 by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *h2 = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");
    h2->SetDirectory(0);

    //Canvas and plot SETTINGS
    TCanvas *c = new TCanvas("c_cent", "Centrality classification", 1200, 900);
    c->SetLeftMargin(0.09);
    c->SetRightMargin(0.14);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
    c->SetLogz();
    c->SetLogy();

    h2->SetTitle("");
    h2->GetXaxis()->SetTitleOffset(1.);
    h2->GetYaxis()->SetTitleOffset(0.9);
    h2->GetZaxis()->SetTitleOffset(0.9);
    h2->GetXaxis()->CenterTitle(true);
    h2->GetYaxis()->CenterTitle(true);
    h2->GetZaxis()->CenterTitle(true);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetYaxis()->SetLabelFont(42);
    h2->GetZaxis()->SetTitleFont(42);
    h2->GetZaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleSize(0.044);
    h2->GetYaxis()->SetTitleSize(0.044);
    h2->GetZaxis()->SetTitleSize(0.044);
    h2->GetXaxis()->SetLabelSize(0.036);
    h2->GetYaxis()->SetLabelSize(0.036);
    h2->GetZaxis()->SetLabelSize(0.036);
    h2->GetXaxis()->SetTickLength(0.012);
    h2->GetYaxis()->SetTickLength(0.012);
    h2->GetZaxis()->SetTickLength(0.015);
    h2->SetStats(0);
    h2->SetMinimum(0.5);
    h2->GetXaxis()->SetTitle("E_{T,sum}^{HF} [GeV]");
    h2->GetYaxis()->SetTitle("p_{T} [GeV]");
    h2->GetZaxis()->SetTitle("N of tracks");
    h2->GetYaxis()->SetRangeUser(0.4, h2->GetYaxis()->GetXmax());
    h2->GetXaxis()->SetRangeUser(0., 150.);
    h2->Draw("COLZ");

    double HF_30_80_low = 2.5;
    double HF_30_80_up  = 14.5;

    double HF_1_30_low  = 14.5;
    double HF_1_30_up   = 44.0;

    double HF_0_1_low   = 44.0;

    double yMin = h2->GetYaxis()->GetXmin();
    double yMax = h2->GetYaxis()->GetXmax();

    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();

// Helper lambda to draw a clipped vertical line
auto DrawVLine = [&](double x, int color, int style = 2){
    if (x <= xMin || x >= xMax) return;
    TLine *l = new TLine(x, yMin, x, yMax);
    l->SetLineColor(color);
    l->SetLineStyle(style);
    l->SetLineWidth(2);
    l->Draw("same");
};

    // Draw centrality boundaries
    // 30–80%
    DrawVLine(HF_30_80_low, kBlack);
    DrawVLine(HF_30_80_up,  kBlack);

    // 1–30%
    DrawVLine(HF_1_30_low,  kGreen+2);
    DrawVLine(HF_1_30_up,   kGreen+2);

    // 0–1%
    DrawVLine(HF_0_1_low,   kRed+1);

    TLatex latex;
    latex.SetTextFont(42);
    latex.SetTextSize(0.038);
    latex.SetTextAlign(22);
    latex.SetTextAngle(90);

    // Adjust y-position slightly below top.
    double yLabel = yMax * 0.45;
    // Draw vertical lines legends.
    latex.DrawLatex( HF_30_80_low + 3.5, yLabel, "30-80%");
    latex.DrawLatex( HF_1_30_low + 3.5, yLabel, "1-30%");
    latex.DrawLatex( HF_0_1_low + 3.5, yLabel, "0-1%");

TLatex latex2;
latex2.SetNDC();              // use normalized coordinates
latex2.SetTextSize(0.04);     // text size
latex2.SetTextFont(42);       // Helvetica-like
latex2.SetTextAlign(11);      // left-aligned, top

TString cmsText = "#bf{CMS} #it{Work in Progress}";
latex2.DrawLatex(0.1, 0.93, cmsText);

TLatex latex3;
latex3.SetNDC();
latex3.SetTextSize(0.042);
latex3.SetTextFont(42);
latex3.SetTextAlign(11);
latex3.DrawLatex(0.7, 0.8, "|#eta| < 1.5");

TLatex latex4;
latex4.SetNDC();
latex4.SetTextSize(0.042);
latex4.SetTextFont(42);
latex4.SetTextAlign(11);
latex4.DrawLatex(0.52, 0.93, "pPb (186.0 nb^{#minus1}) 8.16 TeV");

    c->SaveAs("../../../../mnt/c/Users/lucas/Documents/centrality_classification_HF_vs_pT.pdf");
}