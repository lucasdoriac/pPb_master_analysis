// Makes several kinematic eta cuts and compare their magnitude as a sanity check.

/*
--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.

---------------Lucas Carvalho---------------
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

// ### Macro settings ###
int n_centralities = 3; // Number of centrality classes defined for the dataset.
double delta = 1e-6; // GeV
double EHFmin = 44.0; // GeV
double EHFmax = 250.0; // GeV
TString plot_extension = ".pdf"; // .png to test .pdf to final result.
//TString output_name = "kinCuts_Eta_5TeV"; // Output name of Fig. 1 reproduction. Plots <pT>(Nch) for n_centralities for both collision energies.
TString output_name = "kinCuts_Eta_8TeV"; // Output name of Fig. 1 reproduction. Plots <pT>(Nch) for n_centralities for both collision energies.
TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

// ##############################################################################
// ##############################################################################

// --- Energy cutoff values [GeV] for 5TeV data
std::vector<std::pair<double, double>> CBins_5 = {
        {2.5, 11.5},
        {11.5, 35.0},
        {35.0, 250.0}
    };

// --- Energy cutoff values [GeV] for 8TeV data
std::vector<std::pair<double, double>> CBins_8 = {
        {2.5, 14.5},
        {14.5, 44.0},
        {44.0, 250.0}
    };

// --- Structs ---
struct DataStruct {
    std::string path;
    std::string name;
    std::string label;
    std::string sufix;
};

// --- struct to 5TeV dataset
const DataStruct dataFile_5TeV = {
    "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root",
    "dataFile_5TeV",
    "5TeV",
    "_5TeV"
};

// --- struct to 8TeV dataset
const DataStruct dataFile_8TeV = {
    "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root",
    "dataFile_8TeV",
    "8TeV",
    "_8TeV"
};

// --- Function headers ---
TH1D* get_pt_histogram(const DataStruct& dataFile, double low_eta = -1.5, double high_eta = 1.5);
void draw_CMS_Header(TString latex_text = "#bf{CMS} #it{Work in Progress}", double x = 0.12, double y = 0.93, double text_size = 0.04, int align = 11);
void general_settings(TH1D* hist = nullptr, Color_t color = kBlack);


// --- main() ---
void kinCuts_Eta(){

    gROOT->SetBatch(kTRUE); // This tells ROOT to run in batch mode, i.e. no GUI or pop-ups.

    TH1D* hist_1 = nullptr;
    TH1D* hist_2 = nullptr;
    TH1D* hist_3 = nullptr;
    TH1D* hist_4 = nullptr;
    TH1D* hist_5 = nullptr;

	// One for each cut
	hist_1 = get_pt_histogram(dataFile_8TeV, -0.5, 0.5);
	hist_2 = get_pt_histogram(dataFile_8TeV, -0.8, 0.8);
	hist_3 = get_pt_histogram(dataFile_8TeV, -1.5, 1.5);
	hist_4 = get_pt_histogram(dataFile_8TeV, -2.0, 2.0);
	hist_5 = get_pt_histogram(dataFile_8TeV, -2.4, 2.4);

    TCanvas *c = new TCanvas("c", "kinematic eta cuts", 800, 600);
    gStyle->SetOptStat(0);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);

    general_settings(hist_5, kBlack);
    general_settings(hist_4, kBlue);
    general_settings(hist_3, kRed);
    general_settings(hist_2, kGreen + 2);
    general_settings(hist_1, kMagenta + 1);

    // first draw settings.
    hist_5->SetTitle("");
    hist_5->GetXaxis()->CenterTitle(true);
    hist_5->GetYaxis()->CenterTitle(true);
    hist_5->GetXaxis()->SetTitleOffset(1.);
    hist_5->GetYaxis()->SetTitleOffset(1.);
    hist_5->GetXaxis()->SetTitleFont(42);
    hist_5->GetYaxis()->SetTitleFont(42);
    hist_5->GetXaxis()->SetLabelFont(42);
    hist_5->GetYaxis()->SetLabelFont(42);
    hist_5->GetXaxis()->SetTitleSize(0.05);
    hist_5->GetYaxis()->SetTitleSize(0.044);
    hist_5->GetXaxis()->SetLabelSize(0.036);
    hist_5->GetYaxis()->SetLabelSize(0.036);

    // Draw on canvas c
    hist_5->GetXaxis()->SetTitle("p_{T} [GeV]");
    hist_5->GetYaxis()->SetTitle("N of tracks");
    hist_5->Draw("E1");
    hist_4->Draw("E1 SAME");
    hist_3->Draw("E1 SAME");
    hist_2->Draw("E1 SAME");
    hist_1->Draw("E1 SAME");

    auto leg = new TLegend(0.66, 0.65, 0.86, 0.88);
    leg->AddEntry(hist_5, "|#eta| < 2.4", "lep");
    leg->AddEntry(hist_4, "|#eta| < 2.0", "lep");
    leg->AddEntry(hist_3, "|#eta| < 1.5", "lep");
    leg->AddEntry(hist_2, "|#eta| < 0.8", "lep");
    leg->AddEntry(hist_1, "|#eta| < 0.5", "lep");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->Draw();

    draw_CMS_Header();
    draw_CMS_Header("pPb (186.0 nb^{#minus1}) 8.16 TeV", 0.93, 0.93, 0.038, 31);
    //draw_CMS_Header("pPb (0.509 nb^{#minus1}) 5.02 TeV", 0.93, 0.93, 0.038, 31);

    c->SetLogy();
    c->Update();
    TString full_path = base_output_path + output_name + plot_extension;
    c->SaveAs(full_path);
}


// --- Function definitions ---
TH1D* get_pt_histogram(const DataStruct& dataFile, double low_eta, double high_eta){

	// Return TH1D histogram of reconstructed pT tracks on the range p_T > 0.3 [GeV]. 
    // Projected from the TH3D histogram located at 'Analysis_histograms'.
    // For selected pseudorapidity window and centrality bin/class.
    TFile *file = TFile::Open(dataFile.path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "From " << __func__ << ": error while opening target file." << std::endl;
        return nullptr;
    }

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    if (!dir) {
        std::cerr << "From " << __func__ << ": directory Analysis_histograms not found." << std::endl;
        file->Close();
        return nullptr;
    }
    
    // Load TH3D histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    if (!hist_HFSumPb_vs_pt_eta){
        std::cerr << "From " << __func__ << ": error while loading TH3D histogram." << std::endl;
        file->Close();
        return nullptr;
    }
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    file->Close();
    std::cout << "\n\nRunning on " << dataFile.label << " dataset." << " Centrality class: [" << EHFmin << "," << EHFmax << "]GeV " << std::endl; 

    // Setting pseudorapidity window [low_eta, high_eta].
    int z_min = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(low_eta + delta);
    int z_max = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(high_eta - delta);
    // Makes pseudorapidity cut on [z_min, z_max].
    hist_HFSumPb_vs_pt_eta->GetZaxis()->SetRange(z_min, z_max);

    // Check pseudorapidity window.
    double lowEdge = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinLowEdge(z_min);
    double highEdge = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(z_max);
    printf("\n-> From %s: Integrating over seudorapidity window [%.1f,%.1f] \n", __func__, lowEdge, highEdge);
    //

    // YX-Projection to TH2D by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *hist_HFSumPb_vs_pt = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");
    if (!hist_HFSumPb_vs_pt){
        std::cerr << "From " << __func__ << ": error while loading TH2D histogram." << std::endl;
        return nullptr;
    }
    hist_HFSumPb_vs_pt->SetDirectory(0);

    // Setting centrality class defined by HF energy cutoffs EHFmin, EHFmax.
    int bin_min = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmax - delta);

    // Check selected bins for centrality class.
    lowEdge = hist_HFSumPb_vs_pt->GetXaxis()->GetBinLowEdge(bin_min);
    highEdge = hist_HFSumPb_vs_pt->GetXaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From %s: Integrating over the centrality interval [%.1f,%.1f]GeV \n", __func__, lowEdge, highEdge);
    //

    // Projects TH2D on TH1D for the defined centrality class.
    TH1D *hist_pT = (TH1D*)hist_HFSumPb_vs_pt->ProjectionY("", bin_min, bin_max);
    if (!hist_pT){
        std::cerr << "From " << __func__ << ": error while loading TH1D histogram." << std::endl;
        return nullptr;
    }
    hist_pT->SetDirectory(0);
    hist_pT->GetXaxis()->SetRangeUser(0.,3.);
    hist_pT->SetStats(0);

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

void general_settings(TH1D* hist, Color_t color){

    // Style histogram general settings.
    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(color);
    hist->SetMarkerSize(0.8);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
    hist->SetTitle("");
}