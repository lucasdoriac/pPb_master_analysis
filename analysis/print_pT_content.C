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

// ### Macro settings ###
int SAMPLE = 1e+8; // Number of generated entries with TF1::GetRandom() following the Hagedorn probability distribution function.
int n_centralities = 3; // Number of centrality classes defined for the dataset.
double low_eta = -1.5;
double high_eta = 1.5;
double EHFmin = 0.;
double EHFmax = 250.;
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-6; // GeV
TString plot_extension = ".pdf";
TString output_name = "kinematic_cuts_pT"; // Output name of Fig. 1 reproduction. Plots <pT>(Nch) for n_centralities for both collision energies.
TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

// ##############################################################################
// ##############################################################################


// --- Vector to store the different TGraph's created by the different pT cuts
std::vector<TGraphErrors*> graphs;

// --- pT cuts ---
std::vector<std::pair<double, double>> pT_cuts = {
    {0., 150.}, // Original p_T > 0 GeV.
    {0.3, 150.}, // p_T > 0.3 GeV (no fit).
    {0.5, 150.}, // p_T > 0.5 GeV (no fit).
    {0., 3.}, // p_T < 3.0 GeV.
    {0., 5.}, // p_T < 5.0 GeV.
    {0., 10.} // p_T < 10.0 GeV.
};

// --- Energy cutoff values for 5TeV data [GeV]
std::vector<std::pair<double, double>> arr_5 = {
        {2.5, 11.5},
        {11.5, 35.0},
        {35.0, 250.0}
    };

// --- Energy cutoff values for 8TeV data [GeV]
std::vector<std::pair<double, double>> arr_8 = {
        {2.5, 14.5},
        {14.5, 44.0},
        {44.0, 250.0}
    };

// --- Structs ---
struct DataFile {
    std::string path;
    std::string name;
    std::string label;
    std::string sufix;
};

// --- struct to 5TeV dataset
const DataFile dataFile_5TeV = {
    "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root",
    "dataFile_5TeV",
    "5TeV",
    "_5TeV"
};

// --- struct to 8TeV dataset
const DataFile dataFile_8TeV = {
    "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root",
    "dataFile_8TeV",
    "8TeV",
    "_8TeV"
};

void print_pT_content(){

	// Return TH1 histogram of p_T track distribution on the range p_T > 0.3 [GeV] projected from the TH3 histogram located at 'Analysis_histograms',
    // for selected centrality class and pseudorapidity window.
    TFile *file = TFile::Open(dataFile_5TeV.path.c_str(), "READ");

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    
    // Load TH3 histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    file->Close();

    // Setting pseudorapidity window [low_eta, high_eta].
    int z_min = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(low_eta + delta);
    int z_max = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(high_eta - delta);
    hist_HFSumPb_vs_pt_eta->GetZaxis()->SetRange(z_min, z_max);

    // If one wants to check pseudorapidity window.
    double lowEdge = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinLowEdge(z_min);
    double highEdge = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(z_min);
    double lowEdge_ = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinLowEdge(z_max);
    double highEdge_ = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(z_max);
    printf("\n-> From get_proj_hist: Pseudorapidity window from bin [%.1f,%.1f] to bin [%.1f,%.1f] \n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Projection to TH2 by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *hist_HFSumPb_vs_pt = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");
    hist_HFSumPb_vs_pt->SetDirectory(0);

    // --- Print X-axis binning information for TH2 hist_HFSumPb_vs_pt ---
printf("hist_HFSumPb_vs_pt X-axis has %d bins\n",
       hist_HFSumPb_vs_pt->GetXaxis()->GetNbins());

printf("X-axis range: [%.6f, %.6f]\n",
       hist_HFSumPb_vs_pt->GetXaxis()->GetXmin(),
       hist_HFSumPb_vs_pt->GetXaxis()->GetXmax());

printf("First bin width: %.6f\n",
       hist_HFSumPb_vs_pt->GetXaxis()->GetBinWidth(1));

printf("---- Individual X-bin edges ----\n");
for (int b = 1; b <= hist_HFSumPb_vs_pt->GetXaxis()->GetNbins(); ++b) {
    double low  = hist_HFSumPb_vs_pt->GetXaxis()->GetBinLowEdge(b);
    double high = hist_HFSumPb_vs_pt->GetXaxis()->GetBinUpEdge(b);
    printf("X-bin %3d: [%.6f , %.6f]\n", b, low, high);
}
printf("---------------------------------\n");

    // Setting centrality class defined by HF energy cutoffs EHFmin, EHFmax.
    int bin_min = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmax - delta);

    // If one wants to check selected bins for centrality class.
    lowEdge = hist_HFSumPb_vs_pt->GetXaxis()->GetBinLowEdge(bin_min);
    highEdge = hist_HFSumPb_vs_pt->GetXaxis()->GetBinUpEdge(bin_min);
    lowEdge_ = hist_HFSumPb_vs_pt->GetXaxis()->GetBinLowEdge(bin_max);
    highEdge_ = hist_HFSumPb_vs_pt->GetXaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From get_proj_hist: Integrating from bin [%.1f,%.1f]GeV to bin [%.1f,%.1f]GeV \n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Projects TH2 on TH1 for the defined centrality class.
    TH1D *hist_pT = (TH1D*)hist_HFSumPb_vs_pt->ProjectionY("", bin_min, bin_max);
    hist_pT->SetDirectory(0);


    // --- Print binning information for hist_pT ---
printf("hist_pT has %d bins\n", hist_pT->GetNbinsX());
printf("X-axis range: [%.6f, %.6f]\n",
       hist_pT->GetXaxis()->GetXmin(),
       hist_pT->GetXaxis()->GetXmax());
printf("Bin width: %.6f\n", hist_pT->GetXaxis()->GetBinWidth(1));

printf("---- Individual bin edges ----\n");
for (int b = 1; b <= hist_pT->GetNbinsX(); ++b) {
    double low  = hist_pT->GetXaxis()->GetBinLowEdge(b);
    double high = hist_pT->GetXaxis()->GetBinUpEdge(b);
    printf("Bin %2d: [%.6f , %.6f]\n", b, low, high);
}
printf("--------------------------------\n");

}

