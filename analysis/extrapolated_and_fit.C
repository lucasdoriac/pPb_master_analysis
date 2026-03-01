// Macro to plot extrapolated bins, Hagedorn fit and CMS data on the range p_T < 3.0 GeV for visualization.

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
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-6; // GeV
TString plot_extension = ".pdf";
TString output_name = "extrapolated_and_fit";
TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

// ##############################################################################
// ##############################################################################


// --- Vector to store the different TGraph's created by the different pT cuts
std::vector<TGraphErrors*> graphs;

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


// --- Function headers ---
double get_n_events(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
TH1D* get_pT_TH1Dhistogram(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
TH1D* Hagedorn_extrapolation(TH1D* hist_pT, const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
void plot_graphs();
void draw_CMS_Header(TString latex_text = "#bf{CMS} #it{Work in Progress}", double x = 0.11, double y = 0.93, double text_size = 0.04, int align = 11);


// --- main() ---
void extrapolated_and_fit(){

	gROOT->SetBatch(kTRUE); // This tells ROOT to run in batch mode, i.e. no GUI or pop-ups.

	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-8);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);

	// Track cpu efficiency.
    TStopwatch timer;
    timer.Start();

    for(int i = 0; i < n_centralities; ++i){

    	auto [low5, high5] = CBins_5[i];
    	auto [low8, high8] = CBins_8[i];

        Hagedorn_extrapolation(get_pT_TH1Dhistogram(dataFile_5TeV, low5, high5), dataFile_5TeV, low5, high5);
        Hagedorn_extrapolation(get_pT_TH1Dhistogram(dataFile_8TeV, low8, high8), dataFile_8TeV, low8, high8);

    }// Ending the for loop.

	plot_graphs();

	timer.Stop();
    std::cout << "-> Job finished in "
              << timer.RealTime() << " seconds (wall time), "
              << timer.CpuTime()  << " seconds (CPU time).\n\n";	
}


// --- Function definitions ---
TH1D* get_pT_TH1Dhistogram(const std::string& filename, double EHFmin, double EHFmax){

    // Return TH1 histogram of p_T track distribution on the range p_T > 0.3 [GeV] projected from the TH3 histogram located at 'Analysis_histograms',
    // for selected centrality class and pseudorapidity window.
    TFile *file = TFile::Open(filename.c_str(), "READ");

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    
    // Load TH3 histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    file->Close();

    if(filename == dataFile_5TeV.path){
        std::cout << "\n\nRunning on " << dataFile_5TeV.label << " data set."
        << " Centrality class: [" << EHFmin << "," << EHFmax << "] GeV." << std::endl; 
    }
    else if(filename == dataFile_8TeV.path){
        std::cout << "\n\nRunning on " << dataFile_8TeV.label << " data set."
        << " Centrality class: [" << EHFmin << "," << EHFmax << "] GeV." << std::endl;
    }

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

    return hist_pT;
}

TH1D* Hagedorn_extrapolation(TH1D* hist_pT, const std::string& filename, double EHFmin, double EHFmax){

    TH1D *original_hist = nullptr;
    original_hist = (TH1D*)hist_pT->Clone("original_hist");
    original_hist->SetDirectory(0);
    printf("\n-> From make_hagedorn_extrapolation: original_hist norm = %.3e \n", original_hist->Integral());

    // Get norm only on the target range: 0.3 to 1.5 GeV.
    int bin_min = original_hist->GetXaxis()->FindBin(lower_pt_forFit + delta);
    int bin_max = original_hist->GetXaxis()->FindBin(upper_pt_forFit - delta);
    double original_hist_norm = original_hist->Integral(bin_min, bin_max); // Integral over [0.3,1.5].
    printf("\n-> From make_hagedorn_extrapolation: original_hist norm over [0.3, 1.5] GeV = %.3e \n\n\n", original_hist_norm);

    // Hagedorn TF1. Function declaration section.
    TF1* pT_fit;
    int cc_low = static_cast<int>(EHFmin);
    int cc_up = static_cast<int>(EHFmax);
    static int k_counter = 0;

    // Original.
    pT_fit = new TF1(Form("ptfit_%d_%d_%d",cc_low,cc_up, k_counter++),"[0]*x*pow(1.+1./sqrt(1.-[1]*[1])*(sqrt(x*x+[4]*[4])-x*[1])/[3]/[2],-[3])",0.,upper_pt_forFit);
    pT_fit->SetParameters(7500000000.,0.3,0.1,6.,0.14);//We used these values for initialization        
    pT_fit->FixParameter(4,0.13957);//pion mass
    pT_fit->SetParLimits(2,0.,0.5);//kinetic freeze-out temperature in GeV    
    pT_fit->SetParLimits(3,4.,9.);//n - free parameter no physical meaning
    if(filename == dataFile_5TeV.path){
        pT_fit->FixParameter(1, 0.4034);// related to radial flow velocity - pPb 5TeV
    }
    else if(filename == dataFile_8TeV.path){
        pT_fit->FixParameter(1, 0.5010);//related to radial flow velocity - pPb 8TeV
    }

    // User-defined Hagedorn function fit.
    original_hist->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    original_hist->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    TFitResultPtr fitResult = original_hist->Fit(pT_fit,"NO R EX0 M S","",lower_pt_forFit,upper_pt_forFit);
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    double pValue = TMath::Prob(chi2, ndf);
    std::cout<<"chi2 : "<<chi2<<"; ndf : "<<ndf<<"; pValue : "<<pValue<<std::endl;

    // Create the histogram fit_hist that will be filled with a random generator of tracks following the Hagedorn distribution defined by the fit.
    int n_bins = bin_max; // number of bins between 0 and 1.5 GeV. 15 in this case.
    TH1D *fit_hist = new TH1D("fit_hist", "fit histogram", n_bins, 0., upper_pt_forFit);
    fit_hist->SetDirectory(0);

    int count = 0;
    double fit_hist_norm = 0.0;
    printf("\n\n-> Random generation of tracks with Hagedorn PDF. Sample size = %.1e \n", static_cast<double>(SAMPLE));
        while(count < SAMPLE){
            double x = pT_fit->GetRandom();
            if(x >= 0.3 && x < 1.5) count+=1;
            fit_hist->Fill(x);
        }
    double scale_factor = original_hist_norm / SAMPLE;
    printf("\n\n-> Scale factor = %.3f \n\n", scale_factor);
    fit_hist->Scale(scale_factor);

    fit_hist_norm = fit_hist->Integral(bin_min, bin_max);
    printf("\n-> Hagedorn histogram complete. Final hagedorn histogram norm over [0.3, 1.5] GeV = %.3e \n", fit_hist_norm);

    // Fill extrapolated tracks into histogram extrapolated_hist.
    TH1D *extrapolated_hist = (TH1D*)original_hist->Clone("extrapolated_hist");
    extrapolated_hist->SetDirectory(0);
    
    // Takes first three bin contents from fit_hist. The rest is taken from the real data from original_hist.
    double content, error;
    for(int i = 1; i <= 3; i++){
        content = fit_hist->GetBinContent(i);
        error = fit_hist->GetBinError(i);
        extrapolated_hist->SetBinContent(i, content);
        extrapolated_hist->SetBinError(i, error);
    }
    
    extrapolated_hist->ResetStats();
    return extrapolated_hist;
}

double get_n_events(const std::string& filename, double EHFmin, double EHFmax){
    std::string this_function = __func__;
    TFile *f = TFile::Open(filename.c_str(), "READ");

    TDirectory *dir = (TDirectory*)f->Get("QA_histograms");

    TH1D *hf_hist = (TH1D*)dir->Get("hfSumEtPb");
    hf_hist->SetDirectory(0);
    f->Close();

    int bin_min = hf_hist->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hf_hist->GetXaxis()->FindBin(EHFmax - delta);

    // If one wants to check interval of integration over transversal E_HF energy.
    double lowEdge   = hf_hist->GetXaxis()->GetBinLowEdge(bin_min);
    double highEdge  = hf_hist->GetXaxis()->GetBinUpEdge(bin_min);
    double lowEdge_   = hf_hist->GetXaxis()->GetBinLowEdge(bin_max);
    double highEdge_  = hf_hist->GetXaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From get_n_events: Integrating 'hfSumEtPb' from bin [%.1f,%.1f] GeV to bin [%.1f,%.1f] GeV \n\n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Calculates the number of events n_events detected by the Pb side HF, for the specified centrality class.
    double n_events = hf_hist->Integral(bin_min, bin_max);

    return n_events;
}

void plot_graphs(){

    // Draw them all
    TCanvas* c = new TCanvas("c_finalplot", "Final plot of cs_squared", 800, 600);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.038);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);

    draw_CMS_Header();
    draw_CMS_Header("|#eta| < 1.5", 0.18, 0.22, 0.04, 11);
    //draw_CMS_Header("non-interacting limit", 0.92,0.76, 0.038, 31);
    draw_CMS_Header("pPb (186.0 nb^{#minus1}) 8.16 TeV", 0.16, 0.84, 0.038, 11);
    draw_CMS_Header("pPb (0.509 nb^{#minus1}) 5.02 TeV", 0.16, 0.79, 0.038, 11);

    c->Modified();
    c->Update();
    // Path to final plot and save output plot.
    TString full_path = base_output_path + output_name + plot_extension;
    c->SaveAs(full_path);
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