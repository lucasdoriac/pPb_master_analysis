// Macro developed to reproduce the results for <p_T> vs N_ch for both 5.02 and 8.16 TeV pPb collision energies.
// Unpublished data by the CMS but exposed in the article [https://cds.cern.ch/record/2931093/files/HIN-25-001-pas.pdf].
// This work is part of the CMS Collaboration and uses CMS Preliminary Data.

/*
--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.

---------------Lucas Carvalho---------------*/

// Including function to fit points with same effective volume in order to extract the speed of sound squared.

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
bool hagedorn_fit = true; // Makes fit and extrapolation of CMS data with the Hagedorn function if true. If false, calculates over p_T > 0.3 GeV only.
bool canvas_plot = true; // If true, will produce a plot similar to the Fig. 1 of the CMS reference. If false, will produce only the .dat file.
int SAMPLE = 1e+8; // Number of generated entries with TF1::GetRandom() following the Hagedorn probability distribution function.
int n_centralities = 3; // Number of centrality classes defined for the dataset.
double low_eta = -1.5; // Low edge of pseudorapidity window.
double high_eta = 1.5; // High edge of pseudorapidity window. The p_T distribution will be integrated over (low_eta, high_eta).
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-6; // GeV
const std::string fit_type = "cms";
TString plot_extension = ".pdf";
TString output_name = "mean_pT_vs_Nch_CMS"; // Output name of Fig. 1 reproduction. Plots <pT>(Nch) for n_centralities for both collision energies.
TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

// ##############################################################################
// ##############################################################################


// --- Vectors to store 5TeV data
std::vector<double> mean_pT_data_5TeV;
std::vector<double> mean_pT_errors_5TeV;
std::vector<double> N_ch_values_5TeV;
std::vector<double> N_ch_errors_5TeV;

// --- Vectors to store 8TeV data
std::vector<double> mean_pT_data_8TeV;
std::vector<double> mean_pT_errors_8TeV;
std::vector<double> N_ch_values_8TeV;
std::vector<double> N_ch_errors_8TeV;

// --- Vectors to store cs results
std::vector<double> cs_results;
std::vector<double> cs_errors;


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
    "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root",
    "dataFile_5TeV",
    "5TeV",
    "_5TeV"
};

// --- struct to 8TeV dataset
const DataFile dataFile_8TeV = {
    "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root",
    "dataFile_8TeV",
    "8TeV",
    "_8TeV"
};

// Function to get and return n_events from 'hfSumEtPb' located at QA_histograms.
double get_n_events(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
// Function to calculate mean pT and Nch from a given TH1 histogram.
void calculate_pT_vs_Nch(TH1D* pt_hist, const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
// Function to write mean the results to a datafile.
void printResults_to_datafile();
// Function to plot and save as a figure the results for p_T vs N_ch.
void plot_pT_vs_Nch();
// Function to get and return the projected p_T distribution from the TH3 histogram located at Analysis_histograms, for the pseudorapidity window and centrality class.
TH1D* get_pT_TH1Dhistogram(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
// Function to make the Hagedorn fit and return extrapolated histogram for the p_T > 0 GeV region.
void Hagedorn_extrapolation(TH1D* hist_pT, const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
// Function to keep log of the macro's run.
void macro_log(const std::string& message);
// Function to draw CMS LaTeX header.
void drawCMSHeader(const char* extraText = "#it{Work in Progress}", double x = 0.12, double y = 0.93);
void plot_fig_2();

// --- main() ---
void add_cs_function(){
    gROOT->SetBatch(kTRUE); // This tells ROOT to run in batch mode, i.e. no GUI or pop-ups.

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-8);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);

	// Track cpu efficiency.
    TStopwatch timer;
    timer.Start();

    for(int i = 0; i < n_centralities; ++i){

    auto [low5, high5] = arr_5[i];
    auto [low8, high8] = arr_8[i];

    TH1D *hist = nullptr;
    if(hagedorn_fit){
        // Calculates with hagedorn fitting and p_T extrapolation for p_T > 0 GeV region. Stores the <pT>(Nch) data on global vectors.
        Hagedorn_extrapolation(get_pT_TH1Dhistogram(dataFile_5TeV.path, low5, high5), dataFile_5TeV.path, low5, high5);
        Hagedorn_extrapolation(get_pT_TH1Dhistogram(dataFile_8TeV.path, low8, high8), dataFile_8TeV.path, low8, high8);
    }

    else{
        // Calculates only with detector data for p_T > 0.3 GeV region.
        hist = get_pT_TH1Dhistogram(dataFile_5TeV.path, low5, high5);
    	hist = get_pT_TH1Dhistogram(dataFile_8TeV.path, low8, high8);
    }

    }// Ending the for loop.

    // The functions below only deal with the data vectors.
    if(canvas_plot){
        plot_pT_vs_Nch();
        plot_fig_2();
    }

    printResults_to_datafile(); // This function always needs to be called.

    timer.Stop();
    std::cout << "-> Job finished in "
              << timer.RealTime() << " seconds (wall time), "
              << timer.CpuTime()  << " seconds (CPU time).\n\n";
}


// --- Function definitions ---

void plot_pT_vs_Nch(){
    // Make graphs with TGraphErrors(npoints, x = N_ch, y = <pT>, N_ch errors, <pT> errors).
    // Since TGE doesn't read vectors we use .data() that returns a raw pointer (double*) to the first element of the vector.
    // e.g. mean_pT_data_8TeV.data() -> pointer to the first element (the 30-80% centrality class <pT>).
    TGraphErrors *g5 = new TGraphErrors(n_centralities, N_ch_values_5TeV.data(), mean_pT_data_5TeV.data(), N_ch_errors_5TeV.data(), mean_pT_errors_5TeV.data());
    TGraphErrors *g8 = new TGraphErrors(n_centralities, N_ch_values_8TeV.data(), mean_pT_data_8TeV.data(), N_ch_errors_8TeV.data(), mean_pT_errors_8TeV.data());
    // N_ch errors is assigned nullptr while i don't know if we need to estimate N_ch errors.

    TCanvas *c = new TCanvas("c","canvas", 800, 600);
    // Margins
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Grid & ticks
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis
    c->SetGrid();

    // Background
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);

    // Thicker border/frame
    c->SetFrameLineWidth(2);

    // General settings.
    g8->SetTitle("");
    g8->GetXaxis()->CenterTitle(true);
    g8->GetYaxis()->CenterTitle(true);
    g8->GetXaxis()->SetTitleOffset(1.0);
    g8->GetYaxis()->SetTitleOffset(1.1);
    g8->GetXaxis()->SetTitleFont(42);
    g8->GetYaxis()->SetTitleFont(42);
    g8->GetXaxis()->SetLabelFont(42);
    g8->GetYaxis()->SetLabelFont(42);
    g8->GetXaxis()->SetTitleSize(0.05);
    g8->GetYaxis()->SetTitleSize(0.044);
    g8->GetXaxis()->SetLabelSize(0.036);
    g8->GetYaxis()->SetLabelSize(0.036);

    // 5TeV data settings.
    g5->SetMarkerStyle(20);
    g5->SetMarkerSize(1.2);
    g5->SetMarkerColor(kGreen + 1);

    // 8TeV data settings.
    g8->SetMarkerStyle(20);
    g8->SetMarkerSize(1.2);
    g8->SetMarkerColor(kBlue + 1);

    // Draw on canvas c
    g8->GetXaxis()->SetTitle("N_{ch}");
    g8->GetYaxis()->SetTitle("#LTp_{T}#GT [GeV]");
    g8->GetYaxis()->SetRangeUser(0.48, 0.85);
    g8->GetXaxis()->SetRangeUser(36, 200);
    g8->Draw("APE");
    g5->Draw("PE SAME");

    // Legend
    auto leg = new TLegend(0.24,0.28,0.54,0.50);
    leg->AddEntry(g8,"8.16 TeV","p");
    leg->AddEntry(g5,"5.02 TeV","p");

    // Speed of sound fits
    for(int i = 0; i < n_centralities; ++i){
        double x_sub[2] = {N_ch_values_5TeV[i], N_ch_values_8TeV[i]};
        double y_sub[2] = {mean_pT_data_5TeV[i], mean_pT_data_8TeV[i]};
        double x_err[2] = {N_ch_errors_5TeV[i], N_ch_errors_8TeV[i]};
        double y_err[2] = {mean_pT_errors_5TeV[i], mean_pT_errors_8TeV[i]};
        TGraphErrors *sub_gr = new TGraphErrors(2, x_sub, y_sub, x_err, y_err);

        TF1 *cs_fit = new TF1(Form("fit_%d",i),"[0]*pow(x,[1])", N_ch_values_5TeV[i] - 1.0, N_ch_values_8TeV[i] + 1.0);
        cs_fit->SetParameters(0.25, 0.22);

        cs_fit->SetLineWidth(4);
        cs_fit->SetLineStyle(7);
        cs_fit->SetLineColor(kBlack);
        // Fit section
        sub_gr->Fit(cs_fit,"NO R EX0 Q","",N_ch_values_5TeV[i] - 1.0, N_ch_values_8TeV[i] + 1.0);
        sub_gr->Fit(cs_fit,"NO R EX0 Q","",N_ch_values_5TeV[i] - 1.0, N_ch_values_8TeV[i] + 1.0);
        TFitResultPtr cs_fit_result = sub_gr->Fit(cs_fit,"NO R EX0 M S","",N_ch_values_5TeV[i] - 1.0, N_ch_values_8TeV[i] + 1.0);
        
        cs_fit->Draw("SAME");
        if(i==0) leg->AddEntry(cs_fit,"Fit","l");

        cs_results.push_back(cs_fit->GetParameter(1));
        cs_errors.push_back(cs_fit->GetParError(1));
    }

    // Legend settings
    leg->SetHeader("Data"); // theres a flag "C" that centers the text.
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->Draw();

    // Draw CMS header
    drawCMSHeader();

    c->Modified();
    c->Update();

    // Path to final plot and save output plot.
    TString full_path = base_output_path + output_name + plot_extension;
    c->SaveAs(full_path);
}

TH1D* get_pT_TH1Dhistogram(const std::string& filename, double EHFmin, double EHFmax){

    // Return TH1 histogram of p_T track distribution on the range p_T > 0.3 [GeV] projected from the TH3 histogram located at 'Analysis_histograms',
    // for selected centrality class and pseudorapidity window.
    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        macro_log("Error: from get_pT_TH1Dhistogram function: Could not open file!");
        return nullptr;
    }

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    if (!dir) {
        macro_log("Error: from get_proj_hist function: Directory Analysis_histograms not found!");
        file->Close();
        return nullptr;
    }
    
    // Load TH3 histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    if (!hist_HFSumPb_vs_pt_eta){
        macro_log("Error: from get_proj_hist function: Error while creating TH3 histogram.");
        file->Close();
        return nullptr;
    }
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
    if (!hist_HFSumPb_vs_pt){
    	macro_log("Error: from get_proj_hist function: Error while creating TH2 histogram.");
        return nullptr;
    }

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
    if (!hist_pT){
        macro_log("Error: from get_proj_hist function: Error while creating TH1 histogram.");
        return nullptr;
    }
    hist_pT->SetDirectory(0);

    if(filename == dataFile_5TeV.path){
        std::cout << "\n\nRunning on " << dataFile_5TeV.label << " data set."
        << " Centrality class: [" << EHFmin << "," << EHFmax << "] GeV." << std::endl; 
    }
    else if(filename == dataFile_8TeV.path){
        std::cout << "\n\nRunning on " << dataFile_8TeV.label << " data set."
        << " Centrality class: [" << EHFmin << "," << EHFmax << "] GeV." << std::endl;
    }
    
    printf("\n\n-> From get_pT_TH1Dhistogram: hist_pT norm = %.3e \n", hist_pT->Integral());

    if(hagedorn_fit) return hist_pT; 

    else{
        calculate_pT_vs_Nch(hist_pT, filename, EHFmin, EHFmax);
        return nullptr;
    }

}

void Hagedorn_extrapolation(TH1D* hist_pT, const std::string& filename, double EHFmin, double EHFmax){

	TH1D *original_hist = nullptr;
    original_hist = (TH1D*)hist_pT->Clone("original_hist");
    if(!original_hist){
        macro_log("From make_hagedorn_extrapolation: error while loading hist_pT.");
        return;
    }
    original_hist->SetDirectory(0);
    printf("\n-> From make_hagedorn_extrapolation: original_hist norm = %.3e \n", original_hist->Integral());

    // Get norm only on the target range: 0.3 to 1.5 GeV.
    int bin_min = original_hist->GetXaxis()->FindBin(lower_pt_forFit + delta);
    int bin_max = original_hist->GetXaxis()->FindBin(upper_pt_forFit - delta);
    double original_hist_norm = original_hist->Integral(bin_min, bin_max); // Integral over [0.3,1.5].
    printf("\n-> From make_hagedorn_extrapolation: original_hist norm over [0.3, 1.5] GeV = %.3e \n\n\n", original_hist_norm);

    // #############################################
    // --------THE FIT SECTION STARTS HERE----------
    // #############################################

    // Hagedorn TF1. Function declaration section.
    TF1* pT_fit;
    int cc_low = static_cast<int>(EHFmin);
    int cc_up = static_cast<int>(EHFmax);
    
    // Original.
	pT_fit = new TF1(Form("ptfit_%d_%d",cc_low,cc_up),"[0]*x*pow(1.+1./sqrt(1.-[1]*[1])*(sqrt(x*x+[4]*[4])-x*[1])/[3]/[2],-[3])",0.,upper_pt_forFit);
    pT_fit->SetParameters(7500000000.,0.3,0.1,6.,0.14);//We used these values for initialization        
    pT_fit->FixParameter(4,0.13957);//pion mass
    
    if(filename == dataFile_5TeV.path){
        pT_fit->FixParameter(1, 0.4034);// related to radial flow velocity - pPb 5TeV
    }
    else if(filename == dataFile_8TeV.path){
        pT_fit->FixParameter(1, 0.5010);//related to radial flow velocity - pPb 8TeV
    }

    pT_fit->SetParLimits(2,0.,0.5);//kinetic freeze-out temperature in GeV
    pT_fit->SetParLimits(3,6.,9.);//n - free parameter no physical meaning


    // User-defined Hagedorn function fit.
    original_hist->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    original_hist->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    TFitResultPtr fitResult = original_hist->Fit(pT_fit,"NO R EX0 M S","",lower_pt_forFit,upper_pt_forFit);
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    double pValue = TMath::Prob(chi2, ndf);
    std::cout<<"chi2 : "<<chi2<<"; ndf : "<<ndf<<"; pValue : "<<pValue<<std::endl;

    // #############################################
    // --------THE FIT SECTION ENDS HERE----------
    // #############################################


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
    for(int i = 1; i <= extrapolated_hist->GetNbinsX(); i++){
        if (i <= 3){
            content = fit_hist->GetBinContent(i);
            error = fit_hist->GetBinError(i);
            extrapolated_hist->SetBinContent(i, content);
            extrapolated_hist->SetBinError(i, error);
        }
        else{
            content = original_hist->GetBinContent(i);
            error = original_hist->GetBinError(i);
            extrapolated_hist->SetBinContent(i, content);
            extrapolated_hist->SetBinError(i, error);
        }
    }

    printf("\n-> Extrapolated histogram mean = %f \n", extrapolated_hist->GetMean());
    printf("\n-> Original histogram mean = %f \n", original_hist->GetMean());
extrapolated_hist->ResetStats();
    printf("\n-> Extrapolated histogram mean = %f \n", extrapolated_hist->GetMean());

/*printf("\n\n----------Histogram diagnosis-----------");
printf("\nThe original_hist has %d bins \n", original_hist->GetNbinsX());
printf("\nThe extrapolated_hist has %d bins \n", extrapolated_hist->GetNbinsX());

printf("\n--- fit_hist: First three bins ---\n");
printf("Bin\tLow edge\tHigh edge\tContent\t\tError\n");

for (int i = 1; i <= 3; i++) {
    double low_edge  = fit_hist->GetXaxis()->GetBinLowEdge(i);
    double high_edge = fit_hist->GetXaxis()->GetBinUpEdge(i);
    double content   = fit_hist->GetBinContent(i);
    double error     = fit_hist->GetBinError(i);

    printf("%2d\t%.3f\t\t%.3f\t\t%.3e ± %.1e\n",
           i, low_edge, high_edge, content, error);
}

printf("\n--- Comparing original_hist and extrapolated_hist ---\n");
printf("Bin\tLow edge\tHigh edge\tOriginal\tExtrapolated\n");

int nbins = extrapolated_hist->GetNbinsX();
for (int i = 1; i <= nbins; i++) {
    double low_edge  = extrapolated_hist->GetXaxis()->GetBinLowEdge(i);
    double high_edge = extrapolated_hist->GetXaxis()->GetBinUpEdge(i);
    double orig_content = original_hist->GetBinContent(i);
    double extr_content = extrapolated_hist->GetBinContent(i);

    printf("%2d\t%.5f\t\t%.5f\t\t%.5e\t%.5e\n",
           i, low_edge, high_edge, orig_content, extr_content);
}
printf("--- End of comparison ---\n\n");*/

    calculate_pT_vs_Nch(extrapolated_hist, filename, EHFmin, EHFmax);
}


void printResults_to_datafile(){

    TString full_path5 = base_output_path + output_name + "_5TeV.dat";
    TString full_path8 = base_output_path + output_name + "_8TeV.dat";

    // Write 5 TeV data
    std::ofstream fout5(full_path5);
    if (!fout5.is_open()){
        std::cerr << "Error: could not open file " << full_path5 << std::endl;
        return;
    }
    fout5 << "# N_ch   mean_pT   mean_pT_error\n";

    for(int i = 0; i < N_ch_values_5TeV.size(); i++){
        fout5 << N_ch_values_5TeV[i] << "  "
              << N_ch_errors_5TeV[i] << "  "
              << mean_pT_data_5TeV[i] << "  "
              << mean_pT_errors_5TeV[i] << "\n";
    }
    fout5.close();

    // ---- Write 8 TeV data ----
    std::ofstream fout8(full_path8);
    if(!fout8.is_open()){
        std::cerr << "Error: could not open file " << full_path8 << std::endl;
        return;
    }
    fout8 << "# N_ch   mean_pT   mean_pT_error\n";
    
    for(int i = 0; i < N_ch_values_8TeV.size(); i++){
        fout8 << N_ch_values_8TeV[i] << "  "
              << N_ch_errors_8TeV[i] << "  "
              << mean_pT_data_8TeV[i] << "  "
              << mean_pT_errors_8TeV[i] << "\n";
    }
    fout8.close();
    //---

    std::cout << "\n\n-> Data written to " 
    << output_name << "_5TeV.dat" << " and " 
    << output_name << "_8TeV.dat" << std::endl;
}

void calculate_pT_vs_Nch(TH1D* pt_hist, const std::string& filename, double EHFmin, double EHFmax){

	// Get <p_T> and <p_T> error.
    double mean_pT = pt_hist->GetMean();
    double mean_pT_error = pt_hist->GetMeanError();

    // Sum number of tracks n_tracks as bin_contents.
    double n_tracks_error = 0.0;
    double n_tracks = pt_hist->IntegralAndError(1, pt_hist->GetNbinsX(), n_tracks_error);

    // Get n_events from QA_histograms respective centrality class and calculate N_ch.
    double n_events = get_n_events(filename, EHFmin, EHFmax);
    double N_ch = (n_tracks/n_events);
    double N_ch_error = (n_tracks_error/n_events);


    if(filename == dataFile_5TeV.path){
        mean_pT_data_5TeV.push_back(mean_pT);
        mean_pT_errors_5TeV.push_back(mean_pT_error);
        N_ch_values_5TeV.push_back(N_ch);
        N_ch_errors_5TeV.push_back(N_ch_error);
    }
    else if(filename == dataFile_8TeV.path){
        mean_pT_data_8TeV.push_back(mean_pT);
        mean_pT_errors_8TeV.push_back(mean_pT_error);
        N_ch_values_8TeV.push_back(N_ch);
        N_ch_errors_8TeV.push_back(N_ch_error);
    }

    printf("--- Results --- \n");
    // We need a way to tell this function which datafile and centrality class it's reading.
    printf("-> <p_T> = %f \n", mean_pT);
    printf("-> <p_T> error = %f \n", mean_pT_error);
    printf("-> N_ch = %f \n", N_ch);
    printf("-> N_ch error = %f \n", N_ch_error);
    printf("--------------- \n\n");
    printf("--------------- \n\n");

}

double get_n_events(const std::string& filename, double EHFmin, double EHFmax){
    std::string this_function = __func__;
    TFile *f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
    	macro_log("From " + this_function + ": error while opening the file " + filename);
    	return -1;
    }

    TDirectory *dir = (TDirectory*)f->Get("QA_histograms");
    if (!dir) {
    	macro_log("From " + this_function + ": error while opening 'QA_histograms' directory.");
        f->Close();
        return -1;
    }

    TH1D *hf_hist = (TH1D*)dir->Get("hfSumEtPb");
    if (!hf_hist) {
    	macro_log("From " + this_function + ": histogram 'hfSumEtPb' not found!");
        f->Close();
        return -1;
    }
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

void macro_log(const std::string& message){
    // Set output to macro_log file.
    std::ofstream _log_("mean_pT_vs_Nch_CMS.log", std::ios::app);

    if(_log_.is_open()) {
        _log_ << message << std::endl;
    }
    else {
        std::cerr << "From " << __func__ << ": could not open 'mean_pT_vs_Nch_CMS.log'" << std::endl;
    }
    
}

void drawCMSHeader(const char* extraText, double x, double y){

TLatex latex;
latex.SetNDC();              // use normalized coordinates
latex.SetTextSize(0.04);     // text size
latex.SetTextFont(42);       // Helvetica-like
latex.SetTextAlign(11);      // left-aligned, top

TString cmsText = "#bf{CMS}";
if (extraText && strlen(extraText) > 0) {
cmsText += " ";
cmsText += extraText;
}

latex.DrawLatex(x, y, cmsText);

TLatex latex2;
latex2.SetNDC();
latex2.SetTextSize(0.035);
latex2.SetTextFont(42);
latex2.SetTextAlign(11);  // right-aligned
latex2.DrawLatex(0.16, 0.84, "p_{T} > 0 GeV, |#eta| < 1.5");

TLatex latex3;
latex3.SetNDC();
latex3.SetTextSize(0.035);
latex3.SetTextFont(42);
latex3.SetTextAlign(31);
latex3.DrawLatex(0.90, 0.84, "pPb (186.0 nb^{#minus1}) 8.16 TeV");

TLatex latex4;
latex4.SetNDC();
latex4.SetTextSize(0.035);
latex4.SetTextFont(42);
latex4.SetTextAlign(31);
latex4.DrawLatex(0.90, 0.79, "pPb (0.509 nb^{#minus1}) 5.02 TeV");
}

void plot_fig_2(){

    std::vector<double> T_eff_5TeV;
    std::vector<double> T_eff_8TeV;
    std::vector<double> T_eff;
    for (double x : mean_pT_data_5TeV) T_eff_5TeV.push_back(x*1000./3.);
    for (double x : mean_pT_data_8TeV) T_eff_8TeV.push_back(x*1000./3.);

    for(int j = 0; j < n_centralities; ++j){
        T_eff.push_back((T_eff_5TeV[j] + T_eff_8TeV[j])/2.0);
    }

    TCanvas *c = new TCanvas();
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Grid & ticks
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis
    c->SetGrid();

    // Background
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);

    // Thicker border/frame
    c->SetFrameLineWidth(2);

    TGraphErrors *fig2 = new TGraphErrors(3, T_eff.data(), cs_results.data(), nullptr, nullptr);

    fig2->SetMarkerStyle(20);
    fig2->SetMarkerColor(kBlack);
    fig2->GetXaxis()->SetLimits(130, 360);
    fig2->GetYaxis()->SetRangeUser(0., 0.4);

    fig2->SetTitle("");
    fig2->GetXaxis()->CenterTitle(true);
    fig2->GetYaxis()->CenterTitle(true);
    fig2->GetXaxis()->SetTitleOffset(1.0);
    fig2->GetYaxis()->SetTitleOffset(1.1);
    fig2->GetXaxis()->SetTitleFont(42);
    fig2->GetYaxis()->SetTitleFont(42);
    fig2->GetXaxis()->SetLabelFont(42);
    fig2->GetYaxis()->SetLabelFont(42);
    fig2->GetXaxis()->SetTitleSize(0.05);
    fig2->GetYaxis()->SetTitleSize(0.044);
    fig2->GetXaxis()->SetLabelSize(0.036);
    fig2->GetYaxis()->SetLabelSize(0.036);

    fig2->GetXaxis()->SetTitle("T_{eff} = #LTp_{T}#GT / 3 [MeV]");
    fig2->GetYaxis()->SetTitle("dln #LTp_{T}#GT / dln N_{ch}");
    fig2->Draw("APE");

TLatex latex;
latex.SetNDC();              // use normalized coordinates
latex.SetTextSize(0.04);     // text size
latex.SetTextFont(42);       // Helvetica-like
latex.SetTextAlign(11);      // left-aligned, top

TString cmsText = "#bf{CMS} Work in Progress";
latex.DrawLatex(0.12, 0.93, cmsText);

TLatex latex2;
latex2.SetNDC();
latex2.SetTextSize(0.035);
latex2.SetTextFont(42);
latex2.SetTextAlign(11);  // right-aligned
latex2.DrawLatex(0.16, 0.84, "p_{T} > 0 GeV, |#eta| < 1.5");

TLatex latex3;
latex3.SetNDC();
latex3.SetTextSize(0.035);
latex3.SetTextFont(42);
latex3.SetTextAlign(31);
latex3.DrawLatex(0.90, 0.84, "pPb (186.0 nb^{#minus1}) 8.16 TeV");

TLatex latex4;
latex4.SetNDC();
latex4.SetTextSize(0.035);
latex4.SetTextFont(42);
latex4.SetTextAlign(31);
latex4.DrawLatex(0.90, 0.79, "pPb (0.509 nb^{#minus1}) 5.02 TeV");

    c->Modified();
    c->Update();
    c->SaveAs("../../../../mnt/c/Users/lucas/Documents/fig2.pdf");
}