// Macro to reproduce the results for <p_T>(N_ch) for both 5.02 and 8.16 TeV pPb collision energies.
// The CMS preliminary results can be found as a public note in [https://cds.cern.ch/record/2931093/files/HIN-25-001-pas.pdf].
// This work is part of the CMS Collaboration and uses CMS Preliminary Data.

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
bool hagedorn_fit = true; // Makes fit and extrapolation of CMS data with the Hagedorn function if true. If false, calculates over p_T > 0.3 GeV only.
bool canvas_plot = true; // If true, will produce a plot similar to the Fig. 1 of the CMS reference. If false, will produce only the .dat file.
int SAMPLE = 1e+8; // Number of generated entries with TF1::GetRandom() following the Hagedorn probability distribution function.
int n_centralities = 3; // Number of centrality classes defined for the dataset.
double low_eta = -2.4; // Low edge of pseudorapidity window.
double high_eta = 2.4; // High edge of pseudorapidity window. The p_T distribution will be integrated over (low_eta, high_eta).
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-6; // GeV
const std::string fit_type = "cms"; // Useful flag if one wants to try other functions.
TString plot_extension = ".pdf"; // .png to test .pdf to final result.
TString output1_name = "mean_pT_vs_Nch"; // Output name of Fig. 1 reproduction. Plots <pT>(Nch) for n_centralities for both collision energies.
TString output2_name = "c_s2_vs_Teff"; // Output name of Fig. 2 reproduction. Plots c_s^2(T_eff) for a total of n_centralities.
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

// --- Vectors to store c_s and T_eff data
std::vector<double> cs_results;
std::vector<double> cs_errors;
std::vector<double> T_eff_5TeV;
std::vector<double> T_eff_8TeV;
std::vector<double> T_eff;

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
// Function to get and return n_events from 'hfSumEtPb' located at 'QA_histograms'.
double get_n_events(const DataStruct& dataFile, double EHFmin = 0.0, double EHFmax = 250.0);
// Function to get and return TH1D pT histogram from 'hist_HFSumPb_vs_pt_eta' located at 'Analysis_histograms'.
TH1D* get_pT_histogram(const DataStruct& dataFile, double EHFmin = 0.0, double EHFmax = 250.0);
// Function to draw latex texts on the plots.
void draw_CMS_Header(TString latex_text = "#bf{CMS} #it{Work in Progress}", double x = 0.12, double y = 0.93, double text_size = 0.04, int align = 11);
// Function to plot Fig. 2 of reference.
void plot_cs2_vs_Teff();
// Function to plot Fig. 1 of reference.
void plot_pT_vs_Nch();
// Function to print results to datafile.
void printResults_to_datafile();
// Function to calculate mean p_T as function of N_ch.
void calculate_pT_vs_Nch(TH1D* pt_hist, const std::string& filename, double EHFmin, double EHFmax);
// Function to fit the pT histogram with the Hagedorn function.
void make_hagedorn_fit(TH1D* hist_pT, const DataStruct& dataFile, double EHFmin, double EHFmax);


// --- main() ---
void aux(){

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

    TH1D *hist = nullptr;
    if(hagedorn_fit){
        // Calculates with hagedorn fitting and p_T extrapolation for p_T > 0 GeV region. Stores the <pT>(Nch) data on global vectors.
        make_hagedorn_fit(get_pT_histogram(dataFile_5TeV, low5, high5), dataFile_5TeV, low5, high5);
        make_hagedorn_fit(get_pT_histogram(dataFile_8TeV, low8, high8), dataFile_8TeV, low8, high8);
    }

    else{
        // Calculates only with detector data for p_T > 0.3 GeV region.
        hist = get_pT_histogram(dataFile_5TeV, low5, high5);
        hist = get_pT_histogram(dataFile_8TeV, low8, high8);
    }

    }// Ending the for loop.

    // The functions below only deal with the data vectors.
    if(canvas_plot){
        plot_pT_vs_Nch();
        plot_cs2_vs_Teff();
    }
    printResults_to_datafile(); // This function always needs to be called.


    timer.Stop();
    std::cout << "-> Job finished in "
              << timer.RealTime() << " seconds (wall time), "
              << timer.CpuTime()  << " seconds (CPU time).\n\n";

}

// --- Function definitions ---
double get_n_events(const std::string& filename, double EHFmin, double EHFmax){

    TFile *f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "From " << __func__ << ": error while opening target file." << std::endl;
        return -1;
    }

    TDirectory *dir = (TDirectory*)f->Get("QA_histograms");
    if (!dir) {
        std::cerr << "From " << __func__ << ": error while opening 'QA_histograms' directory." << std::endl;
        f->Close();
        return -1;
    }

    TH1D *hf_hist = (TH1D*)dir->Get("hfSumEtPb"); // Using only Pb-going side of HF calorimeter to estimate number of events.
    if (!hf_hist) {
        std::cerr << "From " << __func__ << ": histogram 'hfSumEtPb' not found." << std::endl;
        f->Close();
        return -1;
    }
    hf_hist->SetDirectory(0);
    f->Close();

    int bin_min = hf_hist->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hf_hist->GetXaxis()->FindBin(EHFmax - delta);

    // Check interval of integration over transversal E_HF energy.
    double lowEdge   = hf_hist->GetXaxis()->GetBinLowEdge(bin_min);
    double highEdge  = hf_hist->GetXaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From %s: Integrating 'hfSumEtPb' on the interval [%.1f,%.1f]GeV \n\n", __func__, lowEdge, highEdge);

    // Calculates the number of total events n_events, detected by the Pb side HF, in the specified centrality bin/class.
    double n_events = hf_hist->Integral(bin_min, bin_max);
    if(n_events > 0) return n_events;
    else {
        std::cerr << "From " << __func__ << ": calculated n_events equal to or less than zero." << std::endl;
        return -1;
    }
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

TH1D* get_pT_histogram(const DataStruct& dataFile, double EHFmin, double EHFmax){

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
    printf("\n-> From %s: hist_pT norm = %.3e \n", __func__, hist_pT->Integral());

    if(hagedorn_fit) return hist_pT; 

    else{
        calculate_pT_vs_Nch(hist_pT, dataFile.path, EHFmin, EHFmax);
        return nullptr;
    }

}

void make_hagedorn_fit(TH1D* hist_pT, const DataStruct& dataFile, double EHFmin, double EHFmax){

    if(!hist_pT){
        std::cerr << "From " << __func__ << ": error while loading TH1D histogram." << std::endl;
        return;
    }
    printf("\n-> From %s: hist_pT norm = %.3e \n\n", __func__, hist_pT->Integral());

    // Get norm only on the target range: 0.3 to 1.5 GeV.
    int bin_min = hist_pT->GetXaxis()->FindBin(lower_pt_forFit + delta);
    int bin_max = hist_pT->GetXaxis()->FindBin(upper_pt_forFit - delta);
    double hist_pT_norm = hist_pT->Integral(bin_min, bin_max); // Integral over [0.3,1.5].

    // Fit section.
    TF1* pT_fit;
    int cc_low = static_cast<int>(EHFmin);
    int cc_up = static_cast<int>(EHFmax);
    
    // User-defined Hagedorn TF1.
    pT_fit = new TF1(Form("ptfit_%d_%d",cc_low,cc_up),"[0]*x*pow(1.+1./sqrt(1.-[1]*[1])*(sqrt(x*x+[4]*[4])-x*[1])/[3]/[2],-[3])",0.,upper_pt_forFit);
    pT_fit->SetParameters(7500000000.,0.3,0.1,6.,0.14);//We used these values for initialization        
    pT_fit->FixParameter(4,0.13957);//pion mass
    pT_fit->SetParLimits(2,0.,0.5);//kinetic freeze-out temperature in GeV    
    pT_fit->SetParLimits(3,4.,9.);//n - free parameter no physical meaning
    if(dataFile.label == "5TeV"){
        pT_fit->FixParameter(1, 0.4034);// related to radial flow velocity - pPb 5TeV
    }
    else if(dataFile.label == "8TeV"){
        pT_fit->FixParameter(1, 0.5010);//related to radial flow velocity - pPb 8TeV
    }

    // Make fit on hist_pT.
    hist_pT->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    hist_pT->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    TFitResultPtr fitResult = hist_pT->Fit(pT_fit,"NO R EX0 M S","",lower_pt_forFit,upper_pt_forFit);
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    double pValue = TMath::Prob(chi2, ndf);
    std::cout << "\n -----> chi2 : " << chi2 << "; ndf : " << ndf << "; pValue : " << pValue << std::endl;
    // End of fit section.

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
    double scale_factor = (hist_pT_norm/SAMPLE);
    printf("\n-> Scale factor = %.3f \n\n", scale_factor);
    fit_hist->Scale(scale_factor);

    fit_hist_norm = fit_hist->Integral(bin_min, bin_max);
    printf("\n-> Fit histogram complete. Final fit histogram norm over [0.3, 1.5]GeV = %.3e \n", fit_hist_norm);
    printf("\n-> Compare with hist_pT norm over [0.3, 1.5]GeV = %.3e \n", hist_pT_norm);

    // Fill simulated tracks into histogram extrapolated_hist.
    //TH1D *extrapolated_hist = (TH1D*)hist_pT->Clone("extrapolated_hist");    
    double pt_kinetic_cut = 3.0; //GeV
    int pt_kinetic_cut_bin = hist_pT->GetXaxis()->FindBin(pt_kinetic_cut - delta);
    printf("\n-> Kinetic bin = %d \n", pt_kinetic_cut_bin);
    TH1D *extrapolated_hist = new TH1D("extrapolated_hist", "Extrapolated Histogram", pt_kinetic_cut_bin, 0.,pt_kinetic_cut);
    extrapolated_hist->SetDirectory(0);

    // Takes first three bin contents from fit_hist. The rest is taken from the real data from hist_pT.
    double content, error;
    // Before it was extrapolated_hist->GetNbinsX(). Now its pt_kinetic_cut_bin.
    for(int i = 1; i <= pt_kinetic_cut_bin; i++){
        if (i <= 3){
            content = fit_hist->GetBinContent(i);
            error = fit_hist->GetBinError(i);
            extrapolated_hist->SetBinContent(i, content);
            extrapolated_hist->SetBinError(i, error);
        }
        else{
            content = hist_pT->GetBinContent(i);
            error = hist_pT->GetBinError(i);
            extrapolated_hist->SetBinContent(i, content);
            extrapolated_hist->SetBinError(i, error);
        }
    }
    extrapolated_hist->ResetStats();

    calculate_pT_vs_Nch(extrapolated_hist, dataFile.path, EHFmin, EHFmax);
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

    printf("---- Results ---- \n");
    // We need a way to tell this function which datafile and centrality class it's reading.
    printf("-> <p_T> = %f \n", mean_pT);
    printf("-> <p_T> error = %f \n", mean_pT_error);
    printf("-> N_ch = %f \n", N_ch);
    printf("-> N_ch error = %f \n", N_ch_error);
    printf("----------------- \n");
    printf("----------------- \n\n");
}

void printResults_to_datafile(){

    TString full_path5 = base_output_path + output1_name + "_5TeV.dat";
    TString full_path8 = base_output_path + output1_name + "_8TeV.dat";

    // Write 5 TeV data
    std::ofstream fout5(full_path5);
    if (!fout5.is_open()){
        std::cerr << "Error: could not open file " << full_path5 << std::endl;
        return;
    }
    fout5 << "# N_ch    N_ch_error  mean_pT     mean_pT_error\n";

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
    fout8 << "# N_ch    N_ch_error  mean_pT     mean_pT_error\n";
    
    for(int i = 0; i < N_ch_values_8TeV.size(); i++){
        fout8 << N_ch_values_8TeV[i] << "  "
              << N_ch_errors_8TeV[i] << "  "
              << mean_pT_data_8TeV[i] << "  "
              << mean_pT_errors_8TeV[i] << "\n";
    }
    fout8.close();
    //---

    std::cout << "\n\n-> Data written to " 
    << output1_name << "_5TeV.dat" << " and " 
    << output1_name << "_8TeV.dat" << std::endl;
}

void plot_pT_vs_Nch(){
    // Make graphs with TGraphErrors(npoints, x = N_ch, y = <pT>, N_ch errors, <pT> errors).
    // Since TGE doesn't read vectors we use .data() that returns a raw pointer (double*) to the first element of the vector.
    // e.g. mean_pT_data_8TeV.data() -> pointer to the first element (the 30-80% centrality class <pT>).
    TGraphErrors *g5 = new TGraphErrors(n_centralities, N_ch_values_5TeV.data(), mean_pT_data_5TeV.data(), N_ch_errors_5TeV.data(), mean_pT_errors_5TeV.data());
    TGraphErrors *g8 = new TGraphErrors(n_centralities, N_ch_values_8TeV.data(), mean_pT_data_8TeV.data(), N_ch_errors_8TeV.data(), mean_pT_errors_8TeV.data());

    TCanvas *c = new TCanvas("c","canvas", 800, 600);
    // Margins
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Grid & ticks
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis

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
    //g8->GetYaxis()->SetRangeUser(0.48, 0.85);
    //g8->GetXaxis()->SetRangeUser(36, 200);
    g8->Draw("APE");
    g5->Draw("PE SAME");
    g8->GetXaxis()->SetLimits(36, 200);
    g8->GetYaxis()->SetRangeUser(0.48, 0.85);

    // Legend
    auto leg = new TLegend(0.26,0.26,0.56,0.48);
    leg->AddEntry(g8,"8.16 TeV","p");
    leg->AddEntry(g5,"5.02 TeV","p");

    // Speed of sound fits
    for(int i = 0; i < n_centralities; ++i){
        double x_sub[2] = {N_ch_values_5TeV[i], N_ch_values_8TeV[i]};
        double y_sub[2] = {mean_pT_data_5TeV[i], mean_pT_data_8TeV[i]};
        double x_err[2] = {N_ch_errors_5TeV[i], N_ch_errors_8TeV[i]};
        double y_err[2] = {mean_pT_errors_5TeV[i], mean_pT_errors_8TeV[i]};
        TGraphErrors *sub_gr = new TGraphErrors(2, x_sub, y_sub, x_err, y_err);

        TF1 *cs_fit = new TF1(Form("fit_%d",i),"[0]*pow(x,[1])", N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
        cs_fit->SetParameters(0.25, 0.22);

        cs_fit->SetLineWidth(4);
        cs_fit->SetLineStyle(7);
        cs_fit->SetLineColor(kBlack);
        // Fit section
        sub_gr->Fit(cs_fit,"NO R EX0 Q","",N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
        sub_gr->Fit(cs_fit,"NO R EX0 Q","",N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
        TFitResultPtr cs_fit_result = sub_gr->Fit(cs_fit,"NO R EX0 M S","",N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
        
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

    // Draw headers
    draw_CMS_Header();
    draw_CMS_Header("p_{T} > 0 GeV, |#eta| < 1.5", 0.16, 0.84, 0.038, 11);
    draw_CMS_Header("pPb (186.0 nb^{#minus1}) 8.16 TeV", 0.90, 0.84, 0.038, 31);
    draw_CMS_Header("pPb (0.509 nb^{#minus1}) 5.02 TeV", 0.90, 0.79, 0.038, 31);

    c->Modified();
    c->Update();

    // Path to final plot and save output plot.
    TString full_path = base_output_path + output1_name + plot_extension;
    c->SaveAs(full_path);
}

void plot_cs2_vs_Teff(){

    std::vector<double> T_eff_5TeV;
    std::vector<double> T_eff_8TeV;
    std::vector<double> T_eff;
    for (double x : mean_pT_data_5TeV) T_eff_5TeV.push_back(x*1000./3.);
    for (double x : mean_pT_data_8TeV) T_eff_8TeV.push_back(x*1000./3.);

    for(int j = 0; j < n_centralities; ++j){
        T_eff.push_back((T_eff_5TeV[j] + T_eff_8TeV[j])/2.0);
    }

    TCanvas *c = new TCanvas("c2","canvas2", 700, 600);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Grid & ticks
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis

    // Background
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);

    // Thicker border/frame
    c->SetFrameLineWidth(2);

    TGraphErrors *fig2 = new TGraphErrors(T_eff.size(), T_eff.data(), cs_results.data(), nullptr, nullptr);
    // Get Hijing and Trajectum results.
    TFile *trajectum_file = TFile::Open("../../cs2_Trajectum_FCALCent_eta1_2.root", "READ");
    TFile *hijing_file = TFile::Open("../../cs2_MC-MB-HIJING_HF4eta5_trks1p0_Fit0p0-2p0_BoostInvariant.root", "READ");

    TGraphErrors *trajectum_graph = (TGraphErrors*) trajectum_file->Get("Graph;1");
    TGraphErrors *hijing_graph = (TGraphErrors*) hijing_file->Get("Graph;2");
    trajectum_file->Close();
    hijing_file->Close();

    // Style settings
    fig2->SetMarkerStyle(20);
    fig2->SetMarkerSize(1.0);
    fig2->SetMarkerColor(kBlack);
    fig2->SetLineColor(kBlack);
    fig2->SetLineWidth(3);

    trajectum_graph->SetMarkerStyle(25);
    trajectum_graph->SetMarkerSize(1.0);
    trajectum_graph->SetMarkerColor(kGreen + 2);
    trajectum_graph->SetLineColor(kGreen + 2);
    trajectum_graph->SetLineWidth(3);
    trajectum_graph->SetTitle("");

    hijing_graph->SetMarkerStyle(25);
    hijing_graph->SetMarkerSize(1.0);
    hijing_graph->SetMarkerColor(kMagenta);
    hijing_graph->SetLineColor(kMagenta);
    hijing_graph->SetLineWidth(3);
    hijing_graph->SetTitle("");

    // Add non-interacting limit dashed line.
    TLine *line = new TLine(280, 0.33, 380, 0.33);
    line->SetLineColor(kGray + 2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);

    // Plot settings - only needed for fig2 (data).
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
    fig2->GetXaxis()->SetLimits(130, 380);
    fig2->GetYaxis()->SetRangeUser(0., 0.4);
    fig2->GetXaxis()->SetTitle("T_{eff} = #LTp_{T}#GT / 3 [MeV]");
    fig2->GetYaxis()->SetTitle("dln #LTp_{T}#GT / dln N_{ch}");
    fig2->Draw("APE");
    trajectum_graph->Draw("PE SAME");
    hijing_graph->Draw("PE SAME");
    line->Draw("SAME");

    auto leg = new TLegend(0.58,0.20,0.80,0.36);
    leg->AddEntry(fig2,"Data","lep");
    leg->AddEntry(trajectum_graph,"pPb Trajectum","lep");
    leg->AddEntry(hijing_graph,"pPb Hijing","lep");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->Draw();

    draw_CMS_Header();
    draw_CMS_Header("non-interacting limit", 0.9,0.73, 0.038,31);
    draw_CMS_Header("p_{T} < 3.0 GeV, |#eta| < 2.4", 0.16, 0.20, 0.038, 11);
    draw_CMS_Header("pPb (186.0 nb^{#minus1}) 8.16 TeV", 0.16, 0.84, 0.038, 11);
    draw_CMS_Header("pPb (0.509 nb^{#minus1}) 5.02 TeV", 0.16, 0.79, 0.038, 11);

    c->Modified();
    c->Update();

    // Path to final plot and save output plot.
    TString full_path = base_output_path + output2_name + plot_extension;
    c->SaveAs(full_path);
}