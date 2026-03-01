// Macro to apply kinematic cuts on the pT range.

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

// --- Function headers ---
double get_n_events(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
TH1D* get_pT_TH1Dhistogram(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0, double pt_low = 0.0, double pt_high = 150.);
TH1D* Hagedorn_extrapolation(TH1D* hist_pT, const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
void plot_graphs();
void draw_CMS_Header(TString latex_text = "#bf{CMS} #it{Work in Progress}", double x = 0.11, double y = 0.93, double text_size = 0.04, int align = 11);


// --- main () ---
void kinCuts_pT(){

    gROOT->SetBatch(kTRUE); // This tells ROOT to run in batch mode, i.e. no GUI or pop-ups.

	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-8);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);

	// Track cpu efficiency.
    TStopwatch timer;
    timer.Start();

    for(size_t k = 0; k < pT_cuts.size(); ++k){

    	double pt_low = pT_cuts[k].first;
        double pt_high = pT_cuts[k].second;
        printf("\n\n-> pT kinematic cut: [%.1f,%.1f] \n", pt_low, pt_high);
        
        // Deterministic seed for this eta window for reproducibility.
        gRandom->SetSeed(12345 + k);
        printf("Seed for this pT cut = %u \n", gRandom->GetSeed());

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

        bool hagedorn_fit;
        if(k == 1 || k == 2) hagedorn_fit = false;
        else hagedorn_fit = true;

        for(size_t i = 0; i < n_centralities; ++i){

        	auto [low5, high5] = arr_5[i];
            auto [low8, high8] = arr_8[i];

            TH1D *hist_5 = nullptr;
            TH1D *hist_8 = nullptr;
            double mean_pT, mean_pT_error, n_tracks, n_tracks_error;
            double n_events, N_ch, N_ch_error;

            if(!hagedorn_fit){
                hist_5 = get_pT_TH1Dhistogram(dataFile_5TeV.path, low5, high5, pt_low, pt_high);
            }
            else hist_5 = Hagedorn_extrapolation(get_pT_TH1Dhistogram(dataFile_5TeV.path, low5, high5, pt_low, pt_high), dataFile_5TeV.path, low5, high5);

            // Get <p_T> and <p_T> error.
            mean_pT = hist_5->GetMean();
            mean_pT_error = hist_5->GetMeanError();

            // Sum number of tracks n_tracks as bin_contents.
            n_tracks_error = 0.0;
            n_tracks = hist_5->IntegralAndError(1, hist_5->GetNbinsX(), n_tracks_error);

            // Get n_events from QA_histograms respective centrality class and calculate N_ch.
            n_events = get_n_events(dataFile_5TeV.path, low5, high5);
            N_ch = (n_tracks/n_events);
            N_ch_error = (n_tracks_error/n_events);

            mean_pT_data_5TeV.push_back(mean_pT);
            mean_pT_errors_5TeV.push_back(mean_pT_error);
            N_ch_values_5TeV.push_back(N_ch);
            N_ch_errors_5TeV.push_back(N_ch_error);

            printf("--- Results --- \n");
            printf("-> <p_T> = %f \n", mean_pT);
            printf("-> <p_T> error = %f \n", mean_pT_error);
            printf("-> N_ch = %f \n", N_ch);
            printf("-> N_ch error = %f \n", N_ch_error);
            printf("--------------- \n\n");
            printf("--------------- \n\n");
            
            if(!hagedorn_fit){
                hist_8 = get_pT_TH1Dhistogram(dataFile_8TeV.path, low8, high8, pt_low, pt_high);
            }
            else hist_8 = Hagedorn_extrapolation(get_pT_TH1Dhistogram(dataFile_8TeV.path, low8, high8, pt_low, pt_high), dataFile_8TeV.path, low8, high8);

            // Get <p_T> and <p_T> error.
            mean_pT = hist_8->GetMean();
            mean_pT_error = hist_8->GetMeanError();

            // Sum number of tracks n_tracks as bin_contents.
            n_tracks_error = 0.0;
            n_tracks = hist_8->IntegralAndError(1, hist_8->GetNbinsX(), n_tracks_error);

            // Get n_events from QA_histograms respective centrality class and calculate N_ch.
            n_events = get_n_events(dataFile_8TeV.path, low8, high8);
            N_ch = (n_tracks/n_events);
            N_ch_error = (n_tracks_error/n_events);

            mean_pT_data_8TeV.push_back(mean_pT);
            mean_pT_errors_8TeV.push_back(mean_pT_error);
            N_ch_values_8TeV.push_back(N_ch);
            N_ch_errors_8TeV.push_back(N_ch_error);

            printf("--- Results --- \n");
            printf("-> <p_T> = %f \n", mean_pT);
            printf("-> <p_T> error = %f \n", mean_pT_error);
            printf("-> N_ch = %f \n", N_ch);
            printf("-> N_ch error = %f \n", N_ch_error);
            printf("--------------- \n\n");
            printf("--------------- \n\n");
        }//Ending the i-loop.

        // --- Vectors to store c_s data
        std::vector<double> cs_results;
        std::vector<double> cs_errors;

        // --- Speed of sound fits ---
        for(int i = 0; i < n_centralities; ++i){
            double x_sub[2] = {N_ch_values_5TeV[i], N_ch_values_8TeV[i]};
            double y_sub[2] = {mean_pT_data_5TeV[i], mean_pT_data_8TeV[i]};
            double x_err[2] = {N_ch_errors_5TeV[i], N_ch_errors_8TeV[i]};
            double y_err[2] = {mean_pT_errors_5TeV[i], mean_pT_errors_8TeV[i]};
            TGraphErrors *sub_gr = new TGraphErrors(2, x_sub, y_sub, x_err, y_err);

            TF1 *cs_fit = new TF1(Form("fit_%d",i),"[0]*pow(x,[1])", N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
            cs_fit->SetParameters(0.25, 0.22);

            // Fit section
            sub_gr->Fit(cs_fit,"NO R EX0 Q","",N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
            sub_gr->Fit(cs_fit,"NO R EX0 Q","",N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
            TFitResultPtr cs_fit_result = sub_gr->Fit(cs_fit,"NO R EX0 M S","",N_ch_values_5TeV[i], N_ch_values_8TeV[i]);
        
            cs_results.push_back(cs_fit->GetParameter(1));
            cs_errors.push_back(cs_fit->GetParError(1));
        }
    
        // --- Vectors to store T_eff data
        std::vector<double> T_eff_5TeV;
        std::vector<double> T_eff_8TeV;
        std::vector<double> T_eff;
        for (double x : mean_pT_data_5TeV) T_eff_5TeV.push_back(x*1000./3.);
        for (double x : mean_pT_data_8TeV) T_eff_8TeV.push_back(x*1000./3.);

        for(int j = 0; j < n_centralities; ++j){
            T_eff.push_back((T_eff_5TeV[j] + T_eff_8TeV[j])/2.0);
        }

        // Store the data on a TGraphError.
        TGraphErrors *graph = new TGraphErrors(T_eff.size(), T_eff.data(), cs_results.data(), nullptr, nullptr); // Graph for this eta window.
        // Store the TGE on a vector of TGE's.
        graphs.push_back(graph);

    }//Ending the k-loop.

    // Make final plot.
    plot_graphs();

    timer.Stop();
    std::cout << "-> Job finished in "
              << timer.RealTime() << " seconds (wall time), "
              << timer.CpuTime()  << " seconds (CPU time).\n\n";	
}

// --- Function definitions ---
TH1D* get_pT_TH1Dhistogram(const std::string& filename, double EHFmin, double EHFmax, double pt_low, double pt_high){

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

    // Make kinematic cut on the specified pT range.
    TH1D* hist_pT_cut = (TH1D*)hist_pT->Clone("hist_pT_cut");
    hist_pT_cut->SetDirectory(0);

    for(int j = 1; j <= hist_pT_cut->GetNbinsX(); ++j) {

        double xlow  = hist_pT_cut->GetXaxis()->GetBinLowEdge(j);
        double xhigh = hist_pT_cut->GetXaxis()->GetBinUpEdge(j);

        if(xhigh > pt_high || xlow < pt_low) hist_pT_cut->SetBinContent(j, 0.);
    }

    // Print binning information for hist_pT_cut just to be sure.
    printf("hist_pT_cut has %d bins\n", hist_pT_cut->GetNbinsX());
    printf("X-axis range: [%.6f, %.6f]\n",
       hist_pT_cut->GetXaxis()->GetXmin(),
       hist_pT_cut->GetXaxis()->GetXmax());
    printf("Bin width: %.6f\n", hist_pT_cut->GetXaxis()->GetBinWidth(1));

    printf("---- Individual bin edges ----\n");
    for (int b = 1; b <= hist_pT_cut->GetNbinsX(); ++b) {
        double low  = hist_pT_cut->GetXaxis()->GetBinLowEdge(b);
        double high = hist_pT_cut->GetXaxis()->GetBinUpEdge(b);
        printf("Bin %2d: [%.6f , %.6f]\n", b, low, high);
    }
    printf("--------------------------------\n");

    return hist_pT_cut;
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
    //calculate_pT_vs_Nch(extrapolated_hist, filename, EHFmin, EHFmax);
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
    
    int colors[6] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+7};
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerSize(1.);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetTitle("");
    }

    // Settings only needed for first plot.
    graphs[0]->GetXaxis()->CenterTitle(true);
    graphs[0]->GetYaxis()->CenterTitle(true);
    graphs[0]->GetXaxis()->SetTitleOffset(1.1);
    graphs[0]->GetYaxis()->SetTitleOffset(1.1);
    graphs[0]->GetXaxis()->SetTitleFont(42);
    graphs[0]->GetYaxis()->SetTitleFont(42);
    graphs[0]->GetXaxis()->SetLabelFont(42);
    graphs[0]->GetYaxis()->SetLabelFont(42);
    graphs[0]->GetXaxis()->SetTitleSize(0.044);
    graphs[0]->GetYaxis()->SetTitleSize(0.044);
    graphs[0]->GetXaxis()->SetLabelSize(0.036);
    graphs[0]->GetYaxis()->SetLabelSize(0.036);

    // Axis range and title.
    graphs[0]->GetXaxis()->SetNdivisions(505, kTRUE);  // 5 major divisions = 50/5 = 10 units each
    graphs[0]->GetXaxis()->SetLimits(140, 400);
    graphs[0]->GetYaxis()->SetRangeUser(0.08, 0.4);
    graphs[0]->GetXaxis()->SetTitle("T_{eff} = #LTp_{T}#GT / 3 [MeV]");
    graphs[0]->GetYaxis()->SetTitle("dln #LTp_{T}#GT / dln N_{ch}");

    // Get Hijing and Trajectum results.
    TFile *trajectum_file = TFile::Open("../../cs2_Trajectum_FCALCent_eta1_2.root", "READ");
    TFile *hijing_file = TFile::Open("../../cs2_MC-MB-HIJING_HF4eta5_trks1p0_Fit0p0-2p0_BoostInvariant.root", "READ");

    TGraphErrors *trajectum_graph = (TGraphErrors*) trajectum_file->Get("Graph;1");
    TGraphErrors *hijing_graph = (TGraphErrors*) hijing_file->Get("Graph;2");
    trajectum_file->Close();
    hijing_file->Close();

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

    graphs[0]->Draw("APL");
    for (size_t i = 1; i < graphs.size(); ++i){
        graphs[i]->Draw("PL SAME");
    }
    //line->Draw("SAME");
    trajectum_graph->Draw("PE SAME");
    hijing_graph->Draw("PE SAME");

    auto leg = new TLegend(0.70,0.56,0.90,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->AddEntry(graphs[0],"p_{T} > 0 GeV","lep");
    leg->AddEntry(graphs[1],"p_{T} > 0.3 GeV","lep");
    leg->AddEntry(graphs[2],"p_{T} > 0.5 GeV","lep");
    leg->AddEntry(graphs[3],"0 < p_{T} < 3 GeV","lep");
    leg->AddEntry(graphs[4],"0 < p_{T} < 5 GeV","lep");
    leg->AddEntry(graphs[5],"0 < p_{T} < 10 GeV","lep");
    leg->AddEntry(trajectum_graph,"pPb Trajectum","lep");
    leg->AddEntry(hijing_graph,"pPb Hijing","lep");
    leg->Draw();

    draw_CMS_Header();
    draw_CMS_Header("|#eta| < 1.5", 0.16, 0.74, 0.038, 11);
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