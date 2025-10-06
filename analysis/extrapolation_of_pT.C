// Macro developed to reproduce the results for <p_T> vs N_ch for both 5.02 and 8.16 TeV pPb collision energies.
// Unpublished data by the CMS but exposed in the article [https://cds.cern.ch/record/2931093/files/HIN-25-001-pas.pdf].
// This work is part of the CMS Collaboration and uses CMS Preliminary Data.

/*
The macro currently correctly estimates the <p_T>(N_ch) value for the 0-1% centrality class for the 5TeV dataset.
Now we change the macro to test with other centralities and datasets.
This is currently a copy of the working macro 'hagedorn_fit_2.C'.



The main change of the macros hagedorn_fit.C to the study_fit.C ones are the change in the way we fill the fit_hist.
In the study_fits we were using a reasonable sample (10^8) to randomly generate entries in fit_hist according to the Hagedorn PDF.
However, that method was shown to overestimate both <p_T> and underestimate N_ch.
The reason for that had to be a poor simulation of tracks on the 0.0 and 0.3 GeV region.
In other words, the method where we first sample random entries and later scale the fit_hist with Scale(), with the correct scale_factor,
was insufficient to produce the desired results in the [0, 0.3] GeV region.
The solution I found was to change the sample value and still use TF1::GetRandom() to match the desired normalization.
The initial sample size in which fit_hist is filled in a faster way is defined now by
    long SAMPLE = static_cast<long>(original_hist_norm);
Where original_hist_norm is the Integral of the specific class p_T distribution from 0.3 to 1.5 GeV (the desired fit interval).
After the filling, the Integral of the fit_hist on the same interval will always be less than original_hist_norm because some tracks will have been generated on the interval 0 to 0.3 GeV.
That way, we always complement the random track generation with the Hagedorn PDF using a while loop
    while (fit_hist_norm <= original_hist_norm){
        fit_hist->Fill(pT_fit->GetRandom());
        fit_hist_norm = fit_hist->Integral(bin_min, bin_max);
    }
That will always be activated for the argument explained above.
The macro will stay in the while loop until the Integral of both histograms match on the interval 0.3 to 1.5 GeV.
When the macro exits the loop, it correctly calculates the mean p_T, along with its error, and the corresponding N_ch.

// --- To do's ---
// - 

--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => > 35 GeV, > 44 GeV.

---------------Lucas Carvalho---------------*/

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
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

// --- Macro settings --- # Change according to what the user wants.
int SAMPLE = 1e+8;
double low_eta = -1.5;
double high_eta = 1.5;
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-6; // GeV
const std::string output_extension = ".png";
const std::string output_name = "fit";
const std::string base_output_path = "../../../../mnt/c/Users/lucas/Documents/";


// --- Vectors to store 5TeV data
std::vector<double> mean_pT_data_5TeV;
std::vector<double> mean_pT_errors_5TeV;
std::vector<double> N_ch_values_5TeV;

// --- Vectors to store 8TeV data
std::vector<double> mean_pT_data_8TeV;
std::vector<double> mean_pT_errors_8TeV;
std::vector<double> N_ch_values_8TeV;

// --- Structs
struct DataFile {
    std::string path;
    std::string label;
    std::string sufix;
};

// --- struct to 5TeV dataset
const DataFile data_File_5TeV = {
    "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root",
    "5TeV",
    "_5TeV"
};

// --- struct to 8TeV dataset
const DataFile data_File_8TeV = {
    "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root",
    "8TeV",
    "_8TeV"
};

// --- Function prototypes ---
double get_n_events(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
TH1D* get_proj_hist(const std::string& filename, double EHFmin = 0.0, double EHFmax = 250.0);
void fit_Hagedorn(TH1D* hist_pT, double EHFmin = 0.0, double EHFmax = 250.0);

//void print_pT_vs_Nch();

// --- main() ---
void extrapolation_of_pT(){
    //gROOT->SetBatch(kTRUE); // This tells ROOT to run in Batch mode, i. e. no GUI or pop-ups.

    //fit_Hagedorn(get_proj_hist(data_File_5TeV.path, 2.5, 11.5), 2.5, 11.5);
    //fit_Hagedorn(get_proj_hist(data_File_5TeV.path, 11.5, 35.0), 11.5, 35.0);
    fit_Hagedorn(get_proj_hist(data_File_5TeV.path, 35.0), 35.0);

    //fit_Hagedorn(get_proj_hist(data_File_8TeV.path, 2.5, 14.5), 2.5, 14.5);
    //fit_Hagedorn(get_proj_hist(data_File_8TeV.path, 14.5, 44.0), 14.5, 44.0);
    //fit_Hagedorn(get_proj_hist(data_File_8TeV.path, 44.0), 44.0);


}

// --- Function definitions ---
double get_n_events(const std::string& filename, double EHFmin, double EHFmax){
    
    // Return number of events n_events in a determined E_HF range [EHFmin, EHFmax], or centrality class, 
    // by integrating the 'hfSumEtPb' TH1 histogram in "QA_histograms".

    TFile *f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: from get_n_events function: Could not open file!" << std::endl;
        exit(1);
    }

    TDirectory *dir = (TDirectory*)f->Get("QA_histograms");
    if (!dir) {
        std::cerr << "Error: from get_n_events function: Directory QA_histograms not found!" << std::endl;
        f->Close();
        exit(1);
    }

    TH1D *hf_hist = (TH1D*)dir->Get("hfSumEtPb");
    if (!hf_hist) {
        std::cerr << "Histogram hfSumEtPb not found!" << std::endl;
        f->Close();
        exit(1);
    }
    hf_hist->SetDirectory(0);
    f->Close();

    int bin_min = hf_hist->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hf_hist->GetXaxis()->FindBin(EHFmax - delta);

    // Just a quick check.
    double lowEdge   = hf_hist->GetXaxis()->GetBinLowEdge(bin_min);
    double highEdge  = hf_hist->GetXaxis()->GetBinUpEdge(bin_min);
    double lowEdge_   = hf_hist->GetXaxis()->GetBinLowEdge(bin_max);
    double highEdge_  = hf_hist->GetXaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From get_n_events: Integrating 'hfSumEtPb' from bin [%.1f,%.1f]GeV to bin [%.1f,%.1f]GeV \n\n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    double n_events = hf_hist->Integral(bin_min, bin_max);
    
    return n_events;
}


TH1D* get_proj_hist(const std::string& filename, double EHFmin, double EHFmax){

    // Return projected TH1D histogram of p_T track distribution on the range p_T > 0.3 [GeV]
    // for selected centrality class and pseudorapidity window.

    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: from get_proj_hist function: Could not open file!" << std::endl;
        exit(1);
    }

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    if (!dir) {
        std::cerr << "Error: from get_proj_hist function: Directory Analysis_histograms not found!" << std::endl;
        file->Close();
        exit(1);
    }
    
    // Load TH3 histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    if (!hist_HFSumPb_vs_pt_eta){
        std::cerr << "Error: from get_proj_hist function: Error while creating TH3 histogram. " << std::endl;
        file->Close();
        exit(1);
    }
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    file->Close();

    // Setting pseudorapidity window [low_eta, high_eta].
    int zmin = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(low_eta + delta);
    int zmax = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(high_eta - delta);
    hist_HFSumPb_vs_pt_eta->GetZaxis()->SetRange(zmin, zmax);

    // Just a quick check.
    double lowEdge = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinLowEdge(zmin);
    double highEdge = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(zmin);
    double lowEdge_ = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinLowEdge(zmax);
    double highEdge_ = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(zmax);
    printf("\n-> From get_proj_hist: Pseudorapidity window from bin [%.1f,%.1f] to bin [%.1f,%.1f] \n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Projection to TH2 by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *hist_HFSumPb_vs_pt = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");
    if (!hist_HFSumPb_vs_pt){
        std::cerr << "Error: from get_proj_hist function: Error while creating TH2 histogram. " << std::endl;
        exit(1);
    }

    // Setting centrality class defined by HF energy cutoffs EHFmin, EHFmax.
    int binXmin = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmin + delta);
    int binXmax = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmax - delta);

    // Just a quick check.
    lowEdge = hist_HFSumPb_vs_pt->GetXaxis()->GetBinLowEdge(binXmin);
    highEdge = hist_HFSumPb_vs_pt->GetXaxis()->GetBinUpEdge(binXmin);
    lowEdge_ = hist_HFSumPb_vs_pt->GetXaxis()->GetBinLowEdge(binXmax);
    highEdge_ = hist_HFSumPb_vs_pt->GetXaxis()->GetBinUpEdge(binXmax);
    printf("\n-> From get_proj_hist: Integrating from bin [%.1f,%.1f]GeV to bin [%.1f,%.1f]GeV \n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Projects TH2 on TH1 for the defined centrality class.
    TH1D *hist_pT = (TH1D*)hist_HFSumPb_vs_pt->ProjectionY("", binXmin, binXmax);
    if (!hist_pT){
        std::cerr << "Error: from get_proj_hist function: Error while creating TH1 histogram. " << std::endl;
        exit(1);
    }
    hist_pT->SetDirectory(0);

    printf("\n -> From get_proj_hist function: hist_pT norm = %.3e \n", hist_pT->Integral());

    return hist_pT;
}

void fit_Hagedorn(TH1D *hist_pT, double EHFmin, double EHFmax){
    
    TH1D *original_hist = nullptr;

    original_hist = (TH1D*)hist_pT->Clone("original_hist");
    if(!original_hist){
        std::cerr << "Error: from fit_Hagedorn: error while loading hist_pT. " << std::endl;
        return;
    }
    original_hist->SetDirectory(0);
    printf("\n\n -> From fit_Hagedorn function: original_hist norm = %.3e \n", original_hist->Integral());

    // Get norm only on the target range: 0.3 to 1.5 GeV.
    int bin_min = original_hist->GetXaxis()->FindBin(lower_pt_forFit + delta);
    int bin_max = original_hist->GetXaxis()->FindBin(upper_pt_forFit - delta);
    double original_hist_norm = original_hist->Integral(bin_min, bin_max); // Integral over [0.3,1.5].
    printf("\n -> From fit_Hagedorn function: original_hist norm over [0.3, 1.5] GeV = %.3e \n\n", original_hist_norm);

    // Hagedorn TF1. Function declaration section.
    TF1* pT_fit;
    
    // Original.
    pT_fit = new TF1(Form("ptfit_%f_%f",EHFmin,EHFmax),"[0]*x*pow(1.+1./sqrt(1.-[1]*[1])*(sqrt(x*x+[4]*[4])-x*[1])/[3]/[2],-[3])",0.,upper_pt_forFit);
    pT_fit->SetParameters(7500000000.,0.3,0.1,6.,0.14);//We used these values for initialization

    // Edit by Lucas.
    //double AMPL = original_hist->GetMaximum();
    //pT_fit = new TF1("ptfit","[0]*x*pow(1.+1./sqrt(1.-[1]*[1])*(sqrt(x*x+[4]*[4])-x*[1])/[3]/[2],-[3])",0.,upper_pt_forFit);
    //pT_fit->SetParameters(AMPL,0.3,0.1,6.,0.14);//We used these values for initialization
    
    pT_fit->FixParameter(4,0.13957);//pion mass
    pT_fit->FixParameter(1,0.5010);//related to radial flow velocity
    pT_fit->SetParLimits(2,0.,0.5);//kinetic freeze-out temperature in GeV
    pT_fit->SetParLimits(3,6.,9.);//n - free parameter no physical meaning


    // Track cpu efficiency.
    TStopwatch timer;
    timer.Start();
    
    /*TCanvas *c = new TCanvas("c","Fit Canvas",800,600);
    original_hist->GetXaxis()->SetRangeUser(0.0, 3.0);  // for example 0–3 GeV
    original_hist->SetMarkerStyle(20);
    original_hist->SetMarkerColor(kBlack);
    original_hist->Draw("E");
    pT_fit->SetLineColor(kRed);*/

    // Hagedorn fit.
    original_hist->Fit(pT_fit,"R","",lower_pt_forFit,upper_pt_forFit);
    // Usando a flag "M" (improved Minuit2 minimizer) abaixo, os resultados sao visualmente melhores. A ver.
    //original_hist->Fit(pT_fit,"RM","",lower_pt_forFit,upper_pt_forFit);
    // A flag "V" retorna um log detalhado sobre o warning que estamos recebendo.
    //original_hist->Fit(pT_fit,"RV","",lower_pt_forFit,upper_pt_forFit);
    // A flag "L" utiliza outro metodo de fitting, ao inves do Chi2. Resultados tambem aparentemente melhores.
    //original_hist->Fit(pT_fit,"RL","",lower_pt_forFit,upper_pt_forFit);


    // Create histogram that will be filled with a Hagedorn distribution.
    int n_bins = bin_max; // number of bins between 0 and 1.5 GeV.
    TH1D *fit_hist = new TH1D("fit_hist", "fit histogram", n_bins, 0.0, 1.5);
    fit_hist->SetDirectory(0);
    double fit_hist_norm = 0.0;

    int count = 0;
    printf("-> Random generation of tracks with Hagedorn PDF. Sample size = %d \n", SAMPLE);    
        while(count < SAMPLE){
            double x = pT_fit->GetRandom();
            if(x >= 0.3 && x<= 1.5) count+=1;
            fit_hist->Fill(x);
        }
    double scale_factor = original_hist_norm / SAMPLE;
    printf("\n -> Scale factor = %.3f \n", scale_factor);
    fit_hist->Scale(scale_factor);

    fit_hist_norm = fit_hist->Integral(bin_min, bin_max);
    printf("-> Fit histogram complete. Final fit_hist norm over [0.3, 1.5] GeV = %.3e \n", fit_hist_norm);
    if(original_hist_norm == fit_hist_norm) printf("-> Normalization over [0.3, 1.5] succesfull. \n");
    else {
        printf("-> From fit_Hagedorn: Normalization over fit interval unsuccesfull! \n");
        exit(1);
    }

    // Create histogram extrapolated_hist. Results will be calculated over THIS histogram.
    TH1D *extrapolated_hist = (TH1D*)original_hist->Clone("extrapolated_hist");
    extrapolated_hist->SetDirectory(0);
    
    // Takes first three bin contents from fit_hist. The rest is taken from the real data from original_hist.
    double content, error;
    for(int i = 1; i <= extrapolated_hist->GetNbinsX(); i++) {
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
    
    // Get <p_T> and <p_T> error.
    double mean_pT = extrapolated_hist->GetMean();
    double mean_pT_error = extrapolated_hist->GetMeanError();

    // Sum number of tracks n_tracks as bin_contents.
    double n_tracks = extrapolated_hist->Integral();

    // Get n_events from QA_histograms respective file and calculate N_ch.
    double n_events = get_n_events(data_File_5TeV.path, EHFmin, EHFmax);
    double N_ch = (n_tracks/n_events);

    printf("--- Results --- \n");
    printf("--- Centrality class 0-1%% for 5TeV dataset --- \n");
    printf("-> Estimated <p_T> = %f \n", mean_pT);
    printf("-> Estimated <p_T> error = %f \n", mean_pT_error);
    printf("-> Estimated N_ch = %f \n", N_ch);
    printf("---\n\n");

    timer.Stop();
    std::cout << "-> Job finished in "
              << timer.RealTime() << " seconds (wall time), "
              << timer.CpuTime()  << " seconds (CPU time).\n\n";

}