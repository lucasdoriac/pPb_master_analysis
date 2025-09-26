// Plots 'hfSumEtPb' and 'hfSumEtp' TH1D histograms on the same TCanvas for each respective collision energy dataset.
// 'hfSumEtPb' contains the TH1D histogram of the sum of the transverse component of the total collision energy deposited on the Pb-going HF detector.
// 'hfSumEtp' same observable but for the p-going HF detector.
//
// --- lucasdoriadecarvalho@gmail.com

// --- Headers
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstring>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

// --- Global settings (variables)
bool set_Normalization = false;
bool set_LogScale = true;
bool set_range = true;
const std::string output_extension = ".pdf";
const std::string output_name = "hfSumEt_sNN_"
const std::string base_output_path = "../../../../mnt/c/Users/lucas/Documents/";

// --- Structs
struct DataFile {
    std::string path;
    std::string label;
};

// --- struct to 5TeV dataset
const DataFile data_File_5TeV = {
    "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root",
    "5TeV"
};

// --- struct to 8TeV dataset
const DataFile data_File_8TeV = {
    "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root",
    "8TeV"
};


// --- Function prototypes
TCanvas* makeCanvas(const char *name = "c", const char *title = "Canvas", int width = 800, int height = 600);
TH1D* get_th1d(const std::string& filename, const std::string& histogram);
void plot_th1(const std::string& coll_energy, TH1D *h1, TH1D *h2);
void refine_hist(const std::string& coll_energy, TH1D *h1, Color_t color_h1, TH1D *h2, Color_t color_h2);
void logError(const std::string& message);


// --- main () ---
void plot_hfSumEt_sNN() {

    gROOT->SetBatch(kTRUE); // This tells ROOT to run in Batch mode, i. e. no GUI or pop-ups.
    TH1D *h1, *h2;

    // Gets histograms from "/QA_histograms/"
    h1 = get_th1d(data_File_5TeV.path, "hfSumEtPb"); //h1-> Pb-going side.
    h2 = get_th1d(data_File_5TeV.path, "hfSumEtp"); //h2-> p-going side.
    // Plots histograms to files.
    if(h1 && h2) plot_th1(data_File_5TeV.label, h1, h2);
    else logError("Skipping 5TeV plot due to missing histogram.");

    h1 = get_th1d(data_File_8TeV.path, "hfSumEtPb");
    h2 = get_th1d(data_File_8TeV.path, "hfSumEtp");
    // Plots histograms to files.
    if(h1 && h2) plot_th1(data_File_8TeV.label, h1, h2);
    else logError("Skipping 8TeV plot due to missing histogram.");
}

// --- Function definitions

TCanvas* makeCanvas(const char* name, const char *title, int width, int height) {

    TCanvas* c = new TCanvas(name, title, width, height);

    // Margins
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Grid & ticks
    c->SetTickx(1);     // sets ticks on top x-axis
    c->SetTicky(1);     // sets ticks on right y-axis

    // Background
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);

    // Thicker border/frame
    c->SetFrameLineWidth(2);  // default is 1, increase for thicker axes box

    // Log scale option
    if(set_LogScale) c->SetLogy();

    return c;
}

TH1D* get_th1d(const std::string& filename, const std::string& histogram) {

    // Gets TH1D histogram located on "QA_histograms" from filename.
    TFile *file = TFile::Open(filename.c_str(), "READ");
    if(!file || file->IsZombie()) {
        logError("Error: from get_th1d: Could not open file " + filename);
        return nullptr;
    }

    TDirectory *dir = (TDirectory*)file->Get("QA_histograms");
    if(!dir) {
        logError("Error: from get_th1d: Could not load directory 'QA_histograms'.");
        file->Close();
        return nullptr;
    }
    
    TH1D *hist = (TH1D*)dir->Get(histogram.c_str());
    if (!hist){
        logError("Error: from get_th1d: Could not load TH1 histogram hist.");
        return nullptr;
    }
    hist->SetDirectory(0);
    file->Close();

    return hist;
}

void plot_th1(const std::string& coll_energy, TH1D *h1, TH1D *h2) {

    // Plots TH1D histograms. Pb vs p going side HF calorimeter measured energies.
    TCanvas* c = makeCanvas();

    // Prepares histograms for plot.
    refine_hist(coll_energy, h1, kGreen, h2, kCyan);

    // Currently plotting as filling. To plot as markers use the commented code.
    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    // Create TLegend.
    TLegend* leg = new TLegend(0.68, 0.70, 0.88, 0.85); // x1, y1, x2, y2.

    // Add entries
    leg->AddEntry(h1, "Pb HF", "f");
    leg->AddEntry(h2, "p HF", "f");

    // Legend box settings.
    leg->SetHeader("Data"); // There's a "C" flag that centers this.
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->Draw();

    // Draw CMS header.
    double x = 0.13;
    double y = 0.94;
    // Draw CMS preliminary.
    TString cmsText = "#bf{CMS} #it{Preliminary}";
    TLatex latex;
    latex.SetNDC();              // use normalized coordinates
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);       // 42 always Helvetica for what i see.
    latex.SetTextAlign(11);      // left-aligned, top
    latex.DrawLatex(x, y, cmsText);

    // Draw luminosity/energy info on the right. This changes depending on what you want.
    TString lumiText;
    if(coll_energy == data_File_5TeV.label) lumiText = "pPb (0.509 nb^{-1})  5.02 TeV";
    else if(coll_energy == data_File_8TeV.label) lumiText = "pPb (186.0 nb^{-1})  8.16 TeV";
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.04);
    latex2.SetTextFont(42);
    latex2.SetTextAlign(31);  // right-aligned
    latex2.DrawLatex(0.94, y, lumiText);


    gPad->RedrawAxis();
    c->Modified();
    c->Update();

    // Output string.
    std::string full_path = base_output_path + output_name + coll_energy + output_extension;

    c->SaveAs(TString(full_path));
    delete c;
}

void refine_hist(const std::string& coll_energy, TH1D* h1, Color_t color_h1, TH1D *h2, Color_t color_h2) {
    
    if(!h1 || !h2) {
        logError("From refine_hist: histograms were not passed correctly as variables.");
        return;
    }
    // General settings - Only needed for h1.
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle("E_{T,sum}^{HF} [GeV]");
    if(set_Normalization) h1->GetYaxis()->SetTitle("Density of events");
    else if (!set_Normalization && !set_LogScale) h1->GetYaxis()->SetTitle("Number of Events (#times 10^{6})");
    else if (!set_Normalization && set_LogScale) h1->GetYaxis()->SetTitle("Number of Events");

    h1->GetXaxis()->CenterTitle(true);
    h1->GetYaxis()->CenterTitle(true);

    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetYaxis()->SetTitleOffset(1.1);
    
    h1->GetXaxis()->SetTitleFont(42); // 42 for verdana. Standard. 62 for bold standard. I liked 132 as well.
    h1->GetYaxis()->SetTitleFont(42);

    h1->GetXaxis()->SetLabelFont(42);
    h1->GetYaxis()->SetLabelFont(42);

    h1->GetXaxis()->SetTitleSize(0.044);
    h1->GetYaxis()->SetTitleSize(0.044);

    h1->GetXaxis()->SetLabelSize(0.036);
    h1->GetYaxis()->SetLabelSize(0.036);


    // Pb-going side histogram SETTINGS (Currently plotting with filling style.)
    h1->SetStats(0);
    h1->SetLineWidth(1);
    h1->SetLineStyle(1);
    h1->SetFillColorAlpha(color_h1 + 1, 0.5);
    h1->SetLineColor(color_h1 + 2);
    h1->SetFillStyle(1001);

    // p-going side histogram SETTINGS.
    h2->SetStats(0);
    h2->SetLineWidth(1);
    h2->SetLineStyle(1);
    h2->SetFillColorAlpha(color_h2 + 1, 0.5);
    h2->SetLineColor(color_h2 + 2);
    h2->SetFillStyle(1001);

    // Normalization settings.
    double Int;
    if(set_Normalization) {
        // By default, the hist->Integral() method doesn't include underflow and overflow bins.
        // That means the interval of integration is implicitly Integral(1, nbins);
        // If we want to include underflow and overflow bins we need to explicitly set Integral(0, nbins+1);
        Int = h1->Integral();
        if(Int > 0) h1->Scale(1.0 / Int);
        Int = h2->Integral();
        if(Int > 0) h2->Scale(1.0 / Int);
    }

    // Settings for better plot of raw counts. No normalization. 
    if(!set_Normalization && !set_LogScale) {
        h1->GetYaxis()->SetNoExponent(kTRUE);
        h1->Scale(1.0 / 1e+6);
        h2->Scale(1.0 / 1e+6);
    }

    // Adjusting axes range. There are helpers to do this automatically in the end of the code if wanted.
    // I generally use them first to see the data range and then adjust the axes range manually for better plots.
    // Change as wanted.
    if(coll_energy == data_File_5TeV.label && set_range){
        h1->GetXaxis()->SetRangeUser(0.0, 120.0);
        h1->GetYaxis()->SetRangeUser(1.0, 1e+8);
    }
    else if(coll_energy == data_File_8TeV.label && set_range){
        h1->GetXaxis()->SetRangeUser(0.0, 160.0);
        h1->GetYaxis()->SetRangeUser(1.0, 1e+9);
    }
}

void logError(const std::string& message) {
    // Set error output to err_log file.
    std::ofstream err_log("plot_hfSumEt_sNN.log", std::ios::app);

    if(err_log.is_open()) {
        err_log << message << std::endl;
    }
    else {
        std::cerr << "From logError: Could not open 'plot_hfSumEt_sNN.log' for writing." << std::endl;
    }
    
}