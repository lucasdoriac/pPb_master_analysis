// Plots 'hfSumEtPb' and 'hfSumEtp' TH1D histograms on the same TCanvas for each respective collision energy dataset.
// 'hfSumEtPb' contains the TH1D histogram of the sum of the transverse component of the total collision energy deposited on the Pb-going HF detector.
// 'hfSumEtp' same observable but for the p-going HF detector.
//
// --- lucasdoriadecarvalho@gmail.com

// To-do's
// - Check Bins functions.
// - Add struct solution to current messy string passing.


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
bool set_LogScale = false;
const std::string output_extension = ".pdf";
const std::string base_output_path = "../../../../mnt/c/Users/lucas/Documents/";

// --- Path to data files
const std::string data_File_5TeV = "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root";
const std::string data_File_8TeV = "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root";

// --- Function prototypes
TCanvas* makeCanvas(const char *name = "c", const char *title = "Canvas", int width = 800, int height = 600);
TH1D* get_th1d(const std::string& filename, const std::string& histogram);
void plot_th1(const std::string& coll_energy, TH1D *h1, TH1D *h2);
void refine_hist(TH1D *h1, Color_t color_h1, TH1D *h2, Color_t color_h2);
int GetLastNonZeroBin(const TH1D *hist);
double GetLastNonZeroX(const TH1D *hist);
void logError(const std::string& message);


// --- main () ---
void plot_hfSumEt_sNN(){

    gROOT->SetBatch(kTRUE);
    TH1D *h1, *h2;

    // Gets histograms from "/QA_histograms/"
    h1 = get_th1d(data_File_5TeV, "hfSumEtPb"); //h1-> Pb-going side.
    h2 = get_th1d(data_File_5TeV, "hfSumEtp"); //h2-> p-going side.
    // Plots histograms to files.
    plot_th1(data_File_5TeV, h1, h2);

    h1 = get_th1d(data_File_8TeV, "hfSumEtPb");
    h2 = get_th1d(data_File_8TeV, "hfSumEtp");
    // Plots histograms to files.
    plot_th1(data_File_8TeV, h1, h2);

}

// --- Function definitions

TCanvas* makeCanvas(const char* name, const char *title, int width, int height) 
{
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

TH1D* get_th1d(const std::string& filename, const std::string& histogram){

    // Gets TH1D histogram located on "QA_histograms" from filename.
    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()){
        logError("Error: from get_th1d: Could not open file " + filename);
        exit(1);
    }

    TDirectory *dir = (TDirectory*)file->Get("QA_histograms");
    if (!dir) {
        logError("Error: from get_th1d: Could not load directory 'QA_histograms'.");
        file->Close();
        exit(1);
    }
    
    TH1D *hist = (TH1D*)dir->Get(histogram.c_str());
    if (!hist){
        logError("Error: from get_th1d: Could not load TH1 histogram hist.");
        exit(1);
    }
    hist->SetDirectory(0);
    file->Close();

    return hist;
}

void plot_th1(const std::string& coll_energy, TH1D *h1, TH1D *h2){

    // Plots TH1D histograms. Pb vs p going side HF calorimeter measured energies.
    TCanvas* c = makeCanvas();

    // Prepares histograms for plot.
    refine_hist(h1, kGreen, h2, kCyan);

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
    if(coll_energy == data_File_5TeV) lumiText = "pPb (0.509 nb^{-1})  5.02 TeV";
    else if(coll_energy == data_File_8TeV) lumiText = "pPb (186.0 nb^{-1})  8.16 TeV";
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.04);
    latex2.SetTextFont(42);
    latex2.SetTextAlign(31);  // right-aligned
    latex2.DrawLatex(0.94, y, lumiText);


    gPad->RedrawAxis();
    c->Modified();
    c->Update();

    // Output configurations. This will be better once we create the structs.
    std::string energy_sufix;
    if (coll_energy == data_File_5TeV) energy_sufix = "_5TeV";
    else if (coll_energy == data_File_8TeV) energy_sufix = "_8TeV";

    std::string base_output_name = "hfSumEt_sNN";
    std::string full_path = base_output_path + base_output_name + energy_sufix + output_extension;

    c->SaveAs(TString(full_path));
    delete c;
}

void refine_hist(TH1D* h1, Color_t color_h1, TH1D *h2, Color_t color_h2){
    
    if(!h1 || !h2){
        logError("From refine_hist: histograms were not passed correctly as variables.");
        exit(1);
    }
    // General settings - Only needed for h1.
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle("E_{T,sum}^{HF} [GeV]");
    if(set_Normalization) h1->GetYaxis()->SetTitle("Density of events");
    else h1->GetYaxis()->SetTitle("Number of Events (#times 10^{6})");

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
    if(set_Normalization){
        Int = h1->Integral();
        if (Int > 0) h1->Scale(1.0 / Int);
        Int = h2->Integral();
        if (Int > 0) h2->Scale(1.0 / Int);
    }

    // Settings for better plot of raw counts. No normalization. 
    if(!set_Normalization && !set_LogScale){
        h1->GetYaxis()->SetNoExponent(kTRUE);
        h1->Scale(1.0 / 1e+6);
        h2->Scale(1.0 / 1e+6);
    }

    // Adjusting X-axis range.
    int lastBin = GetLastNonZeroBin(h1);
    double xmax = GetLastNonZeroX(h1);
    h1->GetXaxis()->SetRangeUser(0.0, xmax + 10.0);
}


int GetLastNonZeroBin(const TH1D *hist) {
    if (!hist) return -1;  // protect against null
    for (int i = hist->GetNbinsX(); i >= 1; --i) {
        if (hist->GetBinContent(i) > 0) {
            return i; // return index of last nonzero bin
        }
    }
    return -1; // no nonzero bins
}

double GetLastNonZeroX(const TH1D *hist) {
    int bin = GetLastNonZeroBin(hist);
    if (bin < 0) return hist->GetXaxis()->GetXmin(); // fallback
    return hist->GetXaxis()->GetBinUpEdge(bin); // upper edge of that bin
}


void logError(const std::string& message){
    // Set off stream to error.log file.
    std::ofstream err_log("plot_hfSumEt_sNN.log", std::ios::app);

    if (err_log.is_open()){
        err_log << message << std::endl;
    }
    else {
        std::cerr << "From logError: Could not open 'plot_hfSumEt_sNN.log' for writing." << std::endl;
    }
    
}