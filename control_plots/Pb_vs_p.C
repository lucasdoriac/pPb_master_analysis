// Macro used to compare energy deposited in the HF calorimeter on Pb and p side for the same collision energy.
// Works for both 8 TeV and 5 TeV.

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TKey.h>
#include <cstring>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>

int GetLastNonZeroBin(const TH1D *hist);
double GetLastNonZeroX(const TH1D *hist);
void style_histogram(const char* target_histogram, TH1D *h1, Color_t color_h1, TH1D *h2, Color_t color_h2, bool set_x_range = false, int sparse_axis = 0);
TLegend* makeLegend(TH1D *h1, const char *label1, TH1D *h2, const char *label2, double x1 = 0.75, double y1 = 0.7, double x2 = 0.95, double y2 = 0.82); 
TCanvas* makeCanvas(const char *name = "c", const char *title = "Canvas", bool logy = false, int width = 1280, int height = 720);
TString return_Xaxis_title(const char* variable, int sparse_axis = 0);
void drawCMSHeader(const char* extraText = "#it{Work in Progress}", const char* lumiText  = "pPb (0.509 nb^{-1}) 5.02 TeV", double x = 0.11, double y = 0.94);

// if 8 -> (186.0 nb^{-1}) 8.16 TeV
// if 5 -> (0.509 nb^{-1}) 5.02 TeV


void Pb_vs_p(int which = 0){

	TFile *f1 = nullptr;

	if(which == 8) {
		f1 = TFile::Open("../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root", "READ");
	}
	else if (which == 5) f1 = TFile::Open("../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root", "READ");

	TDirectory *dir1 = (TDirectory*)f1->Get("QA_histograms");

	TH1D *h1 = (TH1D*)dir1->Get("hfSumEtPb");
	h1->SetDirectory(0);

	TH1D *h2 = (TH1D*)dir1->Get("hfSumEtp");
	h2->SetDirectory(0);

	// Histogram plotting section.
	style_histogram("Pb_vs_p", h1, kOrange, h2, kGreen, true);
	TCanvas* c = makeCanvas("c", "canvas", true);
	h1->Draw("hist f ][");
	h2->Draw("same hist ][");
	gPad->RedrawAxis();
	makeLegend(h1,"HF - Pb side",h2,"HF - p side")->Draw();

	// Draw CMS header.
	drawCMSHeader();

	c->Modified();
	c->Update();
	TString outname = Form("../../../mnt/c/Users/lucas/Documents/Pb_vs_p_%d.png", which);
	c->SaveAs(outname);

}

// --- Function definitions ---

// Histogram settings.
void style_histogram(const char* target_histogram, TH1D *h1, Color_t color_h1, TH1D *h2, Color_t color_h2, bool set_x_range = false, int sparse_axis = 0){
	if(!h1 || !h2) return; // Avoid null pointer crash.

	// General settings - Only needed for h1.
    h1->GetXaxis()->SetTitle("E_{T,sum}^{HF} [GeV]");
	h1->GetYaxis()->SetTitle("Number of events");

	h1->SetTitle("");
    h1->GetXaxis()->SetTitleOffset(1.0);
    h1->GetYaxis()->SetTitleOffset(1.0);
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->SetLabelFont(42);
    h1->GetYaxis()->SetLabelFont(42);
    h1->GetYaxis()->SetTitleFont(42);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.036);
    h1->GetYaxis()->SetTitleSize(0.044);
    h1->GetYaxis()->SetLabelSize(0.036);
	
	//h1->GetXaxis()->SetTitleSize(0.05);
	//h1->GetXaxis()->SetLabelSize(0.04);

	// Adjusting X-axis range.
	if(set_x_range){
		int lastBin = GetLastNonZeroBin(h1);
		double xmax = GetLastNonZeroX(h1);
		double xmin = h1->GetXaxis()->GetXmin();
		if( strcmp(target_histogram, "ptresolution") == 0 ) h1->GetXaxis()->SetRangeUser(xmin, xmax + 0.01);
		else h1->GetXaxis()->SetRangeUser(xmin, xmax + 5.0); 
	}

	// 8 TeV histogram SETTINGS
	h1->SetStats(0);
	h1->SetLineWidth(1);
    h1->SetLineStyle(1);
	h1->SetFillColor(color_h1 - 9);
    h1->SetLineColor(color_h1 + 2);
    h1->SetFillStyle(1001);

    double int1 = h1->Integral("width");
	if (int1 > 0) h1->Scale(1.0 / int1);
	
	// 5 TeV histogram SETTINGS.
	h2->SetStats(0);
	h2->SetLineWidth(1);
    h2->SetLineStyle(1);
	h2->SetFillColor(color_h2 - 9);
    h2->SetLineColor(color_h2 + 2);
    h2->SetFillStyle(3001);

    double int2 = h2->Integral("width");
	if (int2 > 0) h2->Scale(1.0 / int2);
}

// Legend settings function.
TLegend* makeLegend(TH1D *h1, const char *label1, 
                    TH1D *h2, const char *label2,
                    double x1, double y1,
                    double x2, double y2) 
{
    // Create legend
    TLegend* leg = new TLegend(x1, y1, x2, y2);

    // Add entries
    if (h1) leg->AddEntry(h1, label1, "f");
    if (h2) leg->AddEntry(h2, label2, "f");

    // Style legend
    leg->SetBorderSize(0);   // no border
    leg->SetFillStyle(0);    // transparent background
    leg->SetTextSize(0.038);  // readable font size
    leg->SetTextFont(42);    // Text font
    leg->SetMargin(0.2);     // spacing inside
    leg->SetEntrySeparation(0.04);  // smaller separation means shorter marker box

    return leg;
}

// Canvas settings.
TCanvas* makeCanvas(const char* name = "c", const char *title = "Canvas", bool logy = false, 
                    int width = 1280, int height = 720) 
{
    // Create canvas
    TCanvas* c = new TCanvas(name, title, width, height);

    // Margins
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Grid & ticks
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis

    // Background
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);

    // Thicker border/frame
    c->SetFrameLineWidth(2);  // default is 1, increase for thicker axes box

    // Log scale option
    if (logy) c->SetLogy();

    return c;
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

TString return_Xaxis_title(const char* variable, int sparse_axis = 0) {
    std::string var(variable);

    if (var == "hfSumEtPb" || var == "hfSumEtp")
        return "E_{T,sum}^{HF} [GeV]";
    else if (var == "vzhist")
        return "#it{z} coordinate of primary vertex [cm]";
    else if (var == "dxyoversigmadxy")
        return "dxy/#sigmadxy";
    else if (var == "dzoversigmadz")
        return "dz/#sigmadz";
    else if (var == "ptresolution")
        return "p_{T}/#Delta p_{T}";
    else if (var == "hist_reco_trk_corr") {
        if (sparse_axis == 0) return "p_{T} [GeV]";
        if (sparse_axis == 1) return "#eta";
        if (sparse_axis == 2) return "#phi";
        if (sparse_axis == 3) return "Charge [e]";
    }

    return "Unknown variable"; // No known variable was given.
}

void drawCMSHeader(const char* extraText = "#it{Work in Progress}",
                   const char* lumiText  = "pPb (186.0 nb^{-1}) 8.16 TeV, (0.509 nb^{-1}) 5.02 TeV",
                   double x = 0.11, double y = 0.94)
{
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

    // Draw "CMS [Preliminary]" on the left
    latex.DrawLatex(x, y, cmsText);

    // Draw luminosity/energy info on the right (optional)
    if (lumiText && strlen(lumiText) > 0) {
        TLatex latex2;
        latex2.SetNDC();
        latex2.SetTextSize(0.035);
        latex2.SetTextFont(42);
        latex2.SetTextAlign(31);  // right-aligned
        latex2.DrawLatex(0.95, y, lumiText);
    }
}