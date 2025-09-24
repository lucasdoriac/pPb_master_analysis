/* Macro deisgned to open ROOT files from CERN.
Data collected from the CMS detector in 2016.
Data composed by physical measurements from pPb collisions at 5.02 and 8.16 TeV.

This macro will address the following problems of the old version:
-> Generalize code in order to plot whatever histogram object it is passed as parameter. ok
-> Include a flag to check if the histogram does exist in the .root file. ok
-> Control log scale on y-axis via boolean input parameter. ok
-> Assign each histogram name with a x-axis title so we don`t have to manually write all the time. ok
-> Include CMS header function as TLatex object. ok
-> Instead of defining new range for x-axis as only non-zero bins, increase a bit further so we can see the data approaching zero. ok
-> SetDirectory(0) should be immediately after assigning value to TH1D *h. ok
-> The block of code where the histogram is plotted need better format.

Developments so far:
-> Removed check_binWidth function;
-> Included boolean input parameter to control log scale on the histogram y-axis.
-> Included a flag to check if the histogram was found in the .root file.
-> Set new x-axis range to draw a little further from first zero bin.
-> SetDirectory(0) is now declared immediately after assigning value to histogram pointer.
-> Included a helper function to associate each histogram name to correct x-axis title.
-> Included a helper function to draw CMS Header with a TLatex object.
-> Code now plots TH1Ds even if they are stored in a THnSparse object.
-> Included a boolean that controls x-axis scaling. If true, scales range to dataset.
-> Only task remaining is to find the optimal plotting style and code format.

v3.2 -> Included better normalization technique to account for underflow and overflow. 

X-axis range is set AFTER normalization of the histogram relative to all the data range.

---------------Lucas Carvalho---------------*/
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TKey.h>
#include <cstring>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>


// --- Function declarations (prototypes) ---
int GetLastNonZeroBin(const TH1D *hist);
double GetLastNonZeroX(const TH1D *hist);
void style_histogram(const char* target_histogram, TH1D *h1, Color_t color_h1, TH1D *h2, Color_t color_h2, bool set_x_range = false, int sparse_axis = 0);
TLegend* makeLegend(TH1D *h1, const char *label1, TH1D *h2, const char *label2, double x1 = 0.75, double y1 = 0.7, double x2 = 0.95, double y2 = 0.82); 
TCanvas* makeCanvas(const char *name = "c", const char *title = "Canvas", bool logy = false, int width = 1280, int height = 720);
TString return_Xaxis_title(const char* variable, int sparse_axis = 0);
void drawCMSHeader(const char* extraText = "#it{Work in Progress}", const char* lumiText  = "pPb (186.0 nb^{#minus1}) 8.16 TeV, (0.509 nb^{#minus1}) 5.02 TeV", double x = 0.11, double y = 0.94);


// --- main() ---
void control_plots_v3(const char* target_histogram, bool log_y = false, bool set_x_range = false, int sparse_axis = 0){

	gROOT->SetBatch(kTRUE);

	bool found = false; // Boolean to check if target_histogram was found.

	TFile *f1 = TFile::Open("../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root", "READ");
	TFile *f2 = TFile::Open("../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root", "READ");
	if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
    	    std::cerr << "Error while opening files. \n" << std::endl;
    	    exit(1);
    	}

    TIter nextDirKey(f1->GetListOfKeys());
	TKey *dirKey;

	// Loop over top-level directories using their keys.
	while((dirKey = (TKey*) nextDirKey())){

		if (strcmp(dirKey->GetClassName(), "TDirectoryFile") != 0) continue; // If object class is NOT TDirectoryFile, move to the next Key.

		TDirectory *dir1 = (TDirectory*)dirKey->ReadObj(); // dir1 is a pointer to the TDirectoryFile.
        TDirectory *dir2 = (TDirectory*)f2->Get(dir1->GetName()); // dir2 is a pointer to the TDirectoryFile of same name in file2.

        if(!dir2){
        	std::cerr << "Directory " << dir1->GetName() << "not found in " << f2->GetName() << "\n" << std::endl;
        	exit(1);
        }

        std::cout << "Getting directory: " << dir1->GetName() << "\n" << std::endl;

        // Loop over sub-level structure to find target_histogram.
        TIter nextKey(dir1->GetListOfKeys());
        TKey *key;

        while((key = (TKey*)nextKey())){

		std::string hname = key->GetName();// Saves the name of the actual key (file) on hname.
			
		if (strcmp(hname.c_str(), target_histogram) != 0) continue; // If file name isn`t target_histogram, go to the next Key.

			else{ 	found = true;

				TObject *obj1 = key->ReadObj();
				if(!obj1){
					std::cout << "Error while reading TObject. \n" << std::endl;
					exit(1);
				}
				
				std::string classname = obj1->ClassName();
    			std::cout << "Found " << hname << " class: " << classname << std::endl;

    			if(obj1->InheritsFrom(TH1::Class())){

    				// Created the if but still didnt look at inside the brackets.
					TH1D *h1 = (TH1D*) dir1->Get(hname.c_str());
					h1->SetDirectory(0);
					TH1D *h2 = (TH1D*) dir2->Get(hname.c_str());
					h2->SetDirectory(0);

					if (!h1 || !h2){
						std::cout << "Error: Histogram found but could not load " << hname << "from TH1 \n" << std::endl;
						exit(1);
					}
		
					// Histogram plotting section.
					style_histogram(target_histogram, h1, kRed, h2, kBlue, set_x_range);
					TCanvas* c = makeCanvas(hname.c_str(), hname.c_str(), log_y);
					h1->Draw("hist f ][");
					h2->Draw("same hist ][");
					gPad->RedrawAxis();
					makeLegend(h1,"8 TeV",h2,"5 TeV")->Draw();

					// Draw CMS header.
					drawCMSHeader();

					c->Modified();
					c->Update();
					TString outname = "../../../mnt/c/Users/lucas/Documents/" + hname + ".png";
					c->SaveAs(outname);

					std::cout << "Integral h1 = " << h1->Integral() << std::endl;
					std::cout << "Integral h2 = " << h2->Integral() << std::endl;

					break; // Histogram found and plotted. Break the while loop of sub-level structure.
				}

				else if (obj1->InheritsFrom(THnSparse::Class())) {

        		THnSparseD *h1 = (THnSparseD*) dir1->Get(hname.c_str());
				THnSparseD *h2 = (THnSparseD*) dir2->Get(hname.c_str());
				if (!h1 || !h2){
					std::cout << "Error: Histogram found but could not load " << hname << "from THnSparse \n" << std::endl;
					exit(1);
					}

        			TH1D *proj1 = (TH1D*)h1->Projection(sparse_axis);  // project axis
        			proj1->SetDirectory(0);
        			TH1D *proj2 = (TH1D*)h2->Projection(sparse_axis);  // project axis
        			proj2->SetDirectory(0);
	        		if (proj1 && proj2) std::cout << "Projections from THnSparse succesfull\n" << std::endl; 
    			
    				// Histogram plotting section.
					style_histogram(target_histogram, proj1, kRed, proj2, kBlue, set_x_range, sparse_axis);
					TCanvas* c = makeCanvas(hname.c_str(), hname.c_str(), log_y);
					proj1->Draw("hist f ][");
					proj2->Draw("same hist ][");
					gPad->RedrawAxis();
					makeLegend(proj1,"8 TeV",proj2,"5 TeV")->Draw();

					// Draw CMS header.
					drawCMSHeader();
					
					c->Modified();
					c->Update();
					TString outname = "../../../mnt/c/Users/lucas/Documents/" + hname + "_axis_" + sparse_axis + ".png";
					c->SaveAs(outname);

					std::cout << "Integral h1 = " << proj1->Integral() << std::endl;
					std::cout << "Integral h2 = " << proj2->Integral() << std::endl;

					break; // Histogram found and plotted. Break the while loop of sub-level structure.
	    		
	    		}

    			else {
        				std::cerr << "Unsupported class: " << classname << std::endl;
        				exit(1);
    			}

			}	

		}// Switching sub-level file.

	if(found){
		std::cout << "Histogram created for " << target_histogram << "\n" << std::endl;
		break; // Histogram found. Breaking while loop of top-level structure.
	}

	}// Switching top-level directory.

	// End of top-level while loop.

	if (!found) std::cerr << "Error: histogram \"" << target_histogram << "\" not found in either directory." << std::endl;
	
	f1->Close();
	f2->Close();
}

// --- Function definitions ---

// Histogram settings.
void style_histogram(const char* target_histogram, TH1D *h1, Color_t color_h1, TH1D *h2, Color_t color_h2, bool set_x_range = false, int sparse_axis = 0){
	if(!h1 || !h2) return; // Avoid null pointer crash.

	// General settings - Only needed for h1.
    h1->GetXaxis()->SetTitle(return_Xaxis_title(target_histogram, sparse_axis));
	
    if (strcmp(target_histogram, "hfSumEtPb") == 0 || strcmp(target_histogram, "hfSumEtp") == 0 || strcmp(target_histogram, "vzhist") == 0){
    	h1->GetYaxis()->SetTitle("Probability density of events");
	}
    else h1->GetYaxis()->SetTitle("Probability density of tracks");

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

	// 8 TeV histogram SETTINGS
	h1->SetStats(0);
	h1->SetLineWidth(1);
    h1->SetLineStyle(1);
	h1->SetFillColor(color_h1 - 9);
    h1->SetLineColor(color_h1 + 2);
    h1->SetFillStyle(1001);

    double int1 = h1->Integral(0, h1->GetNbinsX() + 1);
	if (int1 > 0) h1->Scale(1.0 / int1);
	
	// 5 TeV histogram SETTINGS.
	h2->SetStats(0);
	h2->SetLineWidth(1);
    h2->SetLineStyle(1);
	h2->SetFillColor(color_h2 - 9);
    h2->SetLineColor(color_h2 + 2);
    h2->SetFillStyle(3001);

    double int2 = h2->Integral(0, h2->GetNbinsX() + 1);
	if (int2 > 0) h2->Scale(1.0 / int2);

	// Adjusting X-axis plotting range.
	if(set_x_range){
		int lastBin = GetLastNonZeroBin(h1);
		double xmax = GetLastNonZeroX(h1);
		double xmin = h1->GetXaxis()->GetXmin();
		if( strcmp(target_histogram, "ptresolution") == 0 ) h1->GetXaxis()->SetRangeUser(xmin, xmax + 0.01);
		else h1->GetXaxis()->SetRangeUser(xmin, xmax + 5.0);
	}
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
        return "dxy/#sigma_{dxy}";
    else if (var == "dzoversigmadz")
        return "dz/#sigma_{dz}";
    else if (var == "ptresolution")
        return "#Delta p_{T}/p_{T}";
    else if (var == "hist_reco_trk_corr") {
        if (sparse_axis == 0) return "p_{T} [GeV]";
        if (sparse_axis == 1) return "#eta";
        if (sparse_axis == 2) return "#phi";
        if (sparse_axis == 3) return "Charge [e]";
    }

    return "Unknown variable"; // No known variable was given.
}

void drawCMSHeader(const char* extraText = "#it{Work in Progress}",
                   const char* lumiText  = "pPb (186.0 nb^{#minus1}) 8.16 TeV, (0.509 nb^{#minus1}) 5.02 TeV",
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