// Macro developed to plot the HF energy vs eta distribution.
// This work is part of the CMS Collaboration and uses CMS Preliminary Data.

/*
--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.*/

double low_pT = 0.3;
double high_pT = 5.0;
double delta = 1e-6;
double EHFmin = 0.0;
double EHFmax = 250.0;

TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

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

void hfSumEtPb_vs_eta(){

    gROOT->SetBatch(kTRUE); // This tells ROOT to run in batch mode, i.e. no GUI or pop-ups.

    TFile *file = TFile::Open(dataFile_8TeV.path.c_str(), "READ");

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    
    // Load TH3 histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    file->Close();

    // Setting pT window [low_pT, high_pT].
    int y_min = hist_HFSumPb_vs_pt_eta->GetYaxis()->FindBin(low_pT + delta);
    int y_max = hist_HFSumPb_vs_pt_eta->GetYaxis()->FindBin(high_pT - delta);
    hist_HFSumPb_vs_pt_eta->GetYaxis()->SetRange(y_min, y_max);

    // Projection to TH2 by integrating on pT window [low_pT, high_pT].
    TH2D *h2 = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("xz");
    h2->SetDirectory(0);

    //Canvas and plot SETTINGS
    TCanvas *c1 = new TCanvas("c_cent", "Et vs eta", 1200, 900);
    c1->SetLeftMargin(0.09);
    c1->SetRightMargin(0.14);
    c1->SetBottomMargin(0.12);
    c1->SetTopMargin(0.08);
    c1->SetFillColor(0);   // white/transparent
    c1->SetFrameFillColor(0);
    c1->SetFrameLineWidth(2);
    c1->SetLogz();
    //c1->SetLogy();

    h2->SetTitle("");
    h2->GetXaxis()->SetTitleOffset(1.);
    h2->GetYaxis()->SetTitleOffset(0.9);
    h2->GetZaxis()->SetTitleOffset(0.9);
    h2->GetXaxis()->CenterTitle(true);
    h2->GetYaxis()->CenterTitle(true);
    h2->GetZaxis()->CenterTitle(true);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetYaxis()->SetLabelFont(42);
    h2->GetZaxis()->SetTitleFont(42);
    h2->GetZaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleSize(0.044);
    h2->GetYaxis()->SetTitleSize(0.044);
    h2->GetZaxis()->SetTitleSize(0.044);
    h2->GetXaxis()->SetLabelSize(0.036);
    h2->GetYaxis()->SetLabelSize(0.036);
    h2->GetZaxis()->SetLabelSize(0.036);
    h2->GetXaxis()->SetTickLength(0.012);
    h2->GetYaxis()->SetTickLength(0.012);
    h2->GetZaxis()->SetTickLength(0.015);
    h2->SetStats(0);
//    h2->SetMinimum(0.5);
    h2->GetYaxis()->SetTitle("E_{HF} [GeV]");
    h2->GetXaxis()->SetTitle("#eta [GeV]");
    h2->GetZaxis()->SetTitle("N of tracks");
    //h2->GetYaxis()->SetRangeUser(0.4, h2->GetYaxis()->GetXmax());
    h2->GetXaxis()->SetRangeUser(-2.6, 2.6);
    h2->Draw("COLZ");

    c1->SaveAs("../../../../mnt/c/Users/lucas/Documents/Et_vs_eta.png");

    // Setting centrality class defined by HF energy cutoffs EHFmin, EHFmax.
    int bin_min = h2->GetYaxis()->FindBin(EHFmin + delta);
    int bin_max = h2->GetYaxis()->FindBin(EHFmax - delta);

    // If one wants to check selected bins for centrality class.
    double lowEdge = h2->GetYaxis()->GetBinLowEdge(bin_min);
    double highEdge = h2->GetYaxis()->GetBinUpEdge(bin_min);
    double lowEdge_ = h2->GetYaxis()->GetBinLowEdge(bin_max);
    double highEdge_ = h2->GetYaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From get_proj_hist: Integrating from bin [%.1f,%.1f]GeV to bin [%.1f,%.1f]GeV \n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Projects TH2 on TH1 for the defined centrality class.
    TH1D *hist_eta = (TH1D*)h2->ProjectionX("", bin_min, bin_max);
    hist_eta->SetDirectory(0);

    TCanvas *c2 = new TCanvas("C", "eta", 800, 600);
    c2->SetLeftMargin(0.09);
    c2->SetRightMargin(0.14);
    c2->SetBottomMargin(0.12);
    c2->SetTopMargin(0.08);
    c2->SetFillColor(0);   // white/transparent
    c2->SetFrameFillColor(0);
    c2->SetFrameLineWidth(2);

    hist_eta->SetTitle("");
    hist_eta->GetXaxis()->SetTitleOffset(1.);
    hist_eta->GetYaxis()->SetTitleOffset(0.9);
    hist_eta->GetXaxis()->CenterTitle(true);
    hist_eta->GetYaxis()->CenterTitle(true);
    hist_eta->GetXaxis()->SetTitleFont(42);
    hist_eta->GetXaxis()->SetLabelFont(42);
    hist_eta->GetYaxis()->SetTitleFont(42);
    hist_eta->GetYaxis()->SetLabelFont(42);
    hist_eta->GetXaxis()->SetTitleSize(0.044);
    hist_eta->GetYaxis()->SetTitleSize(0.044);
    hist_eta->GetXaxis()->SetLabelSize(0.036);
    hist_eta->GetYaxis()->SetLabelSize(0.036);
    hist_eta->GetXaxis()->SetTickLength(0.012);
    hist_eta->GetYaxis()->SetTickLength(0.012);
    hist_eta->SetStats(0);
    hist_eta->GetXaxis()->SetTitle("#eta [GeV]");
    hist_eta->GetYaxis()->SetTitle("N of tracks");
    hist_eta->GetXaxis()->SetRangeUser(-2.6, 2.6);
 
    hist_eta->SetStats(0);
	hist_eta->SetLineWidth(1);
    hist_eta->SetLineStyle(1);
    hist_eta->SetFillStyle(1001);

    hist_eta->Draw("hist");

    c2->SaveAs("../../../../mnt/c/Users/lucas/Documents/eta.png");
}