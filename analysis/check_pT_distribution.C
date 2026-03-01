/*
This macro can plot the pT distribution used by cesar and my data to check consistency.
Can plot hagedorn fits.
Used to plot hagedorn fits vs. pT distribution for 0-1% centrality class for both 5TeV and 8TeV datasets, in order to show fit quality.
Final plot also has extrapolated data.
*/

/*--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.*/

double EHFmin = 44.0;
double EHFmax = 250.0;

double delta = 1e-6; // GeV
double low_eta = -1.5; // Low edge of pseudorapidity window.
double high_eta = 1.5; // High edge of pseudorapidity window. The p_T distribution will be integrated over (low_eta, high_eta).
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
int SAMPLE = 1e+8; // Number of generated entries with TF1::GetRandom() following the Hagedorn probability distribution function.

//const std::string myFile = "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root";
const std::string myFile = "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v13-10-02-25_tot.root";
//const std::string cesarFile = "../../fout_5TeV-cs2_HF4eta5_trks1p0_Fit0p3-1p5.root";
const std::string cesarFile = "../../fout_8TeV-cs2_HF4eta5_trks1p0_Fit0p3-1p5.root";

TH1D* plot_pT_cesarData();
TF1* plot_pT_cesarFit();
void plot_pT_myData();

void check_pT_distribution(){

	plot_pT_myData();
}

TF1* plot_pT_cesarFit(){
    TFile *f = TFile::Open(cesarFile.c_str(), "READ");
//    TF1 *fit = (TF1*)f->Get("ptfit_1_500");
//    TF1 *fit = (TF1*)f->Get("ptfit_6_29");
//    TF1 *fit = (TF1*)f->Get("ptfit_71_500"); //5tev
    TF1 *fit = (TF1*)f->Get("ptfit_89_500"); //8tev

    return fit;
}

TH1D* plot_pT_cesarData(){
    TFile *f = TFile::Open(cesarFile.c_str(), "READ");
//    TH1D *hist = (TH1D*)f->Get("hist_pT_1_500");
//    TH1D *hist = (TH1D*)f->Get("clone_hist_1D_pT_fullRange_1_500");
//    TH1D *hist = (TH1D*)f->Get("hist_pT_6_29"); // Supondo que 6-29 corresponda a classe de centralidade 80%.
//    TH1D *hist = (TH1D*)f->Get("hist_pT_6_23"); // Supondo que 6-29 corresponda a classe de centralidade 80%. 5tev
//    TH1D *hist = (TH1D*)f->Get("hist_pT_30_88"); // Supondo que 30-88 corresponda a classe de centralidade 30%.
//    TH1D *hist = (TH1D*)f->Get("hist_pT_71_500"); // Supondo que 71-500 corresponda a classe de centralidade 1%. 5tev
    TH1D *hist = (TH1D*)f->Get("hist_pT_89_500"); // Supondo que 89-500 corresponda a classe de centralidade 1%. 8tev

    return hist;
}

void plot_pT_myData(){
    TFile *f = TFile::Open(myFile.c_str(), "READ");

    TDirectory *dir = (TDirectory*)f->Get("Analysis_histograms");
    
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
    hist_HFSumPb_vs_pt_eta->SetDirectory(0);
    f->Close();

    // Setting pseudorapidity window [low_eta, high_eta].
    int z_min = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(low_eta + delta);
    int z_max = hist_HFSumPb_vs_pt_eta->GetZaxis()->FindBin(high_eta - delta);
    hist_HFSumPb_vs_pt_eta->GetZaxis()->SetRange(z_min, z_max);

    // Projection to TH2 by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *hist_HFSumPb_vs_pt = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");

    // Setting centrality class defined by HF energy cutoffs EHFmin, EHFmax.
    int bin_min = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hist_HFSumPb_vs_pt->GetXaxis()->FindBin(EHFmax - delta);

    // Projects TH2 on TH1 for the defined centrality class.
    TH1D *hist_pT = (TH1D*)hist_HFSumPb_vs_pt->ProjectionY("", bin_min, bin_max);
    hist_pT->SetDirectory(0);

    hist_pT->GetXaxis()->SetRangeUser(0.0,2.5);
    hist_pT->SetStats(0);

    TF1 *cesar_fit = plot_pT_cesarFit();
    TH1D *cesar_hist = plot_pT_cesarData();
    cesar_hist->SetStats(0);

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
    //pT_fit->FixParameter(1, 0.4034);//related to radial flow velocity - pPb 5TeV
    pT_fit->FixParameter(1, 0.5010);//related to radial flow velocity - pPb 8TeV
    pT_fit->SetParLimits(2,0.,0.5);//kinetic freeze-out temperature in GeV
    pT_fit->SetParLimits(3,4.,9.);
    //pT_fit->SetParLimits(3,6.,9.);//n - free parameter no physical meaning - old

    // User-defined Hagedorn function fit.
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-8);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    hist_pT->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    hist_pT->Fit(pT_fit,"NO R EX0 Q","",lower_pt_forFit,upper_pt_forFit);
    TFitResultPtr fitResult = hist_pT->Fit(pT_fit,"NO R EX0 M S","",lower_pt_forFit,upper_pt_forFit);
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    double pValue = TMath::Prob(chi2, ndf);
    std::cout<<"chi2 : "<<chi2<<"; ndf : "<<ndf<<"; pValue : "<<pValue<<std::endl;

    // #############################################
    // --------THE FIT SECTION ENDS HERE----------
    // #############################################

    // Get norm only on the target range: 0.3 to 1.5 GeV.
    int Bin_min = hist_pT->GetXaxis()->FindBin(lower_pt_forFit + delta);
    int Bin_max = hist_pT->GetXaxis()->FindBin(upper_pt_forFit - delta);
    double hist_pT_norm = hist_pT->Integral(Bin_min, Bin_max); // Integral over [0.3,1.5].
    printf("\n-> From make_hagedorn_extrapolation: original_hist norm over [0.3, 1.5] GeV = %.3e \n\n", hist_pT_norm);

    int n_bins = Bin_max; // number of bins between 0 and 1.5 GeV. 15 in this case.
    TH1D *fit_hist = new TH1D("fit_hist", "fit histogram", n_bins, 0., upper_pt_forFit);
    fit_hist->SetDirectory(0);

    int count = 0;
    double fit_hist_norm = 0.0;
        while(count < SAMPLE){
            double x = pT_fit->GetRandom();
            if(x >= 0.3 && x < 1.5) count+=1;
            fit_hist->Fill(x);
        }
    double scale_factor = hist_pT_norm/SAMPLE;
    fit_hist->Scale(scale_factor);

    fit_hist_norm = fit_hist->Integral(Bin_min, Bin_max);
    printf("\n-> Fit histogram complete. Final fit histogram norm over [0.3, 1.5] GeV = %.3e \n", fit_hist_norm);

    // Keep only the first three bins of fit_hist
int n_bins_to_keep = 3;
TH1D *fit_hist_first3 = (TH1D*)fit_hist->Clone("fit_hist_first3");
for (int i = n_bins_to_keep + 1; i <= fit_hist_first3->GetNbinsX(); ++i) {
    fit_hist_first3->SetBinContent(i, 0);
    fit_hist_first3->SetBinError(i, 0);
}

// --------DRAWING SECTION------

TCanvas *c = new TCanvas("c", "pT comparison", 800, 600);
c->SetLeftMargin(0.12);
c->SetRightMargin(0.035);
c->SetBottomMargin(0.12);
c->SetTopMargin(0.08);

c->SetFillColor(0);   // white/transparent
c->SetFrameFillColor(0);
c->SetFrameLineWidth(2);

//My hist and fit
hist_pT->SetMarkerStyle(20);
hist_pT->SetMarkerSize(0.8);
hist_pT->SetMarkerColor(kBlack);
hist_pT->SetLineColor(kBlack);
hist_pT->SetStats(0);

fit_hist_first3->SetMarkerStyle(20);
fit_hist_first3->SetMarkerSize(0.8);
fit_hist_first3->SetMarkerColor(kGreen + 1);
fit_hist_first3->SetLineColor(kGreen + 1);
fit_hist_first3->SetStats(0);

pT_fit->SetLineColor(kRed);
pT_fit->SetLineWidth(2);
pT_fit->SetLineStyle(8); // dashed line
//pT_fit->SetLineStyle(1); // dashed line
pT_fit->SetTitle("");

pT_fit->GetXaxis()->CenterTitle(true);
pT_fit->GetYaxis()->CenterTitle(true);
pT_fit->GetXaxis()->SetTitleOffset(1.0);
pT_fit->GetYaxis()->SetTitleOffset(1.1);
pT_fit->GetXaxis()->SetTitleFont(42);
pT_fit->GetYaxis()->SetTitleFont(42);
pT_fit->GetXaxis()->SetLabelFont(42);
pT_fit->GetYaxis()->SetLabelFont(42);
pT_fit->GetXaxis()->SetTitleSize(0.05);
pT_fit->GetYaxis()->SetTitleSize(0.044);
pT_fit->GetXaxis()->SetLabelSize(0.036);
pT_fit->GetYaxis()->SetLabelSize(0.036);

pT_fit->GetXaxis()->SetTitle("p_{T} [GeV]");
pT_fit->GetYaxis()->SetTitle("Number of tracks");

//Cesar hist and fit.
cesar_hist->SetMarkerStyle(4);
cesar_hist->SetMarkerSize(0.8);
cesar_hist->SetMarkerColor(kGreen);
cesar_hist->SetLineColor(kGreen);
cesar_hist->SetStats(0);

cesar_fit->SetLineColor(kBlue);
cesar_fit->SetLineWidth(2);
cesar_fit->SetLineStyle(2); // solid line


pT_fit->Draw();     // Your fitted function
hist_pT->Draw("E1 SAME");      // Draw histogram with error bars
//cesar_fit->Draw("SAME");   // Cesar’s TF1
//cesar_hist->Draw("E1 SAME");      // Draw histogram with error bars
fit_hist_first3->Draw("E1 SAME");

// --- Add legend ---
TLegend *leg = new TLegend(0.65, 0.65, 0.93, 0.81);
leg->AddEntry(hist_pT, "CMS data", "lep");
leg->AddEntry(fit_hist_first3, "p_{T} < 0.3", "lep");
leg->AddEntry(pT_fit, "Fit", "l");
//leg->AddEntry(cesar_hist, "Cesar data", "lep");
//leg->AddEntry(cesar_fit, "Cesar Fit", "l");
leg->SetBorderSize(0);
leg->SetTextSize(0.035);
leg->Draw();

TLatex latex;
latex.SetNDC();              // use normalized coordinates
latex.SetTextSize(0.04);     // text size
latex.SetTextFont(42);       // Helvetica-like
latex.SetTextAlign(11);      // left-aligned, top

TString cmsText = "#bf{CMS} #it{Work in Progress}";
latex.DrawLatex(0.19, 0.93, cmsText);

TLatex latex2;
latex2.SetNDC();
latex2.SetTextSize(0.038);
latex2.SetTextFont(42);
latex2.SetTextAlign(11);
latex2.DrawLatex(0.66, 0.58, "|#eta| < 1.5, 0-1%");

TLatex latex3;
latex3.SetNDC();
latex3.SetTextSize(0.035);
latex3.SetTextFont(42);
latex3.SetTextAlign(31);
//latex3.DrawLatex(0.92, 0.93, "pPb (0.509 nb^{#minus1}) 5.02 TeV");
latex3.DrawLatex(0.92, 0.93, "pPb (186.0 nb^{#minus1}) 8.16 TeV");

//c->SaveAs("../../../../mnt/c/Users/lucas/Documents/myFull_fit_comparison.pdf");
//c->SaveAs("../../../../mnt/c/Users/lucas/Documents/observed_plus_fit_plus_extrap_5TeV.pdf");
c->SaveAs("../../../../mnt/c/Users/lucas/Documents/observed_plus_fit_plus_extrap_8TeV.pdf");
//c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_5TeV_6_29_fit_comparison.pdf");
//c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_30_88_fit_comparison.pdf");
//c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_89_500_fit_comparison_5TeV.pdf");
//c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_89_500_fit_comparison_8TeV.pdf");

    /*// --- Style settings ---
    hist_pT->SetTitle("6_29 centrality class");
    //hist_pT->SetTitle("30_88 centrality class");
    //hist_pT->SetTitle("89_500 centrality class");
    hist_pT->SetLineColor(kBlue);
    hist_pT->SetMarkerStyle(20);
    hist_pT->SetMarkerColor(kBlue);
    hist_pT->SetMarkerSize(0.8);

    hist->SetLineColor(kRed);
    hist->SetMarkerStyle(4);
    hist->SetMarkerColor(kRed);
    hist->SetMarkerSize(0.8);

    // --- Draw on the same canvas ---
    TCanvas *c = new TCanvas("c", "pT comparison", 800, 600);
    //c->SetLogy();
    //c->SetLogx();

    hist_pT->Draw("E1");      // Draw with error bars
    hist->Draw("E1 SAME");    // Draw over it

    // --- Add legend for clarity ---
    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.85);
    leg->AddEntry(hist_pT, "My data", "lep");
    leg->AddEntry(hist, "Cesar data", "lep");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->Draw();

    // --- Axis titles ---
    hist_pT->GetXaxis()->SetTitle("p_{T}");
    hist_pT->GetYaxis()->SetTitle("Counts");
    hist_pT->GetXaxis()->SetTitleOffset(1.4);*/

    //c->SaveAs("../../../../mnt/c/Users/lucas/Documents/myFull_pT_hist.pdf");
    //c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_pT_hist.pdf");
    //c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_6_29_pT_hist.pdf");
    //c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_30_88_pT_hist.pdf");
    //c->SaveAs("../../../../mnt/c/Users/lucas/Documents/my_89_500_pT_hist.pdf");

}