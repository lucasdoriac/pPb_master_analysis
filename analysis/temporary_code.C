// Test code to calculate mean pT and Nch at the particular case of 30-80% on 8TeV dataset.

// ### Macro settings ###
int SAMPLE = 1e+8; // Number of generated entries with TF1::GetRandom() following the Hagedorn probability distribution function.
double low_eta = -1.5; // Low edge of pseudorapidity window.
double high_eta = 1.5; // High edge of pseudorapidity window. The p_T distribution will be integrated over (low_eta, high_eta).
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-5; // GeV
double EHFmin = 2.5; // GeV
double EHFmax = 14.5; // GeV


double get_n_events();
void get_pT_histogram();

void temporary_code(){

    get_pT_histogram();

}


double get_n_events(){

	TFile* f = TFile::Open("../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root", "READ");

    TDirectory *dir = (TDirectory*)f->Get("QA_histograms");

    TH1D *hf_hist = (TH1D*)dir->Get("hfSumEtPb");
    hf_hist->SetDirectory(0);
    f->Close();

    int bin_min = hf_hist->GetXaxis()->FindBin(EHFmin + delta); // EHFmin = 2.5
    int bin_max = hf_hist->GetXaxis()->FindBin(EHFmax - delta); // EHFmax = 14.5

    for (int i = 1; i <= hf_hist->GetNbinsX(); ++i) {
    double xlow  = hf_hist->GetXaxis()->GetBinLowEdge(i);
    double xhigh = hf_hist->GetXaxis()->GetBinUpEdge(i);
    std::cout << " Bin " << i << ": [" << xlow << ", " << xhigh << "]" << std::endl;
}

    // Optionally, print the upper edge of the last bin explicitly
    std::cout << "-> Last upper edge: " << hf_hist->GetXaxis()->GetBinUpEdge(hf_hist->GetNbinsX()) << std::endl;

    std::cout << "-> Evaluating n_events over [" << EHFmin << "," << EHFmax << "]" << std::endl; 

    // Calculates the number of events n_events detected by the Pb side HF, for the specified centrality class.
    double n_events = hf_hist->Integral(bin_min, bin_max);

    return n_events;
}

void get_pT_histogram(){

    TFile *file = TFile::Open("../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root", "READ");

    TDirectory *dir = (TDirectory*)file->Get("Analysis_histograms");
    
    // Load TH3 histogram from "Analysis_histograms/hist_HFSumPb_vs_pt_eta".
    TH3D *hist_HFSumPb_vs_pt_eta = (TH3D*)dir->Get("hist_HFSumPb_vs_pt_eta");
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

    // Assume you already have a TH3D* named h3
    int nBinsZ = hist_HFSumPb_vs_pt_eta->GetNbinsZ();

    // Print lower edges of each Z bin
    for (int iz = 1; iz <= nBinsZ; ++iz) {
        double zlow  = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinLowEdge(iz);
        double zhigh = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(iz);
        std::cout << "Bin " << iz << ": [" << zlow << ", " << zhigh << "]" << std::endl;
    }

// If you also want the upper edge of the *last* bin:
double zmax = hist_HFSumPb_vs_pt_eta->GetZaxis()->GetBinUpEdge(nBinsZ);
std::cout << "Last upper edge: " << zmax << std::endl;

    // Projection to TH2 by integrating on pseudorapidity window [low_eta, high_eta].
    TH2D *hist_HFSumPb_vs_pt = (TH2D*) hist_HFSumPb_vs_pt_eta->Project3D("yx");

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


    int b_min = hist_pT->GetXaxis()->FindBin(lower_pt_forFit + delta);
    int b_max = hist_pT->GetXaxis()->FindBin(upper_pt_forFit - delta);
    double hist_pT_norm = hist_pT->Integral(b_min, b_max); // Integral over [0.3,1.5].
    printf("\n-> From make_hagedorn_extrapolation: original_hist norm over [0.3, 1.5] GeV = %.3e \n\n\n", hist_pT_norm);


    // #############################################
    // --------THE FIT SECTION STARTS HERE----------
    // #############################################

    // Hagedorn TF1. Function declaration section.
    TF1* pT_fit;
    
    // Original.
	pT_fit = new TF1("pt_fit","[0]*x*pow(1.+1./sqrt(1.-[1]*[1])*(sqrt(x*x+[4]*[4])-x*[1])/[3]/[2],-[3])",0.,upper_pt_forFit);
    pT_fit->SetParameters(7500000000.,0.3,0.1,6.,0.14);//We used these values for initialization        
    pT_fit->FixParameter(4,0.13957);//pion mass
    pT_fit->FixParameter(1,0.5010);//related to radial flow velocity - pPb 8TeV
    pT_fit->SetParLimits(2,0.,0.5);//kinetic freeze-out temperature in GeV
    pT_fit->SetParLimits(3,6.,9.);//n - free parameter no physical meaning


    // User-defined Hagedorn function fit.
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


    // Create the histogram fit_hist that will be filled with a random generator of tracks following the Hagedorn distribution defined by the fit.
    int n_bins = b_max; // number of bins between 0 and 1.5 GeV. 15 in this case.
    printf("\n-> number of bins in fit_hist must be 15. now is = %d \n", n_bins);
    TH1D *fit_hist = new TH1D("fit_hist", "fit histogram", n_bins, 0., upper_pt_forFit);
    fit_hist->SetDirectory(0);

    printf("\n-> Constant bin width in fit_hist = %.3f \n", fit_hist->GetXaxis()->GetBinWidth(1));
    for (int i = 1; i <= fit_hist->GetNbinsX(); ++i) {
    double xlow  = fit_hist->GetXaxis()->GetBinLowEdge(i);
    double xhigh = fit_hist->GetXaxis()->GetBinUpEdge(i);
    std::cout << " Bin " << i << ": [" << xlow << ", " << xhigh << "]" << std::endl;
}


    int count = 0;
    double fit_hist_norm = 0.0;
    printf("\n\n-> Random generation of tracks with Hagedorn PDF. Sample size = %.1e \n", static_cast<double>(SAMPLE));
        while(count < SAMPLE){
            double x = pT_fit->GetRandom();
            if(x >= 0.3 && x < 1.5) count+=1;
            fit_hist->Fill(x);
        }
    double scale_factor = hist_pT_norm / SAMPLE;
    printf("\n\n-> Scale factor = %.3f \n\n", scale_factor);
    fit_hist->Scale(scale_factor);

    fit_hist_norm = fit_hist->Integral(b_min, b_max);
    printf("\n-> Hagedorn histogram complete. Final hagedorn histogram norm over [0.3, 1.5] GeV = %.3e \n", fit_hist_norm);

    // Fill extrapolated tracks into histogram extrapolated_hist.
    TH1D *extrapolated_hist = (TH1D*)hist_pT->Clone("extrapolated_hist");
    extrapolated_hist->SetDirectory(0);
    extrapolated_hist->ResetStats();
    
    // Takes first three bin contents from fit_hist. The rest is taken from the real data from original_hist.
    double content, error;
    for(int i = 1; i <= 3; ++i){
            content = fit_hist->GetBinContent(i);
            error = fit_hist->GetBinError(i);
            extrapolated_hist->SetBinContent(i, content);
            extrapolated_hist->SetBinError(i, error);
    }

    double mean_pT = extrapolated_hist-> GetMean();
    double mean_pT_error = extrapolated_hist->GetMeanError();

    double n_tracks_error = 0.0;
    double n_tracks = extrapolated_hist->IntegralAndError(1, extrapolated_hist->GetNbinsX(), n_tracks_error);

    double n_events = get_n_events();
    double N_ch = (n_tracks/n_events);
    double N_ch_error = (n_tracks_error/n_events);

    printf("\n-> --- Results for 30-80%% centrality class for 8 TeV dataset! --- \n");
    printf("\n-> <pT> = %.5f \n", mean_pT);
    printf("\n-> <pT> error = %.5f \n", mean_pT_error);
    printf("\n-> N_ch = %.5f \n", N_ch);
    printf("\n-> N_ch error = %.5f \n", N_ch_error);

    return;
}
