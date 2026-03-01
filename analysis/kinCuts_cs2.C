/*
--- Energy cutoff values used by the CMS ref.
(Centrality class)(5.02 TeV)(8.16 TeV)
30 - 80% => 2.5–11.5 GeV, 2.5–14.5 GeV.
1 - 30%  => 11.5–35 GeV, 14.5–44 GeV.
0 - 1%   => >35 GeV, >44 GeV.

---------------Lucas Carvalho---------------
*/

// ### Macro settings ###
int SAMPLE = 1e+8; // Number of generated entries with TF1::GetRandom() following the Hagedorn probability distribution function.
int n_centralities = 3; // Number of centrality classes defined for the dataset.
double lower_pt_forFit = 0.3; // GeV
double upper_pt_forFit = 1.5; // GeV
double delta = 1e-6; // GeV
TString plot_extension = ".pdf"; // .png to test .pdf to final result.
TString final_output_name = "c_s2_vs_Teff_kCuts";
TString base_output_path = "../../../../mnt/c/Users/lucas/Documents/"; // Base path where the outputs will be saved.

// ##############################################################################
// ##############################################################################

// --- Vector to store the different TGraph's created by the different pseudorapidity cuts
std::vector<TGraphErrors*> graphs;

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

// --- Pseudorapidity cuts ---
std::vector<std::pair<double, double>> Eta_cuts = {
	{-2.4, 2.4},
	{-2.0, 2.0},
	{-1.5, 1.5},
	{-0.8, 0.8},
	{-0.5, 0.5}
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
double get_n_events(const std::string& filename, double EHFmin, double EHFmax);
TH1D* get_pT_histogram(const DataStruct& dataFile, double EHFmin = 0.0, double EHFmax = 250.0, double low_eta = -2.5, double high_eta = 2.5);
TH1D* make_hagedorn_fit(TH1D* hist_pT, const DataStruct& dataFile, double EHFmin, double EHFmax);
void plot_graphs();
void draw_CMS_Header(TString latex_text = "#bf{CMS} #it{Work in Progress}", double x = 0.11, double y = 0.93, double text_size = 0.04, int align = 11);

// --- main() ---
void kinCuts_cs2(){

    gROOT->SetBatch(kTRUE); // This tells ROOT to run in batch mode, i.e. no GUI or pop-ups.
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-8);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
	
    TStopwatch timer;
    timer.Start();

    for(size_t k = 0; k < Eta_cuts.size(); ++k){
    	double low_eta = Eta_cuts[k].first;
    	double high_eta = Eta_cuts[k].second;
    	printf("\n\n-> Pseudorapidity window: [%.1f,%.1f] \n", low_eta, high_eta);

        // Deterministic seed for this eta window for reproducibility.
        gRandom->SetSeed(12345 + k);
        printf("Seed for this eta window = %u\n", gRandom->GetSeed());

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

    	for(size_t i = 0; i < n_centralities; ++i){

    		auto [low5, high5] = CBins_5[i];
    		auto [low8, high8] = CBins_8[i];

        	TH1D *hist_5 = nullptr;
        	TH1D *hist_8 = nullptr;
        	double mean_pT, mean_pT_error, n_tracks, n_tracks_error;
        	double n_events, N_ch, N_ch_error;

        	hist_5 = make_hagedorn_fit(get_pT_histogram(dataFile_5TeV, low5, high5, low_eta, high_eta), dataFile_5TeV, low5, high5);

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

            if(low_eta == -1.5){
                printf("\n-> <pT>_5 = %.5f \n", mean_pT);
                printf("\n-> Nch_5 = %.5f \n", N_ch);
            }

        	mean_pT_data_5TeV.push_back(mean_pT);
        	mean_pT_errors_5TeV.push_back(mean_pT_error);
        	N_ch_values_5TeV.push_back(N_ch);
        	N_ch_errors_5TeV.push_back(N_ch_error);
        		
        	hist_8 = make_hagedorn_fit(get_pT_histogram(dataFile_8TeV, low8, high8, low_eta, high_eta), dataFile_8TeV, low8, high8);

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

            if(low_eta == -1.5){
                printf("\n-> <pT>_8 = %.5f \n", mean_pT);
                printf("\n-> Nch_8 = %.5f \n", N_ch);
            }

        	mean_pT_data_8TeV.push_back(mean_pT);
        	mean_pT_errors_8TeV.push_back(mean_pT_error);
        	N_ch_values_8TeV.push_back(N_ch);
        	N_ch_errors_8TeV.push_back(N_ch_error);

    	}// Ending the for loop on i.

    	// --- Vectors to store c_s data
		std::vector<double> cs_results;
		std::vector<double> cs_errors;

    	// Speed of sound fits
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

    }// Ending the for loop on k.

    // Make final plot.
	plot_graphs();

	// Write data on .dat file.
    //print_data();

    timer.Stop();
    std::cout << "-> Job finished in "
              << timer.RealTime() << " seconds (wall time), "
              << timer.CpuTime()  << " seconds (CPU time).\n\n";
}


// --- Function definitions ---
TH1D* make_hagedorn_fit(TH1D* hist_pT, const DataStruct& dataFile, double EHFmin, double EHFmax){

	if(!hist_pT){
        std::cerr << "From " << __func__ << ": error while loading TH1D histogram." << std::endl;
        return nullptr;
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
    	printf("\n-> Fit parameter [1] = 0.4034 \n");
    }
    else if(dataFile.label == "8TeV"){
        pT_fit->FixParameter(1, 0.5010);//related to radial flow velocity - pPb 8TeV
    	printf("\n-> Fit parameter [1] = 0.5010 \n");
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
    //double pt_kinetic_cut = 3.0; //GeV
    //int pt_kinetic_cut_bin = hist_pT->GetXaxis()->FindBin(pt_kinetic_cut - delta);
    //printf("\n-> Kinetic bin = %d \n", pt_kinetic_cut_bin);
    //TH1D *extrapolated_hist = new TH1D("extrapolated_hist", "Extrapolated Histogram", pt_kinetic_cut_bin, 0.,pt_kinetic_cut);
    TH1D *extrapolated_hist = (TH1D*)hist_pT->Clone("extrapolated_hist");    
    extrapolated_hist->SetDirectory(0);

    // Takes first three bin contents from fit_hist. The rest is taken from the real data from hist_pT.
    double content, error;
    for(int i = 1; i <= 3; i++){
        content = fit_hist->GetBinContent(i);
        error = fit_hist->GetBinError(i);
        extrapolated_hist->SetBinContent(i, content);
        extrapolated_hist->SetBinError(i, error);
    }
    extrapolated_hist->ResetStats();

    return extrapolated_hist;
}

TH1D* get_pT_histogram(const DataStruct& dataFile, double EHFmin, double EHFmax, double low_eta, double high_eta){

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

    return hist_pT;
}

void plot_graphs(){

	// Draw them all
    TCanvas* c = new TCanvas("c_finalplot", "Final plot of cs_squared", 800, 600);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.038);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);     // ticks on top x-axis
    c->SetTicky(1);     // ticks on right y-axis
    c->SetFillColor(0);   // white/transparent
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
    
    int colors[5] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1};
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerStyle(21);
        graphs[i]->SetMarkerSize(1.);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetTitle("");
    }

    // Add non-interacting limit dashed line.
    //TLine *line = new TLine(260, 0.33, 280, 0.33);
    //line->SetLineColor(kGray + 2);
    //line->SetLineWidth(3);
    //line->SetLineStyle(9);

    // Settings only needed for first plot.
    graphs[0]->GetXaxis()->CenterTitle(true);
    graphs[0]->GetYaxis()->CenterTitle(true);
    graphs[0]->GetXaxis()->SetTitleOffset(1.);
    graphs[0]->GetYaxis()->SetTitleOffset(1.);
    graphs[0]->GetXaxis()->SetTitleFont(42);
    graphs[0]->GetYaxis()->SetTitleFont(42);
    graphs[0]->GetXaxis()->SetLabelFont(42);
    graphs[0]->GetYaxis()->SetLabelFont(42);
    graphs[0]->GetXaxis()->SetTitleSize(0.044);
    graphs[0]->GetYaxis()->SetTitleSize(0.044);
    graphs[0]->GetXaxis()->SetLabelSize(0.036);
    graphs[0]->GetYaxis()->SetLabelSize(0.036);

    // Axis range and title.
    graphs[0]->GetXaxis()->SetLimits(210, 260);
    graphs[0]->GetYaxis()->SetRangeUser(0.15, 0.38);
    graphs[0]->GetXaxis()->SetTitle("T_{eff} = #LTp_{T}#GT / 3 [MeV]");
    graphs[0]->GetYaxis()->SetTitle("dln #LTp_{T}#GT / dln N_{ch}");

    graphs[0]->Draw("AP");
    for (size_t i = 1; i < graphs.size(); ++i){
        graphs[i]->Draw("PL SAME");
    }
    //line->Draw("SAME");

    auto leg = new TLegend(0.70,0.66,0.90,0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->AddEntry(graphs[0],"|#eta| < 2.4","lep");
    leg->AddEntry(graphs[1],"|#eta| < 2.0","lep");
    leg->AddEntry(graphs[2],"|#eta| < 1.5","lep");
    leg->AddEntry(graphs[3],"|#eta| < 0.8","lep");
    leg->AddEntry(graphs[4],"|#eta| < 0.5","lep");
    leg->Draw();

    draw_CMS_Header();
    //draw_CMS_Header("non-interacting limit", 0.92,0.76, 0.038, 31);
    draw_CMS_Header("p_{T} > 0 GeV", 0.16, 0.21, 0.04, 11);
    draw_CMS_Header("pPb (186.0 nb^{#minus1}) 8.16 TeV", 0.16, 0.84, 0.038, 11);
    draw_CMS_Header("pPb (0.509 nb^{#minus1}) 5.02 TeV", 0.16, 0.79, 0.038, 11);

    c->Modified();
    c->Update();
    // Path to final plot and save output plot.
    TString full_path = base_output_path + final_output_name + plot_extension;
    c->SaveAs(full_path);
}

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