// Macro to test getbin helpers and root functions.

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

// --- Global variables
const char* target_hist = "hfSumEtPb";

void logError(const std::string& message);

void test(){

    gROOT->SetBatch(kTRUE); // This tells ROOT to run in Batch mode, i. e. no GUI or pop-ups.

    TFile *file = TFile::Open(data_File_5TeV.path.c_str(), "READ");
    if (!file || file->IsZombie()){
        logError("Error: from get_th1d: Could not open file " + data_File_5TeV.path);
        return;
    }

    TDirectory *dir = (TDirectory*)file->Get("QA_histograms");
    if (!dir) {
        logError("Error: from get_th1d: Could not load directory 'QA_histograms'.");
        file->Close();
        return;
    }

    TH1D* hist = (TH1D*)dir->Get(target_hist);
    if (!hist){
        logError("Error: from get_th1d: Could not load target histogram hist.");
        return;
    }
    hist->SetDirectory(0);
    file->Close();

    int n_bins = hist->GetNbinsX();

    double underflow_content = hist->GetBinContent(0);
    double overflow_content = hist->GetBinContent(n_bins+1);
    
    printf("Underflow content = %.f \n", underflow_content);
    printf("Overflow content = %.f \n", overflow_content);

    int i = 5;

    double bin_content = hist->GetBinContent(i);
    double low_edge = hist->GetXaxis()->GetBinLowEdge(i);
    double up_edge = hist->GetXaxis()->GetBinUpEdge(i);

    printf("Bin %d content = %f \n", bin_content);
    printf("Bin %d lower edge = %f \n", low_edge);
    printf("Bin %d upper edge = %f \n", up_edge);

//    std::cout << "low edge = " << low_edge << ", up edge = " << up_edge << std::endl;
//    std::cout << "bin content = " << bin_content << std::endl;

//    std::cout << "integral without anything" << hist->Integral() << "Integral with explicit bins" << hist->Integral(0,250) << std::endl;


}


void logError(const std::string& message) {
    // Set off stream to error.log file.
    std::ofstream err_log("test.log", std::ios::app);

    if (err_log.is_open()){
        err_log << message << std::endl;
    }
    else {
        std::cerr << "From logError: Could not open 'test.log' for writing." << std::endl;
    }
    
}