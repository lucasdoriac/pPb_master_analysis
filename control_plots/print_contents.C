// Prints the content of both datafiles.

// --- Headers
#include <iostream> // std::cout, std::cerr
#include <fstream> // std::ofstream
#include <TFile.h> // TFile
#include <TDirectory.h> // TDirectory
#include <TClass.h> // TObject::InheritsFrom()
#include <TKey.h> // TKey
#include <string> // std::string


// --- Path to data files
const std::string data_File_5TeV = "../../pPb_meanpT_vs_Nch_histos_5TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root";
const std::string data_File_8TeV = "../../pPb_meanpT_vs_Nch_histos_8TeV_MBonly_PUGPlus_HFSumEtEta4_TrkEta2p4_v12-09-01-25_tot.root";

// --- Function prototypes
void print_content(const std::string& datafile);
void logError(const std::string& message);

// --- main() ---
void print_contents(){

    std::cout << "\n 5 TeV dataset: " << std::endl;
    print_content(data_File_5TeV);

    std::cout << "\n";

    std::cout << "\n 8 TeV dataset: " << std::endl;
    print_content(data_File_8TeV);
}

// --- Function definitions
void print_content(const std::string& datafile) {

    TFile *file = TFile::Open(datafile.c_str(), "READ");
    if(!file || file->IsZombie()) {
        logError("Error: Could not open " + datafile);
        return;
    }

    std::cout << "--- Contents of file ---" << std::endl;

    TIter nextKey(file->GetListOfKeys());
    TKey *key;
    while ( (key = (TKey*)nextKey()) ){

        TObject *obj = key->ReadObj();
        if(!obj){
            logError("From print_content: Error getting TObject inside " + datafile);
            return;
        }
        // Print current object's Name and Class.
        std::cout << "Name: " << obj->GetName()
                  << " | Class: " << obj->ClassName() << std::endl;

        // If it's a directory, also list its contents.
        if (obj->InheritsFrom(TDirectory::Class())){

            TDirectory *dir = (TDirectory*)obj;
            if(!dir){
                logError("From print_content: Error while getting TDirectory class in" + datafile);
                return;
            }
            
            TIter nextDirKey(dir->GetListOfKeys());
            TKey *subKey;
            while ((subKey = (TKey*)nextDirKey())){

                // Sub-object located in the directory.
                TObject *sub_obj = subKey->ReadObj();
                if(!sub_obj){
                    logError("From print_content: Error while read sub-object inside " + datafile);
                    continue;
                }
                std::cout << "   -> " << sub_obj->GetName()
                          << " | Class: " << sub_obj->ClassName() << std::endl;
            }
        }
    }

    file->Close();
}

void logError(const std::string& message){
    // Set off stream to error.log file.
    std::ofstream err_log("print_contents.log", std::ios::app);

    if (err_log.is_open()){
        err_log << message << std::endl;
    }
    else {
        std::cerr << "From logError: Could not open 'print_contents.log' for writing." << std::endl;
    }
    
}