#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>


int reduce_root_file(const char* input_file,
                     const char* output_file,
                     const char* treename = "DecayTree",
                     double min_mass = 50,
                     double max_mass = 600)
{
    
    TFile* infile = TFile::Open(input_file, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error opening the file " << input_file << std::endl;
        return 1;
    }

    
    TTree* tree = (TTree*)infile->Get(treename);
    if (!tree) {
        std::cerr << "Error finding the tree '" << treename << "'" << std::endl;
        return 2;
    }

    
    tree->SetBranchStatus("*", 0); 
    tree->SetBranchStatus("rho_M", 1);  
    tree->SetBranchStatus("FillNumber", 1);
    tree->SetBranchStatus("RUNNUMBER", 1);

    
    TFile* outfile = new TFile(output_file, "RECREATE");
    if (!outfile || outfile->IsZombie()) {
        std::cerr << "Error creating output file " << output_file << std::endl;
        return 3;
    }

    
    TString cut = Form("rho_M > %f && rho_M < %f", min_mass, max_mass);

    
    TTree* newtree = tree->CopyTree(cut);

    
    newtree->Write();
    outfile->Close();
    infile->Close();

    std::cout << "✅ Saved in: " << output_file << std::endl;
    std::cout << "Events copied: " << newtree->GetEntries() << std::endl;

    return 0;
}


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Use: " << argv[0] << " input_file.root output_file.root" << std::endl;
        return 1;
    }

    const char* input_file = argv[1];
    const char* output_file = argv[2];

    return reduce_root_file(input_file, output_file);
}
