#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <vector>
#include <filesystem>

void mergeTTrees() {
    namespace fs = std::filesystem;
    TFile *outputFile = new TFile("merged.root", "RECREATE");

    std::vector<std::string> treeNames = {"SumProballTracks", "ThisTrack", "OtherTracks",
                                          "HighChargeClusters", "ClusterCandidates", "McTruth"};

    // Maps to hold original trees and their clones in the output file
    std::map<std::string, TTree*> originalTrees;
    std::map<std::string, TTree*> clonedTrees;

    // Iterate over each tree name to initialize cloned trees
    for (const auto& name : treeNames) {
        originalTrees[name] = new TChain(name.c_str());
        clonedTrees[name] = nullptr;  // Initialize with nullptr to be cloned later
    }

    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
        std::string filename = entry.path().string();
        if (filename.find("sigma") != std::string::npos && entry.path().extension() == ".root") {
            for (auto& name : treeNames) {
                ((TChain*)originalTrees[name])->Add(filename.c_str());
            }
        }
    }

    // Clone trees and copy entries
    for (const auto& name : treeNames) {
        if (originalTrees[name]->GetEntries() > 0) {
            outputFile->cd();                                        // Ensure we are in the right directory of the output file
            clonedTrees[name] = originalTrees[name]->CloneTree(0);   // Clone structure
            clonedTrees[name]->CopyEntries(originalTrees[name]);     // Copy all entries
            clonedTrees[name]->Write();                              // Write the cloned tree to the output file
        }
    }

    outputFile->Close();
    std::cout << "All trees have been merged and saved into merged.root." << std::endl;

    // Cleanup: delete original trees
    for (auto& pair : originalTrees) {
        delete pair.second;
    }
}
