// OutputManager.h - ROOT output and histogram management
#ifndef RAINIER_OUTPUT_MANAGER_H
#define RAINIER_OUTPUT_MANAGER_H

#include "Config.h"
#include <string>
#include <map>
#include <memory>

// Forward declarations for ROOT
class TFile;
class TH1D;
class TH2D;
class TTree;

namespace rainier {

// Forward declarations
class DecaySimulator;
class CascadeEvent;

/**
 * @brief Manages ROOT output files and histograms
 */
class OutputManager {
public:
    OutputManager(const Config& config, int runNumber);
    ~OutputManager();
    
    /**
     * @brief Save results from one realization
     */
    void saveRealization(int realization, const DecaySimulator& simulator);
    
    /**
     * @brief Finalize and close output file
     */
    void finalize();

private:
    void createHistograms(int realization);
    void fillHistograms(int realization, const CascadeEvent& event);
    void saveParameters();
    
    Config config_;
    int runNumber_;
    
    TFile* outputFile_;
    TTree* cascadeTree_;
    
    // Histograms per realization
    std::map<std::string, TH1D*> histograms1D_;
    std::map<std::string, TH2D*> histograms2D_;
    
    // Tree branches
    std::vector<double> treeGammaEnergies_;
    std::vector<double> treeFinalEnergies_;
    std::vector<double> treeTimes_;
    double treeInitialEx_;
    int treeInitialSpin_;
    int treeInitialParity_;
};

} // namespace rainier

#endif // RAINIER_OUTPUT_MANAGER_H