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
class Nucleus;

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
     * @brief Fill level spectrum histograms from generated continuum levels
     * @param realization Realization number
     * @param nucleus Nucleus containing the generated levels
     */
    void fillLevelSpectra(int realization, const Nucleus& nucleus);
    
    /**
     * @brief Create level density histograms showing ρ(E) in levels/MeV
     * @param realization Realization number
     */
    void createLevelDensityHistograms(int realization);

    /**
     * @brief Fill level density histograms from the level density model and discrete levels
     * @param realization Realization number
     * @param simulator DecaySimulator containing the level density model
     */
    void fillLevelDensityHistograms(int realization, const DecaySimulator& simulator);
    /**
     * @brief Finalize and close output file
     */
    void finalize();

private:
    void createHistograms(int realization);
    void createLevelSpectraHistograms(int realization);
    void fillHistograms(int realization, const CascadeEvent& event);
    void saveParameters();
    
    Config config_;
    int runNumber_;
    
    TFile* outputFile_;
    TTree* cascadeTree_;
    
    // Histograms per realization
    std::map<std::string, TH1D*> histograms1D_;
    std::map<std::string, TH2D*> histograms2D_;
    
    // Level spectrum histograms (one per spin value)
    std::map<std::string, TH1D*> levelSpectraHistograms_;
    
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
