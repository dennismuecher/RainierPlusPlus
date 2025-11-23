// main.cpp - Fixed version with discrete level loading
#include "Config.h"
#include "core/Nucleus.h"
#include "simulation/DecaySimulator.h"
#include "io/OutputManager.h"
#include <iostream>
#include <exception>

using namespace rainier;

void printHeader() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                        ║\n";
    std::cout << "║   RAINIER++ - Modern C++ Edition                      ║\n";
    std::cout << "║                                                        ║\n";
    std::cout << "║   Randomizer of Assorted Initial Nuclear              ║\n";
    std::cout << "║   Intensities and Emissions of Radiation              ║\n";
    std::cout << "║                                                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    try {
        printHeader();
        
        // Load configuration
        std::string configFile = (argc > 1) ? argv[1] : "config.dat";
        std::cout << "═══ Loading Configuration ═══\n";
        std::cout << "Configuration from: " << configFile << "\n";
        Config config = Config::loadFromFile(configFile);
        config.print();
        
        // Initialize nucleus
        std::cout << "\n═══ Initializing Nucleus ═══\n";
        Nucleus nucleus(config.nucleus.Z, config.nucleus.A);
        nucleus.setSn(config.nucleus.Sn);
        
        // Load discrete levels from file
        std::cout << "Loading discrete levels from: " << config.nucleus.levelsFile << "\n";
        std::cout << "Maximum levels to read: " << config.continuum.maxDiscreteLevel << "\n";
        nucleus.loadDiscreteLevels(config.nucleus.levelsFile, 
                                   config.continuum.maxDiscreteLevel);
        
        std::cout << "Loaded " << nucleus.getNumDiscreteLevels() 
                  << " discrete levels\n";
        std::cout << "Critical energy (Ecrit): " 
                  << nucleus.getCriticalEnergy() << " MeV\n";
        
        // Create output manager
        int runNumber = 1; // Could be passed as argument
        OutputManager outputMgr(config, runNumber);
        
        // Run simulation
        std::cout << "\n═══ Starting Simulations ═══\n";
        std::cout << "Population mode: ";
        
        switch(config.initialExcitation.mode) {
            case Config::InitialExcitationConfig::Mode::SINGLE:
                std::cout << "SINGLE (E=" << config.initialExcitation.excitationEnergy 
                         << " MeV, J=" << config.initialExcitation.spin << ")\n";
                break;
            case Config::InitialExcitationConfig::Mode::SELECT:
                std::cout << "SELECT (Beta-decay-like)\n";
                break;
            case Config::InitialExcitationConfig::Mode::SPREAD:
                std::cout << "SPREAD (E=" << config.initialExcitation.meanEnergy 
                         << " ± " << config.initialExcitation.energySpread << " MeV)\n";
                break;
            case Config::InitialExcitationConfig::Mode::FULL_REACTION:
                std::cout << "FULL_REACTION (TALYS population)\n";
                break;
        }

        for (int real = 0; real < config.simulation.numRealizations; ++real) {
            std::cout << "Realization " << real + 1 << " / " 
                      << config.simulation.numRealizations << "\n";

            nucleus.buildContinuumLevels(config, real);
            outputMgr.fillLevelSpectra(real, nucleus);
            
            DecaySimulator simulator(nucleus, config, real);

            simulator.run();

            outputMgr.saveRealization(real, simulator);
            std::cout << "  Completed.\n\n";
        }
        
        // Finalize output
        std::cout << "\n═══ Finalizing ═══\n";
        outputMgr.finalize();
        
        std::cout << "\n╔════════════════════════════════════════════════════════╗\n";
        std::cout << "║   Simulation Complete                                 ║\n";
        std::cout << "╚════════════════════════════════════════════════════════╝\n";
        std::cout << "\nOutput files:\n";
        std::cout << "  - " << config.output.outputFile << " (ROOT file)\n";
        std::cout << "  - " << config.output.paramFile << " (parameters)\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n╔════════════════════════════════════════════════════════╗\n";
        std::cerr << "║   ERROR                                               ║\n";
        std::cerr << "╚════════════════════════════════════════════════════════╝\n";
        std::cerr << "\nError: " << e.what() << "\n\n";
        return 1;
    }
}
