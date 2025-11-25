// main.cpp - RAINIER++ main entry point with YAML config support
#include "Config.h"
#include "core/Nucleus.h"
#include "simulation/DecaySimulator.h"
#include "io/OutputManager.h"
#include <iostream>
#include <string>
#include <stdexcept>

using namespace rainier;

void printBanner() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║   RAINIER++ - Modern C++ Edition                      ║\n";
    std::cout << "║   Nuclear Decay Cascade Monte Carlo Simulation        ║\n";
    std::cout << "║                                                        ║\n";
    std::cout << "║   Now with YAML configuration support!                ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [config_file]\n";
    std::cout << "\n";
    std::cout << "Arguments:\n";
    std::cout << "  config_file  Path to YAML configuration file\n";
    std::cout << "               (default: config/example.yaml)\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << programName << "                      # Use default config\n";
    std::cout << "  " << programName << " my_config.yaml      # Use custom config\n";
    std::cout << "  " << programName << " ../config/test.yaml # Use config from parent dir\n";
    std::cout << "\n";
    std::cout << "Note: YAML format provides:\n";
    std::cout << "  - Native comment support with #\n";
    std::cout << "  - Cleaner, more readable syntax\n";
    std::cout << "  - Better organization of parameters\n";
    std::cout << "\n";
}

int main(int argc, char** argv) {
    printBanner();
    
    // Determine config file
    std::string configFile = "../config/example.yaml";
    if (argc > 1) {
        if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        configFile = argv[1];
    }
    
    try {
        // Load configuration
        std::cout << "Loading configuration from: " << configFile << "\n";
        Config config = Config::loadFromFile(configFile);
        config.print();
        
        // Initialize nucleus
        std::cout << "\n═══ Initializing Nucleus ═══\n";
        Nucleus nucleus(config.nucleus.Z, config.nucleus.A, config);
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
                std::cout << "SELECT (Beta-decay-like with " 
                         << config.initialExcitation.selectStates.size() << " states)\n";
                break;
            case Config::InitialExcitationConfig::Mode::SPREAD:
                std::cout << "SPREAD (E=" << config.initialExcitation.meanEnergy
                         << " ± " << config.initialExcitation.energySpread << " MeV)\n";
                break;
            case Config::InitialExcitationConfig::Mode::FULL_REACTION:
                std::cout << "FULL_REACTION (from " 
                         << config.initialExcitation.populationFile << ")\n";
                break;
        }
        
        std::cout << "\nRealizations: " << config.simulation.numRealizations << "\n";
        std::cout << "Events per realization: " << config.simulation.eventsPerRealization << "\n";
        
        // Run all realizations
        for (int i = 0; i < config.simulation.numRealizations; ++i) {
            std::cout << "\n─── Realization " << (i+1) << "/" 
                     << config.simulation.numRealizations << " ───\n";
            
            // Build continuum levels for this realization
            nucleus.buildContinuumLevels(config, i);
            
            // Create simulator with realization number
            DecaySimulator simulator(nucleus, config, i);
                        
            // Run simulation
            simulator.run();
                        
            // Save output
            outputMgr.saveRealization(i, simulator);
         
        }
        
        // Finalize output
        outputMgr.finalize();
        
        std::cout << "\n═══ Simulation Complete ═══\n";
        std::cout << "Output saved to: " << config.output.outputFile << "\n";
        std::cout << "Parameters saved to: " << config.output.paramFile << "\n";
        std::cout << "\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << "\n\n";
        return 1;
    }
}
