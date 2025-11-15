// main.cpp - RAINIER++ main entry point with JSON config support
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
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [config_file]\n";
    std::cout << "\n";
    std::cout << "Arguments:\n";
    std::cout << "  config_file  Path to JSON configuration file (default: config/example.json)\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << programName << "                    # Use default config\n";
    std::cout << "  " << programName << " my_config.json    # Use custom config\n";
    std::cout << "\n";
}

int main(int argc, char** argv) {
    printBanner();
    
    // Determine config file
    std::string configFile = "config/example.json";
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
        std::cout << "═══ Initializing Nucleus ═══\n";
        Nucleus nucleus(config.nucleus.Z, config.nucleus.A, config);
        
        std::cout << "Creating " << nucleus.getNumDiscreteLevels() 
                  << " discrete levels...\n";
        std::cout << "Critical energy (Ecrit): " 
                  << nucleus.getCriticalEnergy() << " MeV\n";
        
        // Create output manager
        int runNumber = 1; // Could be passed as argument
        OutputManager outputMgr(config, runNumber);
        outputMgr.initialize();
        
        // Create simulator
        DecaySimulator simulator(config, nucleus, outputMgr);
        
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
        
        std::cout << "Realizations: " << config.simulation.numRealizations << "\n";
        std::cout << "Events per realization: " << config.simulation.eventsPerRealization << "\n";
        std::cout << "\n";
        
        simulator.run();
        
        // Finalize output
        std::cout << "\n═══ Finalizing ═══\n";
        outputMgr.finalize();
        
        std::cout << "\n╔════════════════════════════════════════════════════════╗\n";
        std::cout << "║   Simulation Complete                                 ║\n";
        std::cout << "╚════════════════════════════════════════════════════════╝\n";
        std::cout << "\n";
        std::cout << "Output files:\n";
        std::cout << "  - " << config.output.outputFile << " (ROOT file)\n";
        std::cout << "  - " << config.output.paramFile << " (parameters)\n";
        std::cout << "\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n╔════════════════════════════════════════════════════════╗\n";
        std::cerr << "║   ERROR                                               ║\n";
        std::cerr << "╚════════════════════════════════════════════════════════╝\n";
        std::cerr << "\nError: " << e.what() << "\n\n";
        return 1;
    }
}
