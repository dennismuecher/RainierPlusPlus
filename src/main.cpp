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
    std::cout << "║   RAINIER 2.0 - Modern C++ Edition                    ║\n";
    std::cout << "║                                                        ║\n";
    std::cout << "║   Randomizer of Assorted Initial Nuclear              ║\n";
    std::cout << "║   Intensities and Emissions of Radiation              ║\n";
    std::cout << "║                                                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    printHeader();

    std::string configFile = "config.json";
    int runNumber = 1;

    if (argc > 1) configFile = argv[1];
    if (argc > 2) runNumber = std::stoi(argv[2]);

    try {
        std::cout << "Loading configuration from: " << configFile << "\n";
        Config config;
        try {
            config = Config::loadFromFile(configFile);
        } catch (const std::exception&) {
            std::cout << "Config file not found. Using defaults.\n";
            config = Config::createDefault(60, 144);
        }

        std::string errorMsg;
        if (!config.validate(errorMsg)) {
            std::cerr << "Config validation failed: " << errorMsg << "\n";
            return 1;
        }

        config.print();

        std::cout << "=== Initializing Nucleus ===\n";
        Nucleus nucleus(config.nucleus.Z, config.nucleus.A);
        nucleus.setSn(config.nucleus.Sn);

        nucleus.loadDiscreteLevels(config.nucleus.levelsFile, 20);
        nucleus.printSummary();

        OutputManager output(config, runNumber);

        std::cout << "\n=== Starting Simulations ===\n";
        std::cout << "Realizations: " << config.simulation.numRealizations << "\n";
        std::cout << "Events per realization: " << config.simulation.eventsPerRealization << "\n\n";

        for (int real = 0; real < config.simulation.numRealizations; ++real) {
            std::cout << "Realization " << real + 1 << " / " 
                      << config.simulation.numRealizations << "\n";

            nucleus.buildContinuumLevels(config, real);

            DecaySimulator simulator(nucleus, config, real);
            simulator.run();

            output.saveRealization(real, simulator);
            std::cout << "  Completed.\n\n";
        }

        output.finalize();

        std::cout << "=== Simulation Complete ===\n";
        std::cout << "This is a STUB implementation - physics models need completion\n\n";

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
