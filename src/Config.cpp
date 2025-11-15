#include "Config.h"
#include <iostream>

namespace rainier {

Config Config::loadFromFile(const std::string& filename) {
    std::cout << "Note: JSON parsing not yet implemented\n";
    return createDefault(60, 144);
}

void Config::saveToFile(const std::string& filename) const {
    std::cout << "Note: Config saving not yet implemented\n";
}

Config Config::createDefault(int Z, int A) {
    Config config;
    config.nucleus.Z = Z;
    config.nucleus.A = A;
    config.nucleus.Sn = 7.0;
    config.simulation.numRealizations = 1;
    config.simulation.eventsPerRealization = 100;
    return config;
}

bool Config::validate(std::string& errorMessage) const {
    if (nucleus.Z <= 0 || nucleus.A <= 0) {
        errorMessage = "Invalid Z or A";
        return false;
    }
    return true;
}

void Config::print() const {
    std::cout << "\n=== Configuration ===\n";
    std::cout << "Nucleus: Z=" << nucleus.Z << " A=" << nucleus.A << "\n";
    std::cout << "Realizations: " << simulation.numRealizations << "\n";
    std::cout << "Events: " << simulation.eventsPerRealization << "\n\n";
}

} // namespace rainier
