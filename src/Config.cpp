// Config.cpp - Enhanced configuration with JSON support
#include "Config.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace rainier {

// Simple JSON parser (avoids external dependencies)
// For production use, consider nlohmann/json or similar
namespace {
    std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(" \t\n\r");
        if (first == std::string::npos) return "";
        size_t last = str.find_last_not_of(" \t\n\r");
        return str.substr(first, last - first + 1);
    }
    
    std::pair<std::string, std::string> parseLine(const std::string& line) {
        size_t colon = line.find(':');
        if (colon == std::string::npos) return {"", ""};
        
        std::string key = trim(line.substr(0, colon));
        std::string value = trim(line.substr(colon + 1));
        
        // Remove quotes and commas
        if (!key.empty() && key.front() == '"') key = key.substr(1, key.length() - 2);
        if (!value.empty() && value.back() == ',') value = value.substr(0, value.length() - 1);
        if (!value.empty() && value.front() == '"' && value.back() == '"') 
            value = value.substr(1, value.length() - 2);
            
        return {key, value};
    }
}

Config Config::loadFromFile(const std::string& filename) {
    Config config;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open config file: " << filename << std::endl;
        std::cerr << "Using default configuration.\n";
        return createDefault(60, 144);
    }
    
    std::cout << "Loading configuration from: " << filename << std::endl;
    
    std::string line;
    std::string currentSection;
    
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#' || line[0] == '/' || line == "{" || line == "}") continue;
        
        // Check for section headers
        if (line.find("nucleus") != std::string::npos && line.find(":") != std::string::npos) {
            currentSection = "nucleus";
            continue;
        }
        if (line.find("levelDensity") != std::string::npos) {
            currentSection = "levelDensity";
            continue;
        }
        if (line.find("spinCutoff") != std::string::npos) {
            currentSection = "spinCutoff";
            continue;
        }
        if (line.find("gammaStrength") != std::string::npos) {
            currentSection = "gammaStrength";
            continue;
        }
        if (line.find("initialExcitation") != std::string::npos) {
            currentSection = "initialExcitation";
            continue;
        }
        if (line.find("simulation") != std::string::npos) {
            currentSection = "simulation";
            continue;
        }
        if (line.find("output") != std::string::npos) {
            currentSection = "output";
            continue;
        }
        if (line.find("continuum") != std::string::npos) {
            currentSection = "continuum";
            continue;
        }
        if (line.find("internalConversion") != std::string::npos) {
            currentSection = "internalConversion";
            continue;
        }
        if (line.find("parity") != std::string::npos) {
            currentSection = "parity";
            continue;
        }
        
        auto [key, value] = parseLine(line);
        if (key.empty()) continue;
        
        try {
            // Nucleus section
            if (currentSection == "nucleus") {
                if (key == "Z") config.nucleus.Z = std::stoi(value);
                else if (key == "A") config.nucleus.A = std::stoi(value);
                else if (key == "Sn") config.nucleus.Sn = std::stod(value);
                else if (key == "levelsFile") config.nucleus.levelsFile = value;
            }
            // Level Density section
            else if (currentSection == "levelDensity") {
                if (key == "model") {
                    if (value == "BSFG") config.levelDensity.model = Config::LevelDensityConfig::Model::BSFG;
                    else if (value == "CTM") config.levelDensity.model = Config::LevelDensityConfig::Model::CTM;
                    else if (value == "TABLE") config.levelDensity.model = Config::LevelDensityConfig::Model::TABLE;
                    else if (value == "USER_DEFINED") config.levelDensity.model = Config::LevelDensityConfig::Model::USER_DEFINED;
                }
                else if (key == "a") config.levelDensity.a = std::stod(value);
                else if (key == "E1") config.levelDensity.E1 = std::stod(value);
                else if (key == "T") config.levelDensity.T = std::stod(value);
                else if (key == "E0") config.levelDensity.E0 = std::stod(value);
                else if (key == "useEnergyDependentA") config.levelDensity.useEnergyDependentA = (value == "true");
                else if (key == "tableFile") config.levelDensity.tableFile = value;
            }
            // Spin Cutoff section
            else if (currentSection == "spinCutoff") {
                if (key == "model") {
                    if (value == "VON_EGIDY_05") config.spinCutoff.model = Config::SpinCutoffConfig::Model::VON_EGIDY_05;
                    else if (value == "TALYS") config.spinCutoff.model = Config::SpinCutoffConfig::Model::TALYS;
                    else if (value == "SINGLE_PARTICLE") config.spinCutoff.model = Config::SpinCutoffConfig::Model::SINGLE_PARTICLE;
                }
                else if (key == "spinCutoffD") config.spinCutoff.spinCutoffD = std::stod(value);
            }
            // Gamma Strength section
            else if (currentSection == "gammaStrength") {
                if (key == "e1Model") {
                    if (value == "GEN_LORENTZ") config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::GEN_LORENTZ;
					else if (value == "EGLO") config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::EGLO;
                    else if (value == "KMF") config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::KMF;
                }
            }
            // Initial Excitation section
            
            else if (currentSection == "initialExcitation") {
                if (key == "mode") {
                    if (value == "SINGLE") config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::SINGLE;
                    else if (value == "SELECT") config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::SELECT;
                    else if (value == "SPREAD") config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::SPREAD;
                    else if (value == "FULL_REACTION") config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::FULL_REACTION;
                }
                
                else if (key == "SINGLE_excitationEnergy") config.initialExcitation.excitationEnergy = std::stod(value);
                else if (key == "SINGLE_spin") config.initialExcitation.spin = std::stod(value);
                else if (key == "SINGLE_parity") config.initialExcitation.parity = std::stoi(value);
                else if (key == "SPREAD_meanEnergy") config.initialExcitation.meanEnergy = std::stod(value);
                else if (key == "SPREAD_energySpread") config.initialExcitation.energySpread = std::stod(value);
                
                else if (key.find("SELECT_states") != std::string::npos) {
                    // This is the start of the SELECT_states array
                    // Read until we find the closing bracket
                    std::vector<InitialExcitationConfig::SelectState> states;
                    std::cout <<"Alive" <<std::endl;
                    while (std::getline(file, line)) {
                        line = trim(line);
                        if (line.empty() || line[0] == '/' || line[0] == '#') continue;
                        std::cout <<"Alive1" <<std::endl;
                        // Check for end of array
                        if (line.find("]") != std::string::npos) break;
                        std::cout <<"Alive2" <<std::endl;
                        // Skip array opening and comments
                        if (line.find("[") != std::string::npos) continue;
                        if (line.find("_SELECT") != std::string::npos) continue;
                        std::cout <<"Alive3" <<std::endl;
                        // Parse state object
                        if (line.find("{") != std::string::npos) {
                            InitialExcitationConfig::SelectState state;
                            std::cout <<"Alive4" <<std::endl;
                            // Read state properties
                            while (std::getline(file, line)) {
                                line = trim(line);
                                if (line.empty()) continue;
                                if (line.find("}") != std::string::npos) break;
                                
                                auto [stateKey, stateValue] = parseLine(line);
                                if (stateKey == "energy") {
                                    state.energy = std::stod(stateValue);
                                } else if (stateKey == "spin") {
                                    state.spin = std::stod(stateValue);
                                } else if (stateKey == "parity") {
                                    state.parity = std::stoi(stateValue);
                                } else if (stateKey == "branchingRatio") {
                                    state.branchingRatio = std::stod(stateValue);
                                }
                            }
                            states.push_back(state);
                        }
                    }
                    
                    config.initialExcitation.selectStates = states;
                    
                    // Validate that branching ratios sum to 1.0
                    double brSum = 0.0;
                    for (const auto& state : states) {
                        brSum += state.branchingRatio;
                    }
                    if (std::abs(brSum - 1.0) > 0.01) {
                        std::cerr << "Warning: SELECT branching ratios sum to " << brSum
                                  << " (should be 1.0)\n";
                    }
                }
                
            }
            
            // Simulation section
            else if (currentSection == "simulation") {
                if (key == "numRealizations") config.simulation.numRealizations = std::stoi(value);
                else if (key == "eventsPerRealization") config.simulation.eventsPerRealization = std::stoi(value);
                else if (key == "randomSeed") config.simulation.randomSeed = std::stoul(value);
                else if (key == "useParallel") config.simulation.useParallel = (value == "true");
            }
            // Output section
            else if (currentSection == "output") {
                if (key == "outputFile") config.output.outputFile = value;
                else if (key == "paramFile") config.output.paramFile = value;
                else if (key == "saveTree") config.output.saveTree = (value == "true");
                else if (key == "gammaSpectrumBins") config.output.gammaSpectrumBins = std::stoi(value);
                else if (key == "maxGammaEnergy") config.output.maxGammaEnergy = std::stod(value);
            }
            // Continuum section
            else if (currentSection == "continuum") {
                if (key == "energySpacing") config.continuum.energySpacing = std::stod(value);
                else if (key == "distribution") {
                    if (value == "POISSON") config.continuum.distribution = Config::ContinuumConfig::Distribution::POISSON;
                    else if (value == "WIGNER") config.continuum.distribution = Config::ContinuumConfig::Distribution::WIGNER;
                }
            }
            // Internal Conversion section
            else if (currentSection == "internalConversion") {
                if (key == "enabled") config.internalConversion.enabled = (value == "true");
                else if (key == "model") {
                    if (value == "BRICC") config.internalConversion.model = Config::InternalConversionConfig::Model::BRICC;
                    else if (value == "TABLE") config.internalConversion.model = Config::InternalConversionConfig::Model::TABLE;
                }
            }
            // Parity section
            else if (currentSection == "parity") {
                if (key == "model") {
                    if (value == "EQUIPARTITION") config.parity.model = Config::ParityConfig::Model::EQUIPARTITION;
                    else if (value == "ENERGY_DEPENDENT") config.parity.model = Config::ParityConfig::Model::ENERGY_DEPENDENT;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Error parsing config line: " << line << std::endl;
            std::cerr << "  Error: " << e.what() << std::endl;
        }
    }
    
    file.close();
    
    std::string errorMsg;
    if (!config.validate(errorMsg)) {
        std::cerr << "Configuration validation failed: " << errorMsg << std::endl;
        throw std::runtime_error("Invalid configuration");
    }
    
    std::cout << "Configuration loaded successfully." << std::endl;
    return config;
}

void Config::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not create config file: " << filename << std::endl;
        return;
    }
    
    file << "{\n";
    file << "  \"nucleus\": {\n";
    file << "    \"Z\": " << nucleus.Z << ",\n";
    file << "    \"A\": " << nucleus.A << ",\n";
    file << "    \"Sn\": " << nucleus.Sn << ",\n";
    file << "    \"levelsFile\": \"" << nucleus.levelsFile << "\"\n";
    file << "  },\n";
    
    file << "  \"levelDensity\": {\n";
    file << "    \"model\": \"";
    switch(levelDensity.model) {
        case LevelDensityConfig::Model::BSFG: file << "BSFG"; break;
        case LevelDensityConfig::Model::CTM: file << "CTM"; break;
        case LevelDensityConfig::Model::TABLE: file << "TABLE"; break;
        case LevelDensityConfig::Model::USER_DEFINED: file << "USER_DEFINED"; break;
    }
    file << "\",\n";
    file << "    \"a\": " << levelDensity.a << ",\n";
    file << "    \"E1\": " << levelDensity.E1 << "\n";
    file << "  },\n";
    
    file << "  \"initialExcitation\": {\n";
    file << "    \"mode\": \"";
    switch(initialExcitation.mode) {
        case InitialExcitationConfig::Mode::SINGLE: file << "SINGLE"; break;
        case InitialExcitationConfig::Mode::SELECT: file << "SELECT"; break;
        case InitialExcitationConfig::Mode::SPREAD: file << "SPREAD"; break;
        case InitialExcitationConfig::Mode::FULL_REACTION: file << "FULL_REACTION"; break;
    }
    file << "\",\n";
    file << "    \"excitationEnergy\": " << initialExcitation.excitationEnergy << ",\n";
    file << "    \"spin\": " << initialExcitation.spin << ",\n";
    file << "    \"parity\": " << initialExcitation.parity << "\n";
    file << "  },\n";
    
    file << "  \"simulation\": {\n";
    file << "    \"numRealizations\": " << simulation.numRealizations << ",\n";
    file << "    \"eventsPerRealization\": " << simulation.eventsPerRealization << "\n";
    file << "  },\n";
    
    file << "  \"output\": {\n";
    file << "    \"outputFile\": \"" << output.outputFile << "\",\n";
    file << "    \"gammaSpectrumBins\": " << output.gammaSpectrumBins << "\n";
    file << "  }\n";
    file << "}\n";
    
    file.close();
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
    if (simulation.numRealizations <= 0) {
        errorMessage = "Number of realizations must be positive";
        return false;
    }
    if (simulation.eventsPerRealization <= 0) {
        errorMessage = "Events per realization must be positive";
        return false;
    }
    return true;
}

void Config::print() const {
    std::cout << "\n╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║          RAINIER++ Configuration Summary              ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Nucleus:\n";
    std::cout << "  Z = " << nucleus.Z << ", A = " << nucleus.A << "\n";
    std::cout << "  Sn = " << nucleus.Sn << " MeV\n";
    std::cout << "  Levels file: " << nucleus.levelsFile << "\n\n";
    
    std::cout << "Level Density Model: ";
    switch(levelDensity.model) {
        case LevelDensityConfig::Model::BSFG: 
            std::cout << "Back-Shifted Fermi Gas\n";
            std::cout << "  a = " << levelDensity.a << " MeV^-1\n";
            std::cout << "  E1 = " << levelDensity.E1 << " MeV\n";
            break;
        case LevelDensityConfig::Model::CTM:
            std::cout << "Constant Temperature Model\n";
            std::cout << "  T = " << levelDensity.T << " MeV\n";
            break;
        case LevelDensityConfig::Model::TABLE:
            std::cout << "Table-based\n";
            break;
        case LevelDensityConfig::Model::USER_DEFINED:
            std::cout << "User-defined\n";
            break;
    }
    std::cout << "\n";
    
    std::cout << "Spin Cutoff Model: ";
    switch(spinCutoff.model) {
        case SpinCutoffConfig::Model::VON_EGIDY_05: std::cout << "Von Egidy 2005\n"; break;
        case SpinCutoffConfig::Model::TALYS: std::cout << "TALYS\n"; break;
        case SpinCutoffConfig::Model::SINGLE_PARTICLE: std::cout << "Single Particle\n"; break;
        default: std::cout << "Other\n"; break;
    }
    std::cout << "\n";
    
    std::cout << "Initial Excitation Mode: ";
    switch(initialExcitation.mode) {
        case InitialExcitationConfig::Mode::SINGLE:
            std::cout << "SINGLE (DICEBOX-like)\n";
            std::cout << "  E = " << initialExcitation.excitationEnergy << " MeV\n";
            std::cout << "  J = " << initialExcitation.spin << "\n";
            std::cout << "  π = " << (initialExcitation.parity == 1 ? "+" : "-") << "\n";
            break;
        case InitialExcitationConfig::Mode::SELECT:
            std::cout << "SELECT (Beta-decay-like)\n";
            std::cout << "  Number of states: " << initialExcitation.selectStates.size() << "\n";
            for (size_t i = 0; i < initialExcitation.selectStates.size(); ++i) {
                const auto& state = initialExcitation.selectStates[i];
                std::cout << "    State " << i+1 << ": "
                          << "E=" << state.energy << " MeV, "
                          << "J=" << state.spin << ", "
                          << "π=" << (state.parity == 1 ? "+" : "-") << ", "
                          << "BR=" << state.branchingRatio << "\n";
            }
            break;
        case InitialExcitationConfig::Mode::SPREAD:
            std::cout << "SPREAD (Energy spread with spin distribution)\n";
            std::cout << "  Mean E = " << initialExcitation.meanEnergy << " MeV\n";
            std::cout << "  Spread = " << initialExcitation.energySpread << " MeV\n";
            break;
        case InitialExcitationConfig::Mode::FULL_REACTION:
            std::cout << "FULL_REACTION (TALYS population)\n";
            break;
    }
    std::cout << "\n";
    
    std::cout << "Simulation:\n";
    std::cout << "  Realizations: " << simulation.numRealizations << "\n";
    std::cout << "  Events per realization: " << simulation.eventsPerRealization << "\n";
    std::cout << "  Random seed: " << simulation.randomSeed << "\n";
    std::cout << "  Parallel: " << (simulation.useParallel ? "Yes" : "No") << "\n\n";
    
    std::cout << "Output:\n";
    std::cout << "  File: " << output.outputFile << "\n";
    std::cout << "  Save tree: " << (output.saveTree ? "Yes" : "No") << "\n";
    std::cout << "  Gamma spectrum bins: " << output.gammaSpectrumBins << "\n\n";
    
    std::cout << "Internal Conversion: " << (internalConversion.enabled ? "Enabled" : "Disabled") << "\n\n";
}

} // namespace rainier
