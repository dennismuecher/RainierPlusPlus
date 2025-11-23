// Nucleus.cpp - Complete implementation of nuclear level scheme
#include "core/Nucleus.h"
#include "Config.h"
#include "io/LevelFileReader.h"
#include "utils/PhysicsConstants.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <TRandom2.h>
#include <fstream>
#include <sstream>

namespace rainier {

Nucleus::Nucleus(int Z, int A) 
    : Z_(Z), A_(A), Sn_(0.0), criticalEnergy_(0.0),
      totalContinuumLevels_(0), maxSpinBin_(20),  // Reasonable default
      maxEnergy_(10.0), energySpacing_(0.1), numEnergyBins_(0) {
    
    if (Z <= 0 || A <= 0) {
        throw std::invalid_argument("Z and A must be positive");
    }
    if (Z > A) {
        throw std::invalid_argument("Z cannot exceed A");
    }
}

void Nucleus::loadDiscreteLevels(const std::string& filename, int maxLevels) {
    std::cout << "Loading discrete levels from: " << filename << std::endl;
    
    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open levels file: " + filename);
        }
        
        std::string line;
        bool foundNucleus = false;
        int totalLevelsInFile = 0;
        
        // Clear both vectors
        discreteLevels_.clear();
        allDiscreteLevels_.clear();
        
        while (std::getline(file, line)) {
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#') continue;
            
            // Check if this is a nucleus header line
            // Format: "144Nd  144   60  202  294  101   15    7.817030    7.968790"
            std::istringstream iss(line);
            std::string nuclideName;
            int fileA, fileZ, numLevels;
            
            // Try to parse as header line
            if (iss >> nuclideName >> fileA >> fileZ >> numLevels) {
                // Check if this matches our nucleus
                if (fileA == A_ && fileZ == Z_) {
                    foundNucleus = true;
                    totalLevelsInFile = numLevels;
                    std::cout << "  Found nucleus " << nuclideName
                              << " (A=" << A_ << ", Z=" << Z_
                              << ") with " << numLevels << " levels" << std::endl;
                    
                    // Now read exactly numLevels lines
                    int levelIndex = 0;
                    int linesRead = 0;
                    
                    while (linesRead < numLevels && std::getline(file, line)) {
                        // Skip empty lines
                        if (line.empty()) continue;
                        
                        // Parse level line
                        // Format: "  1   0.000000   0.0  1  7.227E+22  ..."
                        std::istringstream levelStream(line);
                        int levelNum;
                        double energy, spin, lifetime;
                        int parity;
                        
                        if (!(levelStream >> levelNum >> energy >> spin >> parity >> lifetime)) {
                            std::cerr << "Warning: Skipping malformed level line: " << line << std::endl;
                            linesRead++;
                            continue;
                        }
                        if (parity !=0 && parity !=1)
                        {
                            std::cout <<"Error in parity in line: " << linesRead <<": parity is" <<parity<<std::endl;
                        }
                        // Create level
                        auto level = std::make_shared<DiscreteLevel>(energy, spin, parity, lifetime);
                        level->setLevelIndex(levelIndex);
                        
                        // Store in allDiscreteLevels_ (ALL levels for this nucleus)
                        allDiscreteLevels_.push_back(level);
                        
                        // Only add to discreteLevels_ if within maxLevels limit (for simulation)
                        if (maxLevels <= 0 || static_cast<int>(discreteLevels_.size()) < maxLevels) {
                            discreteLevels_.push_back(level);
                        }
                        
                        levelIndex++;
                        linesRead++;
                    }
                    
                    // We've read all levels for our nucleus, we can stop
                    break;
                    
                } else if (foundNucleus) {
                    // We already found our nucleus and now hit another one, stop
                    break;
                }
            }
        }
        
        file.close();
        
        if (!foundNucleus) {
            throw std::runtime_error("Nucleus A=" + std::to_string(A_) +
                                   ", Z=" + std::to_string(Z_) + " not found in file");
        }
        
        if (discreteLevels_.empty()) {
            throw std::runtime_error("No discrete levels loaded for A=" + std::to_string(A_));
        }
        
        // Critical energy is automatically the highest energy in the truncated set
        criticalEnergy_ = discreteLevels_.back()->getEnergy();
        
        std::cout << "Loaded " << allDiscreteLevels_.size()
                  << " total discrete levels from file for A=" << A_ << std::endl;
        std::cout << "Using " << discreteLevels_.size()
                  << " discrete levels up to " << criticalEnergy_
                  << " MeV for simulation" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading discrete levels: " << e.what() << std::endl;
        std::cerr << "Creating minimal level scheme for testing..." << std::endl;
        
        // Create minimal scheme if file reading fails
        auto gs = std::make_shared<DiscreteLevel>(0.0, 0.0, 1, 1e20);
        gs->setLevelIndex(0);
        discreteLevels_.push_back(gs);
        allDiscreteLevels_.push_back(gs);
        criticalEnergy_ = 0.0;
    }
}

void Nucleus::buildContinuumLevels(const Config& config, int realization) {
    std::cout << "Building continuum levels (realization " << realization << ")..." << std::endl;
    
    // Clear previous continuum levels
    continuumLevels_.clear();
    totalContinuumLevels_ = 0;
    
    // Calculate energy binning
    maxEnergy_ = config.initialExcitation.excitationEnergy + 2.0;  // Add buffer
    if (config.continuum.forceBinNumber) {
        numEnergyBins_ = config.continuum.numBins;
        energySpacing_ = (maxEnergy_ - criticalEnergy_) / numEnergyBins_;
    } else {
        energySpacing_ = config.continuum.energySpacing;
        numEnergyBins_ = static_cast<int>((maxEnergy_ - criticalEnergy_) / energySpacing_) + 1;
    }
    
    std::cout << "  Energy range: " << criticalEnergy_ << " - " << maxEnergy_ << " MeV" << std::endl;
    std::cout << "  Energy spacing: " << energySpacing_ << " MeV" << std::endl;
    std::cout << "  Number of bins: " << numEnergyBins_ << std::endl;
    
    // Build levels according to distribution type
    if (config.continuum.distribution == Config::ContinuumConfig::Distribution::POISSON) {
        buildPoissonLevels(config, realization);
    } else {
        buildWignerLevels(config, realization);
    }
    
    std::cout << "  Total continuum levels: " << totalContinuumLevels_ << std::endl;
}

void Nucleus::buildPoissonLevels(const Config& config, int realization) {
    // Poisson distribution for level spacings
    // Each E-J-π bin gets a random number of levels from Poisson distribution
    // based on the average level density
    
    (void)config; // Unused for now - will use when level density models integrated
    
    TRandom2 rng(1 + realization);  // Seed with realization number
    
    // We need level density and spin cutoff models to calculate average number of levels
    // For now, use a simple approximation
    // TODO: Use actual level density model once implemented
    
    int maxLevelsInBin = 0;
    
    for (int ex = 0; ex < numEnergyBins_; ++ex) {
        double binCenterEnergy = criticalEnergy_ + (ex + 0.5) * energySpacing_;
        
        for (int spb = 0; spb < maxSpinBin_; ++spb) {
            for (int par = 0; par < 2; ++par) {
                double spin = Level::binToSpin(spb, isEvenA());
                
                // Simple level density estimate (will be replaced with proper model)
                // ρ(E,J,π) ≈ ρ_total(E) × P(J) × P(π)
                // For now, use a crude approximation
                double avgLevels = this->estimateLevelDensity(binCenterEnergy, spin, par) 
                                 * energySpacing_;
                
                if (avgLevels < 0.01) avgLevels = 0.01;  // Minimum for Poisson
                
                // Generate random number of levels
                int numLevels = rng.Poisson(avgLevels);
                
                if (numLevels > maxLevelsInBin) {
                    maxLevelsInBin = numLevels;
                }
                
                // Create the levels
                std::vector<std::shared_ptr<ContinuumLevel>> binLevels;
                for (int lvl = 0; lvl < numLevels; ++lvl) {
                    // Random energy within bin
                    double energy = binCenterEnergy + 
                                  rng.Uniform(-energySpacing_/2.0, energySpacing_/2.0);
                    
                    auto level = std::make_shared<ContinuumLevel>(
                        energy, spin, par, ex, lvl);
                    
                    binLevels.push_back(level);
                    totalContinuumLevels_++;
                }
                
                if (!binLevels.empty()) {
                    continuumLevels_[std::make_tuple(ex, spb, par)] = binLevels;
                }
            }
        }
    }
    
    std::cout << "  Maximum levels in any bin: " << maxLevelsInBin << std::endl;
}

void Nucleus::buildWignerLevels(const Config& config, int realization) {
    // Wigner distribution for level spacings
    // Levels are distributed according to Wigner spacing distribution:
    // P(s) = (π/2) s exp(-πs²/4)
    // where s is the spacing in units of mean spacing
    
    (void)config; // Unused for now - will use when level density models integrated
    
    TRandom2 rng(1 + realization);
    
    // Build levels independently for each J-π combination
    for (int spb = 0; spb < maxSpinBin_; ++spb) {
        for (int par = 0; par < 2; ++par) {
            double spin = Level::binToSpin(spb, isEvenA());
            
            // Calculate expected cumulative number of levels vs energy
            std::vector<double> expectedCumulative(numEnergyBins_);
            for (int ex = 0; ex < numEnergyBins_; ++ex) {
                double binCenterEnergy = criticalEnergy_ + (ex + 0.5) * energySpacing_;
                double density = this->estimateLevelDensity(binCenterEnergy, spin, par);
                
                if (ex == 0) {
                    expectedCumulative[ex] = density * energySpacing_;
                } else {
                    expectedCumulative[ex] = expectedCumulative[ex-1] + 
                                           density * energySpacing_;
                }
            }
            
            // Generate Wigner-distributed spacings
            // Start with a random initial spacing
            double wignerSum = 2.0 / std::sqrt(constants::PI) * 
                             std::sqrt(-std::log(rng.Uniform(1.0)));
            
            // Place levels according to Wigner distribution
            for (int ex = 0; ex < numEnergyBins_; ++ex) {
                std::vector<std::shared_ptr<ContinuumLevel>> binLevels;
                int levelInBin = 0;
                
                while (wignerSum < expectedCumulative[ex]) {
                    // Place a level here
                    double binCenterEnergy = criticalEnergy_ + (ex + 0.5) * energySpacing_;
                    double energy = binCenterEnergy + 
                                  rng.Uniform(-energySpacing_/2.0, energySpacing_/2.0);
                    
                    auto level = std::make_shared<ContinuumLevel>(
                        energy, spin, par, ex, levelInBin);
                    
                    binLevels.push_back(level);
                    levelInBin++;
                    totalContinuumLevels_++;
                    
                    // Generate next Wigner spacing
                    wignerSum += 2.0 / std::sqrt(constants::PI) * 
                               std::sqrt(-std::log(rng.Uniform(1.0)));
                }
                
                if (!binLevels.empty()) {
                    continuumLevels_[std::make_tuple(ex, spb, par)] = binLevels;
                }
            }
        }
    }
}

double Nucleus::estimateLevelDensity(double energy, double spin, int parity) const {
    // Simple estimate for level density used during continuum construction
    // This is a simplified BSFG model - the full simulation will use
    // the proper LevelDensityModel classes
    
    if (energy < criticalEnergy_) {
        return 0.0;
    }
    
    // Simple BSFG approximation matching original RAINIER physics
    double U = energy - 1.0;  // Effective excitation energy (E - E1)
    if (U < 0.1) U = 0.1;
    
    double a = A_ / 8.0;  // Level density parameter (A/8 is typical)
    
    // Spin cutoff parameter (Von Egidy formula)
    double sigma2 = 0.0146 * std::pow(A_, 5.0/3.0) * 
                   (1.0 + std::sqrt(1.0 + 4.0 * a * U)) / (2.0 * a);
    
    // Total density (BSFG formula)
    double rhoTotal = std::exp(2.0 * std::sqrt(a * U)) / 
                     (12.0 * std::sqrt(2.0 * sigma2) * std::pow(a, 0.25) * std::pow(U, 1.25));
    
    // Spin distribution
    double spinFactor = (2.0 * spin + 1.0) * 
                       std::exp(-(spin + 0.5) * (spin + 0.5) / (2.0 * sigma2)) / 
                       sigma2;
    
    // Parity distribution (equipartition)
    double parityFactor = 0.5;
    
    return rhoTotal * spinFactor * parityFactor;
}

std::shared_ptr<DiscreteLevel> Nucleus::getDiscreteLevel(int index) const {
    if (index < 0 || index >= static_cast<int>(discreteLevels_.size())) {
        throw std::out_of_range("Discrete level index out of range");
    }
    return discreteLevels_[index];
}

const std::vector<std::shared_ptr<ContinuumLevel>>& 
Nucleus::getContinuumLevels(int energyBin, int spinBin, int parity) const {
    static std::vector<std::shared_ptr<ContinuumLevel>> empty;
    
    auto key = std::make_tuple(energyBin, spinBin, parity);
    auto it = continuumLevels_.find(key);
    
    if (it != continuumLevels_.end()) {
        return it->second;
    }
    return empty;
}

int Nucleus::getNumLevelsInBin(int energyBin, int spinBin, int parity) const {
    auto key = std::make_tuple(energyBin, spinBin, parity);
    auto it = continuumLevels_.find(key);
    
    if (it != continuumLevels_.end()) {
        return it->second.size();
    }
    return 0;
}

double Nucleus::getBinCentralEnergy(int bin) const {
    return criticalEnergy_ + (bin + 0.5) * energySpacing_;
}

int Nucleus::getEnergyBin(double energy) const {
    if (energy < criticalEnergy_) {
        return -1;  // In discrete region
    }
    
    int bin = static_cast<int>((energy - criticalEnergy_) / energySpacing_);
    
    if (bin >= numEnergyBins_) {
        return numEnergyBins_ - 1;  // Cap at maximum
    }
    
    return bin;
}

void Nucleus::printSummary() const {
    std::cout << "\n=== Nucleus Summary ===" << std::endl;
    std::cout << "Z = " << Z_ << ", A = " << A_ << ", N = " << getN() << std::endl;
    std::cout << "Even-A: " << (isEvenA() ? "yes" : "no") << std::endl;
    std::cout << "Sn = " << Sn_ << " MeV" << std::endl;
    std::cout << "\nDiscrete region:" << std::endl;
    std::cout << "  Levels: " << discreteLevels_.size() << std::endl;
    std::cout << "  Critical energy: " << criticalEnergy_ << " MeV" << std::endl;
    
    if (totalContinuumLevels_ > 0) {
        std::cout << "\nContinuum region:" << std::endl;
        std::cout << "  Energy range: " << criticalEnergy_ << " - " << maxEnergy_ << " MeV" << std::endl;
        std::cout << "  Energy bins: " << numEnergyBins_ << std::endl;
        std::cout << "  Bin spacing: " << energySpacing_ << " MeV" << std::endl;
        std::cout << "  Total levels: " << totalContinuumLevels_ << std::endl;
        std::cout << "  Max spin bin: " << maxSpinBin_ << std::endl;
    }
    
    // Print first few discrete levels
    std::cout << "\nFirst discrete levels:" << std::endl;
    int numToPrint = std::min(5, static_cast<int>(discreteLevels_.size()));
    for (int i = 0; i < numToPrint; ++i) {
        auto level = discreteLevels_[i];
        std::cout << "  " << i << ": E=" << level->getEnergy() 
                  << " MeV, J=" << level->getSpin()
                  << (level->getParity() == 1 ? "+" : "-")
                  << ", gammas=" << level->getTransitions().size()
                  << std::endl;
    }
    std::cout << std::endl;
}

} // namespace rainier
