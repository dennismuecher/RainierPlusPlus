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

Nucleus::Nucleus(int Z, int A, const Config& config)
    : Z_(Z), A_(A), Sn_(0.0), criticalEnergy_(0.0),
      totalContinuumLevels_(0), maxSpinBin_(20),  // Reasonable default
      maxEnergy_(10.0), energySpacing_(0.1), numEnergyBins_(0),
    levelDensity_(nullptr), spinCutoff_(nullptr) {
    
    if (Z <= 0 || A <= 0) {
        throw std::invalid_argument("Z and A must be positive");
    }
    if (Z > A) {
        throw std::invalid_argument("Z cannot exceed A");
    }
        // Create level density model based on config
            if (config.levelDensity.model == Config::LevelDensityConfig::Model::BSFG) {
                levelDensity_ = std::make_unique<BackShiftedFermiGas>(
                    config.levelDensity.a, config.levelDensity.E1,
                    config.levelDensity.useEnergyDependentA,
                    config.levelDensity.aAsymptotic,
                    config.levelDensity.shellCorrectionW,
                    config.levelDensity.dampingGamma,
                    A_
                );
            }
        
            else if (config.levelDensity.model == Config::LevelDensityConfig::Model::CTM) {
                levelDensity_ = std::make_unique<ConstantTemperature>(
                    config.levelDensity.T,
                    config.levelDensity.E0,
                     A_,
                    Z
                );
            }
            else {
                throw std::runtime_error("Unsupported level density model");
            }
        
        // Create spin cutoff model based on config
        if (config.spinCutoff.model == Config::SpinCutoffConfig::Model::VON_EGIDY_05) {
            spinCutoff_ = std::make_unique<VonEgidy05>(
                levelDensity_,
                A_,
                config.spinCutoff.useOsloShift ? config.spinCutoff.osloShift : 0.0
            );
        }
        else if (config.spinCutoff.model == Config::SpinCutoffConfig::Model::SINGLE_PARTICLE) {
            spinCutoff_ = std::make_unique<SingleParticle>(
                levelDensity_,
                A_,
                config.spinCutoff.useOsloShift ? config.spinCutoff.osloShift : 0.0
            );
        }
        else if (config.spinCutoff.model == Config::SpinCutoffConfig::Model::RIGID_SPHERE) {
            spinCutoff_ = std::make_unique<RigidSphere>(
                levelDensity_,
                A_,
                config.spinCutoff.useOsloShift ? config.spinCutoff.osloShift : 0.0
            );
        }
       
        else if (config.spinCutoff.model == Config::SpinCutoffConfig::Model::TALYS) {
            spinCutoff_ = std::make_unique<TALYSSpinCutoff>(
                levelDensity_,
                A_
                Sn_
                config.spinCutoff.spinCutoffD,
                config.spinCutoff.Ed,
                config.levelDensity.aAsymptotic,
                config.spinCutoff.useOsloShift ? config.spinCutoff.osloShift : 0.0
            );
        }
        
}

void Nucleus::loadDiscreteLevels(const std::string& filename, int maxLevels) {
    std::cout << "Loading discrete levels from: " << filename << std::endl;
    
    // First read ALL levels for allDiscreteLevels_
    auto allLevelData = LevelFileReader::readLevelFile(filename, Z_, A_, 0);  // 0 = all
    allDiscreteLevels_ = LevelFileReader::createDiscreteLevels(allLevelData);
    
    // Then read limited levels for discreteLevels_ (for simulation)
    auto limitedLevelData = LevelFileReader::readLevelFile(filename, Z_, A_, maxLevels);
    discreteLevels_ = LevelFileReader::createDiscreteLevels(limitedLevelData);
    
    if (!discreteLevels_.empty()) {
        criticalEnergy_ = discreteLevels_.back()->getEnergy();
    }
    
    std::cout << "Loaded " << allDiscreteLevels_.size() << " total levels" << std::endl;
    std::cout << "Using " << discreteLevels_.size() << " levels for simulation" << std::endl;
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

void Nucleus::buildPoissonLevels(int realization) {
    // Poisson distribution for level spacings
    // Each E-J-π bin gets a random number of levels from Poisson distribution
    // based on the average level density
        
    TRandom2 rng(1 + realization);  // Seed with realization number
    
    int maxLevelsInBin = 0;
    
    for (int ex = 0; ex < numEnergyBins_; ++ex) {
        double binCenterEnergy = criticalEnergy_ + (ex + 0.5) * energySpacing_;
        
        for (int spb = 0; spb < maxSpinBin_; ++spb) {
            for (int par = 0; par < 2; ++par) {
                double spin = Level::binToSpin(spb, isEvenA());
                
            
                double density = levelDensity_->getDensity(binCenterEnergy, spin, par);
                double avgLevels = density * energySpacing_;
                
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

void Nucleus::buildWignerLevels(int realization) {
    // Wigner distribution for level spacings
    // Levels are distributed according to Wigner spacing distribution:
    // P(s) = (π/2) s exp(-πs²/4)
    // where s is the spacing in units of mean spacing
    
    
    TRandom2 rng(1 + realization);
    
    // Build levels independently for each J-π combination
    for (int spb = 0; spb < maxSpinBin_; ++spb) {
        for (int par = 0; par < 2; ++par) {
            double spin = Level::binToSpin(spb, isEvenA());
            
            // Calculate expected cumulative number of levels vs energy
            std::vector<double> expectedCumulative(numEnergyBins_);
            for (int ex = 0; ex < numEnergyBins_; ++ex) {
                double binCenterEnergy = criticalEnergy_ + (ex + 0.5) * energySpacing_;
                double density = levelDensity_->getDensity(binCenterEnergy, spin, par);
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
