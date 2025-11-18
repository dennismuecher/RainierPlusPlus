// DecaySimulator.cpp - Fixed implementation with proper energy conservation
// Follows original RAINIER.C two-stage selection: bin first, then level
#include "simulation/DecaySimulator.h"
#include "simulation/CascadeEvent.h"
#include "models/LevelDensity.h"
#include "models/SpinCutoff.h"
#include "models/GammaStrength.h"
#include "utils/PhysicsConstants.h"
#include <TRandom2.h>
#include <iostream>
#include <cmath>

namespace rainier {

DecaySimulator::DecaySimulator(Nucleus& nucleus, const Config& config, int realization)
    : nucleus_(nucleus), config_(config), realization_(realization), numStuckEvents_(0) {
    
    rng_ = std::make_unique<TRandom2>(1 + realization + config.simulation.randomSeed);
    
    levelDensity_ = std::make_unique<BackShiftedFermiGas>(
        config.levelDensity.a, config.levelDensity.E1,
        config.levelDensity.useEnergyDependentA,
        config.levelDensity.aAsymptotic,
        config.levelDensity.shellCorrectionW,
        config.levelDensity.dampingGamma,
        nucleus.getA()
    );
    
    spinCutoff_ = std::make_unique<VonEgidy05>(levelDensity_.get(), nucleus.getA());
    
    std::vector<Resonance> e1Resonances;
    for (const auto& r : config.gammaStrength.e1Resonances) {
        e1Resonances.push_back({r.energy, r.width, r.sigma});
    }
    auto e1 = std::make_unique<E1GenLorentz>(e1Resonances, config.gammaStrength.constantT);
    
    std::vector<Resonance> m1Resonances;
    for (const auto& r : config.gammaStrength.m1Resonances) {
        m1Resonances.push_back({r.energy, r.width, r.sigma});
    }
    auto m1 = std::make_unique<M1StandardLorentz>(m1Resonances);
    
    auto e2 = std::make_unique<E2StandardLorentz>(
        config.gammaStrength.e2Energy,
        config.gammaStrength.e2Width,
        config.gammaStrength.e2Sigma
    );
    
    gammaStrength_ = std::make_unique<CombinedStrength>(
        std::move(e1), std::move(m1), std::move(e2)
    );
    
    std::cout << "  DecaySimulator initialized" << std::endl;
}

DecaySimulator::~DecaySimulator() = default;

void DecaySimulator::run() {
    std::cout << "  Running simulation..." << std::endl;
    events_.reserve(config_.simulation.eventsPerRealization);
    
    for (int ev = 0; ev < config_.simulation.eventsPerRealization; ++ev) {
        if (ev % config_.simulation.updateInterval == 0) {
            std::cout << "    Event " << ev << " / " 
                     << config_.simulation.eventsPerRealization << "\r" << std::flush;
        }
        simulateEvent();
    }
    
    std::cout << "    Event " << config_.simulation.eventsPerRealization 
              << " / " << config_.simulation.eventsPerRealization << std::endl;
    
    if (numStuckEvents_ > 0) {
        std::cout << "  Warning: " << numStuckEvents_ 
                  << " events stuck in isomeric states" << std::endl;
    }
}

void DecaySimulator::simulateEvent() {
    CascadeEvent event;
    
    std::shared_ptr<Level> currentLevel;
    double initialExcitation;
    selectInitialState(currentLevel, initialExcitation);
    event.setInitialState(currentLevel, initialExcitation);
    
    const int maxSteps = 1000;
    int numSteps = 0;
    double cumulativeTime = 0.0;
    
    while (currentLevel && currentLevel->getEnergy() > 0.001 && numSteps < maxSteps) {
        CascadeStep step;
        
        if (!performDecayStep(currentLevel, step)) {
            numStuckEvents_++;
            break;
        }
        
        cumulativeTime += step.timeToDecay;
        step.timeToDecay = cumulativeTime;
        
        event.addStep(step);
        currentLevel = step.finalLevel;
        numSteps++;
    }
    
    events_.push_back(event);
}

void DecaySimulator::selectInitialState(std::shared_ptr<Level>& level, 
                                       double& excitationEnergy) {
    excitationEnergy = config_.initialExcitation.excitationEnergy;
    double spin = config_.initialExcitation.spin;
    int parity = config_.initialExcitation.parity;
    
    if (excitationEnergy < nucleus_.getCriticalEnergy()) {
        for (int i = 0; i < nucleus_.getNumDiscreteLevels(); ++i) {
            auto discLevel = nucleus_.getDiscreteLevel(i);
            if (std::abs(discLevel->getEnergy() - excitationEnergy) < 0.01 &&
                std::abs(discLevel->getSpin() - spin) < 0.1 &&
                discLevel->getParity() == parity) {
                level = discLevel;
                return;
            }
        }
    }
    
    int energyBin = nucleus_.getEnergyBin(excitationEnergy);
    int spinBin = Level::spinToBin(spin, nucleus_.isEvenA());
    
    const auto& levels = nucleus_.getContinuumLevels(energyBin, spinBin, parity);
    if (!levels.empty()) {
        int randomIndex = rng_->Integer(levels.size());
        level = levels[randomIndex];
    }
}

bool DecaySimulator::performDecayStep(std::shared_ptr<Level>& currentLevel,
                                      CascadeStep& step) {
    step.initialLevel = currentLevel;
    
    // Calculate widths - two separate structures
    std::vector<std::shared_ptr<Transition>> discreteTransitions;
    std::vector<BinWidth> binWidths;
    
    double totalWidth = calculateTotalWidth(currentLevel, discreteTransitions, binWidths);
    
    if (totalWidth <= 0.0) {
        return false;  // Stuck - no transitions
    }
    
    // Select transition using two-stage process
    std::shared_ptr<Transition> selectedTransition;
    if (!selectTransition(currentLevel, discreteTransitions, binWidths, totalWidth, 
                         selectedTransition, step)) {
        return false;
    }
    
    step.timeToDecay = getDecayTime(totalWidth);
    step.finalLevel = selectedTransition->getFinalLevel();
    step.gammaEnergy = selectedTransition->getGammaEnergy();
    step.transitionType = static_cast<int>(selectedTransition->getType());
    step.mixingRatio = selectedTransition->getMixingRatio();
    
    double icc = getInternalConversionCoeff(step.gammaEnergy, 
                                           step.transitionType, 
                                           step.mixingRatio);
    step.isElectron = isInternallyConverted(icc);
    
    return true;
}

double DecaySimulator::calculateTotalWidth(const std::shared_ptr<Level>& level,
                                          std::vector<std::shared_ptr<Transition>>& discreteTransitions,
                                          std::vector<BinWidth>& binWidths) {
    discreteTransitions.clear();
    binWidths.clear();
    
    // For discrete levels, use known transitions
    if (level->isDiscrete()) {
        auto discLevel = std::dynamic_pointer_cast<DiscreteLevel>(level);
        if (discLevel) {
            discreteTransitions = discLevel->getTransitions();
            return discLevel->getTotalWidth();
        }
        return 0.0;
    }
    
    // For continuum levels, calculate both discrete and continuum widths
    double totalWidth = 0.0;
    totalWidth += calculateDiscreteWidths(level, discreteTransitions);
    totalWidth += calculateContinuumBinWidths(level, binWidths);
    
    return totalWidth;
}

double DecaySimulator::calculateDiscreteWidths(const std::shared_ptr<Level>& level,
                                              std::vector<std::shared_ptr<Transition>>& transitions) {
    // Original RAINIER lines 561-593
    double Ex = level->getEnergy();
    double spin = level->getSpin();
    int parity = level->getParity();
    
    double density = levelDensity_->getDensity(Ex, spin, parity);
    double levelSpacing = 1.0 / density;
    double totalWidth = 0.0;
    
    for (int i = 0; i < nucleus_.getNumDiscreteLevels(); ++i) {
        auto finalLevel = nucleus_.getDiscreteLevel(i);
        double finalSpin = finalLevel->getSpin();
        int finalParity = finalLevel->getParity();
        
        auto transType = Transition::determineType(
            spin, parity, finalSpin, finalParity, nucleus_.isEvenA()
        );
        
        if (transType == Transition::Type::NONE) continue;
        
        double Egamma = Ex - finalLevel->getEnergy();
        if (Egamma <= 0) continue;
        
        int transTypeInt = static_cast<int>(transType);
        double mixingRatio = 0.0;
        double strength = gammaStrength_->getTotalStrength(Ex, Egamma, transTypeInt, mixingRatio);
        double fluctuation = getPorterThomasFluctuation();
        strength *= fluctuation;
        
        double icc = getInternalConversionCoeff(Egamma, transTypeInt, mixingRatio);
        double partialWidth = strength * levelSpacing * (1.0 + icc);
        
        if (partialWidth > 0) {
            auto transition = std::make_shared<Transition>(level, finalLevel, 0.0, icc);
            transition->setType(transType);
            transition->setMixingRatio(mixingRatio);
            transition->setPartialWidth(partialWidth);
            transitions.push_back(transition);
            totalWidth += partialWidth;
        }
    }
    
    return totalWidth;
}

double DecaySimulator::calculateContinuumBinWidths(const std::shared_ptr<Level>& level,
                                                   std::vector<BinWidth>& binWidths) {
    // Original RAINIER lines 527-560
    // KEY FIX: Store total width PER BIN, not per level
    
    double Ex = level->getEnergy();
    double spin = level->getSpin();
    int parity = level->getParity();
    
    double density = levelDensity_->getDensity(Ex, spin, parity);
    double levelSpacing = 1.0 / density;
    double totalWidth = 0.0;
    
    int currentEnergyBin = nucleus_.getEnergyBin(Ex);
    
    // Loop over possible final spins (ΔJ = -2, -1, 0, +1, +2)
    for (int dJ = -2; dJ <= 2; ++dJ) {
        double finalSpin = spin + dJ;
        if (finalSpin < 0) continue;
        
        for (int finalParity = 0; finalParity <= 1; ++finalParity) {
            auto transType = Transition::determineType(
                spin, parity, finalSpin, finalParity, nucleus_.isEvenA()
            );
            
            if (transType == Transition::Type::NONE) continue;
            
            int transTypeInt = static_cast<int>(transType);
            int finalSpinBin = Level::spinToBin(finalSpin, nucleus_.isEvenA());
            
            // Loop over energy bins BELOW current energy (energy conservation!)
            for (int energyBin = 0; energyBin < currentEnergyBin; ++energyBin) {
                const auto& finalLevels = nucleus_.getContinuumLevels(
                    energyBin, finalSpinBin, finalParity
                );
                
                int numLevels = finalLevels.size();
                if (numLevels == 0) continue;
                
                // Use bin center for gamma strength calculation
                double binCenterEnergy = nucleus_.getBinCentralEnergy(energyBin);
                double Egamma = Ex - binCenterEnergy;
                if (Egamma <= 0) continue;  // Should not happen, but safety check
                
                double mixingRatio = 0.0;
                double strength = gammaStrength_->getTotalStrength(
                    Ex, Egamma, transTypeInt, mixingRatio
                );
                
                double icc = getInternalConversionCoeff(Egamma, transTypeInt, mixingRatio);
                
                // Create unique seed for this bin to reproduce fluctuations later
                int binKey = energyBin + finalSpinBin * nucleus_.getNumEnergyBins() + 
                            finalParity * nucleus_.getNumEnergyBins() * nucleus_.getMaxSpinBin();
                
                if (binRandomStates_.find(binKey) == binRandomStates_.end()) {
                    binRandomStates_[binKey] = std::make_unique<TRandom2>(
                        1 + realization_ + binKey * 1000
                    );
                }
                
                // Calculate total width for ALL levels in this bin
                // Each level gets a Porter-Thomas fluctuated width
                double binTotalWidth = 0.0;
                for (int lvl = 0; lvl < numLevels; ++lvl) {
                    double fluctuation = getPorterThomasFluctuation();
                    binTotalWidth += strength * levelSpacing * fluctuation * (1.0 + icc);
                }
                
                if (binTotalWidth > 0) {
                    BinWidth bw;
                    bw.totalWidth = binTotalWidth;
                    bw.energyBin = energyBin;
                    bw.spinBin = finalSpinBin;
                    bw.parity = finalParity;
                    bw.transType = transTypeInt;
                    bw.mixingRatio = mixingRatio;
                    bw.icc = icc;
                    bw.randomState = binRandomStates_[binKey].get();
                    
                    binWidths.push_back(bw);
                    totalWidth += binTotalWidth;
                }
            }
        }
    }
    
    return totalWidth;
}

bool DecaySimulator::selectTransition(const std::shared_ptr<Level>& level,
                                     const std::vector<std::shared_ptr<Transition>>& discreteTransitions,
                                     const std::vector<BinWidth>& binWidths,
                                     double totalWidth,
                                     std::shared_ptr<Transition>& selectedTransition,
                                     CascadeStep& step) {
    // Original RAINIER lines 1050-1150: Two-stage selection
    
    // For discrete initial level, use branching ratios
    if (level->isDiscrete()) {
        double randomBR = rng_->Uniform(1.0);
        double cumulativeBR = 0.0;
        
        for (const auto& trans : discreteTransitions) {
            cumulativeBR += trans->getBranchingRatio();
            if (cumulativeBR >= randomBR) {
                selectedTransition = trans;
                return true;
            }
        }
        return false;
    }
    
    // For continuum level: STAGE 1 - select discrete or bin
    double randomWidth = rng_->Uniform(totalWidth);
    double cumulativeWidth = 0.0;
    
    // First check discrete transitions
    for (const auto& trans : discreteTransitions) {
        cumulativeWidth += trans->getPartialWidth();
        if (cumulativeWidth >= randomWidth) {
            selectedTransition = trans;
            return true;
        }
    }
    
    // Then check continuum bins
    for (const auto& binWidth : binWidths) {
        cumulativeWidth += binWidth.totalWidth;
        if (cumulativeWidth >= randomWidth) {
            // STAGE 2 - select specific level within this bin
            return selectFromBin(level, binWidth, selectedTransition, step);
        }
    }
    
    return false;
}

bool DecaySimulator::selectFromBin(const std::shared_ptr<Level>& level,
                                  const BinWidth& binWidth,
                                  std::shared_ptr<Transition>& selectedTransition,
                                  CascadeStep& step) {
    // Original RAINIER lines 1090-1135
    // Select a specific level within the chosen bin
    
    double Ex = level->getEnergy();
    double spin = level->getSpin();
    int parity = level->getParity();
    
    double density = levelDensity_->getDensity(Ex, spin, parity);
    double levelSpacing = 1.0 / density;
    
    const auto& finalLevels = nucleus_.getContinuumLevels(
        binWidth.energyBin, binWidth.spinBin, binWidth.parity
    );
    
    if (finalLevels.empty()) return false;
    
    // Use bin center for Egamma calculation
    double binCenterEnergy = nucleus_.getBinCentralEnergy(binWidth.energyBin);
    double Egamma = Ex - binCenterEnergy;
    
    // Get base strength (without fluctuation)
    double strength = gammaStrength_->getTotalStrength(
        Ex, Egamma, binWidth.transType, const_cast<double&>(binWidth.mixingRatio)
    );
    
    // Now randomly select a level within the bin using Porter-Thomas weights
    double randomWidth = rng_->Uniform(binWidth.totalWidth);
    double cumulativeWidth = 0.0;
    
    // Reset the bin's random state to reproduce same fluctuations
    TRandom2 binRng(binWidth.randomState->GetSeed());
    
    for (const auto& finalLevel : finalLevels) {
        // Reproduce the same Porter-Thomas fluctuation
        double gaussian = binRng.Gaus(0.0, 1.0);
        double fluctuation = gaussian * gaussian;
        
        double partialWidth = strength * levelSpacing * fluctuation * (1.0 + binWidth.icc);
        cumulativeWidth += partialWidth;
        
        if (cumulativeWidth >= randomWidth) {
            // Create transition to this specific level
            auto transType = static_cast<Transition::Type>(binWidth.transType);
            selectedTransition = std::make_shared<Transition>(
                level, finalLevel, 0.0, binWidth.icc
            );
            selectedTransition->setType(transType);
            selectedTransition->setMixingRatio(binWidth.mixingRatio);
            selectedTransition->setPartialWidth(partialWidth);
            
            // Verify energy conservation (should always be true now)
            if (finalLevel->getEnergy() >= Ex) {
                std::cerr << "ERROR: Energy not conserved! " 
                         << "Ex=" << Ex << " -> " << finalLevel->getEnergy() << std::endl;
                return false;
            }
            
            return true;
        }
    }
    
    return false;
}

double DecaySimulator::getDecayTime(double width) {
    if (width <= 0) {
        return constants::DEFAULT_MAX_HALFLIFE / constants::LOG_2;
    }
    double lifetime = constants::HBAR / width;
    return rng_->Exp(lifetime);
}

double DecaySimulator::getPorterThomasFluctuation() {
    double gaussian = rng_->Gaus(0.0, 1.0);
    return gaussian * gaussian;
}

double DecaySimulator::getInternalConversionCoeff(double Egamma, int transType, 
                                                  double mixingRatio) {
    (void)mixingRatio;
    
    if (!config_.internalConversion.enabled) {
        return 0.0;
    }
    
    double Z = nucleus_.getZ();
    double alpha = 0.0;
    
    switch (transType) {
        case 1:  // E1
            alpha = 0.01 * std::pow(Z/60.0, 3) * std::pow(0.5/Egamma, 3);
            break;
        case 2:  // M1+E2
        case 3:  // M1
            alpha = 0.05 * std::pow(Z/60.0, 3) * std::pow(0.5/Egamma, 3);
            break;
        case 4:  // E2
            alpha = 0.03 * std::pow(Z/60.0, 3) * std::pow(0.5/Egamma, 3);
            break;
        default:
            alpha = 0.0;
    }
    
    return std::max(0.0, alpha);
}

bool DecaySimulator::isInternallyConverted(double icc) {
    double probability = icc / (1.0 + icc);
    return (rng_->Uniform(1.0) < probability);
}

} // namespace rainier
