// DecaySimulator.cpp - Complete Monte Carlo simulation engine
// OPTIMIZED VERSION - Matches original RAINIER performance
// Physics from original RAINIER.C lines 503-1450
#include "simulation/DecaySimulator.h"
#include "simulation/CascadeEvent.h"
#include "models/LevelDensity.h"
#include "models/SpinCutoff.h"
#include "models/GammaStrength.h"
#include "utils/PhysicsConstants.h"
#include <TRandom2.h>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace rainier {

DecaySimulator::DecaySimulator(Nucleus& nucleus, const Config& config, int realization)
    : nucleus_(nucleus), 
      config_(config), 
      realization_(realization), 
      numStuckEvents_(0) {
    
    // Calculate array sizes
    maxDiscreteLevels_ = std::max(100, nucleus.getNumDiscreteLevels());
    maxContinuumBins_ = nucleus.getNumEnergyBins() * 50 * 2;  // E * J * P
    
    // Pre-allocate arrays ONCE (like original RAINIER lines 1292-1298)
    // These will be reused for every decay step
    discreteWidths_.resize(maxDiscreteLevels_, 0.0);
    continuumWidths_.resize(maxContinuumBins_, 0.0);
    
    int maxTransitions = maxDiscreteLevels_ + maxContinuumBins_;
    transitionTypes_.resize(maxTransitions, 0);
    finalLevelIndices_.resize(maxTransitions, -1);
    cumulativeWidths_.resize(maxTransitions, 0.0);
    
    // Initialize random number generator
    rng_ = std::make_unique<TRandom2>(1 + realization + config.simulation.randomSeed);
    
    // Create level density model
    levelDensity_ = std::make_unique<BackShiftedFermiGas>(
        config.levelDensity.a,
        config.levelDensity.E1,
        config.levelDensity.useEnergyDependentA,
        config.levelDensity.aAsymptotic,
        config.levelDensity.shellCorrectionW,
        config.levelDensity.dampingGamma,
        nucleus.getA()
    );
    
    // Create spin cutoff model
    spinCutoff_ = std::make_unique<VonEgidy05>(levelDensity_.get(), nucleus.getA());
    
    // Create gamma strength functions
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
    
    std::cout << "  DecaySimulator initialized (optimized array-based)" << std::endl;
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
                  << " events got stuck in isomeric states" << std::endl;
    }
}

void DecaySimulator::simulateEvent() {
    CascadeEvent event;
    
    // Select initial state
    std::shared_ptr<Level> currentLevel;
    double initialExcitation;
    selectInitialState(currentLevel, initialExcitation);
    
    event.setInitialState(currentLevel, initialExcitation);
    
    // Cascade until ground state or stuck
    const int maxSteps = 1000;
    int numSteps = 0;
    double cumulativeTime = 0.0;
    
    while (currentLevel && currentLevel->getEnergy() > 0.001 && numSteps < maxSteps) {
        CascadeStep step;
        
        if (!performDecayStep(currentLevel, step)) {
            // Stuck in isomeric state
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
    // For now, use simple single-state initialization
    excitationEnergy = config_.initialExcitation.excitationEnergy;
    double spin = config_.initialExcitation.spin;
    int parity = config_.initialExcitation.parity;
    
    // Check if in discrete region
    if (excitationEnergy < nucleus_.getCriticalEnergy()) {
        // Find matching discrete level
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
    
    // In continuum - find a level in the appropriate bin
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
    
    // Zero out arrays for reuse (faster than reallocating)
    std::fill(discreteWidths_.begin(), discreteWidths_.end(), 0.0);
    std::fill(continuumWidths_.begin(), continuumWidths_.end(), 0.0);
    
    // Calculate widths directly into pre-allocated arrays
    double totalWidth = calculateWidthsOptimized(currentLevel);
    
    if (totalWidth <= 0.0) {
        return false;
    }
    
    // Select transition using fast array-based search
    int selectedIndex = selectTransitionOptimized(totalWidth);
    if (selectedIndex < 0) {
        return false;
    }
    
    // Fill step from selected index
    fillStepFromIndex(selectedIndex, step, currentLevel);
    
    // Get decay time
    step.timeToDecay = getDecayTime(totalWidth);
    
    // Update current level
    currentLevel = step.finalLevel;
    
    return true;
}

double DecaySimulator::calculateWidthsOptimized(const std::shared_ptr<Level>& level) {
    
    // For discrete levels, use known transitions
    if (level->isDiscrete()) {
        auto discLevel = std::dynamic_pointer_cast<DiscreteLevel>(level);
        const auto& transitions = discLevel->getTransitions();
        
        double totalWidth = 0.0;
        for (size_t i = 0; i < transitions.size() && i < discreteWidths_.size(); ++i) {
            discreteWidths_[i] = transitions[i]->getPartialWidth();
            transitionTypes_[i] = static_cast<int>(transitions[i]->getType());
            finalLevelIndices_[i] = -static_cast<int>(i) - 1;  // Negative = discrete
            cumulativeWidths_[i] = totalWidth + discreteWidths_[i];
            totalWidth = cumulativeWidths_[i];
        }
        return totalWidth;
    }
    
    // For continuum levels - calculate widths
    double Ex = level->getEnergy();
    double Ji = level->getSpin();
    int pari = level->getParity();
    
    double density = levelDensity_->getDensity(Ex, Ji, pari);
    if (density <= 0) {
        return 0.0;
    }
    double levelSpacing = 1.0 / density;
    
    // Seed for Porter-Thomas per level (like original line 519)
    TRandom2 ranStr(1 + realization_ + static_cast<int>(Ex * 1000));
    
    double totalWidth = 0.0;
    int transitionIndex = 0;
    
    // Loop over spin changes (ΔJ = -2, -1, 0, +1, +2 for dipole and quadrupole)
    for (int dJ = -2; dJ <= 2; ++dJ) {
        double Jf = Ji + dJ;
        if (Jf < 0) continue;
        
        for (int parf = 0; parf <= 1; ++parf) {
            
            auto transType = Transition::determineType(Ji, pari, Jf, parf, nucleus_.isEvenA());
            if (transType == Transition::Type::NONE) continue;
            
            int transTypeInt = static_cast<int>(transType);
            int finalSpinBin = Level::spinToBin(Jf, nucleus_.isEvenA());
            
            // ================================================================
            // TRANSITIONS TO DISCRETE LEVELS (original lines 567-586)
            // ================================================================
            for (int i = 0; i < nucleus_.getNumDiscreteLevels(); ++i) {
                auto finalLevel = nucleus_.getDiscreteLevel(i);
                
                if (std::abs(finalLevel->getSpin() - Jf) > 0.1 ||
                    finalLevel->getParity() != parf) continue;
                
                double Egamma = Ex - finalLevel->getEnergy();
                if (Egamma <= 0) continue;
                
                // Calculate strength with Porter-Thomas fluctuation
                double mixingRatio = 0.0;
                double strength = gammaStrength_->getTotalStrength(Ex, Egamma, transTypeInt, mixingRatio);
                
                // Porter-Thomas: χ²(1) distribution (original line 577)
                double g = ranStr.Gaus(0.0, 1.0);
                double ptFactor = g * g;
                
                double icc = getInternalConversionCoeff(Egamma, transTypeInt, mixingRatio);
                double width = strength * levelSpacing * ptFactor * (1.0 + icc);
                
                if (width > 0 && transitionIndex < static_cast<int>(discreteWidths_.size())) {
                    discreteWidths_[transitionIndex] = width;
                    transitionTypes_[transitionIndex] = transTypeInt;
                    finalLevelIndices_[transitionIndex] = -i - 1;  // Negative for discrete
                    cumulativeWidths_[transitionIndex] = totalWidth + width;
                    totalWidth = cumulativeWidths_[transitionIndex];
                    transitionIndex++;
                }
            }
            
            // ================================================================
            // TRANSITIONS TO CONTINUUM (original lines 537-562)
            // KEY OPTIMIZATION: Calculate strength ONCE per bin, not per level!
            // ================================================================
            int currentEBin = nucleus_.getEnergyBin(Ex);
            
            for (int eBin = 0; eBin < currentEBin; ++eBin) {
                const auto& finalLevels = nucleus_.getContinuumLevels(eBin, finalSpinBin, parf);
                int nLevelsInBin = finalLevels.size();
                if (nLevelsInBin == 0) continue;
                
                double binCenterE = nucleus_.getBinCentralEnergy(eBin);
                double Egamma = Ex - binCenterE;
                if (Egamma <= 0) continue;
                
                // Calculate strength ONCE for the bin (original line 548)
                double mixingRatio = 0.0;
                double strength = gammaStrength_->getTotalStrength(Ex, Egamma, transTypeInt, mixingRatio);
                double icc = getInternalConversionCoeff(Egamma, transTypeInt, mixingRatio);
                
                // Pack bin index
                int binIndex = eBin + finalSpinBin * nucleus_.getNumEnergyBins() + 
                              parf * nucleus_.getNumEnergyBins() * 50;
                
                if (binIndex >= static_cast<int>(continuumWidths_.size())) continue;
                
                // Calculate Porter-Thomas ONCE for the entire bin (not per level!)
                double g = ranStr.Gaus(0.0, 1.0);
                double ptFactor = g * g;
                
                // Width for entire bin (strength * spacing * PT * ICC * number_of_levels)
                double widthForBin = strength * levelSpacing * ptFactor * (1.0 + icc) * nLevelsInBin;
                
                if (widthForBin > 0 && binIndex < static_cast<int>(continuumWidths_.size())) {
                    continuumWidths_[binIndex] = widthForBin;
                    int idx = maxDiscreteLevels_ + binIndex;
                    if (idx < static_cast<int>(transitionTypes_.size())) {
                        transitionTypes_[idx] = transTypeInt;
                        finalLevelIndices_[idx] = binIndex;
                        cumulativeWidths_[idx] = totalWidth + widthForBin;
                        totalWidth = cumulativeWidths_[idx];
                    }
                }
            }
        }
    }
    
    return totalWidth;
}

int DecaySimulator::selectTransitionOptimized(double totalWidth) {
    double randomWidth = rng_->Uniform(totalWidth);
    
    // Linear search through cumulative widths (could use binary search if needed)
    int maxIndex = maxDiscreteLevels_ + maxContinuumBins_;
    for (int i = 0; i < maxIndex; ++i) {
        if (cumulativeWidths_[i] >= randomWidth && cumulativeWidths_[i] > 0) {
            return i;
        }
    }
    
    return -1;  // Failed to select
}

void DecaySimulator::fillStepFromIndex(int index, CascadeStep& step,
                                       const std::shared_ptr<Level>& initialLevel) {
    if (index < maxDiscreteLevels_) {
        // Discrete transition
        int finalIndex = -finalLevelIndices_[index] - 1;
        if (finalIndex >= 0 && finalIndex < nucleus_.getNumDiscreteLevels()) {
            step.finalLevel = nucleus_.getDiscreteLevel(finalIndex);
        }
    } else {
        // Continuum transition - decode bin index
        int binIndex = finalLevelIndices_[index];
        int eBin = binIndex % nucleus_.getNumEnergyBins();
        int temp = binIndex / nucleus_.getNumEnergyBins();
        int spinBin = temp % 50;
        int parf = temp / 50;
        
        // Select random level from bin
        const auto& levels = nucleus_.getContinuumLevels(eBin, spinBin, parf);
        if (!levels.empty()) {
            int levelInBin = rng_->Integer(levels.size());
            step.finalLevel = levels[levelInBin];
        }
    }
    
    if (step.finalLevel) {
        step.gammaEnergy = initialLevel->getEnergy() - step.finalLevel->getEnergy();
        step.transitionType = transitionTypes_[index];
        step.mixingRatio = 0.0;
        
        double icc = getInternalConversionCoeff(step.gammaEnergy, step.transitionType, 0.0);
        step.isElectron = isInternallyConverted(icc);
    }
}

double DecaySimulator::getDecayTime(double width) {
    // Original RAINIER lines 601-605: τ = ℏ/Γ, then sample exponentially
    if (width <= 0) {
        return constants::DEFAULT_MAX_HALFLIFE / constants::LOG_2;
    }
    
    double lifetime = constants::HBAR / width;
    return rng_->Exp(lifetime);
}

double DecaySimulator::getPorterThomasFluctuation() {
    // Porter-Thomas distribution: χ²(ν=1) → Gaussian squared
    double gaussian = rng_->Gaus(0.0, 1.0);
    return gaussian * gaussian;
}

double DecaySimulator::getInternalConversionCoeff(double Egamma, int transType, 
                                                  double mixingRatio) {
    (void)mixingRatio;  // Will be used when full BrIcc implemented
    
    if (!config_.internalConversion.enabled) {
        return 0.0;
    }
    
    // Simplified ICC - rough approximation
    // Real implementation would call BrIcc or use tables
    double Z = nucleus_.getZ();
    double alpha = 0.0;
    
    // Rough scaling: α ∝ Z³/E³
    switch (transType) {
        case 1:  // E1
            alpha = 0.01 * std::pow(Z/60.0, 3) * std::pow(0.5/Egamma, 3);
            break;
        case 2:  // M1+E2 (use M1)
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
