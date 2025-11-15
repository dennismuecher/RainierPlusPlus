// DecaySimulator.cpp - Complete Monte Carlo simulation engine
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
    : nucleus_(nucleus), config_(config), realization_(realization), numStuckEvents_(0) {
    
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
    // Will be expanded for other modes in future
    
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
    
    // Calculate decay widths and build transition list
    std::vector<std::shared_ptr<Transition>> transitions;
    double totalWidth = calculateTotalWidth(currentLevel, transitions);
    
    if (totalWidth <= 0.0 || transitions.empty()) {
        return false;  // Stuck - no allowed transitions
    }
    
    // Select transition
    std::shared_ptr<Transition> selectedTransition;
    if (!selectTransition(currentLevel, transitions, totalWidth, selectedTransition)) {
        return false;
    }
    
    // Get decay time
    step.timeToDecay = getDecayTime(totalWidth);
    
    // Fill step information
    step.finalLevel = selectedTransition->getFinalLevel();
    step.gammaEnergy = selectedTransition->getGammaEnergy();
    step.transitionType = static_cast<int>(selectedTransition->getType());
    step.mixingRatio = selectedTransition->getMixingRatio();
    
    // Check for internal conversion
    double icc = getInternalConversionCoeff(step.gammaEnergy, 
                                           step.transitionType, 
                                           step.mixingRatio);
    step.isElectron = isInternallyConverted(icc);
    
    return true;
}

double DecaySimulator::calculateTotalWidth(const std::shared_ptr<Level>& level,
                                          std::vector<std::shared_ptr<Transition>>& transitions) {
    transitions.clear();
    
    // For discrete levels, use known transitions
    if (level->isDiscrete()) {
        auto discLevel = std::dynamic_pointer_cast<DiscreteLevel>(level);
        if (discLevel) {
            transitions = discLevel->getTransitions();
            return discLevel->getTotalWidth();
        }
        return 0.0;
    }
    
    // For continuum levels, calculate widths
    return calculateContinuumWidth(level, transitions);
}

double DecaySimulator::calculateContinuumWidth(const std::shared_ptr<Level>& level,
                                              std::vector<std::shared_ptr<Transition>>& transitions) {
    // Original RAINIER lines 503-600
    
    double Ex = level->getEnergy();
    double spin = level->getSpin();
    int parity = level->getParity();
    
    // Calculate level spacing
    double density = levelDensity_->getDensity(Ex, spin, parity);
    double levelSpacing = 1.0 / density;
    
    double totalWidth = 0.0;
    
    // Loop over possible final states (ΔJ = -2, -1, 0, +1, +2)
    for (int dJ = -2; dJ <= 2; ++dJ) {
        double finalSpin = spin + dJ;
        if (finalSpin < 0) continue;
        
        for (int finalParity = 0; finalParity <= 1; ++finalParity) {
            // Determine transition type
            auto transType = Transition::determineType(
                spin, parity, finalSpin, finalParity, nucleus_.isEvenA()
            );
            
            if (transType == Transition::Type::NONE) continue;
            
            int transTypeInt = static_cast<int>(transType);
            int finalSpinBin = Level::spinToBin(finalSpin, nucleus_.isEvenA());
            
            // Transitions to discrete levels
            for (int i = 0; i < nucleus_.getNumDiscreteLevels(); ++i) {
                auto finalLevel = nucleus_.getDiscreteLevel(i);
                
                if (std::abs(finalLevel->getSpin() - finalSpin) < 0.1 &&
                    finalLevel->getParity() == finalParity) {
                    
                    double Egamma = Ex - finalLevel->getEnergy();
                    if (Egamma <= 0) continue;
                    
                    // Get gamma strength with fluctuation
                    double mixingRatio = 0.0;
                    double strength = gammaStrength_->getTotalStrength(
                        Ex, Egamma, transTypeInt, mixingRatio
                    );
                    
                    // Apply Porter-Thomas fluctuation
                    double fluctuation = getPorterThomasFluctuation();
                    strength *= fluctuation;
                    
                    // Calculate partial width
                    double icc = getInternalConversionCoeff(Egamma, transTypeInt, mixingRatio);
                    double partialWidth = strength * levelSpacing * (1.0 + icc);
                    
                    if (partialWidth > 0) {
                        auto transition = std::make_shared<Transition>(
                            level, finalLevel, 0.0, icc
                        );
                        transition->setType(transType);
                        transition->setMixingRatio(mixingRatio);
                        transition->setPartialWidth(partialWidth);
                        transitions.push_back(transition);
                        totalWidth += partialWidth;
                    }
                }
            }
            
            // Transitions to continuum levels
            int currentEnergyBin = nucleus_.getEnergyBin(Ex);
            for (int energyBin = 0; energyBin < currentEnergyBin; ++energyBin) {
                const auto& finalLevels = nucleus_.getContinuumLevels(
                    energyBin, finalSpinBin, finalParity
                );
                
                if (finalLevels.empty()) continue;
                
                double binCenterEnergy = nucleus_.getBinCentralEnergy(energyBin);
                double Egamma = Ex - binCenterEnergy;
                if (Egamma <= 0) continue;
                
                // Calculate strength once per bin
                double mixingRatio = 0.0;
                double strength = gammaStrength_->getTotalStrength(
                    Ex, Egamma, transTypeInt, mixingRatio
                );
                
                double icc = getInternalConversionCoeff(Egamma, transTypeInt, mixingRatio);
                
                // Add contribution for each level in bin
                for (const auto& finalLevel : finalLevels) {
                    double fluctuation = getPorterThomasFluctuation();
                    double partialWidth = strength * levelSpacing * fluctuation * (1.0 + icc);
                    
                    if (partialWidth > 0) {
                        auto transition = std::make_shared<Transition>(
                            level, finalLevel, 0.0, icc
                        );
                        transition->setType(transType);
                        transition->setMixingRatio(mixingRatio);
                        transition->setPartialWidth(partialWidth);
                        transitions.push_back(transition);
                        totalWidth += partialWidth;
                    }
                }
            }
        }
    }
    
    return totalWidth;
}

bool DecaySimulator::selectTransition(const std::shared_ptr<Level>& level,
                                     const std::vector<std::shared_ptr<Transition>>& transitions,
                                     double totalWidth,
                                     std::shared_ptr<Transition>& selectedTransition) {
    // For discrete levels, use branching ratios
    if (level->isDiscrete()) {
        double randomBR = rng_->Uniform(1.0);
        double cumulativeBR = 0.0;
        
        for (const auto& trans : transitions) {
            cumulativeBR += trans->getBranchingRatio();
            if (cumulativeBR >= randomBR) {
                selectedTransition = trans;
                return true;
            }
        }
        return false;
    }
    
    // For continuum levels, use widths
    double randomWidth = rng_->Uniform(totalWidth);
    double cumulativeWidth = 0.0;
    
    for (const auto& trans : transitions) {
        cumulativeWidth += trans->getPartialWidth();
        if (cumulativeWidth >= randomWidth) {
            selectedTransition = trans;
            return true;
        }
    }
    
    return false;
}

double DecaySimulator::getDecayTime(double width) {
    // Original RAINIER lines 601-605
    // τ = ℏ/Γ, then sample exponentially
    if (width <= 0) {
        return constants::DEFAULT_MAX_HALFLIFE / constants::LOG_2;
    }
    
    double lifetime = constants::HBAR / width;
    return rng_->Exp(lifetime);
}

double DecaySimulator::getPorterThomasFluctuation() {
    // Original RAINIER lines 411-422
    // Porter-Thomas distribution: χ²(ν=1) → Gaussian squared
    
    double gaussian = rng_->Gaus(0.0, 1.0);
    return gaussian * gaussian;
}

double DecaySimulator::getInternalConversionCoeff(double Egamma, int transType, 
                                                  double mixingRatio) {
    // Simplified ICC - full BrIcc integration would go here
    // For now, use simple approximations
    
    (void)mixingRatio;  // Will be used when full ICC implemented
    
    if (!config_.internalConversion.enabled) {
        return 0.0;
    }
    
    // Very crude approximation
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
