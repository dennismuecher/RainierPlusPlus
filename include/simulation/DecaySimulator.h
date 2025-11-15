// DecaySimulator.h - Main Monte Carlo decay simulation engine
#ifndef RAINIER_DECAY_SIMULATOR_H
#define RAINIER_DECAY_SIMULATOR_H

#include "core/Nucleus.h"
#include "simulation/CascadeEvent.h"
#include "Config.h"
#include <memory>
#include <vector>

// Forward declarations
class TRandom2;

namespace rainier {

// Forward declarations
class LevelDensityModel;
class SpinCutoffModel;
class GammaStrengthFunction;
class CombinedStrength;

/**
 * @brief Main Monte Carlo simulation engine for gamma cascades
 * 
 * Simulates complete gamma-ray cascades from initial excitation to ground state.
 * Uses statistical models for continuum region and known data for discrete levels.
 */
class DecaySimulator {
public:
    /**
     * @brief Construct simulator
     * @param nucleus Nuclear level scheme
     * @param config Simulation configuration
     * @param realization Realization number (for random seed)
     */
    DecaySimulator(Nucleus& nucleus, const Config& config, int realization);
    
    ~DecaySimulator();
    
    /**
     * @brief Run the complete simulation
     */
    void run();
    
    /**
     * @brief Get all simulated events
     */
    const std::vector<CascadeEvent>& getEvents() const { return events_; }
    
    /**
     * @brief Get a specific event
     */
    const CascadeEvent& getEvent(int index) const { return events_[index]; }
    
    /**
     * @brief Get number of events simulated
     */
    int getNumEvents() const { return events_.size(); }

private:
    // Core simulation methods
    void simulateEvent();
    void selectInitialState(std::shared_ptr<Level>& level, double& excitationEnergy);
    bool performDecayStep(std::shared_ptr<Level>& currentLevel, 
                         CascadeStep& step);
    
    // Width calculations
    double calculateTotalWidth(const std::shared_ptr<Level>& level,
                              std::vector<std::shared_ptr<Transition>>& transitions);
    
    double calculateContinuumWidth(const std::shared_ptr<Level>& level,
                                  std::vector<std::shared_ptr<Transition>>& transitions);
    
    // Transition selection
    bool selectTransition(const std::shared_ptr<Level>& level,
                         const std::vector<std::shared_ptr<Transition>>& transitions,
                         double totalWidth,
                         std::shared_ptr<Transition>& selectedTransition);
    
    // Helper methods
    double getDecayTime(double width);
    double getPorterThomasFluctuation();
    double getInternalConversionCoeff(double Egamma, int transType, double mixingRatio);
    bool isInternallyConverted(double icc);
    
    // Member variables
    Nucleus& nucleus_;
    Config config_;
    int realization_;
    
    std::unique_ptr<TRandom2> rng_;
    std::unique_ptr<LevelDensityModel> levelDensity_;
    std::unique_ptr<SpinCutoffModel> spinCutoff_;
    std::unique_ptr<CombinedStrength> gammaStrength_;
    
    std::vector<CascadeEvent> events_;
    
    // Statistics
    int numStuckEvents_;
};

} // namespace rainier

#endif // RAINIER_DECAY_SIMULATOR_H