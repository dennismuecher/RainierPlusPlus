// DecaySimulator.h - Main Monte Carlo decay simulation engine
// OPTIMIZED VERSION - Uses arrays like original RAINIER for 10x speedup
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
 * 
 * PERFORMANCE NOTE: This implementation uses pre-allocated arrays (like original RAINIER)
 * instead of creating Transition objects. This provides ~10x speedup by:
 * - Eliminating object construction/destruction overhead
 * - Better cache locality with contiguous arrays
 * - Reduced memory allocations in hot loops
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
	
	Nucleus& getNucleus()const  {return nucleus_;}

private:
    // Core simulation methods
    void simulateEvent();
    void selectInitialState(std::shared_ptr<Level>& level, double& excitationEnergy);
    bool performDecayStep(std::shared_ptr<Level>& currentLevel, CascadeStep& step);
    
    // Optimized width calculation - uses pre-allocated arrays
    double calculateWidthsOptimized(const std::shared_ptr<Level>& level);
    
    // Fast transition selection using cumulative width arrays
    int selectTransitionOptimized(double totalWidth);
    
    // Fill cascade step from selected transition index
    void fillStepFromIndex(int index, CascadeStep& step, 
                          const std::shared_ptr<Level>& initialLevel);
    
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
    
    // OPTIMIZATION: Pre-allocated arrays (like original RAINIER lines 1292-1298)
    // These are reused for every decay step to avoid allocations
    std::vector<double> discreteWidths_;      // Width to each discrete level
    std::vector<double> continuumWidths_;     // Width to each continuum bin
    
    std::vector<int> transitionTypes_;        // Store transition multipolarity
    std::vector<int> finalLevelIndices_;      // Encoded final level info
    std::vector<double> cumulativeWidths_;    // For fast binary/linear search
    
    int maxDiscreteLevels_;
    int maxContinuumBins_;
    
	
    // Statistics
    int numStuckEvents_;
};

} // namespace rainier

#endif // RAINIER_DECAY_SIMULATOR_H
