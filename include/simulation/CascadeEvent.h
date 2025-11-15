// CascadeEvent.h - Container for a single cascade decay event
#ifndef RAINIER_CASCADE_EVENT_H
#define RAINIER_CASCADE_EVENT_H

#include "core/Level.h"
#include "core/Transition.h"
#include <vector>
#include <memory>

namespace rainier {

/**
 * @brief Represents a single step in a cascade
 */
struct CascadeStep {
    std::shared_ptr<Level> initialLevel;
    std::shared_ptr<Level> finalLevel;
    double gammaEnergy;           // MeV
    double timeToDecay;           // fs (from creation of initial level)
    bool isElectron;              // True if internal conversion
    int transitionType;           // 0=none, 1=E1, 2=M1+E2, 3=M1, 4=E2
    double mixingRatio;           // δ² for M1+E2
    
    CascadeStep() : gammaEnergy(0), timeToDecay(0), isElectron(false),
                    transitionType(0), mixingRatio(0) {}
};

/**
 * @brief Container for a complete cascade event from initial excitation to ground state
 */
class CascadeEvent {
public:
    CascadeEvent();
    
    /**
     * @brief Set initial state
     */
    void setInitialState(std::shared_ptr<Level> level, double excitationEnergy);
    
    /**
     * @brief Add a decay step
     */
    void addStep(const CascadeStep& step);
    
    /**
     * @brief Add a decay step with explicit parameters
     */
    void addStep(std::shared_ptr<Level> initialLevel,
                 std::shared_ptr<Level> finalLevel,
                 double gammaEnergy,
                 double timeToDecay,
                 bool isElectron,
                 int transitionType,
                 double mixingRatio);
    
    /**
     * @brief Get all steps in the cascade
     */
    const std::vector<CascadeStep>& getSteps() const { return steps_; }
    
    /**
     * @brief Get number of steps
     */
    int getNumSteps() const { return steps_.size(); }
    
    /**
     * @brief Get initial excitation energy
     */
    double getInitialExcitation() const { return initialExcitation_; }
    
    /**
     * @brief Get initial level
     */
    std::shared_ptr<Level> getInitialLevel() const { return initialLevel_; }
    
    /**
     * @brief Get total cascade time
     */
    double getTotalTime() const;
    
    /**
     * @brief Check if cascade reached ground state
     */
    bool reachedGroundState() const;
    
    /**
     * @brief Get number of gammas emitted (not electrons)
     */
    int getNumGammas() const;
    
    /**
     * @brief Get number of electrons from internal conversion
     */
    int getNumElectrons() const;
    
    /**
     * @brief Clear the event
     */
    void clear();

private:
    std::shared_ptr<Level> initialLevel_;
    double initialExcitation_;
    std::vector<CascadeStep> steps_;
};

} // namespace rainier

#endif // RAINIER_CASCADE_EVENT_H
