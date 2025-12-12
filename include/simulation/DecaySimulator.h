// DecaySimulator.h - Monte Carlo cascade decay simulator
// Fixed version with proper energy conservation (RAINIER.C lines 503-1450)
#ifndef RAINIER_DECAY_SIMULATOR_H
#define RAINIER_DECAY_SIMULATOR_H

#include "core/Nucleus.h"
#include "core/Level.h"
#include "Config.h"
#include "simulation/CascadeEvent.h"
#include <vector>
#include <memory>
#include <map>

class TRandom2;

namespace rainier {

class LevelDensityModel;
class SpinCutoffModel;
class CombinedStrength;

// Structure to hold width information for a specific E-J-π bin
struct BinWidth {
    double totalWidth;              // Total width for all levels in this bin
    int energyBin;
    int spinBin;
    int parity;
    int transType;
    double mixingRatio;
    double icc;
    TRandom2* randomState;         // To reproduce Porter-Thomas fluctuations
    
    BinWidth() : totalWidth(0), energyBin(0), spinBin(0), parity(0), 
                 transType(0), mixingRatio(0), icc(0), randomState(nullptr) {}
};

class DecaySimulator {
public:
    DecaySimulator(Nucleus& nucleus, const Config& config, int realization);
    ~DecaySimulator();
    
    void run();
    
    const std::vector<CascadeEvent>& getEvents() const { return events_; }
    const CascadeEvent& getEvent(int index) const { return events_[index]; }
    int getNumEvents() const { return events_.size(); }
	
	Nucleus& getNucleus()const  {return nucleus_;}
    
    const LevelDensityModel& getLevelDensityModel() const { return *levelDensity_; }
    const SpinCutoffModel& getSpinCutoffModel() const { return *spinCutoff_; }

private:
    void simulateEvent();
    void selectInitialState(std::shared_ptr<Level>& level, double& excitationEnergy);
    bool performDecayStep(std::shared_ptr<Level>& currentLevel, CascadeStep& step);
    
    // Two-stage width calculation like original RAINIER
    double calculateTotalWidth(const std::shared_ptr<Level>& level,
                              std::vector<std::shared_ptr<Transition>>& discreteTransitions,
                              std::vector<BinWidth>& binWidths);
    
    double calculateDiscreteWidths(const std::shared_ptr<Level>& level,
                                  std::vector<std::shared_ptr<Transition>>& transitions);
    
    double calculateContinuumBinWidths(const std::shared_ptr<Level>& level,
                                       std::vector<BinWidth>& binWidths);
    
    // Two-stage selection: bin first, then level
    bool selectTransition(const std::shared_ptr<Level>& level,
                         const std::vector<std::shared_ptr<Transition>>& discreteTransitions,
                         const std::vector<BinWidth>& binWidths,
                         double totalWidth,
                         std::shared_ptr<Transition>& selectedTransition,
                         CascadeStep& step);
    
    bool selectFromBin(const std::shared_ptr<Level>& level,
                      const BinWidth& binWidth,
                      std::shared_ptr<Transition>& selectedTransition,
                      CascadeStep& step);
    
    double getDecayTime(double width);
    double getPorterThomasFluctuation();
    double getInternalConversionCoeff(double Egamma, int transType, double mixingRatio);
    bool isInternallyConverted(double icc);
    
    void selectBetaDecayState(std::shared_ptr<Level>& level, double& excitationEnergy);
    std::vector<double> getAllowedSpins(double parentSpin) const;
    
    Nucleus& nucleus_;
    Config config_;
    int realization_;
    
    std::unique_ptr<TRandom2> rng_;
    
    const LevelDensityModel* levelDensity_;
    const SpinCutoffModel* spinCutoff_;

    std::unique_ptr<CombinedStrength> gammaStrength_;
    
    std::vector<CascadeEvent> events_;
    std::map<int, std::unique_ptr<TRandom2>> binRandomStates_;
    
    int numStuckEvents_;
};

} // namespace rainier

#endif
