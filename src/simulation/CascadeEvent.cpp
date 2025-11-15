// CascadeEvent.cpp - Implementation of cascade event container
#include "simulation/CascadeEvent.h"

namespace rainier {

CascadeEvent::CascadeEvent() 
    : initialLevel_(nullptr), initialExcitation_(0.0) {
}

void CascadeEvent::setInitialState(std::shared_ptr<Level> level, double excitationEnergy) {
    initialLevel_ = level;
    initialExcitation_ = excitationEnergy;
}

void CascadeEvent::addStep(const CascadeStep& step) {
    steps_.push_back(step);
}

void CascadeEvent::addStep(std::shared_ptr<Level> initialLevel,
                           std::shared_ptr<Level> finalLevel,
                           double gammaEnergy,
                           double timeToDecay,
                           bool isElectron,
                           int transitionType,
                           double mixingRatio) {
    CascadeStep step;
    step.initialLevel = initialLevel;
    step.finalLevel = finalLevel;
    step.gammaEnergy = gammaEnergy;
    step.timeToDecay = timeToDecay;
    step.isElectron = isElectron;
    step.transitionType = transitionType;
    step.mixingRatio = mixingRatio;
    
    steps_.push_back(step);
}

double CascadeEvent::getTotalTime() const {
    if (steps_.empty()) {
        return 0.0;
    }
    return steps_.back().timeToDecay;
}

bool CascadeEvent::reachedGroundState() const {
    if (steps_.empty()) {
        return false;
    }
    
    auto finalLevel = steps_.back().finalLevel;
    return (finalLevel && finalLevel->getEnergy() < 0.001);  // Ground state
}

int CascadeEvent::getNumGammas() const {
    int count = 0;
    for (const auto& step : steps_) {
        if (!step.isElectron) {
            count++;
        }
    }
    return count;
}

int CascadeEvent::getNumElectrons() const {
    int count = 0;
    for (const auto& step : steps_) {
        if (step.isElectron) {
            count++;
        }
    }
    return count;
}

void CascadeEvent::clear() {
    initialLevel_ = nullptr;
    initialExcitation_ = 0.0;
    steps_.clear();
}

} // namespace rainier
