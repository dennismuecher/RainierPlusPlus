#include "core/Level.h"
#include "core/DiscreteLevel.h"
#include "core/ContinuumLevel.h"
#include "core/Transition.h"
#include "utils/PhysicsConstants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace rainier {

// ===== Level Base Class =====
Level::Level(double energy, double spin, int parity)
    : energy_(energy), spin_(spin), parity_(parity) {
    
    if (energy < 0.0) {
        throw std::invalid_argument("Energy must be non-negative");
    }
    if (spin < 0.0) {
        throw std::invalid_argument("Spin must be non-negative");
    }
    if (parity != 0 && parity != 1) {
        throw std::invalid_argument("Parity must be 0 or 1");
    }
    
    spinBin_ = static_cast<int>(spin);
}

double Level::getLifetime() const {
    double width = getTotalWidth();
    if (width <= 0.0) {
        return constants::DEFAULT_MAX_HALFLIFE / constants::LOG_2;
    }
    return constants::HBAR / width;
}

int Level::spinToBin(double spin, bool isEvenA) {
    if (isEvenA) {
        return static_cast<int>(spin + 0.001);
    } else {
        return static_cast<int>(spin - 0.499);
    }
}

double Level::binToSpin(int bin, bool isEvenA) {
    if (isEvenA) {
        return static_cast<double>(bin);
    } else {
        return static_cast<double>(bin) + 0.5;
    }
}

// ===== DiscreteLevel =====
DiscreteLevel::DiscreteLevel(double energy, double spin, int parity, double halfLife)
    : Level(energy, spin, parity), halfLife_(halfLife), levelIndex_(-1) {
    
    if (halfLife <= 0.0) {
        halfLife_ = constants::DEFAULT_MAX_HALFLIFE;
    }
}

double DiscreteLevel::getTotalWidth() const {
    if (halfLife_ >= constants::DEFAULT_MAX_HALFLIFE) {
        return 0.0;
    }
    return constants::HBAR * constants::LOG_2 / halfLife_;
}

void DiscreteLevel::addTransition(std::shared_ptr<Transition> transition) {
    if (!transition) {
        throw std::invalid_argument("Cannot add null transition");
    }
    transitions_.push_back(transition);
}

void DiscreteLevel::normalizeBranchingRatios() {
    double totalBR = 0.0;
    for (const auto& trans : transitions_) {
        totalBR += trans->getBranchingRatio();
    }
    if (totalBR <= 0.0) {
        throw std::runtime_error("Total branching ratio is zero");
    }
    for (auto& trans : transitions_) {
        trans->setBranchingRatio(trans->getBranchingRatio() / totalBR);
    }
}

// ===== ContinuumLevel =====
ContinuumLevel::ContinuumLevel(double energy, double spin, int parity,
                               int energyBin, int levelInBin)
    : Level(energy, spin, parity), 
      totalWidth_(0.0), energyBin_(energyBin), levelInBin_(levelInBin) {
}

double ContinuumLevel::getHalfLife() const {
    if (totalWidth_ <= 0.0) {
        return constants::DEFAULT_MAX_HALFLIFE;
    }
    return constants::HBAR * constants::LOG_2 / totalWidth_;
}

void ContinuumLevel::addTransition(std::shared_ptr<Transition> transition) {
    if (!transition) {
        throw std::invalid_argument("Cannot add null transition");
    }
    transitions_.push_back(transition);
}

// ===== Transition =====
Transition::Transition(std::shared_ptr<Level> initialLevel,
                       std::shared_ptr<Level> finalLevel,
                       double branchingRatio, double icc)
    : initialLevel_(initialLevel), finalLevel_(finalLevel),
      branchingRatio_(branchingRatio), icc_(icc),
      type_(Type::NONE), mixingRatio_(0.0), partialWidth_(0.0) {
    
    if (!initialLevel || !finalLevel) {
        throw std::invalid_argument("Transition requires valid levels");
    }
    if (branchingRatio < 0.0 || branchingRatio > 1.0) {
        throw std::invalid_argument("BR must be 0-1");
    }
    if (icc < 0.0) {
        throw std::invalid_argument("ICC must be non-negative");
    }
}

double Transition::getGammaEnergy() const {
    return initialLevel_->getEnergy() - finalLevel_->getEnergy();
}

double Transition::conversionProbability() const {
    return icc_ / (1.0 + icc_);
}

Transition::Type Transition::determineType(double spinI, int parityI,
                                          double spinF, int parityF,
                                          bool isEvenA) {
    int spinBinI = Level::spinToBin(spinI, isEvenA);
    int spinBinF = Level::spinToBin(spinF, isEvenA);
    int dSpinBin = std::abs(spinBinF - spinBinI);
    int dParity = std::abs(parityF - parityI);
    
    if (spinBinF < 0 || spinBinI < 0) return Type::NONE;
    if (spinBinI == 0 && spinBinF == 0) return Type::NONE;
    
    if (dSpinBin == 0) {
        if (spinBinI > 0 && spinBinF > 0) {
            return (dParity == 0) ? Type::M1_E2 : Type::E1;
        }
        return Type::NONE;
    } else if (dSpinBin == 1) {
        if (spinBinI > 0 && spinBinF > 0) {
            return (dParity == 0) ? Type::M1_E2 : Type::E1;
        } else {
            return (dParity == 0) ? Type::M1_PURE : Type::E1;
        }
    } else if (dSpinBin == 2) {
        return (dParity == 0) ? Type::E2_PURE : Type::NONE;
    }
    
    return Type::NONE;
}

} // namespace rainier
