#ifndef RAINIER_CONTINUUM_LEVEL_H
#define RAINIER_CONTINUUM_LEVEL_H

#include "Level.h"
#include "Transition.h"
#include <vector>
#include <memory>

namespace rainier {

class ContinuumLevel : public Level {
public:
    
    ContinuumLevel(double energy, double spin, int parity,
                   int energyBin, int spinBin, int levelInBin);
    
    double getTotalWidth() const override { return totalWidth_; }
    const std::vector<std::shared_ptr<Transition>>& getTransitions() const override {
        return transitions_;
    }
    bool isDiscrete() const override { return false; }
    double getHalfLife() const override;
    
    void setTotalWidth(double width) { totalWidth_ = width; }
    void addTransition(std::shared_ptr<Transition> transition);
    
    int getEnergyBin() const { return energyBin_; }
    int getLevelInBin() const { return levelInBin_; }
    int getSpinBin() const { return spinBin_; }

 
    void clearTransitions() { transitions_.clear(); totalWidth_ = 0.0; }

private:
    double totalWidth_;
    int energyBin_;
    int spinBin_;
    int levelInBin_;
    
    
    std::vector<std::shared_ptr<Transition>> transitions_;
};

} // namespace rainier

#endif
