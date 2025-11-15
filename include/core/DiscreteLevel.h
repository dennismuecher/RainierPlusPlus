#ifndef RAINIER_DISCRETE_LEVEL_H
#define RAINIER_DISCRETE_LEVEL_H

#include "Level.h"
#include "Transition.h"
#include <vector>
#include <memory>

namespace rainier {

class DiscreteLevel : public Level {
public:
    DiscreteLevel(double energy, double spin, int parity, double halfLife);
    
    double getTotalWidth() const override;
    const std::vector<std::shared_ptr<Transition>>& getTransitions() const override {
        return transitions_;
    }
    bool isDiscrete() const override { return true; }
    double getHalfLife() const override { return halfLife_; }
    
    void addTransition(std::shared_ptr<Transition> transition);
    void setHalfLife(double halfLife) { halfLife_ = halfLife; }
    
    int getLevelIndex() const { return levelIndex_; }
    void setLevelIndex(int index) { levelIndex_ = index; }
    
    void normalizeBranchingRatios();

private:
    double halfLife_;
    int levelIndex_;
    std::vector<std::shared_ptr<Transition>> transitions_;
};

} // namespace rainier

#endif
