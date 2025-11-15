#ifndef RAINIER_LEVEL_H
#define RAINIER_LEVEL_H

#include <vector>
#include <memory>

namespace rainier {

class Transition;

class Level {
public:
    Level(double energy, double spin, int parity);
    virtual ~Level() = default;

    double getEnergy() const { return energy_; }
    double getSpin() const { return spin_; }
    int getParity() const { return parity_; }
    int getSpinBin() const { return spinBin_; }
    
    virtual double getTotalWidth() const = 0;
    virtual const std::vector<std::shared_ptr<Transition>>& getTransitions() const = 0;
    virtual bool isDiscrete() const = 0;
    virtual double getHalfLife() const = 0;

    double getLifetime() const;
    
    static int spinToBin(double spin, bool isEvenA);
    static double binToSpin(int bin, bool isEvenA);

protected:
    double energy_;
    double spin_;
    int parity_;
    int spinBin_;
};

} // namespace rainier

#endif
