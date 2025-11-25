// LevelDensity.h - Nuclear level density models
#ifndef RAINIER_LEVEL_DENSITY_H
#define RAINIER_LEVEL_DENSITY_H

#include <memory>

namespace rainier {

// Forward declarations
class Config;
class SpinCutoffModel;

/**
 * @brief Abstract base class for level density models
 */
class LevelDensityModel {
public:
    virtual ~LevelDensityModel() = default;

    virtual double getLevelDensityParameter(double Ex) const = 0;
    virtual double getDensity(double Ex) const = 0;
    virtual double getDensity(double Ex, double spin, int parity) const = 0;
    virtual double getEffectiveEnergy(double Ex) const = 0;
    
    /**
     * @brief Set spin cutoff model (called after construction)
     */
    virtual void setSpinCutoffModel(const SpinCutoffModel* spinCutoff) = 0;
};

/**
 * @brief Back-Shifted Fermi Gas (BSFG) level density model
 */
class BackShiftedFermiGas : public LevelDensityModel {
public:
    BackShiftedFermiGas(double a, double E1, bool useEnergyDependent,
                        double aAsymptotic, double shellCorrection,
                        double dampingGamma, int A);

    double getDensity(double Ex) const override;
    double getDensity(double Ex, double spin, int parity) const override;
    double getEffectiveEnergy(double Ex) const override;
    double getLevelDensityParameter(double Ex) const override;
    void setSpinCutoffModel(const SpinCutoffModel* spinCutoff) override;

private:
    double a_;
    double E1_;
    bool useEnergyDependent_;
    double aAsymptotic_;
    double shellCorrection_;
    double dampingGamma_;
    int A_;
    const SpinCutoffModel* spinCutoff_;  // Added
};

/**
 * @brief Constant Temperature Model (CTM)
 */
class ConstantTemperature : public LevelDensityModel {
public:
    ConstantTemperature(double T, double E0, int A, int Z);

    double getDensity(double Ex) const override;
    double getDensity(double Ex, double spin, int parity) const override;
    double getEffectiveEnergy(double Ex) const override;
    double getLevelDensityParameter(double Ex) const override;
    void setSpinCutoffModel(const SpinCutoffModel* spinCutoff) override;
    
private:
    double T_;
    double E0_;
    const SpinCutoffModel* spinCutoff_;  // Added
};

/**
 * @brief Helper class for parity distribution
 */
class ParityDistribution {
public:
    enum class Model { EQUIPARTITION, ENERGY_DEPENDENT };
    
    static double getFactor(double Ex, int parity, bool isEvenA,
                          Model model, double parC = 0.5, double parD = 2.0);
};

} // namespace rainier

#endif // RAINIER_LEVEL_DENSITY_H
