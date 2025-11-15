// LevelDensity.h - Nuclear level density models
#ifndef RAINIER_LEVEL_DENSITY_H
#define RAINIER_LEVEL_DENSITY_H

#include <memory>

namespace rainier {

// Forward declaration
class Config;

/**
 * @brief Abstract base class for level density models
 * 
 * Level density ρ(E, J, π) gives the number of levels per MeV
 * at excitation energy E with spin J and parity π.
 */
class LevelDensityModel {
public:
    virtual ~LevelDensityModel() = default;

    /**
     * @brief Calculate total level density at excitation energy
     * @param Ex Excitation energy (MeV)
     * @return Total density (MeV^-1)
     */
    virtual double getDensity(double Ex) const = 0;

    /**
     * @brief Calculate level density for specific spin and parity
     * @param Ex Excitation energy (MeV)
     * @param spin Angular momentum
     * @param parity Parity (0=negative, 1=positive)
     * @return Density (MeV^-1)
     */
    virtual double getDensity(double Ex, double spin, int parity) const = 0;

    /**
     * @brief Get effective excitation energy (E - shift)
     */
    virtual double getEffectiveEnergy(double Ex) const = 0;
};

/**
 * @brief Back-Shifted Fermi Gas (BSFG) level density model
 * 
 * ρ(E) = 1/(12√2σ) * exp(2√(aU)) / (a^(1/4) U^(5/4))
 * where U = E - E1 (effective energy)
 */
class BackShiftedFermiGas : public LevelDensityModel {
public:
    BackShiftedFermiGas(double a, double E1, bool useEnergyDependent,
                        double aAsymptotic, double shellCorrection,
                        double dampingGamma, int A);

    double getDensity(double Ex) const override;
    double getDensity(double Ex, double spin, int parity) const override;
    double getEffectiveEnergy(double Ex) const override;

    /**
     * @brief Get level density parameter a(E)
     */
    double getLevelDensityParameter(double Ex) const;

private:
    double a_;              // Level density parameter (MeV^-1)
    double E1_;             // Back-shift energy (MeV)
    bool useEnergyDependent_; // Use a(E) instead of constant a
    double aAsymptotic_;    // Asymptotic a value
    double shellCorrection_; // Shell correction W (MeV)
    double dampingGamma_;   // Damping parameter
    int A_;                 // Mass number
};

/**
 * @brief Constant Temperature Model (CTM)
 * 
 * ρ(E) = exp(U/T) / T
 * where U = E - E0
 */
class ConstantTemperature : public LevelDensityModel {
public:
    ConstantTemperature(double T, double E0);

    double getDensity(double Ex) const override;
    double getDensity(double Ex, double spin, int parity) const override;
    double getEffectiveEnergy(double Ex) const override;

private:
    double T_;   // Nuclear temperature (MeV)
    double E0_;  // Energy shift (MeV)
};

/**
 * @brief Helper class for parity distribution
 */
class ParityDistribution {
public:
    enum class Model { EQUIPARTITION, ENERGY_DEPENDENT };
    
    /**
     * @brief Get parity distribution factor
     */
    static double getFactor(double Ex, int parity, bool isEvenA,
                          Model model, double parC = 0.5, double parD = 2.0);
};

} // namespace rainier

#endif // RAINIER_LEVEL_DENSITY_H
