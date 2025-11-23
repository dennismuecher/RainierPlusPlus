// SpinCutoff.h - Spin cutoff parameter models
#ifndef RAINIER_SPIN_CUTOFF_H
#define RAINIER_SPIN_CUTOFF_H

#include <memory>

namespace rainier {

// Forward declarations
class LevelDensityModel;

/**
 * @brief Spin cutoff parameter σ²
 * 
 * The spin distribution is approximately:
 * P(J) ∝ (2J+1) exp(-(J+1/2)² / 2σ²) / σ²
 */
class SpinCutoffModel {
public:
    virtual ~SpinCutoffModel() = default;

    /**
     * @brief Calculate spin cutoff parameter squared
     * @param Ex Excitation energy (MeV)
     * @return σ² (dimensionless)
     */
    virtual double getSigmaSquared(double Ex) const = 0;

    /**
     * @brief Get spin distribution factor
     * @param Ex Excitation energy
     * @param spin Angular momentum
     * @return Distribution factor (normalized)
     */
    double getSpinDistribution(double Ex, double spin) const;
};

/**
 * @brief Von Egidy & Bucurescu (2005) model
 * PRC 72, 044311 (2005)
 * σ² = 0.0146 A^(5/3) [1 + √(1 + 4aU)] / (2a)
 */
class VonEgidy05 : public SpinCutoffModel {
public:
    VonEgidy05(const LevelDensityModel* densityModel, int A, double userE1Shift = 0.0);
    double getSigmaSquared(double Ex) const override;

private:
    const LevelDensityModel* densityModel_;
    int A_;
    double userE1Shift_;  // Oslo energy shift

};

/**
 * @brief TALYS default spin cutoff model
 * Interpolates between discrete and asymptotic regions
 */

class TALYSSpinCutoff : public SpinCutoffModel {
public:
    TALYSSpinCutoff(const LevelDensityModel* densityModel, int A, double Sn,
                    double spinCutoffD, double Ed, double aAsymptotic,
                    double osloShift = 0.0);
    double getSigmaSquared(double Ex) const override;

private:
    const LevelDensityModel* densityModel_;
    int A_;
    double spinCutoffD_;
    double Ed_;
    double Sn_;
    double aAsymptotic_;
    double osloShift_;
};

/**
 * @brief Single-particle model
 * Gholami et al., PRC 75, 044308 (2007)
 * σ² = 0.1461 √(aU) A^(2/3)
 */
class SingleParticle : public SpinCutoffModel {
public:
    SingleParticle(const LevelDensityModel* densityModel, int A);
    double getSigmaSquared(double Ex) const override;

private:
    const LevelDensityModel* densityModel_;
    int A_;
};

/**
 * @brief Rigid sphere model
 * Grimes et al., PRC 10 (1974) 2373
 * σ² = 0.0145 √(U/a) A^(5/3)
 */
class RigidSphere : public SpinCutoffModel {
public:
    RigidSphere(const LevelDensityModel* densityModel, int A);
    double getSigmaSquared(double Ex) const override;

private:
    const LevelDensityModel* densityModel_;
    int A_;
};

/**
 * @brief Von Egidy & Bucurescu (2009) empirical formula
 * PRC 80, 054310 (2009)
 * σ² = 0.391 A^0.675 (E-Δ/2)^0.312
 */
class VonEgidy09 : public SpinCutoffModel {
public:
    VonEgidy09(int A, double pairingEnergy);
    double getSigmaSquared(double Ex) const override;

private:
    int A_;
    double pairingEnergy_;
};

} // namespace rainier

#endif // RAINIER_SPIN_CUTOFF_H
