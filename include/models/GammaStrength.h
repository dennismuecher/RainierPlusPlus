// GammaStrength.h - Gamma-ray strength function models
#ifndef RAINIER_GAMMA_STRENGTH_H
#define RAINIER_GAMMA_STRENGTH_H

#include <vector>
#include <memory>

namespace rainier {

/**
 * @brief Resonance parameters for Lorentzian models
 */
struct Resonance {
    double energy;  // Resonance energy (MeV)
    double width;   // Resonance width (MeV)
    double sigma;   // Cross section (mb)
};

/**
 * @brief Base class for gamma strength functions
 * 
 * Returns f_XL(Eγ) × Eγ^(2L+1) where L is the multipolarity
 */
class GammaStrengthFunction {
public:
    virtual ~GammaStrengthFunction() = default;
    
    /**
     * @brief Get strength × Eγ^(2L+1)
     * @param Ex Excitation energy of initial state (MeV)
     * @param Egamma Gamma-ray energy (MeV)
     * @return Strength × Eγ^(2L+1) (MeV^-(2L+2))
     */
    virtual double getStrength(double Ex, double Egamma) const = 0;
};

// ============================================================================
// E1 Strength Functions
// ============================================================================

/**
 * @brief Generalized Lorentzian (GenLor) E1 strength
 * Kopecky & Uhl, PRC 41, 1941 (1990)
 * Default model for RAINIER
 */
class E1GenLorentz : public GammaStrengthFunction {
public:
    E1GenLorentz(const std::vector<Resonance>& resonances, double constantT);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    std::vector<Resonance> resonances_;
    double constantT_;  // Constant temperature (MeV) or 0 for temperature from Ex
    
    double getTemperature(double Ex) const;
};

/**
 * @brief Enhanced Generalized Lorentzian (EGLO) for A >= 148
 * Kopecky et al., NPA 567, 576 (1994)
 */
class E1EGLO : public GammaStrengthFunction {
public:
    E1EGLO(const std::vector<Resonance>& resonances, double constantT, int A);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    std::vector<Resonance> resonances_;
    double constantT_;
    int A_;
    
    double getTemperature(double Ex) const;
    double getEnhancementFactor(double Egamma) const;
};

/**
 * @brief Kadmenskij, Markushev, Furman (KMF) model
 */
class E1KMF : public GammaStrengthFunction {
public:
    E1KMF(const std::vector<Resonance>& resonances);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    std::vector<Resonance> resonances_;
};

/**
 * @brief Standard Lorentzian E1
 */
class E1StandardLorentz : public GammaStrengthFunction {
public:
    E1StandardLorentz(const std::vector<Resonance>& resonances);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    std::vector<Resonance> resonances_;
};

// ============================================================================
// M1 Strength Functions
// ============================================================================

/**
 * @brief Standard Lorentzian M1 strength
 */
class M1StandardLorentz : public GammaStrengthFunction {
public:
    M1StandardLorentz(const std::vector<Resonance>& resonances,
                      bool useUpbend = false,
                      double upbendConstant = 0.1,
                      double upbendExponent = 0.5);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    std::vector<Resonance> resonances_;
    bool useUpbend_;
    double upbendConstant_;
    double upbendExponent_;
};

/**
 * @brief Single-particle M1 strength
 */
class M1SingleParticle : public GammaStrengthFunction {
public:
    M1SingleParticle(double sigma);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    double sigma_;  // Single-particle strength (mb)
};

// ============================================================================
// E2 Strength Functions
// ============================================================================

/**
 * @brief Standard Lorentzian E2 strength
 */
class E2StandardLorentz : public GammaStrengthFunction {
public:
    E2StandardLorentz(double energy, double width, double sigma);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    Resonance resonance_;
};

/**
 * @brief Single-particle E2 strength
 */
class E2SingleParticle : public GammaStrengthFunction {
public:
    E2SingleParticle(double sigma);
    double getStrength(double Ex, double Egamma) const override;
    
private:
    double sigma_;
};

// ============================================================================
// Combined Strength Function
// ============================================================================

/**
 * @brief Combines E1, M1, E2 strengths for mixed transitions
 */
class CombinedStrength {
public:
    CombinedStrength(std::unique_ptr<GammaStrengthFunction> e1,
                     std::unique_ptr<GammaStrengthFunction> m1,
                     std::unique_ptr<GammaStrengthFunction> e2);
    
    /**
     * @brief Get total strength for a transition
     * @param Ex Excitation energy
     * @param Egamma Gamma energy
     * @param transType Transition type (E1, M1, E2, M1+E2)
     * @param mixingRatio Output: δ² for M1+E2 transitions
     * @return Total strength × Eγ^(2L+1)
     */
    double getTotalStrength(double Ex, double Egamma, int transType, 
                           double& mixingRatio) const;
    
private:
    std::unique_ptr<GammaStrengthFunction> e1_;
    std::unique_ptr<GammaStrengthFunction> m1_;
    std::unique_ptr<GammaStrengthFunction> e2_;
};

} // namespace rainier

#endif // RAINIER_GAMMA_STRENGTH_H
