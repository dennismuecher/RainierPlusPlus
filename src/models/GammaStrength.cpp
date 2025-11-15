// GammaStrength.cpp - Implementation of gamma strength functions
// Physics from original RAINIER.C lines 363-472
#include "models/GammaStrength.h"
#include "utils/PhysicsConstants.h"
#include <cmath>

namespace rainier {

// ============================================================================
// E1 Generalized Lorentzian
// ============================================================================

E1GenLorentz::E1GenLorentz(const std::vector<Resonance>& resonances, double constantT)
    : resonances_(resonances), constantT_(constantT) {
}

double E1GenLorentz::getTemperature(double Ex) const {
    if (constantT_ > 0.0) {
        return constantT_;
    }
    // Temperature from excitation energy (simple approximation)
    // T ≈ √(E/a) where a ≈ A/8
    return std::sqrt(Ex / 15.0);  // Crude estimate
}

double E1GenLorentz::getStrength(double Ex, double Egamma) const {
    // Generalized Lorentzian E1 strength
    // Original RAINIER lines 376-385
    
    double T = getTemperature(Ex - Egamma);  // Temperature at final state
    double strength = 0.0;
    
    for (const auto& res : resonances_) {
        // Energy-dependent width
        // Γ(Eγ) = Γ₀ (Eγ² + 4π²T²) / E₀²
        double Gamma = res.width * (Egamma * Egamma + 
                                    constants::FOUR_PI_SQUARED * T * T) /
                      (res.energy * res.energy);
        
        // Term 1: Lorentzian
        double term1 = Egamma * Gamma / 
                      (std::pow(Egamma * Egamma - res.energy * res.energy, 2) + 
                       std::pow(Egamma * Gamma, 2));
        
        // Term 2: Zero-energy limit
        // 0.7 is the Fermi liquid parameter
        double term2 = 0.7 * res.width * constants::FOUR_PI_SQUARED * T * T /
                      std::pow(res.energy, 5);
        
        strength += constants::K_X1 * res.sigma * res.width * (term1 + term2);
    }
    
    return strength * std::pow(Egamma, 3);  // × Eγ³
}

// ============================================================================
// E1 EGLO (Enhanced Generalized Lorentzian)
// ============================================================================

E1EGLO::E1EGLO(const std::vector<Resonance>& resonances, double constantT, int A)
    : resonances_(resonances), constantT_(constantT), A_(A) {
}

double E1EGLO::getTemperature(double Ex) const {
    if (constantT_ > 0.0) {
        return constantT_;
    }
    return std::sqrt(Ex / 15.0);
}

double E1EGLO::getEnhancementFactor(double Egamma) const {
    // EGLO enhancement for A >= 148
    // Original RAINIER lines 387-393
    
    if (A_ < 148) {
        return 1.0;  // No enhancement
    }
    
    double E0 = 4.5;  // Reference energy (MeV)
    double k = 1.0 + 0.09 * std::pow(A_ - 148, 2) * std::exp(-0.18 * (A_ - 148));
    
    double chi = k + (1.0 - k) * (Egamma - E0) / (resonances_[0].energy - E0);
    
    return chi;
}

double E1EGLO::getStrength(double Ex, double Egamma) const {
    double T = getTemperature(Ex - Egamma);
    double strength = 0.0;
    
    for (const auto& res : resonances_) {
        double Gamma = res.width * (Egamma * Egamma + 
                                    constants::FOUR_PI_SQUARED * T * T) /
                      (res.energy * res.energy);
        
        // Apply EGLO enhancement
        double chi = getEnhancementFactor(Egamma);
        Gamma *= chi;
        
        double term1 = Egamma * Gamma / 
                      (std::pow(Egamma * Egamma - res.energy * res.energy, 2) + 
                       std::pow(Egamma * Gamma, 2));
        
        double term2 = 0.7 * res.width * constants::FOUR_PI_SQUARED * T * T /
                      std::pow(res.energy, 5);
        
        strength += constants::K_X1 * res.sigma * res.width * (term1 + term2);
    }
    
    return strength * std::pow(Egamma, 3);
}

// ============================================================================
// E1 KMF (Kadmenskij, Markushev, Furman)
// ============================================================================

E1KMF::E1KMF(const std::vector<Resonance>& resonances)
    : resonances_(resonances) {
}

double E1KMF::getStrength(double Ex, double Egamma) const {
    // KMF model
    // Original RAINIER lines 395-398
    
    double strength = 0.0;
    
    for (const auto& res : resonances_) {
        double term = 0.7 * res.energy * res.width /
                     std::pow(Egamma * Egamma - res.energy * res.energy, 2);
        
        strength += constants::K_X1 * res.sigma * res.width * term;
    }
    
    return strength * std::pow(Egamma, 3);
}

// ============================================================================
// E1 Standard Lorentzian
// ============================================================================

E1StandardLorentz::E1StandardLorentz(const std::vector<Resonance>& resonances)
    : resonances_(resonances) {
}

double E1StandardLorentz::getStrength(double Ex, double Egamma) const {
    // Standard Lorentzian
    // Original RAINIER lines 404-408
    
    (void)Ex;  // Not used in standard Lorentzian
    
    double strength = 0.0;
    
    for (const auto& res : resonances_) {
        double term = Egamma * res.width /
                     (std::pow(Egamma * Egamma - res.energy * res.energy, 2) +
                      std::pow(Egamma * res.width, 2));
        
        strength += constants::K_X1 * res.sigma * res.width * term;
    }
    
    return strength * std::pow(Egamma, 3);
}

// ============================================================================
// M1 Standard Lorentzian
// ============================================================================

M1StandardLorentz::M1StandardLorentz(const std::vector<Resonance>& resonances,
                                     bool useUpbend,
                                     double upbendConstant,
                                     double upbendExponent)
    : resonances_(resonances),
      useUpbend_(useUpbend),
      upbendConstant_(upbendConstant),
      upbendExponent_(upbendExponent) {
}

double M1StandardLorentz::getStrength(double Ex, double Egamma) const {
    // M1 Standard Lorentzian
    // Original RAINIER lines 425-434
    
    (void)Ex;  // Not used
    
    double strength = 0.0;
    
    for (const auto& res : resonances_) {
        double term = Egamma * res.width * res.width /
                     (std::pow(Egamma * Egamma - res.energy * res.energy, 2) +
                      std::pow(Egamma * res.width, 2));
        
        strength += constants::K_X1 * res.sigma * term;
    }
    
    // Add upbend if requested
    if (useUpbend_) {
        double upbend = upbendConstant_ * std::exp(-upbendExponent_ * Egamma);
        strength += upbend;
    }
    
    return strength * std::pow(Egamma, 3);
}

// ============================================================================
// M1 Single Particle
// ============================================================================

M1SingleParticle::M1SingleParticle(double sigma)
    : sigma_(sigma) {
}

double M1SingleParticle::getStrength(double Ex, double Egamma) const {
    // Constant single-particle strength
    // Original RAINIER lines 442-444
    
    (void)Ex;
    
    return sigma_ * std::pow(Egamma, 3);
}

// ============================================================================
// E2 Standard Lorentzian
// ============================================================================

E2StandardLorentz::E2StandardLorentz(double energy, double width, double sigma)
    : resonance_{energy, width, sigma} {
}

double E2StandardLorentz::getStrength(double Ex, double Egamma) const {
    // E2 Standard Lorentzian
    // Original RAINIER lines 447-451
    // Note: TALYS formula has Eγ in denominator (units correction)
    
    (void)Ex;
    
    double term = resonance_.width * resonance_.width /
                 (Egamma * (std::pow(Egamma * Egamma - resonance_.energy * resonance_.energy, 2) +
                            std::pow(Egamma * resonance_.width, 2)));
    
    double strength = constants::K_X2 * resonance_.sigma * term;
    
    return strength * std::pow(Egamma, 5);  // × Eγ⁵
}

// ============================================================================
// E2 Single Particle
// ============================================================================

E2SingleParticle::E2SingleParticle(double sigma)
    : sigma_(sigma) {
}

double E2SingleParticle::getStrength(double Ex, double Egamma) const {
    // Constant single-particle E2 strength
    // Original RAINIER lines 459-461
    
    (void)Ex;
    
    return sigma_ * std::pow(Egamma, 5);
}

// ============================================================================
// Combined Strength
// ============================================================================

CombinedStrength::CombinedStrength(std::unique_ptr<GammaStrengthFunction> e1,
                                   std::unique_ptr<GammaStrengthFunction> m1,
                                   std::unique_ptr<GammaStrengthFunction> e2)
    : e1_(std::move(e1)), m1_(std::move(m1)), e2_(std::move(e2)) {
}

double CombinedStrength::getTotalStrength(double Ex, double Egamma, int transType,
                                         double& mixingRatio) const {
    // Get total strength based on transition type
    // Original RAINIER lines 464-502
    // transType: 0=None, 1=E1, 2=M1+E2, 3=M1, 4=E2
    
    mixingRatio = 0.0;
    
    switch (transType) {
        case 0:  // No transition
            return 0.0;
            
        case 1:  // Pure E1
            if (e1_) {
                return e1_->getStrength(Ex, Egamma);
            }
            return 0.0;
            
        case 2: {  // Mixed M1+E2
            double strM1 = 0.0;
            double strE2 = 0.0;
            
            if (m1_) strM1 = m1_->getStrength(Ex, Egamma);
            if (e2_) strE2 = e2_->getStrength(Ex, Egamma);
            
            // Mixing ratio δ² = Γ_E2 / Γ_M1 = str_E2 / str_M1
            // (factor of Eγ cancels)
            if (strM1 > 0.0) {
                mixingRatio = strE2 / strM1;
            }
            
            return strM1 + strE2;
        }
            
        case 3:  // Pure M1
            if (m1_) {
                return m1_->getStrength(Ex, Egamma);
            }
            return 0.0;
            
        case 4:  // Pure E2
            if (e2_) {
                return e2_->getStrength(Ex, Egamma);
            }
            return 0.0;
            
        default:
            return 0.0;
    }
}

} // namespace rainier
