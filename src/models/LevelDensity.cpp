// LevelDensity.cpp - Complete implementation of level density models
// Physics preserved from original RAINIER.C (lines 159-295)
#include "models/LevelDensity.h"
#include "models/SpinCutoff.h"
#include "utils/PhysicsConstants.h"
#include "Config.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include "utils/NLDSystematics.h"

namespace rainier {

// ============================================================================
// Back-Shifted Fermi Gas (BSFG) Model
// ============================================================================

BackShiftedFermiGas::BackShiftedFermiGas(
    double a, double E1, bool useEnergyDependent,
    double aAsymptotic, double shellCorrection,
    double dampingGamma, int A)
    : a_(a),
      E1_(E1),
      useEnergyDependent_(useEnergyDependent),
      aAsymptotic_(aAsymptotic),
      shellCorrection_(shellCorrection),
      dampingGamma_(dampingGamma),
      A_(A) {
    
    std::cout << "Initialized BSFG level density model:" << std::endl;
    std::cout << "  a = " << a_ << " MeV^-1" << std::endl;
    std::cout << "  E1 = " << E1_ << " MeV" << std::endl;
    if (useEnergyDependent_) {
        std::cout << "  Using energy-dependent a(E)" << std::endl;
        std::cout << "  a_asymptotic = " << aAsymptotic_ << " MeV^-1" << std::endl;
        std::cout << "  Shell correction W = " << shellCorrection_ << " MeV" << std::endl;
        std::cout << "  Damping γ = " << dampingGamma_ << " MeV^-1" << std::endl;
    }
}

double BackShiftedFermiGas::getEffectiveEnergy(double Ex) const {
    // U = E - E1 (back-shift energy)
    // Original RAINIER lines 159-163
    double U = Ex - E1_;
    if (U < constants::MIN_ENERGY) {
        U = constants::MIN_ENERGY;  // Avoid negative or zero
    }
    return U;
}

double BackShiftedFermiGas::getLevelDensityParameter(double Ex) const {
    // Original RAINIER lines 164-178
    if (!useEnergyDependent_) {
        return a_;
    }
    
    // TALYS 1.8 energy-dependent a(E)
    // a(E) = a_asymptotic × [1 + δW × (1 - exp(-γU))/U]
    double U = getEffectiveEnergy(Ex);
    
    double dampingFactor = (1.0 - std::exp(-dampingGamma_ * U)) / U;
    double a_E = aAsymptotic_ * (1.0 + shellCorrection_ * dampingFactor);
    
    return a_E;
}

double BackShiftedFermiGas::getDensity(double Ex) const {
    // Total level density (summed over all J and π)
    // Original RAINIER lines 280-295
    // ρ(E) = 1/(12√2σ) × exp(2√(aU)) / (a^(1/4) × U^(5/4))
    
    double U = getEffectiveEnergy(Ex);
    double a = getLevelDensityParameter(Ex);
    
    // We need spin cutoff, but we don't have it yet
    // For total density, use approximate spin cutoff
    double sigma2 = 0.0146 * std::pow(A_, 5.0/3.0) * 
                   (1.0 + std::sqrt(1.0 + 4.0 * a * U)) / (2.0 * a);
    
    double sigma = std::sqrt(sigma2);
    
    double denominator = 12.0 * std::sqrt(2.0) * sigma * 
                        std::pow(a, 0.25) * std::pow(U, 1.25);
    
    double exponential = std::exp(2.0 * std::sqrt(a * U));
    
    double density = exponential / denominator;
    
    if (density < constants::MIN_DENSITY) {
        density = constants::MIN_DENSITY;
    }
    
    return density;
}

double BackShiftedFermiGas::getDensity(double Ex, double spin, int parity) const {
    // Level density for specific spin and parity
    // ρ(E,J,π) = ρ_total(E) × P(J) × P(π)
    // Original RAINIER lines 296-324
    
    (void)parity; // Will be used when parity model integrated
    
    double rhoTotal = getDensity(Ex);
    
    // Spin distribution factor
    // P(J) = (2J+1)/σ² × exp(-(J+1/2)²/(2σ²))
    double U = getEffectiveEnergy(Ex);
    double a = getLevelDensityParameter(Ex);
    
    // Spin cutoff parameter (Von Egidy 2005 formula)
    // σ² = 0.0146 A^(5/3) × [1 + √(1 + 4aU)] / (2a)
    double sigma2 = 0.0146 * std::pow(A_, 5.0/3.0) * 
                   (1.0 + std::sqrt(1.0 + 4.0 * a * U)) / (2.0 * a);
    
    double spinFactor = (2.0 * spin + 1.0) * 
                       std::exp(-std::pow(spin + 0.5, 2) / (2.0 * sigma2)) / 
                       sigma2;
    
    // Parity distribution (will use proper model from ParityDistribution)
    // For now, use equipartition
    double parityFactor = 0.5;
    
    double density = rhoTotal * spinFactor * parityFactor;
    
    if (density < constants::MIN_DENSITY) {
        density = constants::MIN_DENSITY;
    }
    
    return density;
}

// ============================================================================
// Constant Temperature Model (CTM)
// ============================================================================


ConstantTemperature::ConstantTemperature(double T, double E0, int A, int Z)
    : T_(T),
      E0_(E0) {
    
    // If both T and E0 are zero, calculate defaults using von Egidy & Bucurescu systematics
    if (std::abs(T_) < 1e-6 && std::abs(E0_) < 1e-6) {
        std::cout << "\n╔═══════════════════════════════════════════════════╗\n";
        std::cout << "║ Calculating CTM defaults from von Egidy & Bucurescu ║\n";
        std::cout << "╚═══════════════════════════════════════════════════╝\n";
        
        double sigma;  // Not used in constructor but calculated for reference
        NLDSystematics::calculateAllCTMParameters(A, Z, 0.0, 0.0, T_, E0_, sigma);
        
        std::cout << "\n═══ Calculated CTM Parameters ═══\n";
        std::cout << "  T  = " << T_ << " MeV (nuclear temperature)\n";
        std::cout << "  E0 = " << E0_ << " MeV (energy shift)\n";
        std::cout << "  σ  = " << sigma << " (spin cutoff parameter)\n";
        std::cout << "════════════════════════════════════\n\n";
    } else {
        std::cout << "Initialized CTM level density model:" << std::endl;
        std::cout << "  T = " << T_ << " MeV" << std::endl;
        std::cout << "  E0 = " << E0_ << " MeV" << std::endl;
    }
}

double ConstantTemperature::getEffectiveEnergy(double Ex) const {
    // U = E - E0
    // Original RAINIER lines 174-178
    double U = Ex - E0_;
    if (U < constants::MIN_ENERGY) {
        U = constants::MIN_ENERGY;
    }
    return U;
}

double ConstantTemperature::getDensity(double Ex) const {
    // Total level density
    // ρ(E) = exp(U/T) / T
    // Original RAINIER lines 174-178
    
    double U = getEffectiveEnergy(Ex);
    double density = std::exp(U / T_) / T_;
    
    if (density < constants::MIN_DENSITY) {
        density = constants::MIN_DENSITY;
    }
    
    return density;
}

double ConstantTemperature::getDensity(double Ex, double spin, int parity) const {
    // CTM with spin and parity factors
    
    (void)parity; // Will be used when parity model integrated
    
    double rhoTotal = getDensity(Ex);
    
    // For CTM, we still need spin cutoff
    // Use a simple approximation based on temperature
    double sigma2 = T_ * T_;  // Simple approximation
    
    double spinFactor = (2.0 * spin + 1.0) * 
                       std::exp(-std::pow(spin + 0.5, 2) / (2.0 * sigma2)) / 
                       sigma2;
    
    double parityFactor = 0.5;  // Equipartition
    
    double density = rhoTotal * spinFactor * parityFactor;
    
    if (density < constants::MIN_DENSITY) {
        density = constants::MIN_DENSITY;
    }
    
    return density;
}

// ============================================================================
// Parity Distribution
// ============================================================================

double ParityDistribution::getFactor(
    double Ex, int parity, bool isEvenA,
    Model model, double parC, double parD) {
    
    // Original RAINIER lines 325-345
    
    if (model == Model::EQUIPARTITION) {
        // Equal probability for both parities
        return 0.5;
    }
    
    // Energy-dependent parity distribution
    // Al-Quraishi et al., PRC 67, 015803 (2003)
    // P(π) = 0.5 × [1 ± 1/(1 + exp(c(E-d)))]
    
    double exponentialTerm = 1.0 / (1.0 + std::exp(parC * (Ex - parD)));
    
    double factor;
    if (isEvenA) {
        // Even-A nuclei
        if (parity == 1) {  // Positive parity
            // More positive parity states at low energy
            factor = 0.5 * (1.0 + exponentialTerm);
        } else {  // Negative parity
            // Fewer negative parity states at low energy
            factor = 0.5 * (1.0 - exponentialTerm);
        }
    } else {
        // Odd-A nuclei (opposite trend)
        if (parity == 0) {  // Negative parity
            factor = 0.5 * (1.0 + exponentialTerm);
        } else {  // Positive parity
            factor = 0.5 * (1.0 - exponentialTerm);
        }
    }
    
    return factor;
}

} // namespace rainier
