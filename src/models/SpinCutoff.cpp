// SpinCutoff.cpp - Implementation of spin cutoff models
// Physics from original RAINIER.C lines 232-295
#include "models/SpinCutoff.h"
#include "models/LevelDensity.h"
#include "utils/PhysicsConstants.h"
#include <cmath>
#include <iostream>

namespace rainier {

// ============================================================================
// Base Class
// ============================================================================

double SpinCutoffModel::getSpinDistribution(double Ex, double spin) const {
    double sigma2 = getSigmaSquared(Ex);
    
    // P(J) = (2J+1)/σ² × exp(-(J+1/2)²/(2σ²))
    double distribution = (2.0 * spin + 1.0) * 
                         std::exp(-std::pow(spin + 0.5, 2) / (2.0 * sigma2)) / 
                         sigma2;
    
    return distribution;
}

// ============================================================================
// Von Egidy 2005 Model
// ============================================================================


VonEgidy05::VonEgidy05(const LevelDensityModel* densityModel, int A, double osloShift)
    : densityModel_(densityModel), A_(A), osloShift_(osloShift) {
}


double VonEgidy05::getSigmaSquared(double Ex) const {
    // Von Egidy & Bucurescu, PRC 72, 044311 (2005)
    // Original RAINIER lines 233-236
    // σ² = 0.0146 A^(5/3) [1 + √(1 + 4aU)] / (2a)
    
    double shiftedEx = Ex - osloShift_;  // Apply Oslo shift
    if (shiftedEx < 0.0) shiftedEx = 1e-8;  // Safety check
    
    double U = densityModel_->getEffectiveEnergy(Ex);
    
    // Get level density parameter - need to cast to access method
    const auto* bsfg = dynamic_cast<const BackShiftedFermiGas*>(densityModel_);
    double a = 15.0; // Default
    if (bsfg) {
        a = bsfg->getLevelDensityParameter(Ex);
    }
    
    double A_factor = 0.0146 * std::pow(A_, 5.0/3.0);
    double sigma2 = A_factor * (1.0 + std::sqrt(1.0 + 4.0 * a * U)) / (2.0 * a);
    
    return sigma2;
}

// ============================================================================
// TALYS Spin Cutoff Model
// ============================================================================

TALYSSpinCutoff::TALYSSpinCutoff(const LevelDensityModel* densityModel, 
                                 int A, double Sn,
                                 double spinCutoffD, double Ed, 
                                 double aAsymptotic)
    : densityModel_(densityModel), A_(A), spinCutoffD_(spinCutoffD),
      Ed_(Ed), Sn_(Sn), aAsymptotic_(aAsymptotic) {
}

TALYSSpinCutoff::TALYSSpinCutoff(const LevelDensityModel* densityModel,
                                 int A, double Sn,
                                 double spinCutoffD, double Ed,
                                 double aAsymptotic, double osloShift)
    : densityModel_(densityModel), A_(A), spinCutoffD_(spinCutoffD),
      Ed_(Ed), Sn_(Sn), aAsymptotic_(aAsymptotic), osloShift_(osloShift) {
}

double TALYSSpinCutoff::getSigmaSquared(double Ex) const {
    // Apply Oslo shift to Ex first (like original RAINIER)
    double shiftedEx = Ex - osloShift_;
    if (shiftedEx < 0.0) shiftedEx = 1e-8;
    
    double spinCutoffD2 = spinCutoffD_ * spinCutoffD_;
    
    // Discrete region - use original Ex for comparison, not shifted
    if (Ex <= Ed_) {
        return spinCutoffD2;
    }
    
    // Calculate dEff and dLDa using shifted Ex (matches original RAINIER)
    double U = densityModel_->getEffectiveEnergy(shiftedEx);
    const auto* bsfg = dynamic_cast<const BackShiftedFermiGas*>(densityModel_);
    double a = 15.0;
    if (bsfg) {
        a = bsfg->getLevelDensityParameter(shiftedEx);
    }
    
    // Calculate sigma at Sn using shifted Sn
    double shiftedSn = Sn_ - osloShift_;
    if (shiftedSn < 0.0) shiftedSn = 1e-8;
    
    double U_Sn = densityModel_->getEffectiveEnergy(shiftedSn);
    double a_Sn = 15.0;
    if (bsfg) {
        a_Sn = bsfg->getLevelDensityParameter(shiftedSn);
    }
    
    double spinCutoffSn2 = 0.01389 * std::pow(A_, 5.0/3.0) / aAsymptotic_ *
                          std::sqrt(a_Sn * U_Sn);
    
    // Interpolation region (use original Ex for comparison)
    if (Ex < Sn_) {
        double fraction = (Ex - Ed_) / (Sn_ - Ed_);
        return spinCutoffD2 + fraction * (spinCutoffSn2 - spinCutoffD2);
    }
    
    // Asymptotic region (Ex >= Sn) - use shifted values
    double sigma2 = 0.01389 * std::pow(A_, 5.0/3.0) / aAsymptotic_ *
                   std::sqrt(a * U);
    
    return sigma2;
}

// ============================================================================
// Single Particle Model
// ============================================================================

SingleParticle::SingleParticle(const LevelDensityModel* densityModel, int A, double osloShift)
    : densityModel_(densityModel), A_(A), osloShift_(osloShift) {
}

double SingleParticle::getSigmaSquared(double Ex) const {
    // Gholami et al., PRC 75, 044308 (2007)
    // Original RAINIER lines 239-241
    // σ² = 0.1461 √(aU) A^(2/3)
    
    double shiftedEx = Ex - osloShift_;
    if (shiftedEx < 0.0) shiftedEx = 1e-8;
    
    double U = densityModel_->getEffectiveEnergy(Ex);
    
    const auto* bsfg = dynamic_cast<const BackShiftedFermiGas*>(densityModel_);
    double a = 15.0;
    if (bsfg) {
        a = bsfg->getLevelDensityParameter(Ex);
    }
    
    double sigma2 = 0.1461 * std::sqrt(a * U) * std::pow(A_, 2.0/3.0);
    
    return sigma2;
}

// ============================================================================
// Rigid Sphere Model
// ============================================================================

RigidSphere::RigidSphere(const LevelDensityModel* densityModel, int A, double osloShift)
    : densityModel_(densityModel), A_(A), osloShift_(osloShift) {
}

double RigidSphere::getSigmaSquared(double Ex) const {
    // Grimes et al., PRC 10 (1974) 2373
    // Original RAINIER lines 243-245
    // σ² = 0.0145 √(U/a) A^(5/3)
    
    double shiftedEx = Ex - osloShift_;
    if (shiftedEx < 0.0) shiftedEx = 1e-8;
    
    double U = densityModel_->getEffectiveEnergy(Ex);
    
    const auto* bsfg = dynamic_cast<const BackShiftedFermiGas*>(densityModel_);
    double a = 15.0;
    if (bsfg) {
        a = bsfg->getLevelDensityParameter(Ex);
    }
    
    double sigma2 = 0.0145 * std::sqrt(U / a) * std::pow(A_, 5.0/3.0);
    
    return sigma2;
}

// ============================================================================
// Von Egidy 2009 Model
// ============================================================================

VonEgidy09::VonEgidy09(int A, double pairingEnergy)
    : A_(A), pairingEnergy_(pairingEnergy) {
}

double VonEgidy09::getSigmaSquared(double Ex) const {
    // Von Egidy & Bucurescu, PRC 80, 054310 (2009)
    // Original RAINIER lines 247-251
    // σ² = 0.391 A^0.675 (E - Δ/2)^0.312
    
    double ExPaired = Ex - 0.5 * pairingEnergy_;
    if (ExPaired < 0.000001) {
        ExPaired = 0.000001;  // Avoid zero/negative
    }
    
    double sigma2 = 0.391 * std::pow(A_, 0.675) * std::pow(ExPaired, 0.312);
    
    return sigma2;
}

} // namespace rainier
