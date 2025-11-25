// NLDSystematics.cpp - Implementation of von Egidy & Bucurescu systematics
#include "utils/NLDSystematics.h"
#include <cmath>
#include <iostream>

namespace rainier {

// ============================================================================
// Main CTM Parameter Calculation (from robin.f lines ~260-280)
// ============================================================================

void NLDSystematics::calculateCTMParameters(int A, int Z, double Pd, double S,
                                           double& T, double& E0) {
    // Calculate Pa' (modified pairing energy)
    // From robin.f: Pa_prime = Pd(i)/1000.
    //               IF(ieo.EQ.0.OR.ieo.EQ.1)Pa_prime = -Pd(i)/1000.
    double Pa_prime = calculatePaPrime(Pd, A, Z);
    
    // Calculate S' (modified shell correction)
    // From robin.f: S_prime = REAL(S) + 0.5 * Pa_prime
    double S_prime = S + 0.5 * Pa_prime;
    
    // Calculate temperature T
    // From robin.f: tt = FLOAT(A0)**(-0.66666)/(0.0597 + 0.00198 * S_prime)
    // Note: A^(-2/3) = A^(-0.66666...)
    double A_pow = std::pow(static_cast<double>(A), -2.0/3.0);
    T = A_pow / (0.0597 + 0.00198 * S_prime);
    
    // Calculate E0 (energy shift)
    // From robin.f: E0 = -1.004 + 0.5 * Pa_prime
    E0 = -1.004 + 0.5 * Pa_prime;
}

// ============================================================================
// Spin Cutoff for CTM (from robin.f lines ~258)
// ============================================================================

double NLDSystematics::calculateCTMSpinCutoff(int A) {
    // From robin.f: sig = (0.98*(FLOAT(A0)**(0.29)))
    // This is σ (not σ²)
    return 0.98 * std::pow(static_cast<double>(A), 0.29);
}

// ============================================================================
// Combined Calculation
// ============================================================================

void NLDSystematics::calculateAllCTMParameters(int A, int Z,
                                              double Pd, double S,
                                              double& T, double& E0, double& sigma) {
    // Estimate parameters if not provided
    if (std::abs(Pd) < 1e-6) {
        Pd = estimateDeuteronPairing(A, Z);
        std::cout << "  Estimated Pd = " << Pd << " MeV" << std::endl;
    }
    
    if (std::abs(S) < 1e-6) {
        S = estimateShellCorrection(A, Z);
        std::cout << "  Estimated S = " << S << " MeV" << std::endl;
    }
    
    // Calculate T and E0
    calculateCTMParameters(A, Z, Pd, S, T, E0);
    
    // Calculate spin cutoff
    sigma = calculateCTMSpinCutoff(A);
}

// ============================================================================
// Helper Functions
// ============================================================================

int NLDSystematics::getPairingType(int A, int Z) {
    // From robin.f:
    // N0=A0-Z0
    // IF((N0/2)*2.EQ.N0.AND.(Z0/2)*2.EQ.Z0)ieo=0 ! even-even nucleus
    // IF((N0/2)*2.NE.N0.AND.(Z0/2)*2.EQ.Z0)ieo=1 ! even-Z odd-N
    // IF((N0/2)*2.EQ.N0.AND.(Z0/2)*2.NE.Z0)ieo=2 ! odd-Z even-N
    // IF((N0/2)*2.NE.N0.AND.(Z0/2)*2.NE.Z0)ieo=3 ! odd-odd nucleus
    
    int N = A - Z;
    bool even_N = (N % 2 == 0);
    bool even_Z = (Z % 2 == 0);
    
    if (even_N && even_Z) return 0;  // even-even
    if (!even_N && even_Z) return 1; // odd-N, even-Z
    if (even_N && !even_Z) return 2; // even-N, odd-Z
    return 3;                         // odd-odd
}

double NLDSystematics::calculatePaPrime(double Pd, int A, int Z) {
    // From robin.f lines ~244-245:
    // Pa_prime = Pd(i)/1000.
    // IF(ieo.EQ.0.OR.ieo.EQ.1)Pa_prime = -Pd(i)/1000.
    
    int ieo = getPairingType(A, Z);
    
    if (ieo == 0 || ieo == 1) {
        // even-even or odd-A with even-Z
        return -Pd;
    } else {
        // odd-A with odd-Z, or odd-odd
        return +Pd;
    }
}

// ============================================================================
// Parameter Estimation (simplified versions)
// ============================================================================

double NLDSystematics::estimateDeuteronPairing(int A, int Z) {
    // Simplified estimation based on nucleus type
    // From typical values observed in robin output and mass tables
    
    int ieo = getPairingType(A, Z);
    int N = A - Z;
    
    // Empirical pairing gap formula (approximate)
    // Δ ~ 12/√A for neutrons and protons separately
    // Pd combines both
    
    double delta_n = 12.0 / std::sqrt(static_cast<double>(A));
    double delta_p = 12.0 / std::sqrt(static_cast<double>(A));
    
    double Pd;
    switch (ieo) {
        case 0: // even-even: both gaps contribute
            Pd = delta_n + delta_p;
            break;
        case 1: // odd-N, even-Z: only proton gap
            Pd = delta_p;
            break;
        case 2: // even-N, odd-Z: only neutron gap
            Pd = delta_n;
            break;
        case 3: // odd-odd: no pairing gaps
            Pd = 0.0;
            break;
        default:
            Pd = 1.0; // fallback
    }
    
    // Typical values from robin:
    // Nd-144 (even-even): Pd ~ 2.698 MeV
    // For A~150: Pd ~ 1.5-3 MeV for even-even
    //            Pd ~ 0.8-1.5 MeV for odd-A
    //            Pd ~ 0-0.5 MeV for odd-odd
    
    return Pd;
}

double NLDSystematics::estimateShellCorrection(int A, int Z) {
    // Simplified shell correction estimation
    // Based on distance from magic numbers: 2, 8, 20, 28, 50, 82, 126
    
    const int magicNumbers[] = {2, 8, 20, 28, 50, 82, 126};
    const int nMagic = 7;
    
    int N = A - Z;
    
    // Find distance to nearest magic number for Z and N
    auto distanceToMagic = [&](int value) -> int {
        int minDist = 999;
        for (int i = 0; i < nMagic; ++i) {
            int dist = std::abs(value - magicNumbers[i]);
            if (dist < minDist) {
                minDist = dist;
            }
        }
        return minDist;
    };
    
    int distZ = distanceToMagic(Z);
    int distN = distanceToMagic(N);
    
    // Shell effects are strongest near magic numbers
    // At magic numbers: S ~ -8 to -5 MeV (more bound)
    // Mid-shell: S ~ 0 to +2 MeV (less bound)
    
    // Simplified formula:
    // S decreases (becomes more negative) near magic numbers
    double S_Z = 3.0 * (1.0 - std::exp(-distZ / 10.0));
    double S_N = 3.0 * (1.0 - std::exp(-distN / 10.0));
    
    // Combine contributions
    double S = -(S_Z + S_N) / 2.0;
    
    // Add some mass-region dependence
    // Heavy nuclei (A > 200) typically have S ~ -2 to 0
    // Medium nuclei (100 < A < 200) have S ~ -4 to +2
    // Light nuclei (A < 100) more variable
    
    if (A > 200) {
        S *= 0.5; // Reduce shell effects in heavy region
    }
    
    // Typical values from robin for Nd-144 (Z=60, A=144):
    // S ~ -1.5 to 0 MeV (mid-shell region)
    
    // For now, use a simple estimate
    // More accurate would require full mass table
    return S;
}

} // namespace rainier
