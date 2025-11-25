// NLDSystematics.h - Calculate nuclear level density parameters from systematics
// Based on von Egidy & Bucurescu formulas (robin.f from Oslo method software)
#ifndef RAINIER_NLD_SYSTEMATICS_H
#define RAINIER_NLD_SYSTEMATICS_H

namespace rainier {

/**
 * @brief Utility class for calculating NLD parameters from global systematics
 * 
 * Implements the von Egidy & Bucurescu (2005, 2009) systematics
 * as used in the Oslo method software (robin.f/rhobin.f).
 * 
 * References:
 * - von Egidy & Bucurescu, PRC 72, 044311 (2005)
 * - von Egidy & Bucurescu, PRC 80, 054310 (2009)
 */
class NLDSystematics {
public:
    
    /**
     * @brief Calculate CTM parameters using von Egidy & Bucurescu (2009) systematics
     * 
     * Formulas from robin.f:
     *   T  = A^(-2/3) / (0.0597 + 0.00198 * S')
     *   E0 = -1.004 + 0.5 * Pa'
     * 
     * where:
     *   S' = S + 0.5 * Pa'  (modified shell correction)
     *   Pa' = -Pd  for even-even and odd-even nuclei
     *   Pa' = +Pd  for even-odd and odd-odd nuclei
     *   Pd = deuteron pairing energy
     * 
     * @param A Mass number
     * @param Z Proton number
     * @param Pd Deuteron pairing energy (MeV) - from mass differences
     * @param S Shell correction (MeV) - calculated from mass formula
     * @param T Output: Nuclear temperature (MeV)
     * @param E0 Output: Energy shift (MeV)
     */
    static void calculateCTMParameters(int A, int Z, double Pd, double S,
                                      double& T, double& E0);
    
    /**
     * @brief Calculate deuteron pairing energy from mass differences
     * 
     * Formula from atomic mass evaluation:
     *   Pd = [BE(A+1,Z+1) + BE(A-1,Z-1) - 2*BE(A,Z)] / 2
     * 
     * where BE = binding energy
     * 
     * Typical values:
     *   Even-even: Pd ~ 2-4 MeV
     *   Odd-A:     Pd ~ 1-2 MeV
     *   Odd-odd:   Pd ~ 0-1 MeV
     * 
     * @param A Mass number
     * @param Z Proton number
     * @return Pd in MeV (estimated if mass data unavailable)
     */
    static double estimateDeuteronPairing(int A, int Z);
    
    /**
     * @brief Calculate shell correction from semi-empirical mass formula
     * 
     * From robin.f ShellCorr subroutine:
     *   S = M_exp - M_theory
     * 
     * where M_theory is from liquid drop model
     * 
     * Typical values:
     *   Near magic numbers: S ~ -8 to +8 MeV
     *   Mid-shell:          S ~ -2 to +2 MeV
     * 
     * @param A Mass number
     * @param Z Proton number
     * @return Shell correction S in MeV
     */
    static double estimateShellCorrection(int A, int Z);
    
    /**
     * @brief Calculate CTM spin cutoff parameter
     * 
     * Formula from robin.f (E&B 2009):
     *   σ = 0.98 * A^0.29
     * 
     * This is a constant (energy-independent) formula for CTM.
     * 
     * @param A Mass number
     * @return Spin cutoff parameter σ (not σ²)
     */
    static double calculateCTMSpinCutoff(int A);
    
    /**
     * @brief Calculate all CTM parameters at once with automatic estimation
     * 
     * If Pd and S are not provided (set to 0.0), they will be estimated
     * using the mass number and proton number.
     * 
     * @param A Mass number
     * @param Z Proton number
     * @param Pd Deuteron pairing (MeV), if 0.0 will be estimated
     * @param S Shell correction (MeV), if 0.0 will be estimated
     * @param T Output: Nuclear temperature (MeV)
     * @param E0 Output: Energy shift (MeV)
     * @param sigma Output: Spin cutoff parameter
     */
    static void calculateAllCTMParameters(int A, int Z,
                                         double Pd, double S,
                                         double& T, double& E0, double& sigma);

private:
    /**
     * @brief Determine nucleus pairing type
     * 
     * @param A Mass number
     * @param Z Proton number
     * @return 0=even-even, 1=odd-A (even-Z), 2=odd-A (odd-Z), 3=odd-odd
     */
    static int getPairingType(int A, int Z);
    
    /**
     * @brief Calculate Pa' (modified pairing energy for formulas)
     * 
     * Pa' = -Pd for even-even and odd-even (ieo=0,1)
     * Pa' = +Pd for even-odd and odd-odd (ieo=2,3)
     * 
     * @param Pd Deuteron pairing energy
     * @param A Mass number
     * @param Z Proton number
     * @return Pa' in MeV
     */
    static double calculatePaPrime(double Pd, int A, int Z);
};

} // namespace rainier

#endif // RAINIER_NLD_SYSTEMATICS_H
