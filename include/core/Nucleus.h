// Nucleus.h - Nuclear level scheme container
#ifndef RAINIER_NUCLEUS_H
#define RAINIER_NUCLEUS_H

#include "DiscreteLevel.h"
#include "ContinuumLevel.h"
#include "models/LevelDensity.h"
#include "models/SpinCutoff.h"

#include <vector>
#include <memory>
#include <map>

namespace rainier {

class Config;

/**
 * @brief Represents a nucleus with its complete level scheme
 * 
 * Contains both discrete (experimentally known) levels and
 * statistically generated continuum levels. Manages the complete
 * nuclear structure for cascade simulations.
 */
class Nucleus {
public:
    /**
     * @brief Construct a nucleus
     * @param Z Atomic number (protons)
     * @param A Mass number (protons + neutrons)
     */
    Nucleus(int Z, int A, const Config& config);

    // Basic properties
    int getZ() const { return Z_; }
    int getA() const { return A_; }
    int getN() const { return A_ - Z_; }
    bool isEvenA() const { return (A_ % 2 == 0); }
    bool isEvenZ() const { return (Z_ % 2 == 0); }
    bool isEvenN() const { return ((A_ - Z_) % 2 == 0); }

    /**
     * @brief Load discrete levels from file
     * @param filename Path to level file
     * @param maxLevels Maximum number of discrete levels to load
     */
    void loadDiscreteLevels(const std::string& filename, int maxLevels);

    /**
     * @brief Build continuum level scheme
     * @param config Configuration parameters
     * @param realization Realization number (for random seed)
     */
    void buildContinuumLevels(const Config& config, int realization);

    /**
     * @brief Get critical energy (top of discrete region)
     */
    double getCriticalEnergy() const { return criticalEnergy_; }

    /**
     * @brief Get maximum excitation energy
     */
    double getMaxEnergy() const { return maxEnergy_; }

    // Access discrete levels
    int getNumDiscreteLevels() const { return discreteLevels_.size(); }
    int getNumAllDiscreteLevels() const { return allDiscreteLevels_.size(); }
        
    std::shared_ptr<DiscreteLevel> getDiscreteLevel(int index) const;
    std::shared_ptr<DiscreteLevel> getAllDiscreteLevel(int index) const {
        return allDiscreteLevels_.at(index);
    }
    const std::vector<std::shared_ptr<DiscreteLevel>>& getDiscreteLevels() const {
        return discreteLevels_;
    }

    /**
    * @brief Get level density model
    */
    const LevelDensityModel* getLevelDensityModel() const {
        return levelDensity_.get();
    }
    
    /**
    * @brief Get spin cutoff model
    */
    const SpinCutoffModel* getSpinCutoffModel() const {
        return spinCutoff_.get();
    }
            
    // Access continuum levels
    int getNumContinuumLevels() const { return totalContinuumLevels_; }
    
    /**
     * @brief Get continuum levels in a specific E-J-π bin
     * @param energyBin Energy bin index
     * @param spinBin Spin bin index
     * @param parity Parity (0 or 1)
     */
    const std::vector<std::shared_ptr<ContinuumLevel>>& 
        getContinuumLevels(int energyBin, int spinBin, int parity) const;
    
    /**
     * @brief Get number of continuum levels in a bin
     */
    int getNumLevelsInBin(int energyBin, int spinBin, int parity) const;

    /**
     * @brief Get central energy of a continuum bin
     */
    double getBinCentralEnergy(int bin) const;

    /**
     * @brief Get energy bin for a given excitation energy
     */
    int getEnergyBin(double energy) const;

    /**
     * @brief Get number of energy bins
     */
    int getNumEnergyBins() const { return numEnergyBins_; }

    /**
     * @brief Get energy spacing between bins
     */
    double getEnergySpacing() const { return energySpacing_; }

    /**
     * @brief Get maximum spin bin
     */
    int getMaxSpinBin() const { return maxSpinBin_; }

    /**
     * @brief Print level scheme summary
     */
    void printSummary() const;

    /**
     * @brief Set neutron separation energy
     */
    void setSn(double Sn) { Sn_ = Sn; }
    double getSn() const { return Sn_; }
    
    /**
         * @brief Set level density model for continuum construction
         */
        void setLevelDensityModel(const LevelDensityModel* model) {
            levelDensity_ = model;
        }
        
        /**
         * @brief Set spin cutoff model for continuum construction
         */
        void setSpinCutoffModel(const SpinCutoffModel* model) {
            spinCutoff_ = model;
        }

private:
    // Index for continuum level storage: maps (energyBin, spinBin, parity) to vector of levels
    using ContinuumKey = std::tuple<int, int, int>;
    
    int Z_;                    // Atomic number
    int A_;                    // Mass number
    double Sn_;                // Neutron separation energy (MeV)
    
    // Discrete levels
    std::vector<std::shared_ptr<DiscreteLevel>> discreteLevels_;
    std::vector<std::shared_ptr<DiscreteLevel>> allDiscreteLevels_;
    double criticalEnergy_;    // Top of discrete region (MeV)
    
    // Continuum levels
    std::map<ContinuumKey, std::vector<std::shared_ptr<ContinuumLevel>>> continuumLevels_;
    int totalContinuumLevels_;
    int maxSpinBin_;
    double maxEnergy_;         // Maximum excitation energy (MeV)
    double energySpacing_;     // Spacing between energy bins (MeV)
    int numEnergyBins_;        // Number of continuum energy bins
    
    
    // Level density models (owning pointers)
       std::unique_ptr<LevelDensityModel> levelDensity_;
       std::unique_ptr<SpinCutoffModel> spinCutoff_;
    
    // Helper methods
    void buildPoissonLevels(int realization);
    void buildWignerLevels(int realization);
};

} // namespace rainier

#endif // RAINIER_NUCLEUS_H
