// LevelFileReader.h - Parse nuclear level data files
#ifndef RAINIER_LEVEL_FILE_READER_H
#define RAINIER_LEVEL_FILE_READER_H

#include "core/DiscreteLevel.h"
#include <string>
#include <vector>
#include <memory>

namespace rainier {

/**
 * @brief Reads nuclear level data from files
 * 
 * File format (compatible with original RAINIER):
 * 
 * Header line:
 *   Element A Z NumLevels
 * 
 * For each level:
 *   LevelNum Energy Spin Parity HalfLife NumGammas
 *   [For each gamma:]
 *     FinalLevelNum GammaEnergy Pg BranchingRatio ICC
 * 
 * Example:
 *   Nd 144 60 20
 *   0   0.000   0.0  1  1e20  0
 *   1   0.697   2.0  1  1e15  1
 *       0   0.697   1.0   1.000   0.05
 */
class LevelFileReader {
public:
    /**
     * @brief Structure to hold gamma transition data
     */
    struct GammaData {
        int finalLevelIndex;      // Index of final level (0 = ground state)
        double gammaEnergy;       // Gamma-ray energy (MeV)
        double pg;                // Total transition probability
        double branchingRatio;    // Branching ratio (fraction)
        double icc;               // Internal conversion coefficient
    };

    /**
     * @brief Structure to hold level data
     */
    struct LevelData {
        int levelIndex;           // Level number (0 = ground state)
        double energy;            // Excitation energy (MeV)
        double spin;              // Angular momentum
        int parity;               // 0 = negative, 1 = positive
        double halfLife;          // Half-life (integer in file, converted to fs)
        std::vector<GammaData> gammas;
    };

    /**
     * @brief Read level file for specific nucleus
     * @param filename Path to level file
     * @param Z Atomic number to search for
     * @param A Mass number to search for
     * @param maxLevels Maximum number of levels to read
     * @return Vector of level data
     * @throws std::runtime_error if file not found or nucleus not found
     */
    static std::vector<LevelData> readLevelFile(
        const std::string& filename,
        int Z, int A,
        int maxLevels);

    /**
     * @brief Convert LevelData to DiscreteLevel objects with transitions
     * @param levelDataVec Vector of level data
     * @param defaultHalfLife Default half-life for unknown values (fs)
     * @return Vector of DiscreteLevel shared pointers with transitions set up
     */
    static std::vector<std::shared_ptr<DiscreteLevel>> 
        createDiscreteLevels(
            const std::vector<LevelData>& levelDataVec,
            double defaultHalfLife = 1e9);

private:
    /**
     * @brief Check if half-life is a placeholder (unknown value)
     */
    static bool isUnknownHalfLife(double halfLife);

    /**
     * @brief Convert half-life from file format to femtoseconds
     */
    static double convertHalfLife(double halfLifeFromFile, int numGammas);
};

} // namespace rainier

#endif // RAINIER_LEVEL_FILE_READER_H
