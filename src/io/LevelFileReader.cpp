// LevelFileReader.cpp - Implementation of level file parser
#include "io/LevelFileReader.h"
#include "utils/PhysicsConstants.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>

namespace rainier {

std::vector<LevelFileReader::LevelData> LevelFileReader::readLevelFile(
    const std::string& filename,
    int targetZ, int targetA,
    int maxLevels) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open level file: " + filename);
    }

    std::cout << "Reading level file: " << filename << std::endl;

    // Search for the nucleus in the file
    std::string line;
    bool foundNucleus = false;
    std::string element;
    int A = 0, Z = 0, numLevelsInFile = 0;

    const int maxSearchLines = 1000000;  // Prevent infinite loop
    int searchCount = 0;

    while (std::getline(file, line) && searchCount < maxSearchLines) {
        searchCount++;
        std::istringstream iss(line);
        
        // Try to parse header line: Element A Z NumLevels
        if (iss >> element >> A >> Z >> numLevelsInFile) {
            if (A == targetA && Z == targetZ) {
                foundNucleus = true;
                std::cout << "Found nucleus: " << element << "-" << A 
                         << " (Z=" << Z << ") with " << numLevelsInFile 
                         << " levels in file" << std::endl;
                break;
            }
        }
    }

    if (!foundNucleus) {
        throw std::runtime_error(
            "Nucleus Z=" + std::to_string(targetZ) + 
            " A=" + std::to_string(targetA) + 
            " not found in file: " + filename);
    }

    if (numLevelsInFile < 2) {
        throw std::runtime_error("No levels found in file");
    }

    // Limit to requested number or available levels
    int levelsToRead = std::min(maxLevels, numLevelsInFile);
    std::vector<LevelData> levels;
    levels.reserve(levelsToRead);

    // Read each level
    for (int lvl = 0; lvl < levelsToRead; ++lvl) {
        if (!std::getline(file, line)) {
            throw std::runtime_error(
                "Unexpected end of file while reading level " + 
                std::to_string(lvl));
        }

        std::istringstream levelStream(line);
        LevelData levelData;

        int levelNum, parity;
        double energy, spin, halfLife;
        int numGammas;

        levelStream >> levelNum >> energy >> spin >> parity >> halfLife >> numGammas;

        if (levelStream.fail()) {
            throw std::runtime_error(
                "Failed to parse level line: " + line);
        }

        // Validate and convert data
        levelNum--;  // Convert to 0-based indexing (file uses 1-based)

        if (levelNum != lvl) {
            std::cerr << "Warning: Level numbering mismatch. Expected " 
                     << lvl << " but got " << levelNum << std::endl;
        }

        // Convert parity: file uses -1 for unknown, we use 0/1
        if (parity == -1) {
            parity = 0;  // Assume negative if unknown
            std::cerr << "Warning: Unknown parity for level " << lvl 
                     << ", assuming negative" << std::endl;
        }

        // Handle special case: half-life might encode number of gammas
        // if actual half-life is unknown (original RAINIER quirk)
        if (isUnknownHalfLife(halfLife)) {
            // Sometimes the half-life field contains the number of gammas
            // when the actual half-life is not measured
            if (static_cast<int>(halfLife) == static_cast<int>(halfLife) && 
                halfLife < 100 && lvl > 0) {
                numGammas = static_cast<int>(halfLife);
                halfLife = constants::DEFAULT_MAX_HALFLIFE;
                std::cerr << "Warning: Level " << lvl 
                         << " has unknown half-life, using default" << std::endl;
            }
        } else {
            // Convert to femtoseconds
            halfLife = halfLife * 1e15;  // Convert to fs
        }

        // Special case: ground state should be stable
        if (lvl == 0) {
            numGammas = 0;
            halfLife = constants::DEFAULT_MAX_HALFLIFE;
        }

        levelData.levelIndex = lvl;
        levelData.energy = energy;
        levelData.spin = spin;
        levelData.parity = parity;
        levelData.halfLife = halfLife;

        // Read gamma transitions
        for (int gam = 0; gam < numGammas; ++gam) {
            if (!std::getline(file, line)) {
                throw std::runtime_error(
                    "Unexpected end of file while reading gamma " + 
                    std::to_string(gam) + " of level " + std::to_string(lvl));
            }

            std::istringstream gammaStream(line);
            GammaData gamma;

            int finalLevel;
            double gammaE, pg, br, icc;

            gammaStream >> finalLevel >> gammaE >> pg >> br >> icc;

            if (gammaStream.fail()) {
                throw std::runtime_error(
                    "Failed to parse gamma line: " + line);
            }

            finalLevel--;  // Convert to 0-based indexing

            gamma.finalLevelIndex = finalLevel;
            gamma.gammaEnergy = gammaE;
            gamma.pg = pg;
            gamma.branchingRatio = br;
            gamma.icc = icc;

            levelData.gammas.push_back(gamma);
        }

        levels.push_back(levelData);
    }

    file.close();

    std::cout << "Successfully read " << levels.size() << " discrete levels" 
              << std::endl;

    return levels;
}

std::vector<std::shared_ptr<DiscreteLevel>> 
LevelFileReader::createDiscreteLevels(
    const std::vector<LevelData>& levelDataVec,
    double defaultHalfLife) {
    
    std::vector<std::shared_ptr<DiscreteLevel>> levels;
    levels.reserve(levelDataVec.size());

    // First pass: create all levels
    for (const auto& data : levelDataVec) {
        auto level = std::make_shared<DiscreteLevel>(
            data.energy,
            data.spin,
            data.parity,
            data.halfLife > 0 ? data.halfLife : defaultHalfLife
        );
        level->setLevelIndex(data.levelIndex);
        levels.push_back(level);
    }

    // Second pass: create transitions
    // We need all levels created first so we can link them
    for (size_t i = 0; i < levelDataVec.size(); ++i) {
        const auto& data = levelDataVec[i];
        auto& initialLevel = levels[i];

        for (const auto& gammaData : data.gammas) {
            // Validate final level index
            if (gammaData.finalLevelIndex < 0 || 
                gammaData.finalLevelIndex >= static_cast<int>(levels.size())) {
                std::cerr << "Warning: Invalid final level index " 
                         << gammaData.finalLevelIndex 
                         << " for transition from level " << i << std::endl;
                continue;
            }

            auto finalLevel = levels[gammaData.finalLevelIndex];

            // Create transition
            auto transition = std::make_shared<Transition>(
                initialLevel,
                finalLevel,
                gammaData.branchingRatio,
                gammaData.icc
            );

            // Determine transition type
            auto transType = Transition::determineType(
                initialLevel->getSpin(),
                initialLevel->getParity(),
                finalLevel->getSpin(),
                finalLevel->getParity(),
                false  // Will be set correctly by Nucleus based on A
            );
            transition->setType(transType);

            initialLevel->addTransition(transition);
        }

        // Normalize branching ratios for this level
        if (!data.gammas.empty()) {
            try {
                initialLevel->normalizeBranchingRatios();
            } catch (const std::exception& e) {
                std::cerr << "Warning: Could not normalize branching ratios for level " 
                         << i << ": " << e.what() << std::endl;
            }
        }
    }

    return levels;
}

bool LevelFileReader::isUnknownHalfLife(double halfLife) {
    // Check if half-life is a placeholder value
    // Original RAINIER uses integers < 100 to indicate unknown half-lives
    if (halfLife == static_cast<int>(halfLife) && halfLife < 100) {
        return true;
    }
    return false;
}

double LevelFileReader::convertHalfLife(double halfLifeFromFile, int numGammas) {
    if (isUnknownHalfLife(halfLifeFromFile)) {
        return constants::DEFAULT_MAX_HALFLIFE;
    }
    return halfLifeFromFile * 1e15;  // Convert to femtoseconds
}

} // namespace rainier
