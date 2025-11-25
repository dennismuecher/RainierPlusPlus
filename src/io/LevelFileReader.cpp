// LevelFileReader.cpp - Fixed to handle missing half-life values
#include "io/LevelFileReader.h"
#include "utils/PhysicsConstants.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <TRandom2.h>

namespace rainier {

std::vector<LevelFileReader::LevelData> LevelFileReader::readLevelFile(
    const std::string& filename,
    int targetZ, int targetA,
    int maxLevels = 0) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open level file: " + filename);
    }

    std::cout << "Reading level file: " << filename << std::endl;
    std::cout << "Searching for Z=" << targetZ << ", A=" << targetA << std::endl;

    // Search for the nucleus in the file
    std::string line;
    bool foundNucleus = false;
    std::string element;
    int A = 0, Z = 0, numLevelsInFile = 0;

    const int maxSearchLines = 1000000;
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
        throw std::runtime_error("No levels found in file (nLvlTot < 2)");
    }

    int levelsToRead = (maxLevels == 0) ? numLevelsInFile : std::min(maxLevels, numLevelsInFile);
    
    std::cout << "Reading " << levelsToRead << " levels (max requested: "
              << maxLevels << ", available: " << numLevelsInFile << ")" << std::endl;
    
    std::vector<LevelData> levels;
    levels.reserve(levelsToRead);

    // Read each level
    for (int lvl = 0; lvl < levelsToRead; ++lvl) {
        if (!std::getline(file, line)) {
            throw std::runtime_error(
                "Unexpected end of file while reading level " + std::to_string(lvl));
        }

        std::istringstream levelStream(line);
        LevelData levelData;

        // Try to parse all 6 fields first
        int levelNum, parity, numGammas;
        double energy, spin, halfLife;
        
        // Read first 4 fields that must be present
        if (!(levelStream >> levelNum >> energy >> spin >> parity)) {
            throw std::runtime_error("Failed to parse level line " + 
                                   std::to_string(lvl) + ": " + line);
        }
        
        // Now try to read the remaining fields (halfLife and numGammas)
        // If halfLife is missing, we'll only get numGammas
        double field5, field6;
        bool hasField5 = static_cast<bool>(levelStream >> field5);
        bool hasField6 = static_cast<bool>(levelStream >> field6);
        
        if (hasField5 && hasField6) {
            // Check if field5 is a small integer (likely numGammas, not halfLife)
            if (field5 == static_cast<int>(field5) && field5 < 100 && lvl > 0) {
                // Small integer in field5: interpret as numGammas, field6 is something else
                halfLife = constants::DEFAULT_MAX_HALFLIFE;
                numGammas = static_cast<int>(field5);
               // std::cout << "  Level " << lvl << ": half-life missing (field5 is small integer), using default\n";
            } else {
                // Normal case: field5=halfLife, field6=numGammas
                halfLife = field5;
                numGammas = static_cast<int>(field6);
            }
        } else if (hasField5 && !hasField6) {
            // Only one field: must be numGammas, halfLife is missing
            halfLife = constants::DEFAULT_MAX_HALFLIFE;
            numGammas = static_cast<int>(field5);
            //std::cout << "  Level " << lvl << ": half-life missing, using default\n";
        }
        else {
            // No additional fields - ground state or bad format
            if (lvl == 0) {
                // Ground state: stable, no gammas
                halfLife = constants::DEFAULT_MAX_HALFLIFE;
                numGammas = 0;
            } else {
                throw std::runtime_error("Failed to parse level line " + 
                                       std::to_string(lvl) + ": insufficient fields");
            }
        }

        // Convert from 1-based to 0-based indexing
        levelNum--;

        if (levelNum != lvl) {
            std::cerr << "Warning: Level numbering mismatch. " 
                     << "Expected " << lvl << " but got " << levelNum << std::endl;
        }
        
        // Convert parity: file format 0=unknown, +1=positive, -1=negative
        //                 internal format 1=positive, 0=negative
        if (parity == 0) {
            // Unknown parity - randomly assign to positive or negative with equal probability
            static TRandom2 parityRNG(42);  // Fixed seed for reproducibility
            parity = parityRNG.Integer(2);  // Returns 0 or 1 with equal probability
            //std::cout << "  Level " << lvl << ": unknown parity, randomly assigned to "
                     // << (parity == 1 ? "positive" : "negative") << std::endl;
        } else if (parity == -1) {
            // Negative parity in file → 0 in internal representation
            parity = 0;
        } else if (parity == 1) {
            // Positive parity in file → stays 1 in internal representation
            // (no change needed, but included for clarity)
        }

        // Convert half-life to femtoseconds
        if (halfLife < constants::DEFAULT_MAX_HALFLIFE) {
            halfLife = halfLife * 1e15;  // Convert to fs
        }

        // Ground state special case
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
        double totalBR = 0.0;
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
                std::cerr << "Failed to parse gamma line: [" << line << "]" << std::endl;
                throw std::runtime_error("Failed to parse gamma transition");
            }

            // Convert from 1-based to 0-based indexing
            finalLevel--;

            gamma.finalLevelIndex = finalLevel;
            gamma.gammaEnergy = gammaE;
            gamma.pg = pg;
            gamma.branchingRatio = br;
            gamma.icc = icc;

            totalBR += br;
            levelData.gammas.push_back(gamma);
        }

        // Normalize branching ratios
        if (totalBR > 0.0 && numGammas > 0) {
            for (auto& gamma : levelData.gammas) {
                gamma.branchingRatio /= totalBR;
            }
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

            // Transition type will be set by Nucleus based on A parity
            transition->setType(Transition::Type::NONE);

            initialLevel->addTransition(transition);
        }
    }

    return levels;
}

} // namespace rainier
