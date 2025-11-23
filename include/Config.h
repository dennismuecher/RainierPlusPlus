// Config.h - Complete configuration management for RAINIER with YAML support
#ifndef RAINIER_CONFIG_H
#define RAINIER_CONFIG_H

#include <string>
#include <vector>

namespace rainier {

/**
 * @brief Configuration class for RAINIER simulation
 * 
 * Supports loading configuration from YAML files using yaml-cpp library.
 * YAML provides cleaner syntax, native comment support, and better readability
 * compared to JSON.
 */
class Config {
public:
    // Nucleus properties
    struct NucleusConfig {
        int Z = 60;
        int A = 144;
        double Sn = 7.0;
        std::string levelsFile = "levels/z060.dat";
    };

    // Level density model configuration
    struct LevelDensityConfig {
        enum class Model { BSFG, CTM, TABLE, USER_DEFINED };
        Model model = Model::BSFG;
        
        // BSFG parameters
        double a = 15.0;              // Level density parameter (MeV^-1)
        double E1 = 1.0;              // Back-shift parameter (MeV)
        
        // CTM parameters
        double T = 0.5;               // Nuclear temperature (MeV)
        double E0 = 0.0;              // Energy shift (MeV)
        
        // Energy-dependent a(E) parameters
        bool useEnergyDependentA = true;
        double aAsymptotic = 15.0;
        double shellCorrectionW = 5.0;
        double dampingGamma = 0.05;
        
        std::string tableFile = "";   // For TABLE model
    };

    struct SpinCutoffConfig {
        enum class Model { VON_EGIDY_05, SINGLE_PARTICLE, RIGID_SPHERE,
                          VON_EGIDY_09, TALYS, USER_DEFINED };
        Model model = Model::TALYS;
        
        double spinCutoffD = 3.0;     // Discrete region spin cutoff
        double Ed = 2.0;               // Matching energy (MeV)
        
        // Oslo energy shift (applies to all models when enabled)
        bool useOsloShift = false;    // Enable Oslo-style energy shift
        double osloShift = 0.12;      // Energy shift (MeV)
    };

    // Gamma strength function configuration
    struct GammaStrengthConfig {
        enum class E1Model { GEN_LORENTZ, EGLO, KMF, KOP_CHRIEN, STD_LORENTZ, USER_DEFINED };
        E1Model e1Model = E1Model::GEN_LORENTZ;
        
        struct E1Resonance {
            double energy = 15.0;
            double width = 5.0;
            double sigma = 300.0;
        };
        std::vector<E1Resonance> e1Resonances;
        double e1ConstantT = 0.5;
        
        enum class M1Model { STD_LORENTZ, SINGLE_PARTICLE, USER_DEFINED };
        M1Model m1Model = M1Model::STD_LORENTZ;
        
        struct M1Resonance {
            double energy = 7.0;
            double width = 4.0;
            double sigma = 1.0;
        };
        std::vector<M1Resonance> m1Resonances;
        
        enum class E2Model { STD_LORENTZ, SINGLE_PARTICLE, USER_DEFINED };
        E2Model e2Model = E2Model::STD_LORENTZ;
        double e2Energy = 12.0;
        double e2Width = 4.0;
        double e2Sigma = 5.0;
        
        // Width Fluctuation Distribution
        enum class WFDModel { PTD, NU, OFF };
        WFDModel wfdModel = WFDModel::PTD;
        double nuParameter = 0.5;  // For NU model
            
        // M1 upbend parameters
        bool m1StrUpbend = false;
        double m1UpbendConst = 5e-8;  // C parameter
        double m1UpbendExp = 1.0;     // A parameter (positive)
        
        // M1 single particle parameters
        double m1SingleParticleSigma = 4e-11;  // MeV^-3

        // E2 single particle parameters
        double e2SingleParticleSigma = 4e-11;  // MeV^-5
    };

    // Initial excitation configuration
    struct InitialExcitationConfig {
        enum class Mode { SINGLE, SELECT, SPREAD, FULL_REACTION };
        Mode mode = Mode::SPREAD;
        
        // SINGLE mode parameters
        double excitationEnergy = 7.0;
        double spin = 0.5;
        int parity = 1;
        
        // SELECT mode parameters
        struct SelectState {
            double energy;
            double spin;
            int parity;
            double branchingRatio;
        };
        std::vector<SelectState> selectStates;
    
        // SPREAD mode parameters
        double meanEnergy = 7.0;
        double energySpread = 0.5;
        
        // FULL_REACTION mode
        std::string populationFile = "";
    };

    // Parity distribution configuration
    struct ParityConfig {
        enum class Model { EQUIPARTITION, ENERGY_DEPENDENT };
        Model model = Model::ENERGY_DEPENDENT;
        
        // Energy-dependent parameters
        double parC = 0.5;
        double parD = 2.0;
    };

    // Continuum level configuration
    struct ContinuumConfig {
        enum class Distribution { POISSON, WIGNER };
        Distribution distribution = Distribution::POISSON;
        
        double energySpacing = 0.025;
        bool forceBinNumber = false;
        int numBins = 100;
        int maxDiscreteLevel = 15;
    };

    // Internal conversion configuration
    struct InternalConversionConfig {
        enum class Model { BRICC, TABLE };
        Model model = Model::BRICC;
        
        bool enabled = true;
        std::string briccModel = "F";  // F, R, or S
        std::string tableFile = "";
    };

    // Simulation parameters
    struct SimulationConfig {
        int numRealizations = 1;
        int eventsPerRealization = 10000;
        int updateInterval = 100;
        int saveInterval = 10000;
        bool useParallel = false;
        int randomSeed = 12345;
    };

    // Output configuration
    struct OutputConfig {
        std::string outputFile = "output.root";
        std::string paramFile = "parameters.dat";
        bool saveTree = true;
        bool saveLevelPopulations = true;
        int gammaSpectrumBins = 1000;
        double maxGammaEnergy = 20.0;
        double maxPlotSpin = 20.0;
    };

    // Configuration sections
    NucleusConfig nucleus;
    LevelDensityConfig levelDensity;
    SpinCutoffConfig spinCutoff;
    GammaStrengthConfig gammaStrength;
    InitialExcitationConfig initialExcitation;
    ParityConfig parity;
    ContinuumConfig continuum;
    InternalConversionConfig internalConversion;
    SimulationConfig simulation;
    OutputConfig output;

    /**
     * @brief Load configuration from YAML file
     * @param filename Path to YAML configuration file
     * @return Loaded configuration
     * @throws std::runtime_error if file cannot be opened or parsed
     */
    static Config loadFromFile(const std::string& filename);

    /**
     * @brief Save configuration to YAML file
     * @param filename Path to output file
     */
    void saveToFile(const std::string& filename) const;

    /**
     * @brief Create default configuration
     * @param Z Atomic number
     * @param A Mass number
     * @return Default configuration for specified nucleus
     */
    static Config createDefault(int Z, int A);

    /**
     * @brief Validate configuration
     * @param errorMsg Output parameter for error message
     * @return true if valid, false otherwise
     */
    bool validate(std::string& errorMsg) const;

    /**
     * @brief Print configuration to console
     */
    void print() const;
};

} // namespace rainier

#endif // RAINIER_CONFIG_H
