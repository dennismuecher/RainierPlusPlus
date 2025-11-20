// Config.h - Complete configuration management for RAINIER
#ifndef RAINIER_CONFIG_H
#define RAINIER_CONFIG_H

#include <string>
#include <vector>

namespace rainier {

/**
 * @brief Configuration class for RAINIER simulation
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

    // Spin cutoff model configuration
    struct SpinCutoffConfig {
        enum class Model { VON_EGIDY_05, SINGLE_PARTICLE, RIGID_SPHERE, 
                          VON_EGIDY_09, TALYS, USER_DEFINED };
        Model model = Model::TALYS;
        
        double spinCutoffD = 3.0;     // Discrete region spin cutoff
        double Ed = 2.0;               // Matching energy (MeV)
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
        std::vector<E1Resonance> e1Resonances = {{15.0, 5.0, 300.0}};
        double constantT = 0.5;
        
        enum class M1Model { STD_LORENTZ, SINGLE_PARTICLE, USER_DEFINED };
        M1Model m1Model = M1Model::STD_LORENTZ;
        
        struct M1Resonance {
            double energy = 7.0;
            double width = 4.0;
            double sigma = 1.0;
        };
        std::vector<M1Resonance> m1Resonances = {{7.0, 4.0, 1.0}};
        
        enum class E2Model { STD_LORENTZ, SINGLE_PARTICLE, USER_DEFINED };
        E2Model e2Model = E2Model::STD_LORENTZ;
        double e2Energy = 12.0;
        double e2Width = 4.0;
        double e2Sigma = 5.0;
    };

	// Internal conversion configuration
	struct InternalConversionConfig {
	    bool enabled = true;
	    enum class Model { BRICC, TABLE };
	    Model model = Model::BRICC;
	    std::string briccModel = "F";  // Frozen orbital
	};

    // Parity distribution configuration
    struct ParityConfig {
        enum class Model { EQUIPARTITION, ENERGY_DEPENDENT };
        Model model = Model::ENERGY_DEPENDENT;
        
        double parC = 0.5;
        double parD = 2.0;
    };

    // Continuum level construction
    struct ContinuumConfig {
        enum class Distribution { POISSON, WIGNER };
        Distribution distribution = Distribution::POISSON;
        
        double energySpacing = 0.1;   // MeV
        bool forceBinNumber = false;
        int numBins = 100;
        
        int maxDiscreteLevel = 20;
    };

    // Initial excitation configuration
    struct InitialExcitationConfig {
        enum class Mode { SINGLE, SELECT, SPREAD, FULL_REACTION };
        Mode mode = Mode::SINGLE;
        
        // SINGLE mode parameters
        double excitationEnergy = 8.0;  // MeV
        double spin = 0.5;
        int parity = 1;
        
        // SELECT mode parameters (beta-decay-like)
        struct SelectState {
            double energy;
            double spin;
            int parity;
            double branchingRatio;
            
            SelectState() : energy(0), spin(0), parity(1), branchingRatio(0) {}
            SelectState(double e, double s, int p, double br)
                : energy(e), spin(s), parity(p), branchingRatio(br) {}
        };
        std::vector<SelectState> selectStates;
        
        // SPREAD mode parameters
        double meanEnergy = 8.0;
        double energySpread = 0.5;
        
        // FULL_REACTION mode parameters
        std::string populationFile;
    };


    // Simulation parameters
    struct SimulationConfig {
        int numRealizations = 1;
        int eventsPerRealization = 100;
        int updateInterval = 100;
        int saveInterval = 10000;
        bool useParallel = false;
        unsigned int randomSeed = 12345;
    };

    // Output configuration
    struct OutputConfig {
        std::string outputFile = "output.root";
        std::string paramFile = "parameters.dat";
        bool saveTree = true;
        bool saveLevelPopulations = true;
        int gammaSpectrumBins = 1000;
        double maxGammaEnergy = 20.0;  // MeV
        double maxPlotSpin = 10.0;
    };

    // Member variables
    NucleusConfig nucleus;
    LevelDensityConfig levelDensity;
    SpinCutoffConfig spinCutoff;
    GammaStrengthConfig gammaStrength;
    ParityConfig parity;
    ContinuumConfig continuum;
    InitialExcitationConfig initialExcitation;
    SimulationConfig simulation;
	InternalConversionConfig internalConversion;
    OutputConfig output;

    // Methods
    Config() = default;
    
    static Config loadFromFile(const std::string& filename);
    void saveToFile(const std::string& filename) const;
    static Config createDefault(int Z, int A);
    bool validate(std::string& errorMessage) const;
    void print() const;
};

} // namespace rainier

#endif // RAINIER_CONFIG_H
