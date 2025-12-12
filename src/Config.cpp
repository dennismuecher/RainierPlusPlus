// Config.cpp - Configuration with YAML support
#include "Config.h"
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace rainier {

Config Config::loadFromFile(const std::string& filename) {
    Config config;
    
    try {
        YAML::Node yaml = YAML::LoadFile(filename);
        
        std::cout << "Loading configuration from: " << filename << std::endl;
        
        // =================================================================
        // NUCLEUS SECTION
        // =================================================================
        if (yaml["nucleus"]) {
            auto nucleus = yaml["nucleus"];
            if (nucleus["Z"]) config.nucleus.Z = nucleus["Z"].as<int>();
            if (nucleus["A"]) config.nucleus.A = nucleus["A"].as<int>();
            if (nucleus["Sn"]) config.nucleus.Sn = nucleus["Sn"].as<double>();
            if (nucleus["levelsFile"]) 
                config.nucleus.levelsFile = nucleus["levelsFile"].as<std::string>();
        }
        
        // =================================================================
        // INITIAL EXCITATION SECTION
        // =================================================================
        if (yaml["initialExcitation"]) {
            auto init = yaml["initialExcitation"];
            
            // Mode selection
            if (init["mode"]) {
                std::string mode = init["mode"].as<std::string>();
                if (mode == "SINGLE") 
                    config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::SINGLE;
                else if (mode == "SELECT") 
                    config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::SELECT;
                else if (mode == "SPREAD") 
                    config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::SPREAD;
                else if (mode == "FULL_REACTION") 
                    config.initialExcitation.mode = Config::InitialExcitationConfig::Mode::FULL_REACTION;
            }
            
            // SINGLE mode parameters
            if (init["SINGLE"]) {
                auto single = init["SINGLE"];
                if (single["excitationEnergy"]) 
                    config.initialExcitation.excitationEnergy = single["excitationEnergy"].as<double>();
                if (single["spin"]) 
                    config.initialExcitation.spin = single["spin"].as<double>();
                if (single["parity"]) 
                    config.initialExcitation.parity = single["parity"].as<int>();
            }
            
            // SELECT mode parameters
            if (init["SELECT"] && init["SELECT"]["states"]) {
                config.initialExcitation.selectStates.clear();
                for (const auto& state : init["SELECT"]["states"]) {
                    Config::InitialExcitationConfig::SelectState s;
                    s.energy = state["energy"].as<double>();
                    s.spin = state["spin"].as<double>();
                    s.parity = state["parity"].as<int>();
                    s.branchingRatio = state["branchingRatio"].as<double>();
                    config.initialExcitation.selectStates.push_back(s);
                }
            }
            
            // SPREAD mode parameters
            if (init["SPREAD"]) {
                auto spread = init["SPREAD"];
                if (spread["meanEnergy"])
                    config.initialExcitation.meanEnergy = spread["meanEnergy"].as<double>();
                if (spread["energySpread"])
                    config.initialExcitation.energySpread = spread["energySpread"].as<double>();
            }
            
            // FULL_REACTION mode parameters
            if (init["FULL_REACTION"] && init["FULL_REACTION"]["populationFile"]) {
                config.initialExcitation.populationFile = 
                    init["FULL_REACTION"]["populationFile"].as<std::string>();
            }
            
            // BETA_DECAY mode parameters
            if (init["BETA_DECAY"]) {
                auto beta = init["BETA_DECAY"];
                if (beta["parentSpin"])
                    config.initialExcitation.parentSpin = beta["parentSpin"].as<double>();
                if (beta["parentParity"])
                    config.initialExcitation.parentParity = beta["parentParity"].as<int>();
                if (beta["Qbeta"])
                    config.initialExcitation.Qbeta = beta["Qbeta"].as<double>();
            }
        }
        
        // =================================================================
        // LEVEL DENSITY SECTION
        // =================================================================
        if (yaml["levelDensity"]) {
            auto ld = yaml["levelDensity"];
            
            // Model selection
            if (ld["model"]) {
                std::string model = ld["model"].as<std::string>();
                if (model == "BSFG") 
                    config.levelDensity.model = Config::LevelDensityConfig::Model::BSFG;
                else if (model == "CTM") 
                    config.levelDensity.model = Config::LevelDensityConfig::Model::CTM;                
                else if (model == "TABLE")
                    config.levelDensity.model = Config::LevelDensityConfig::Model::TABLE;
                else if (model == "USER_DEFINED") 
                    config.levelDensity.model = Config::LevelDensityConfig::Model::USER_DEFINED;
            }
            
            // BSFG parameters
            if (ld["BSFG"]) {
                auto bsfg = ld["BSFG"];
                if (bsfg["a"]) config.levelDensity.a = bsfg["a"].as<double>();
                if (bsfg["E1"]) config.levelDensity.E1 = bsfg["E1"].as<double>();
                if (bsfg["useEnergyDependent"]) 
                    config.levelDensity.useEnergyDependentA = bsfg["useEnergyDependent"].as<bool>();
                if (bsfg["aAsymptotic"]) 
                    config.levelDensity.aAsymptotic = bsfg["aAsymptotic"].as<double>();
                if (bsfg["shellCorrection"]) 
                    config.levelDensity.shellCorrectionW = bsfg["shellCorrection"].as<double>();
                if (bsfg["damping"]) 
                    config.levelDensity.dampingGamma = bsfg["damping"].as<double>();
            }
            
            // CTM parameters
            if (ld["CTM"]) {
                auto ctm = ld["CTM"];
                if (ctm["T"]) config.levelDensity.T = ctm["T"].as<double>();
                if (ctm["E0"]) config.levelDensity.E0 = ctm["E0"].as<double>();
            }
            
            // TABLE parameters
            if (ld["TABLE"] && ld["TABLE"]["file"]) {
                config.levelDensity.tableFile = ld["TABLE"]["file"].as<std::string>();
            }
        }
        
        // =================================================================
        // SPIN CUTOFF SECTION
        // =================================================================
        if (yaml["spinCutoff"]) {
            auto sc = yaml["spinCutoff"];
            
            // Model selection
            if (sc["model"]) {
                std::string model = sc["model"].as<std::string>();
                if (model == "VON_EGIDY_05") 
                    config.spinCutoff.model = Config::SpinCutoffConfig::Model::VON_EGIDY_05;
                else if (model == "VON_EGIDY_09") 
                    config.spinCutoff.model = Config::SpinCutoffConfig::Model::VON_EGIDY_09;
                else if (model == "TALYS") 
                    config.spinCutoff.model = Config::SpinCutoffConfig::Model::TALYS;
                else if (model == "SINGLE_PARTICLE") 
                    config.spinCutoff.model = Config::SpinCutoffConfig::Model::SINGLE_PARTICLE;
                else if (model == "RIGID_SPHERE") 
                    config.spinCutoff.model = Config::SpinCutoffConfig::Model::RIGID_SPHERE;
            }
            
            // TALYS parameters
            if (sc["TALYS"]) {
                auto talys = sc["TALYS"];
                if (talys["spinCutoffD"]) 
                    config.spinCutoff.spinCutoffD = talys["spinCutoffD"].as<double>();
                if (talys["Ed"]) 
                    config.spinCutoff.Ed = talys["Ed"].as<double>();
            }
            
            // Oslo shift (applies to all models)
            if (sc["useOsloShift"])
                config.spinCutoff.useOsloShift = sc["useOsloShift"].as<bool>();
            if (sc["osloShift"])
                config.spinCutoff.osloShift = sc["osloShift"].as<double>();
        }
        
        // =================================================================
        // GAMMA STRENGTH SECTION
        // =================================================================
        if (yaml["gammaStrength"]) {
            auto gs = yaml["gammaStrength"];
            
            // E1 model selection
            if (gs["e1Model"]) {
                std::string model = gs["e1Model"].as<std::string>();
                if (model == "GEN_LORENTZ") 
                    config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::GEN_LORENTZ;
                else if (model == "EGLO") 
                    config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::EGLO;
                else if (model == "KMF") 
                    config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::KMF;
                else if (model == "KOP_CHRIEN") 
                    config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::KOP_CHRIEN;
                else if (model == "STD_LORENTZ") 
                    config.gammaStrength.e1Model = Config::GammaStrengthConfig::E1Model::STD_LORENTZ;
            }
            
            // E1 resonances
            if (gs["e1Resonances"]) {
                config.gammaStrength.e1Resonances.clear();
                for (const auto& res : gs["e1Resonances"]) {
                    Config::GammaStrengthConfig::E1Resonance r;
                    r.energy = res["energy"].as<double>();
                    r.width = res["width"].as<double>();
                    r.sigma = res["sigma"].as<double>();
                    config.gammaStrength.e1Resonances.push_back(r);
                }
            }
            
            if (gs["e1_constantT"]) 
                config.gammaStrength.e1ConstantT = gs["e1_constantT"].as<double>();
            
            // M1 model selection
            if (gs["m1Model"]) {
                std::string model = gs["m1Model"].as<std::string>();
                if (model == "STD_LORENTZ") 
                    config.gammaStrength.m1Model = Config::GammaStrengthConfig::M1Model::STD_LORENTZ;
                else if (model == "SINGLE_PARTICLE") 
                    config.gammaStrength.m1Model = Config::GammaStrengthConfig::M1Model::SINGLE_PARTICLE;
            }
            
            // M1 resonances
            if (gs["m1Resonances"]) {
                config.gammaStrength.m1Resonances.clear();
                for (const auto& res : gs["m1Resonances"]) {
                    Config::GammaStrengthConfig::M1Resonance r;
                    r.energy = res["energy"].as<double>();
                    r.width = res["width"].as<double>();
                    r.sigma = res["sigma"].as<double>();
                    config.gammaStrength.m1Resonances.push_back(r);
                }
            }
            
            // E2 model selection
            if (gs["e2Model"]) {
                std::string model = gs["e2Model"].as<std::string>();
                if (model == "STD_LORENTZ") 
                    config.gammaStrength.e2Model = Config::GammaStrengthConfig::E2Model::STD_LORENTZ;
                else if (model == "SINGLE_PARTICLE") 
                    config.gammaStrength.e2Model = Config::GammaStrengthConfig::E2Model::SINGLE_PARTICLE;
            }
            // M1 upbend
            if (gs["m1StrUpbend"])
                config.gammaStrength.m1StrUpbend = gs["m1StrUpbend"].as<bool>();
            if (gs["m1UpbendConst"])
                config.gammaStrength.m1UpbendConst = gs["m1UpbendConst"].as<double>();
            if (gs["m1UpbendExp"])
                config.gammaStrength.m1UpbendExp = gs["m1UpbendExp"].as<double>();

            //M1 and E2 single particle strength parameters
            if (gs["m1SingleParticleSigma"])
                config.gammaStrength.m1SingleParticleSigma = gs["m1SingleParticleSigma"].as<double>();

            if (gs["e2SingleParticleSigma"])
                config.gammaStrength.e2SingleParticleSigma = gs["e2SingleParticleSigma"].as<double>();
            
            // Width fluctuation
            if (yaml["widthFluctuation"]) {
                auto wfd = yaml["widthFluctuation"];
                if (wfd["model"]) {
                    std::string model = wfd["model"].as<std::string>();
                    if (model == "PTD")
                        config.gammaStrength.wfdModel = Config::GammaStrengthConfig::WFDModel::PTD;
                    else if (model == "NU")
                        config.gammaStrength.wfdModel = Config::GammaStrengthConfig::WFDModel::NU;
                    else if (model == "OFF")
                        config.gammaStrength.wfdModel = Config::GammaStrengthConfig::WFDModel::OFF;
                }
                if (wfd["nuParameter"])
                    config.gammaStrength.nuParameter = wfd["nuParameter"].as<double>();
            }
            
            if (gs["e2Energy"]) config.gammaStrength.e2Energy = gs["e2Energy"].as<double>();
            if (gs["e2Width"]) config.gammaStrength.e2Width = gs["e2Width"].as<double>();
            if (gs["e2Sigma"]) config.gammaStrength.e2Sigma = gs["e2Sigma"].as<double>();
        }
        
        // =================================================================
        // PARITY MODEL SECTION
        // =================================================================
        if (yaml["parityModel"]) {
            auto pm = yaml["parityModel"];
            
            // Model selection
            if (pm["model"]) {
                std::string model = pm["model"].as<std::string>();
                if (model == "EQUIPARTITION") 
                    config.parity.model = Config::ParityConfig::Model::EQUIPARTITION;
                else if (model == "ENERGY_DEPENDENT") 
                    config.parity.model = Config::ParityConfig::Model::ENERGY_DEPENDENT;
            }
            
            // Energy-dependent parameters
            if (pm["ENERGY_DEPENDENT"]) {
                auto ed = pm["ENERGY_DEPENDENT"];
                if (ed["parC"]) config.parity.parC = ed["parC"].as<double>();
                if (ed["parD"]) config.parity.parD = ed["parD"].as<double>();
            }
        }
        
        // =================================================================
        // CONTINUUM SECTION
        // =================================================================
        if (yaml["continuum"]) {
            auto cont = yaml["continuum"];
            
            // Distribution selection
            if (cont["distribution"]) {
                std::string dist = cont["distribution"].as<std::string>();
                if (dist == "POISSON") 
                    config.continuum.distribution = Config::ContinuumConfig::Distribution::POISSON;
                else if (dist == "WIGNER") 
                    config.continuum.distribution = Config::ContinuumConfig::Distribution::WIGNER;
            }
            
            if (cont["energySpacing"]) 
                config.continuum.energySpacing = cont["energySpacing"].as<double>();
            if (cont["forceBinNumber"]) 
                config.continuum.forceBinNumber = cont["forceBinNumber"].as<bool>();
            if (cont["numBins"]) 
                config.continuum.numBins = cont["numBins"].as<int>();
            if (cont["maxDiscreteLevel"]) 
                config.continuum.maxDiscreteLevel = cont["maxDiscreteLevel"].as<int>();
        }
        
        // =================================================================
        // INTERNAL CONVERSION SECTION
        // =================================================================
        if (yaml["internalConversion"]) {
            auto ic = yaml["internalConversion"];
            
            if (ic["enabled"]) 
                config.internalConversion.enabled = ic["enabled"].as<bool>();
            
            // Model selection
            if (ic["model"]) {
                std::string model = ic["model"].as<std::string>();
                if (model == "BRICC") 
                    config.internalConversion.model = Config::InternalConversionConfig::Model::BRICC;
                else if (model == "TABLE") 
                    config.internalConversion.model = Config::InternalConversionConfig::Model::TABLE;
            }
            
            if (ic["briccModel"]) 
                config.internalConversion.briccModel = ic["briccModel"].as<std::string>();
        }
        
        // =================================================================
        // SIMULATION SECTION
        // =================================================================
        if (yaml["simulation"]) {
            auto sim = yaml["simulation"];
            
            if (sim["numRealizations"]) 
                config.simulation.numRealizations = sim["numRealizations"].as<int>();
            if (sim["eventsPerRealization"]) 
                config.simulation.eventsPerRealization = sim["eventsPerRealization"].as<int>();
            if (sim["updateInterval"]) 
                config.simulation.updateInterval = sim["updateInterval"].as<int>();
            if (sim["saveInterval"]) 
                config.simulation.saveInterval = sim["saveInterval"].as<int>();
            if (sim["useParallel"]) 
                config.simulation.useParallel = sim["useParallel"].as<bool>();
            if (sim["randomSeed"]) 
                config.simulation.randomSeed = sim["randomSeed"].as<int>();
        }
        
        // =================================================================
        // OUTPUT SECTION
        // =================================================================
        if (yaml["output"]) {
            auto out = yaml["output"];
            
            if (out["outputFile"]) 
                config.output.outputFile = out["outputFile"].as<std::string>();
            if (out["paramFile"]) 
                config.output.paramFile = out["paramFile"].as<std::string>();
            if (out["saveTree"]) 
                config.output.saveTree = out["saveTree"].as<bool>();
            if (out["saveLevelPopulations"]) 
                config.output.saveLevelPopulations = out["saveLevelPopulations"].as<bool>();
            if (out["gammaSpectrumBins"]) 
                config.output.gammaSpectrumBins = out["gammaSpectrumBins"].as<int>();
            if (out["maxGammaEnergy"]) 
                config.output.maxGammaEnergy = out["maxGammaEnergy"].as<double>();
            if (out["maxPlotSpin"]) 
                config.output.maxPlotSpin = out["maxPlotSpin"].as<double>();
        }
        
        std::cout << "Configuration loaded successfully." << std::endl;
        
        // Validate configuration
        std::string errorMsg;
        if (!config.validate(errorMsg)) {
            throw std::runtime_error("Configuration validation failed: " + errorMsg);
        }
        
        return config;
        
    } catch (const YAML::Exception& e) {
        std::cerr << "YAML parsing error: " << e.what() << std::endl;
        std::cerr << "Using default configuration.\n";
        return createDefault(60, 144);
    } catch (const std::exception& e) {
        std::cerr << "Error loading configuration: " << e.what() << std::endl;
        std::cerr << "Using default configuration.\n";
        return createDefault(60, 144);
    }
}

void Config::saveToFile(const std::string& filename) const {
    YAML::Emitter out;
    out << YAML::BeginMap;
    
    // Nucleus section
    out << YAML::Key << "nucleus";
    out << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "Z" << YAML::Value << nucleus.Z;
    out << YAML::Key << "A" << YAML::Value << nucleus.A;
    out << YAML::Key << "Sn" << YAML::Value << nucleus.Sn;
    out << YAML::Key << "levelsFile" << YAML::Value << nucleus.levelsFile;
    out << YAML::EndMap;
    
    // Initial excitation section
    out << YAML::Key << "initialExcitation";
    out << YAML::Value << YAML::BeginMap;
    
    std::string modeStr;
    switch(initialExcitation.mode) {
        case InitialExcitationConfig::Mode::SINGLE: modeStr = "SINGLE"; break;
        case InitialExcitationConfig::Mode::SELECT: modeStr = "SELECT"; break;
        case InitialExcitationConfig::Mode::SPREAD: modeStr = "SPREAD"; break;
        case InitialExcitationConfig::Mode::FULL_REACTION: modeStr = "FULL_REACTION"; break;
        case InitialExcitationConfig::Mode::BETA_DECAY: modeStr = "BETA_DECAY"; break;

    }
    out << YAML::Key << "mode" << YAML::Value << modeStr;
    
    // Add other sections as needed...
    out << YAML::EndMap;
    
    out << YAML::EndMap;
    
    std::ofstream fout(filename);
    fout << out.c_str();
    fout.close();
    
    std::cout << "Configuration saved to: " << filename << std::endl;
}

Config Config::createDefault(int Z, int A) {
    Config config;
    config.nucleus.Z = Z;
    config.nucleus.A = A;
    config.nucleus.Sn = 7.0;
    
    // Set reasonable defaults for other parameters
    config.levelDensity.a = 15.0;
    config.levelDensity.E1 = 1.0;
    
    config.initialExcitation.excitationEnergy = 7.0;
    config.initialExcitation.spin = 0.5;
    config.initialExcitation.parity = 1;
    
    return config;
}

bool Config::validate(std::string& errorMsg) const {
    // Check nucleus parameters
    if (nucleus.Z <= 0 || nucleus.Z > 120) {
        errorMsg = "Invalid Z (must be 1-120): " + std::to_string(nucleus.Z);
        return false;
    }
    
    if (nucleus.A <= nucleus.Z || nucleus.A > 300) {
        errorMsg = "Invalid A (must be > Z and < 300): " + std::to_string(nucleus.A);
        return false;
    }
    
    if (nucleus.Sn < 0.0) {
        errorMsg = "Negative neutron separation energy";
        return false;
    }
    
    // Check simulation parameters
    if (simulation.numRealizations <= 0) {
        errorMsg = "Number of realizations must be positive";
        return false;
    }
    
    if (simulation.eventsPerRealization <= 0) {
        errorMsg = "Events per realization must be positive";
        return false;
    }
    
    // Check SELECT mode branching ratios
    if (initialExcitation.mode == InitialExcitationConfig::Mode::SELECT) {
        double sum = 0.0;
        for (const auto& state : initialExcitation.selectStates) {
            sum += state.branchingRatio;
        }
        if (std::abs(sum - 1.0) > 1e-6) {
            errorMsg = "SELECT mode branching ratios must sum to 1.0 (sum = " + 
                      std::to_string(sum) + ")";
            return false;
        }
        
    }
    
    if (initialExcitation.mode == InitialExcitationConfig::Mode::BETA_DECAY) {
        if (initialExcitation.Qbeta <= 0.0) {
            errorMsg = "Qbeta must be positive for BETA_DECAY mode";
            return false;
        }
    }
    
    return true;
}

void Config::print() const {
    std::cout << "\n=== RAINIER++ Configuration ===\n";
    std::cout << "Nucleus: Z=" << nucleus.Z << " A=" << nucleus.A 
              << " Sn=" << nucleus.Sn << " MeV\n";
    std::cout << "Levels file: " << nucleus.levelsFile << "\n";
    
    std::cout << "\nInitial Excitation Mode: ";
    switch(initialExcitation.mode) {
        case InitialExcitationConfig::Mode::SINGLE: 
            std::cout << "SINGLE (E=" << initialExcitation.excitationEnergy 
                     << " MeV, J=" << initialExcitation.spin << ")\n";
            break;
        case InitialExcitationConfig::Mode::SELECT: 
            std::cout << "SELECT (" << initialExcitation.selectStates.size() << " states)\n";
            break;
        case InitialExcitationConfig::Mode::SPREAD:
            std::cout << "SPREAD (mean=" << initialExcitation.meanEnergy
                     << " MeV, width=" << initialExcitation.energySpread << " MeV)\n";
            break;
        case InitialExcitationConfig::Mode::FULL_REACTION: 
            std::cout << "FULL_REACTION\n";
            break;
        case InitialExcitationConfig::Mode::BETA_DECAY:
            std::cout << "BETA_DECAY (parent J=" << initialExcitation.parentSpin
                     << (initialExcitation.parentParity > 0 ? "+" : "-")
                     << ", Q=" << initialExcitation.Qbeta << " MeV)\n";
            break;
    }
    
    std::cout << "\nSimulation: " << simulation.numRealizations << " realizations, "
              << simulation.eventsPerRealization << " events each\n";
    std::cout << "Output: " << output.outputFile << "\n";
    std::cout << "===============================\n\n";
}

} // namespace rainier
