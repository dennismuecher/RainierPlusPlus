// OutputManager.cpp
#include "io/OutputManager.h"
#include "simulation/DecaySimulator.h"
#include "simulation/CascadeEvent.h"
#include "core/Nucleus.h"
#include "core/Level.h"
#include "core/DiscreteLevel.h"
#include "core/ContinuumLevel.h"
#include "models/LevelDensity.h"
#include "models/SpinCutoff.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <sstream>
#include <fstream>

namespace rainier {

OutputManager::OutputManager(const Config& config, int runNumber)
    : config_(config), runNumber_(runNumber), outputFile_(nullptr), cascadeTree_(nullptr) {
    
    // Create output file
    std::string filename = "output" + std::to_string(runNumber) + ".root";
    outputFile_ = new TFile(filename.c_str(), "RECREATE");
    
    if (!outputFile_ || outputFile_->IsZombie()) {
        throw std::runtime_error("Failed to create output file: " + filename);
    }
    
    std::cout << "Created output file: " << filename << std::endl;
    
    // Create tree if requested
    if (config_.output.saveTree) {
        cascadeTree_ = new TTree("cascadeTree", "Gamma-ray cascade data");
        cascadeTree_->Branch("gammaEnergies", &treeGammaEnergies_);
        cascadeTree_->Branch("finalEnergies", &treeFinalEnergies_);
        cascadeTree_->Branch("times", &treeTimes_);
        cascadeTree_->Branch("initialEx", &treeInitialEx_);
        cascadeTree_->Branch("initialSpin", &treeInitialSpin_);
        cascadeTree_->Branch("initialParity", &treeInitialParity_);
    }
}

OutputManager::~OutputManager() {
    finalize();
}

void OutputManager::createHistograms(int realization) {
    std::string suffix = "_real" + std::to_string(realization);
    
    // Gamma spectrum
    std::string name = "hGammaSpectrum" + suffix;
    auto hGamma = new TH1D(name.c_str(), 
                           ("Gamma Spectrum, Real " + std::to_string(realization)).c_str(),
                           config_.output.gammaSpectrumBins, 0.0, config_.output.maxGammaEnergy);
    hGamma->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    hGamma->GetYaxis()->SetTitle("Counts");
    histograms1D_[name] = hGamma;
    
    // Level population (E vs J with parity sign)
    name = "h2LevelPopulation" + suffix;
    auto h2Pop = new TH2D(name.c_str(),
                          ("Level Population, Real " + std::to_string(realization)).c_str(),
                          2 * static_cast<int>(config_.output.maxPlotSpin), 
                          -config_.output.maxPlotSpin, 
                          config_.output.maxPlotSpin,
                          200, 0.0, config_.initialExcitation.excitationEnergy);
    h2Pop->GetXaxis()->SetTitle("J#pi");
    h2Pop->GetYaxis()->SetTitle("E_{x} (MeV)");
    h2Pop->GetZaxis()->SetTitle("Counts");
    histograms2D_[name] = h2Pop;
    
    // Ex vs Egamma (for Oslo method analysis)
    name = "h2ExEgamma" + suffix;
    auto h2ExEg = new TH2D(name.c_str(),
                           ("E_{x} vs E_{#gamma}, Real " + std::to_string(realization)).c_str(),
                           300, 0.0, config_.initialExcitation.excitationEnergy * 1000,  // keV
                           300, 0.0, config_.initialExcitation.excitationEnergy * 1000);
    h2ExEg->GetXaxis()->SetTitle("E_{#gamma} (keV)");
    h2ExEg->GetYaxis()->SetTitle("E_{x} (keV)");
    h2ExEg->GetZaxis()->SetTitle("Counts");
    histograms2D_[name] = h2ExEg;
    
    // First generation only
    name = "h2FirstGen" + suffix;
    auto h2First = new TH2D(name.c_str(),
                            ("First Generation, Real " + std::to_string(realization)).c_str(),
                            300, 0.0, config_.initialExcitation.excitationEnergy * 1000,
                            300, 0.0, config_.initialExcitation.excitationEnergy * 1000);
    h2First->GetXaxis()->SetTitle("E_{#gamma} (keV)");
    h2First->GetYaxis()->SetTitle("E_{x} (keV)");
    h2First->GetZaxis()->SetTitle("Counts");
    histograms2D_[name] = h2First;
    
    // Multiplicity
    name = "hMultiplicity" + suffix;
    auto hMult = new TH1D(name.c_str(),
                          ("Cascade Multiplicity, Real " + std::to_string(realization)).c_str(),
                          50, 0, 50);
    hMult->GetXaxis()->SetTitle("Number of Gammas");
    hMult->GetYaxis()->SetTitle("Counts");
    histograms1D_[name] = hMult;
    
    // Cascade time
    name = "hCascadeTime" + suffix;
    auto hTime = new TH1D(name.c_str(),
                          ("Total Cascade Time, Real " + std::to_string(realization)).c_str(),
                          200, 0, 1000);  // fs
    hTime->GetXaxis()->SetTitle("Time (fs)");
    hTime->GetYaxis()->SetTitle("Counts");
    histograms1D_[name] = hTime;
    
    // Initial spin distribution
    name = "hInitialSpin" + suffix;
    auto hSpinI = new TH1D(name.c_str(),
                           ("Initial Spin Distribution, Real " + std::to_string(realization)).c_str(),
                           static_cast<int>(config_.output.maxPlotSpin) + 1, 
                           0, config_.output.maxPlotSpin);
    hSpinI->GetXaxis()->SetTitle("J_{i} (#hbar)");
    hSpinI->GetYaxis()->SetTitle("Counts");
    histograms1D_[name] = hSpinI;
}

void OutputManager::createLevelSpectraHistograms(int realization) {
    std::string suffix = "_real" + std::to_string(realization);
    
    // Determine number of spin bins based on whether A is even or odd
    bool isEvenA = (config_.nucleus.A % 2 == 0);
    int maxSpinBin = static_cast<int>(config_.output.maxPlotSpin);
    std::cout <<"Max Spin Bin in createLevelSpectraHistograms: " <<maxSpinBin<<std::endl;
    // Create one histogram for each spin value
    for (int spinBin = 0; spinBin <= maxSpinBin; ++spinBin) {
        double spin = Level::binToSpin(spinBin, isEvenA);
        
        // Format spin value for histogram name (handle half-integer spins)
        std::ostringstream spinStr;
        if (isEvenA) {
            spinStr << spinBin;
        } else {
            // For odd-A, spins are half-integer
            spinStr << spinBin << "p5";  // e.g., "0p5" for 0.5, "1p5" for 1.5
        }
        
        // Create histogram for this spin (all parities combined)
        std::string name = "hLevelSpectrum_J" + spinStr.str() + suffix;
        
        // Format title with proper spin notation
        std::ostringstream titleStream;
        titleStream << "Level Spectrum J=";
        if (isEvenA) {
            titleStream << spinBin;
        } else {
            titleStream << spinBin << ".5";
        }
        titleStream << ", Real " << realization;
        
        auto hSpectrum = new TH1D(name.c_str(),
                                  titleStream.str().c_str(),
                                  1000,  // Number of bins
                                  0.0,  // Min energy
                                  config_.initialExcitation.excitationEnergy + 2.0);  // Max energy with buffer
        
        hSpectrum->GetXaxis()->SetTitle("E_{x} (MeV)");
        hSpectrum->GetYaxis()->SetTitle("Number of Levels");
        hSpectrum->SetLineColor(spinBin % 9 + 1);  // Different color for each spin
        hSpectrum->SetLineWidth(2);
        
        levelSpectraHistograms_[name] = hSpectrum;
    }
    
    std::cout << "Created " << levelSpectraHistograms_.size() 
              << " level spectrum histograms for realization " << realization << std::endl;
}

void OutputManager::fillLevelSpectra(int realization, const Nucleus& nucleus) {
    std::string suffix = "_real" + std::to_string(realization);
    bool isEvenA = nucleus.isEvenA();
    
    std::cout << "Filling level spectra histograms for realization " << realization << "..." << std::endl;
    
    std::cout << "getNumDiscreteLevels " << nucleus.getNumDiscreteLevels() << "..." << std::endl;
    // First, fill from discrete levels
    for (int i = 0; i < nucleus.getNumDiscreteLevels(); ++i) {
        auto level = nucleus.getDiscreteLevel(i);
        double energy = level->getEnergy();
        double spin = level->getSpin();
        int spinBin = Level::spinToBin(spin, isEvenA);
        
        // Format histogram name
        std::ostringstream spinStr;
        if (isEvenA) {
            spinStr << spinBin;
        } else {
            spinStr << spinBin << "p5";
        }
        
        std::string name = "hLevelSpectrum_J" + spinStr.str() + suffix;
        
        if (levelSpectraHistograms_.count(name)) {
            levelSpectraHistograms_[name]->Fill(energy);
        }
    }
    
    // Now fill from continuum levels
    int totalLevels = 0;
    int maxSpinBin = nucleus.getMaxSpinBin();
    
    for (int ex = 0; ex < nucleus.getNumEnergyBins(); ++ex) {
        for (int spinBin = 0; spinBin <= maxSpinBin; ++spinBin) {
            // Sum over both parities
            for (int parity = 0; parity <= 1; ++parity) {
                const auto& levels = nucleus.getContinuumLevels(ex, spinBin, parity);
                
                // Format histogram name
                std::ostringstream spinStr;
                if (isEvenA) {
                    spinStr << spinBin;
                } else {
                    spinStr << spinBin << "p5";
                }
                
                std::string name = "hLevelSpectrum_J" + spinStr.str() + suffix;
                
                if (levelSpectraHistograms_.count(name)) {
                    for (const auto& level : levels) {
                        levelSpectraHistograms_[name]->Fill(level->getEnergy());
                        totalLevels++;
                    }
                }
            }
        }
    }
    
    std::cout << "  Filled " << totalLevels << " continuum levels into spectra histograms" << std::endl;
}

void OutputManager::fillHistograms(int realization, const CascadeEvent& event) {
    std::string suffix = "_real" + std::to_string(realization);
    
    // Get initial state
    auto initialLevel = event.getInitialLevel();
    if (!initialLevel) return;
    
    double initialEx = event.getInitialExcitation();
    double initialSpin = initialLevel->getSpin();
    int initialParity = initialLevel->getParity();
    
    // Fill initial spin
    std::string name = "hInitialSpin" + suffix;
    if (histograms1D_.count(name)) {
        histograms1D_[name]->Fill(initialSpin);
    }
    
    // Fill multiplicity and cascade time
    int numGammas = event.getNumGammas();
    name = "hMultiplicity" + suffix;
    if (histograms1D_.count(name)) {
        histograms1D_[name]->Fill(numGammas);
    }
    
    name = "hCascadeTime" + suffix;
    if (histograms1D_.count(name)) {
        histograms1D_[name]->Fill(event.getTotalTime());
    }
    
    // Fill level populations and gamma spectra
    bool isFirstStep = true;
    for (const auto& step : event.getSteps()) {
        if (step.isElectron) continue;  // Skip electrons
        
        // Gamma spectrum
        name = "hGammaSpectrum" + suffix;
        if (histograms1D_.count(name)) {
            histograms1D_[name]->Fill(step.gammaEnergy);
        }
        
        // Level population (use parity sign convention)
        name = "h2LevelPopulation" + suffix;
        if (histograms2D_.count(name)) {
            auto initLevel = step.initialLevel;
            if (initLevel) {
                double spinSigned = initLevel->getSpin();
                if (initLevel->getParity() == 0) {
                    spinSigned = -spinSigned - 0.5;  // Negative parity
                } else {
                    spinSigned = spinSigned + 0.5;   // Positive parity
                }
                histograms2D_[name]->Fill(spinSigned, initLevel->getEnergy());
            }
        }
        
        // Ex vs Egamma
        name = "h2ExEgamma" + suffix;
        if (histograms2D_.count(name)) {
            histograms2D_[name]->Fill(step.gammaEnergy * 1000,  // keV
                                     initialEx * 1000);          // keV
        }
        
        // First generation
        if (isFirstStep) {
            name = "h2FirstGen" + suffix;
            if (histograms2D_.count(name)) {
                histograms2D_[name]->Fill(step.gammaEnergy * 1000, initialEx * 1000);
            }
            isFirstStep = false;
        }
    }
}

void OutputManager::createLevelDensityHistograms(int realization) {
    std::string suffix = "_real" + std::to_string(realization);
    
    // Total level density from model (all J, both parities)
    std::string name = "hTotalLevelDensity" + suffix;
    auto hTotal = new TH1D(name.c_str(),
                           ("Total Level Density (Model), Real " + std::to_string(realization)).c_str(),
                           200,  // Number of bins
                           0.0,  // Min energy
                           config_.initialExcitation.excitationEnergy + 2.0);  // Max energy
    hTotal->GetXaxis()->SetTitle("E_{x} (MeV)");
    hTotal->GetYaxis()->SetTitle("#rho(E) (levels/MeV)");
    hTotal->SetLineColor(kBlue);
    hTotal->SetLineWidth(2);
    histograms1D_[name] = hTotal;
    
    // Discrete level density (binned count of discrete levels)
    name = "hDiscreteLevelDensity" + suffix;
    auto hDiscrete = new TH1D(name.c_str(),
                              ("Discrete Level Density, Real " + std::to_string(realization)).c_str(),
                              200,  // Same binning as total
                              0.0,
                              config_.initialExcitation.excitationEnergy + 2.0);
    hDiscrete->GetXaxis()->SetTitle("E_{x} (MeV)");
    hDiscrete->GetYaxis()->SetTitle("#rho(E) (levels/MeV)");
    hDiscrete->SetLineColor(kRed);
    hDiscrete->SetLineWidth(2);
    hDiscrete->SetFillStyle(3004);  // Hatched fill
    hDiscrete->SetFillColor(kRed);
    histograms1D_[name] = hDiscrete;
}
void OutputManager::fillLevelDensityHistograms(int realization, const DecaySimulator& simulator) {
    std::string suffix = "_real" + std::to_string(realization);
    
    // Get the level density model (it returns a reference, not a pointer)
    const auto& levelDensity = simulator.getLevelDensityModel();
    const auto& spinCutoff = simulator.getSpinCutoffModel();
    const auto& nucleus = simulator.getNucleus();
    
    bool isEvenA = nucleus.isEvenA();
    int maxSpinBin = nucleus.getMaxSpinBin();
    
    // Fill total level density from model
    std::string name = "hTotalLevelDensity" + suffix;
    if (histograms1D_.count(name)) {
        auto hTotal = histograms1D_[name];
        
        for (int bin = 1; bin <= hTotal->GetNbinsX(); ++bin) {
            double energy = hTotal->GetBinCenter(bin);
            
            if (energy < 0.001) continue;  // Skip very low energies
            
            // Sum over all spins and both parities
            double totalDensity = 0.0;
            
            for (int spinBin = 0; spinBin <= maxSpinBin; ++spinBin) {
                double spin = Level::binToSpin(spinBin, isEvenA);
                
                // Sum both parities (assuming equal distribution)
                for (int parity = 0; parity <= 1; ++parity) {
                    double density = levelDensity.getDensity(energy, spin, parity);
                    //std::cout <<"Level density for energy" <<energy<< " and spin " << spin << " is " <<density << std::endl;
                    totalDensity += density;
                }
            }
            
            hTotal->SetBinContent(bin, totalDensity);
        }
    }
    
    // Fill discrete level density using a sliding window approach
    // USE ALL LEVELS FROM FILE, not just the truncated set
    name = "hDiscreteLevelDensity" + suffix;
    if (histograms1D_.count(name)) {
        auto hDiscrete = histograms1D_[name];
        
        // Define the energy window for counting levels (e.g., 1 MeV)
        double energyWindow = 1.0;  // MeV - adjustable parameter
        
        for (int bin = 1; bin <= hDiscrete->GetNbinsX(); ++bin) {
            double energy = hDiscrete->GetBinCenter(bin);
            
            // Count levels within ±energyWindow/2 of this energy
            double eMin = energy - energyWindow / 2.0;
            double eMax = energy + energyWindow / 2.0;
            
            int count = 0;
            // Use ALL discrete levels from file (not truncated)
            for (int i = 0; i < nucleus.getNumAllDiscreteLevels(); ++i) {
                auto level = nucleus.getAllDiscreteLevel(i);
                double levelEnergy = level->getEnergy();
                
                if (levelEnergy >= eMin && levelEnergy < eMax) {
                    count++;
                }
            }
            
            // Level density = number of levels / energy window
            double density = static_cast<double>(count) / energyWindow;
            hDiscrete->SetBinContent(bin, density);
        }
    }
    
    std::cout << "  Filled level density histograms (discrete window: 1.0 MeV)" << std::endl;
    std::cout << "  Used " << nucleus.getNumAllDiscreteLevels()
              << " levels from file (vs " << nucleus.getNumDiscreteLevels()
              << " for simulation)" << std::endl;
}

void OutputManager::saveRealization(int realization, const DecaySimulator& simulator) {
    // Create histograms for this realization if not already created
    std::string suffix = "_real" + std::to_string(realization);
    std::string testName = "hGammaSpectrum" + suffix;
    
    if (histograms1D_.count(testName) == 0) {
        createHistograms(realization);
        createLevelSpectraHistograms(realization);
        createLevelDensityHistograms(realization);
    }
    
    // Fill level spectra from the nucleus
    fillLevelSpectra(realization, simulator.getNucleus());
    
    // Fill level density histograms
    fillLevelDensityHistograms(realization, simulator);  // ADD THIS LINE
       
    // Fill event-by-event histograms
    if (config_.output.saveTree) {
        for (const auto& event : simulator.getEvents()) {
            fillHistograms(realization, event);
            
            // Fill tree
            treeGammaEnergies_.clear();
            treeFinalEnergies_.clear();
            treeTimes_.clear();
            
            auto initialLevel = event.getInitialLevel();
            if (initialLevel) {
                treeInitialEx_ = event.getInitialExcitation();
                treeInitialSpin_ = static_cast<int>(initialLevel->getSpin());
                treeInitialParity_ = initialLevel->getParity();
                
                for (const auto& step : event.getSteps()) {
                    if (!step.isElectron) {  // Only gammas
                        treeGammaEnergies_.push_back(step.gammaEnergy);
                        treeFinalEnergies_.push_back(step.finalLevel->getEnergy());
                        treeTimes_.push_back(step.timeToDecay);
                    }
                }
                
                cascadeTree_->Fill();
            }
        }
    }
}

void OutputManager::saveParameters() {
    std::string filename = "parameters" + std::to_string(runNumber_) + ".dat";
    std::ofstream out(filename);
    
    if (!out) {
        std::cerr << "Warning: Could not create parameter file: " << filename << std::endl;
        return;
    }
    
    out << "# RAINIER++ Run Parameters\n";
    out << "# Run number: " << runNumber_ << "\n\n";
    
    out << "[Nucleus]\n";
    out << "Z = " << config_.nucleus.Z << "\n";
    out << "A = " << config_.nucleus.A << "\n";
    out << "Sn = " << config_.nucleus.Sn << "\n";
    out << "levelsFile = " << config_.nucleus.levelsFile << "\n\n";
    
    out << "[Simulation]\n";
    out << "numRealizations = " << config_.simulation.numRealizations << "\n";
    out << "eventsPerRealization = " << config_.simulation.eventsPerRealization << "\n";
    out << "randomSeed = " << config_.simulation.randomSeed << "\n\n";
    
    out << "[InitialExcitation]\n";
    out << "excitationEnergy = " << config_.initialExcitation.excitationEnergy << "\n";
    out << "spin = " << config_.initialExcitation.spin << "\n";
    out << "parity = " << config_.initialExcitation.parity << "\n\n";
    
    out.close();
    std::cout << "Saved parameters to: " << filename << std::endl;
}

void OutputManager::finalize() {
    if (!outputFile_) return;
    
    std::cout << "Finalizing output..." << std::endl;
    
    // Save parameters
    saveParameters();
    
    // Write histograms and tree
    outputFile_->cd();
    
    for (auto& [name, hist] : histograms1D_) {
        hist->Write();
    }
    
    for (auto& [name, hist] : histograms2D_) {
        hist->Write();
    }
    
    for (auto& [name, hist] : levelSpectraHistograms_) {
        hist->Write();
    }
    
    if (cascadeTree_) {
        cascadeTree_->Write();
    }
    
    // Print summary
    std::cout << "Output summary:" << std::endl;
    std::cout << "  1D histograms: " << histograms1D_.size() << std::endl;
    std::cout << "  2D histograms: " << histograms2D_.size() << std::endl;
    std::cout << "  Level spectra histograms: " << levelSpectraHistograms_.size() << std::endl;
    if (cascadeTree_) {
        std::cout << "  Tree entries: " << cascadeTree_->GetEntries() << std::endl;
    }
    
    outputFile_->Close();
    delete outputFile_;
    outputFile_ = nullptr;
    
    std::cout << "Output file closed." << std::endl;
}

} // namespace rainier
