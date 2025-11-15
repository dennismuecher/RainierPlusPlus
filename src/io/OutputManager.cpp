// OutputManager.cpp - Complete ROOT output implementation
// Based on original RAINIER.C lines 1012-1204
#include "io/OutputManager.h"
#include "simulation/DecaySimulator.h"
#include "simulation/CascadeEvent.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <iostream>
#include <fstream>

namespace rainier {

OutputManager::OutputManager(const Config& config, int runNumber)
    : config_(config), runNumber_(runNumber), outputFile_(nullptr), cascadeTree_(nullptr) {
    
    std::string filename = "Run" + std::to_string(runNumber) + ".root";
    outputFile_ = new TFile(filename.c_str(), "RECREATE");
    
    if (!outputFile_ || outputFile_->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << filename << std::endl;
        outputFile_ = nullptr;
        return;
    }
    
    std::cout << "OutputManager created - output file: " << filename << std::endl;
    
    // Create TTree if requested
    if (config_.output.saveTree) {
        cascadeTree_ = new TTree("cascades", "Gamma Cascade Events");
        cascadeTree_->Branch("gammaEnergies", &treeGammaEnergies_);
        cascadeTree_->Branch("finalEnergies", &treeFinalEnergies_);
        cascadeTree_->Branch("times", &treeTimes_);
        cascadeTree_->Branch("initialEx", &treeInitialEx_, "initialEx/D");
        cascadeTree_->Branch("initialSpin", &treeInitialSpin_, "initialSpin/I");
        cascadeTree_->Branch("initialParity", &treeInitialParity_, "initialParity/I");
    }
}

OutputManager::~OutputManager() {
    if (outputFile_) {
        outputFile_->Close();
        delete outputFile_;
    }
}

void OutputManager::saveRealization(int realization, const DecaySimulator& simulator) {
    std::cout << "  Saving realization " << realization << std::endl;
    
    // Create histograms for this realization
    createHistograms(realization);
    
    // Fill histograms with events
    const auto& events = simulator.getEvents();
    for (const auto& event : events) {
        fillHistograms(realization, event);
        
        // Fill tree if enabled
        if (cascadeTree_) {
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
                histograms2D_[name]->Fill(step.gammaEnergy * 1000,
                                         initialEx * 1000);
            }
            isFirstStep = false;
        }
    }
}

void OutputManager::finalize() {
    std::cout << "\nFinalizing output..." << std::endl;
    
    if (outputFile_) {
        outputFile_->cd();
        
        // Write all histograms
        for (auto& pair : histograms1D_) {
            pair.second->Write();
        }
        for (auto& pair : histograms2D_) {
            pair.second->Write();
        }
        
        // Write tree
        if (cascadeTree_) {
            cascadeTree_->Write();
        }
        
        // Save parameters to text file
        saveParameters();
        
        outputFile_->Close();
        std::cout << "Output file closed: " << outputFile_->GetName() << std::endl;
    }
}

void OutputManager::saveParameters() {
    std::string filename = "Param" + std::to_string(runNumber_) + ".dat";
    std::ofstream file(filename);
    
    if (!file) {
        std::cerr << "Warning: Could not create parameter file " << filename << std::endl;
        return;
    }
    
    file << "# RAINIER 2.0 Run Parameters\n";
    file << "# Run number: " << runNumber_ << "\n\n";
    
    file << "# Nucleus\n";
    file << "Z " << config_.nucleus.Z << "\n";
    file << "A " << config_.nucleus.A << "\n";
    file << "Sn " << config_.nucleus.Sn << " MeV\n\n";
    
    file << "# Level Density\n";
    file << "LevelDensityModel BSFG\n";
    file << "a " << config_.levelDensity.a << " MeV^-1\n";
    file << "E1 " << config_.levelDensity.E1 << " MeV\n\n";
    
    file << "# Gamma Strength\n";
    file << "E1Model GenLorentz\n";
    for (const auto& res : config_.gammaStrength.e1Resonances) {
        file << "E1Resonance " << res.energy << " " << res.width << " " << res.sigma << "\n";
    }
    file << "\n";
    
    file << "# Simulation\n";
    file << "Realizations " << config_.simulation.numRealizations << "\n";
    file << "EventsPerRealization " << config_.simulation.eventsPerRealization << "\n";
    file << "InitialExcitation " << config_.initialExcitation.excitationEnergy << " MeV\n";
    
    file.close();
    std::cout << "Parameters saved to: " << filename << std::endl;
}

} // namespace rainier
