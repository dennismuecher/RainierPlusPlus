#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <iomanip>

// Gaussian smoothing (same as before)
TH1D* SmoothGaussian1D(const TH1D* hin, double sigmaBins)
{
    int n = hin->GetNbinsX();
    TH1D* hout = (TH1D*)hin->Clone("hSmooth_tmp"); // temporary name
    hout->Reset();

    int halfWidth = std::max(1, int(3.0 * sigmaBins));
    std::vector<double> kernel(2*halfWidth + 1);

    double norm = 0.0;
    for (int k = -halfWidth; k <= halfWidth; ++k) {
        double w = TMath::Gaus((double)k, 0.0, sigmaBins, true);
        kernel[k + halfWidth] = w;
        norm += w;
    }
    for (double &w : kernel) w /= norm;

    for (int i = 1; i <= n; ++i) {
        double sum = 0.0;
        for (int k = -halfWidth; k <= halfWidth; ++k) {
            int bin = i + k;
            if (bin < 1 || bin > n) continue;
            sum += hin->GetBinContent(bin) * kernel[k + halfWidth];
        }
        hout->SetBinContent(i, sum);
    }

    // ensure histogram will not be owned by any TFile by default
    hout->SetDirectory(0);
    return hout;
}

// Compute normalized fluctuation ACFs (C_norm and C_slide)
void ComputeFluctuationACFs(const TH1D* h, int maxLag, TH1D*& outCnorm, TH1D*& outCslide)
{
    int n = h->GetNbinsX();
    std::vector<double> delta(n);
    int usedBins = 0;
    for (int i = 0; i < n; ++i) {
        double val = h->GetBinContent(i+1);
        delta[i] = val - 1.0;
        ++usedBins;
    }

    double S0 = 0.0;
    for (int i = 0; i < n; ++i) S0 += delta[i]*delta[i];
    double var = (usedBins>0) ? S0 / double(usedBins) : 0.0;

    outCnorm = new TH1D("hAutocorrNorm", "Normalized fluct. ACF;Lag (bins);C_{norm}(lag)", maxLag, 0, maxLag);
    outCslide = new TH1D("hAutocorrSlide", "Slide-style ACF;Lag (bins);C_{slide}(lag)", maxLag, 0, maxLag);

    for (int lag = 0; lag < maxLag; ++lag) {
        double sum = 0.0;
        for (int i = 0; i < n - lag; ++i) {
            sum += delta[i] * delta[i + lag];
        }
        double Cnorm = (S0>0.0) ? (sum / S0) : 0.0;
        double Cslide = 1.0 + var * Cnorm;
        outCnorm->SetBinContent(lag+1, Cnorm);
        outCslide->SetBinContent(lag+1, Cslide);
    }

    // detach from any file
    outCnorm->SetDirectory(0);
    outCslide->SetDirectory(0);
}

// Main fixed macro
void fluctuationStep1_withACF_andD(
        TString infile = "input.root",
        TString outfile = "stationary_output_withACF_fixed.root",
        TString histname = "hLevelSpectrum_J5_real0",
        double smoothSigmaBins = 2.0,
        int maxLag = 200,
        double fwhm_MeV = -1.0
    )
{
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input file " << infile << std::endl;
        return;
    }
    TH1D* hOrig_file = dynamic_cast<TH1D*>(fin->Get(histname));
    if (!hOrig_file) {
        std::cerr << "ERROR: histogram " << histname << " not found in " << infile << std::endl;
        fin->Close();
        return;
    }

    // Clone the original histogram into memory (so it's not owned by the input file)
    TH1D* hOrig = (TH1D*)hOrig_file->Clone("hOrig_inMemory");
    hOrig->SetDirectory(0); // detach from input file
    // we can close input file soon
    fin->Close();

    int nBins = hOrig->GetNbinsX();
    double xMin = hOrig->GetXaxis()->GetXmin();
    double xMax = hOrig->GetXaxis()->GetXmax();
    double binWidth_MeV = (xMax - xMin) / double(nBins);
    std::cout << "Histogram range: " << xMin << " .. " << xMax << " MeV, bins = " << nBins
              << ", bin width = " << binWidth_MeV << " MeV (" << (binWidth_MeV*1000.) << " keV)" << std::endl;

    // smoothing
    TH1D* hSmooth = SmoothGaussian1D(hOrig, smoothSigmaBins);
    hSmooth->SetName("hSmooth"); // user-visible name

    // stationary spectrum
    TH1D* hStat = (TH1D*)hOrig->Clone("hStationary");
    hStat->Reset();
    for (int i = 1; i <= nBins; ++i) {
        double s = hSmooth->GetBinContent(i);
        double val = hOrig->GetBinContent(i);
        if (s > 0.0) hStat->SetBinContent(i, val / s);
        else hStat->SetBinContent(i, 0.0);
    }
    hStat->SetDirectory(0); // detach from any file

    // compute ACFs
    TH1D *hACFnorm = nullptr, *hACFslide = nullptr;
    ComputeFluctuationACFs(hStat, maxLag, hACFnorm, hACFslide);

    // compute S0 and alpha
    double S0 = 0.0;
    for (int i = 1; i <= nBins; ++i) {
        double d = hStat->GetBinContent(i) - 1.0;
        S0 += d*d;
    }
    double alpha = (nBins>0) ? S0 / double(nBins) : 0.0;

    // compute <D> if requested
    double meanD_MeV = -1.0;
    if (fwhm_MeV > 0.0) {
        double sigma_MeV = fwhm_MeV / (2.0 * TMath::Sqrt(2.0 * TMath::Log(2.0)));
        if (alpha > 0.0) meanD_MeV = (2.0 * sigma_MeV * TMath::Sqrt(TMath::Pi())) / alpha;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "FWHM = " << fwhm_MeV << " MeV -> sigma = " << sigma_MeV << " MeV" << std::endl;
    } else {
        std::cout << "FWHM not provided (fwhm_MeV <= 0). Skipping <D> computation." << std::endl;
    }

    std::cout << "---- Results summary ----" << std::endl;
    std::cout << "Used bins (N) = " << nBins << std::endl;
    std::cout << "S0 = " << S0 << std::endl;
    std::cout << "variance var(delta) = alpha = " << alpha << std::endl;
    std::cout << "C_slide(0) = 1 + alpha = " << (1.0 + alpha) << std::endl;
    if (meanD_MeV > 0.0) {
        std::cout << "Extracted mean level spacing <D> = " << meanD_MeV << " MeV"
                  << " = " << (meanD_MeV*1000.0) << " keV" << std::endl;
    }

    // write output file *after* detaching histograms from any input file
    TFile *fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot open output file " << outfile << " for writing." << std::endl;
    } else {
        // Ensure histograms have unique names for file
        hOrig->SetName("hOriginal_inMemory");
        hOrig->Write();
        hSmooth->Write();
        hStat->Write();
        hACFnorm->Write();
        hACFslide->Write();

        TString info;
        info += TString::Format("Nbins=%d;", nBins);
        info += TString::Format("binWidth_MeV=%g;", binWidth_MeV);
        info += TString::Format("smoothSigmaBins=%g;", smoothSigmaBins);
        info += TString::Format("S0=%g;", S0);
        info += TString::Format("alpha=%g;", alpha);
        if (fwhm_MeV > 0.0) {
            info += TString::Format("FWHM_MeV=%g;", fwhm_MeV);
            info += TString::Format("D_MeV=%g;", meanD_MeV);
        } else {
            info += "FWHM_MeV=not_provided;";
        }
        TNamed meta("fluctuation_results", info.Data());
        meta.Write();

        fout->Close();
        std::cout << "Wrote histograms and metadata to " << outfile << std::endl;
    }

    // draw AFTER files are closed and histograms detached from files
    TCanvas *c1 = new TCanvas("c1","Spectra",900,600);
    c1->Divide(2,1);
    c1->cd(1);
    hOrig->SetTitle("Original spectrum");
    hOrig->Draw();
    c1->cd(2);
    hStat->SetTitle("Stationary spectrum (orig / smoothed)");
    hStat->Draw();
    c1->Update();

    TCanvas *c2 = new TCanvas("c2","Autocorrelation",900,450);
    c2->Divide(2,1);
    c2->cd(1);
    hACFnorm->SetTitle("Normalized ACF (C_{norm}), C_norm(0)=1");
    hACFnorm->Draw("HIST");
    c2->cd(2);
    hACFslide->SetTitle("Slide-style ACF (C_{slide} = 1 + alpha * C_norm)");
    hACFslide->Draw("HIST");
    c2->Update();

    // keep objects in memory (they are detached from files)
}
