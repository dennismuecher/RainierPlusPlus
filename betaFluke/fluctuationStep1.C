#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <iostream>
#include <vector>

//---------------------------------------------------------
// Perform Gaussian smoothing manually by convolution
//---------------------------------------------------------
TH1D* SmoothGaussian1D(const TH1D* hin, double sigmaBins)
{
    int n = hin->GetNbinsX();
    TH1D* hout = (TH1D*)hin->Clone("hSmooth");
    hout->Reset();

    // Create kernel (±3 sigma)
    int halfWidth = int(3 * sigmaBins);
    std::vector<double> kernel(2*halfWidth + 1);

    double norm = 0.0;
    for (int k = -halfWidth; k <= halfWidth; k++) {
        double w = TMath::Gaus(k, 0.0, sigmaBins, true); // true = normalized
        kernel[k + halfWidth] = w;
        norm += w;
    }

    // Normalize kernel to 1
    for (double &k : kernel) k /= norm;

    // Convolution
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;
        for (int k = -halfWidth; k <= halfWidth; k++) {
            int bin = i + k;
            if (bin < 1 || bin > n) continue;
            sum += hin->GetBinContent(bin) * kernel[k + halfWidth];
        }
        hout->SetBinContent(i, sum);
    }

    return hout;
}

//---------------------------------------------------------
// Main analysis function
//---------------------------------------------------------
void fluctuationStep1(TString infile = "input.root",
                      TString outfile = "stationary_output.root",
                      TString histname = "hLevelSpectrum_J5_real0",
                      double sigmaBins = 2.0)
{
    // Open file
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening " << infile << std::endl;
        return;
    }

    TH1D *hOrig = (TH1D*)fin->Get(histname);
    if (!hOrig) {
        std::cerr << "Histogram " << histname << " not found!" << std::endl;
        fin->Close();
        return;
    }

    // Gaussian smoothing
    TH1D *hSmooth = SmoothGaussian1D(hOrig, sigmaBins);
    hSmooth->SetName("hSmooth");
    hSmooth->SetTitle("Gaussian-smoothed spectrum");

    // Create stationary spectrum
    TH1D *hStationary = (TH1D*)hOrig->Clone("hStationary");
    hStationary->Reset();
    hStationary->SetTitle("Stationary spectrum");

    for (int i = 1; i <= hOrig->GetNbinsX(); i++) {
        double s = hSmooth->GetBinContent(i);
        if (s > 0)
            hStationary->SetBinContent(i, hOrig->GetBinContent(i) / s);
        else
            hStationary->SetBinContent(i, 0);
    }

    // Write output
    TFile *fout = TFile::Open(outfile, "RECREATE");
    hOrig->Write();
    hSmooth->Write();
    hStationary->Write();
    fout->Close();

    fin->Close();
    std::cout << "Done: stationary spectrum written to " << outfile << std::endl;
}
