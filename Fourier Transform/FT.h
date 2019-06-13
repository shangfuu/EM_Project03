#pragma once
#include <iostream>
#include <complex>
#include <vector>
#define PI 3.141592653589793238460

class FT
{
private:
	
public:
	FT();
	void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void FastFourierTransform(int** InputImage, int** OutputImage,std::complex<double>**FreqRI, int h);
	std::vector<std::complex<double>> FFT(std::vector<std::complex<double>>x);

	void InverseFastFourierTransform(int** InputImage, int** OutputImage, std::complex<double>** FreqRI, int h);
	void InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void LowpassFilter(double** Real, double** Img, int** output, int h, int w);
	void HighpassFilter(double** Real, double** Img, int** filter, int h, int w);

private:

};



