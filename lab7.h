#ifndef WAVELET_H
#define WAVELET_H

#include <vector>
#include <complex>

class Wavelet {

private:

	int N;

	int Stages;

	std::vector<std::complex<double>> u, v;

	std::vector<std::vector<std::complex<double>>> Psi, Phi;

	const double PI;


public:


	Wavelet(int n, int s, int T);

	std::vector<std::complex<double>> U(int l, std::vector<std::complex<double>> vec);

	std::vector<std::complex<double>> convolution(std::vector<std::complex<double>> first, std::vector<std::complex<double>> second);

	~Wavelet();

	void analysis(std::vector<std::complex<double>> z, int stage, std::vector<std::complex<double>>& PHI, std::vector<std::complex<double>>& PSI);

	void synthesis_first(std::vector<std::complex<double>>& z, std::vector<std::complex<double>>& PHI, std::vector<std::complex<double>>& PSI);

	void synthesis_second(std::vector<std::complex<double>>& z, std::vector<std::complex<double>>& P_prev, std::vector<std::complex<double>>& PSI);


	std::complex<double> scalar(std::vector<std::complex<double>> vec, std::vector<std::complex<double>> vec2);

	std::vector<std::complex<double>> shift(std::vector<std::complex<double>> vec, int n);

};


#endif#pragma once
