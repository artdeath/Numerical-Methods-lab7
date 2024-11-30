#ifndef FOURIER_H
#define FOURIER_H

#include <complex>
#include <iostream>
#include <vector>


class Fourier {

private:
	
	int n;

	double a;

	double b;

	double w;

	double f;

	static const double pi;

	std::vector<std::complex<double>> basis;

	std::vector<std::complex<double>> basis_f;

public:

	Fourier();

	Fourier(double A, double B, double W, double F, int N);

	Fourier(const Fourier& Object);

	void DFT(std::vector<std::complex<double>>& Vector);

	void IDFT(std::vector<std::complex<double>>& Vector);

	void FFT(std::vector<std::complex<double>>& Vector);

	void IFFT(std::vector<std::complex<double>>& Vector);

	double phase(std::complex<double> Value) const;

	double amplitude(std::complex<double> Value) const;

	double generate(int j) const;

	double generate_noise(int j) const;

	~Fourier();



};


#endif#pragma once
