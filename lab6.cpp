#include "lab6.h"

const double Fourier::pi = 2 * std::acos(0.);

Fourier::Fourier() : a{ 7. }, b{ 0.6 }, w{ 229, }, f{ pi / 3 }, n{ 512 } {

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			basis.push_back(std::complex<double>(std::cos(2 * pi * i * j / n), -std::sin(2 * pi * i * j / n)));

		}

	}

	if (n % 2 == 0) {

		int m = n / 2;

		for (int i = 0; i < m; i++) {

			for (int j = 0; j < m; j++) {

				basis_f.push_back(std::complex<double>(std::cos(2 * pi * i * j / m), -std::sin(2 * pi * i * j / m)));

			}

		}

	}

}

Fourier::Fourier(double A, double B, double W, double F, int N) : a{ A }, b{ B }, w{ W }, f{ F }, n{ N } {

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			basis.push_back(std::complex<double>(std::cos(2 * pi * i * j / n), -std::sin(2 * pi * i * j / n)));

		}

	}

	if (n % 2 == 0) {

		int m = n / 2;

		for (int i = 0; i < m; i++) {

			for (int j = 0; j < m; j++) {

				basis_f.push_back(std::complex<double>(std::cos(2 * pi * i * j / m), -std::sin(2 * pi * i * j / m)));

			}

		}

	}

}

Fourier::Fourier(const Fourier& Object) : a{ Object.a }, b{ Object.b }, w{ Object.w }, f{ Object.f }, n{ Object.n } {

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			basis.push_back(std::complex<double>(std::cos(2 * pi * i * j / n), -std::sin(2 * pi * i * j / n)));

		}

	}

	if (n % 2 == 0) {

		int m = n / 2;

		for (int i = 0; i < m; i++) {

			for (int j = 0; j < m; j++) {

				basis_f.push_back(std::complex<double>(std::cos(2 * pi * i * j / m), -std::sin(2 * pi * i * j / m)));

			}

		}

	}

}

double Fourier::generate(int j) const {

	double z = 0.;

	z = a + b * std::cos(2 * pi * w * static_cast<double>(j) / static_cast<double>(n) + f);

	return z;

}

void Fourier::DFT(std::vector<std::complex<double>>& Vector) {

	std::vector<std::complex<double>> Temp = Vector;

	for (int i = 0; i < n; i++) {

		Vector.at(i) = 0.;

		for (int j = 0; j < n; j++) {

			Vector.at(i) += Temp.at(j) * basis.at(i * n + j);

		}

	}

	return;

}

void Fourier::IDFT(std::vector<std::complex<double>>& Vector) {

	std::vector<std::complex<double>> Temp = Vector;

	for (int i = 0; i < n; i++) {

		Vector.at(i) = 0.;

		for (int j = 0; j < n; j++) {

			Vector.at(i) += Temp.at(j) * std::conj(basis.at(i * n + j));

		}

		Vector.at(i) /= static_cast<double>(n);

	}

	return;

}

void Fourier::FFT(std::vector<std::complex<double>>& Vector) {

	if (n % 2 != 0) {

		throw std::runtime_error("number of singnals is odd\n");

	}

	std::vector<std::complex<double>> Temp = Vector;

	int m = n / 2;

	for (int i = 0; i < m; i++) {

		Vector.at(2 * i) = 0.;

		Vector.at(2 * i + 1) = 0.;

		for (int j = 0; j < m; j++) {

			Vector.at(2 * i) += Temp.at(2 * j) * basis_f.at(i * m + j);

			Vector.at(2 * i + 1) += Temp.at(2 * j + 1) * basis_f.at(i * m + j);

		}

	}

	Temp = Vector;

	for (int i = 0; i < m; i++) {

		Vector.at(i) = Temp.at(2 * i) + std::complex<double>(std::cos(2 * pi * i / n), -std::sin(2 * pi * i / n)) * Temp.at(2 * i + 1);

		Vector.at(m + i) = Temp.at(2 * i) - std::complex<double>(std::cos(2 * pi * i / n), -std::sin(2 * pi * i / n)) * Temp.at(2 * i + 1);

	}

	return;

}

void Fourier::IFFT(std::vector<std::complex<double>>& Vector) {

	if (n % 2 != 0) {

		throw std::runtime_error("number of singnals is odd\n");

	}

	int m = n / 2;

	FFT(Vector);

	std::vector<std::complex<double>> Temp = Vector;

	for (int j = 1; j < n; j++) {

		Vector.at(j) = Temp.at(n - j) / static_cast<double>(n);

	}

	Vector.at(0) = Temp.at(0) / static_cast<double>(n);

	return;

}

double Fourier::phase(std::complex<double> Value) const {

	return std::atan2(Value.imag(), Value.real());

}

double Fourier::amplitude(std::complex<double> Value) const {

	return std::abs(Value);

}

double Fourier::generate_noise(int j) const {

	double z = 0.;

	z = std::cos(2 * pi * j / n) + 0.01 * std::cos(2 * pi * w * j / n);

	return z;

}

Fourier::~Fourier() {};

