#include <complex>

#include "lab7.h"
#include "lab6.h" //FFT в conv()

Wavelet::Wavelet(int n, int s, int T) : N{ n }, Stages{ s },  PI { 3.1415926535897932 } {

	switch (T) {
	case 1:
	{

		std::complex<double> temp;

		u.push_back(1.0 / sqrt(2.0));

		v.push_back(1.0 / sqrt(2.0));

		for (int i = 1; i < N; i++) {

			temp = std::complex<double>(sqrt(2.0) / N * cos(PI * i / N) * sin(PI * i / 2.0) / sin(PI * i / N), -sqrt(2.0) / N * sin(PI * i / N) * sin(PI * i / 2.0) / sin(PI * i / N));

			u.push_back(temp);

			v.push_back(std::pow(-1, i) * temp);

		}

	}
	break;
	case 2:
	{

		double a = 1. - std::sqrt(10.);

		double b = 1. + std::sqrt(10.);

		double c = std::sqrt(5. + 2. * std::sqrt(10.));

		double d = std::sqrt(2.) / 32.;

		u.push_back((b + c) * d); u.push_back((2. * a + 3. * b + 3. * c) * d); u.push_back((6. * a + 4. * b + 2. * c) * d);
		
		u.push_back((6. * a + 4. * b - 2. * c) * d); u.push_back((2. * a + 3. * b - 3. * c) * d); u.push_back((b - c) * d);

		for (int i = 0; i < N - 6; i++) {

			u.push_back(0);

		}

		v.push_back(-u.at(1)); v.push_back(u.at(0));

		for (int i = 0; i < N - 6; i++) {

			v.push_back(0);

		}

		v.push_back(-u.at(5)); v.push_back(u.at(4));  v.push_back(-u.at(3)); v.push_back(u.at(2));

	}
	break;
	default:
	{
		throw std::runtime_error("Unknown system\n");
	}
	break;
	}



	std::vector<std::vector<std::complex<double>>> u_stages;

	std::vector<std::vector<std::complex<double>>> v_stages;

	u_stages.push_back(u); v_stages.push_back(v);

	for (int i = 1; i < Stages; i++) {

		int filter_n = N / static_cast<int>(std::pow(2, i));

		u_stages.push_back(std::vector<std::complex<double>>());

		v_stages.push_back(std::vector<std::complex<double>>());

		for (int n = 0; n < filter_n; n++) {

			int div = static_cast<int>(std::pow(2, i));

			u_stages.at(i).push_back(0); v_stages.at(i).push_back(0);

			for (int k = 0; k < div; k++) {

				u_stages.at(i).at(n) += u_stages.at(0).at(n + k * N / div);

				v_stages.at(i).at(n) += v_stages.at(0).at(n + k * N / div);

			}
		}
	}

	std::vector<std::vector<std::complex<double>>> f;

	std::vector<std::vector<std::complex<double>>> g;

	std::vector<std::complex<double>> ug; 		//v -> f u -> g

	std::vector<std::complex<double>> vf;

	f.push_back(v_stages.at(0)); g.push_back(u_stages.at(0));

	for (int i = 1; i < Stages; i++) {

		ug = U(i, u_stages.at(i));

		vf = U(i, v_stages.at(i));

		f.push_back(convolution(g.at(i - 1), vf));

		g.push_back(convolution(g.at(i - 1), ug));

	}

	int current = u.size() / static_cast<int>(std::pow(2, Stages));

	for (int i = 0; i < current; i++) {

		Psi.push_back(shift(f.at(Stages - 1), static_cast<int>(std::pow(2, Stages)) * i));

		Phi.push_back(shift(g.at(Stages - 1), static_cast<int>(std::pow(2, Stages))* i));

	}

	ug.clear(); vf.clear();

	for (int i = 0; i < Stages; i++) {

		f.at(i).clear(); g.at(i).clear(); 

		u_stages.at(i).clear(); v_stages.at(i).clear();

	}

	f.clear(); g.clear(); u_stages.clear(); v_stages.clear();

}

std::vector<std::complex<double>> Wavelet::shift(std::vector<std::complex<double>> vec, int n) {

	std::vector<std::complex<double>> result;

	int j = 0; int k = vec.size();

	if (n >= k) {

		n = n % k;

		if (n == 0) {

			return vec;

		}

	}

	result.resize(k);

	for (int i = 0; i < k; i++) {

		result.at(i) = ((i - n) < 0) ? vec.at(i - n + k) : vec.at(i - n);

	}

	return result;

}




std::vector<std::complex<double>> Wavelet::U(int l, std::vector<std::complex<double>> vec) { 

	std::vector<std::complex<double>> result;
	
	int temp = static_cast<int>(std::pow(2, l));
	//если делится на 2^l //U(z)(n)
	int n = vec.size() * temp;

	for (int i = 0; i < n; i++) {
		
		((i % temp) == 0) ? result.push_back(vec.at(i / temp)) : result.push_back(0);

	}

	return result;

}


std::vector<std::complex<double>> Wavelet::convolution(std::vector<std::complex<double>> first, std::vector<std::complex<double>> second) {

	Fourier Object;

	std::vector<std::complex<double>> temp1 = first;

	std::vector<std::complex<double>> temp2 = second;

	Object.FFT(temp1);
	
	Object.FFT(temp2);

	int n = temp1.size();

	for (int i = 0; i < n; i++) {

		temp2.at(i) *= temp1.at(i);

	}

	Object.IFFT(temp2);

	return temp2;

}

std::complex<double> Wavelet::scalar(std::vector<std::complex<double>> vec, std::vector<std::complex<double>> vec2) {

	std::complex<double> result = (0., 0.);

	int n = vec.size();

	for (int i = 0; i < n; i++) {

		result += vec.at(i) * std::complex<double>(vec2.at(i).real(), -vec2.at(i).imag());

	}


	return result;

}

void Wavelet::analysis(std::vector<std::complex<double>> z, int stage, std::vector<std::complex<double>>& PHI, std::vector<std::complex<double>>& PSI) {
//вариант с psi/phi вместо shift(f/g, k)
	int basis_n = u.size() / static_cast<int>(std::pow(2, stage));

	for (int i = 0; i < basis_n; i++) {

		PSI.push_back(scalar(z, Psi.at(i)));

		PHI.push_back(scalar(z, Phi.at(i)));
		
	}

	return;

}

void Wavelet::synthesis_first(std::vector<std::complex<double>>& z, std::vector<std::complex<double>>& PHI, std::vector<std::complex<double>>& PSI) {

	std::vector<std::complex<double>> P;

	std::vector<std::complex<double>> Q;

	int basis_n = Psi.size();

	int z_n = u.size();

	//P-j+1 = P-j + Q-j
	//sum <z, phi>phi + sum <z, psi>psi

	for (int i = 0; i < z_n; i++) {

		std::complex<double> p = std::complex<double>(0, 0);

		std::complex<double> q = std::complex<double>(0, 0);

		for (int j = 0; j < basis_n; j++) {

			p += PHI.at(j) * Phi.at(j).at(i);

			q += PSI.at(j) * Psi.at(j).at(i);

		}

		P.push_back(p); Q.push_back(q);

		z.push_back(P.at(i) + Q.at(i));

	}

}

void Wavelet::synthesis_second(std::vector<std::complex<double>>& z, std::vector<std::complex<double>>& P_prev, std::vector<std::complex<double>>& PSI) {

	//т.к. P-j+1 = P-j + Q-j и P-j - это z_rec
	//чтобы вместо phi_k и phi брала z_rec с предыдущего этапа

	std::vector<std::complex<double>> P = P_prev;

	std::vector<std::complex<double>> Q;

	int basis_n = Psi.size();

	int z_n = u.size();

	for (int i = 0; i < z_n; i++) {

		std::complex<double> q = std::complex<double>(0, 0);

		for (int j = 0; j < basis_n; j++) {

			q += PSI.at(j) * Psi.at(j).at(i);

		}

		Q.push_back(q);

		z.push_back(P.at(i) + Q.at(i));

	}

}

Wavelet::~Wavelet() {

	/*for (auto& i : Psi) {

		i.clear();

	}

	for (auto& i : Phi) {

		i.clear();

	}*/

}