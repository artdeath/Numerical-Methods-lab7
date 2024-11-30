#include "lab7.h"
#include <iostream>
#include <fstream>

int main() {

	std::vector<std::complex<double>> z, phik, psik;

	for (int i = 0; i < 128; i++) {

		z.push_back(0.);

	}
	for (int i = 128; i < 256; i++) {

		z.push_back(std::sin(std::pow(std::abs(i - 128), 1.64) / 128));

	}
	for (int i = 256; i < 384; i++) {

		z.push_back(0.);

	}
	for (int i = 384; i < 448; i++) {

		z.push_back(std::sin(std::pow(std::abs(i - 128), 2.46) / 128));

	}
	for (int i = 448; i < 512; i++) {

		z.push_back(0.);

	}


	Wavelet s1(512, 4, 2);

	std::ofstream ofs("db_z_s4.txt"); std::ofstream ofs2("db_k_s4.txt");

	s1.analysis(z, 4, phik, psik);

	std::vector<std::complex<double>> z_rec;

	s1.synthesis_first(z_rec, phik, psik);

	for (int i = 0; i < z.size(); i++) {
		ofs << z_rec[i].real() << std::endl;
	}
	for (int i = 0; i < phik.size(); i++) {
		ofs2 << psik[i].real() << ' ' << phik[i].real() << std::endl;
	}

	ofs.close(); ofs.clear(); ofs2.close(); ofs2.clear();

	ofs.open("db_z_s3.txt"); ofs2.open("db_k_s3.txt");

	Wavelet s2(512, 3, 2);

	std::vector<std::complex<double>> phik2, psik2;

	s2.analysis(z, 3, phik2, psik2);

	std::vector<std::complex<double>> z_rec2;

	s2.synthesis_second(z_rec2, z_rec, psik2); //z_rec âìåñòî phik2

	for (int i = 0; i < z.size(); i++) {
		ofs << z_rec2[i].real() << std::endl;
	}
	for (int i = 0; i < phik2.size(); i++) {
		ofs2 << psik2[i].real() << ' ' << phik2[i].real() << std::endl;
	}


	ofs.close(); ofs.clear(); ofs2.close(); ofs2.clear();

	ofs.open("db_z_s2.txt"); ofs2.open("db_k_s2.txt");

	Wavelet s3(512, 2, 2);

	std::vector<std::complex<double>> phik3, psik3;

	s3.analysis(z, 2, phik3, psik3);

	std::vector<std::complex<double>> z_rec3;

	s3.synthesis_second(z_rec3, z_rec2, psik3);

	for (int i = 0; i < z.size(); i++) {
		ofs <<  z_rec3[i].real() << std::endl;
	}
	for (int i = 0; i < phik3.size(); i++) {
		ofs2 << psik3[i].real() << ' ' << phik3[i].real() << std::endl;
	}

	ofs.close(); ofs.clear(); ofs2.close(); ofs2.clear();

	ofs.open("db_z_s1.txt"); ofs2.open("db_k_s1.txt");

	Wavelet s4(512, 1, 2);

	std::vector<std::complex<double>> phik4, psik4;

	s4.analysis(z, 1, phik4, psik4);

	std::vector<std::complex<double>> z_rec4;

	s4.synthesis_second(z_rec4, z_rec3, psik4);

	for (int i = 0; i < z.size(); i++) {
		ofs << z_rec4[i].real() << std::endl;
	}
	for (int i = 0; i < phik4.size(); i++) { 
		ofs2 << psik4[i].real() << ' ' << phik4[i].real() << std::endl;
	}

	ofs.close(); ofs.clear(); ofs2.close(); ofs2.clear();

	double Error = 0.;

	for (int i = 0; i < 512; i++) {

		Error += (z.at(i).real() - z_rec4.at(i).real()) * (z.at(i).real() - z_rec4.at(i).real());

	}

	Error /= 512;

	Error = std::sqrt(Error);

	std::cout << Error;


	return 0;

}
