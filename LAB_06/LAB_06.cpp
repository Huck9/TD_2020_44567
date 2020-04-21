#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <bitset>


using namespace std;
int A = 7;
int B = 6;
int C = 5;

template<typename T>
void GenerateData(T* Xtable, T* Ytable, int length, string name)
{
	ofstream outfile;
	outfile.open(name + ".plt");

	for (int i = 0; i < length; i++)
	{
		outfile << Xtable[i] << "\t" << Ytable[i] << endl;
	}
	outfile.close();
}

complex<double>* dft(double* tab, const int N) {

	complex<double>* dtf = new complex<double>[N];
	for (int n = 0; n < N; n++) {
		for (int k = 0; k < N; k++) {
			double m = -2 * M_PI * k * n / N;
			dtf[k] += polar(tab[n], m);
		}
	}
	return dtf;
}

void max(int length, double* Mprim, double* Fk) {
	double max = Mprim[0];
	for (int i = 0; i < length; i++) {
		if (max < Mprim[i]) {
			max = Mprim[i];
		}
	}
	max -= 3;
	int indexMax;
	int indexMin;
	for (int i = 0; i < length / 2; i++)
	{
		if (Mprim[i] >= max)
		{
			indexMin = i;
			break;
		}
	}
	for (int i = length / 2; i >= 0; i--)
	{
		if (Mprim[i] >= max)
		{
			indexMax = i;
			break;
		}
	}
	cout << "Szerokosc pasma:" << Fk[indexMax + 1] - Fk[indexMin] << "HZ" << endl;
}

double* Mk(complex<double>* X, const int N) {
	double* M = new double[N];
	for (int i = 0; i < N; i++)
	{
		M[i] = sqrt(pow(X[i].real(), 2) + pow(X[i].imag(), 2));
	}
	return M;
}

double* Mp(double* M, const int N) {
	double* Mprim = new double[N];
	for (int i = 0; i < N; i++)
	{
		Mprim[i] = 10 * log10(M[i]);
		if (Mprim[i] < 0) {
			Mprim[i] = 0;
		}
	}
	return Mprim;
}

double* FK(const int N, double deltaT) {
	double* Fk = new double[N];
	for (int i = 0; i < N; i++)
	{
		Fk[i] = i * (1 / deltaT) / N;
	}
	return Fk;
}
string reverse(string const& s)
{
	string rev;
	for (int i = s.size() - 1; i >= 0; i--) {
		rev = rev.append(1, s[i]);
	};

	return rev;
}

string tobyte(string a,int sw ) {
	int length = a.length();
	string bit;
	for (int i = 0; i < length; i++)
	{
		bitset<8> b = (bitset<8>)(int(a[i]));
		if (sw == 1)
		{
			//little indian
			string c = b.to_string();
			
			bit += reverse(c);
		}
		else {
			//big indian 
			bit += b.to_string();
		}
	}
	return bit;
}

int main()
{
	cout << tobyte("a", 1) << endl;
	cout << tobyte("a", 2) << endl;
}
