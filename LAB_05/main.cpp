#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <string>

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

double m(double x, int N) {
	double suma = 0;
	for (int i = 1; i <= N; i++)
	{
		suma += (cos(12 * x * pow(i, 2)) + cos(16 * x * i)) / (pow(i, 2));
	}
	return suma;
}


double za(double ka, double t,double fn) {
	return ((ka * m(t, 2) + 1) * cos(2 * M_PI * fn * t));
}
double zp(double kp, double t, double fn) {
	return  cos(2 * M_PI * fn * t + kp * m(t,2));
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

int main() {
	int fn = 800;
	double ka_a = 0.5;
	double kp_a = 1.5;
	double ka_b = 7;
	double kp_b = 2;
	double ka_c = 68;
	double kp_c = 77;
	double deltaT = 1 / 8000.;
	int length = 8000;
	double* tableXm = new double[length];
	double* tableYm = new double[length];
	double xtmp = 0;

	for (int i = 0; i < length; i++)
	{
		tableXm[i] = xtmp;
		tableYm[i] = m(xtmp, 2);
		xtmp += deltaT;
	}
	GenerateData(tableXm, tableYm, length, "m");

	complex<double>* Xm = new complex<double>[length];

	double* Mprimm = new double[length];
	double* Mm = new double[length];
	double* Fkm = new double[length];

	Xm = dft(tableYm, length);
	Mm = Mk(Xm, length);
	Mprimm = Mp(Mm, length);
	Fkm = FK(length, deltaT);
	
	GenerateData(Fkm, Mprimm, length, "dtfm");

	double* tableX = new double[length];
	double* tableY = new double[length];
	double* tableY2 = new double[length];
	 xtmp = 0;
	for (int i = 0; i < length; i++)
	{
		tableX[i] = xtmp;
		tableY[i] = za(ka_a,xtmp,fn);
		xtmp += deltaT;
	}
	GenerateData(tableX, tableY, length, "za_a");

	complex<double>* X = new complex<double>[length];

	double* Mprim = new double[length];
	double* M = new double[length];
	double* Fk = new double[length];

	X = dft(tableY, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, deltaT);
	max(length, Mprim, Fk);
	GenerateData(Fk, Mprim, length, "dtfa_za");

	xtmp = 0;

	for (int i = 0; i < length; i++)
	{
		tableY2[i] = zp(kp_a, xtmp, fn);
		xtmp += deltaT;
	}

	GenerateData(tableX, tableY2, length, "zp_a");

	X = dft(tableY2, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, deltaT);
	max(length, Mprim, Fk);
	GenerateData(Fk, Mprim, length, "dtfa_zp");

	xtmp = 0;
	double* tableYb = new double[length];
	double* tableY2b = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYb[i] = za(ka_b, xtmp, fn);
		xtmp += deltaT;
	}
	GenerateData(tableX, tableYb, length, "za_b");

	X = dft(tableYb, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, deltaT);
	max(length, Mprim, Fk);
	GenerateData(Fk, Mprim, length, "dtfb_za");
	xtmp = 0;

	for (int i = 0; i < length; i++)
	{
		tableY2b[i] = zp(kp_b, xtmp, fn);
		xtmp += deltaT;
	}

	GenerateData(tableX, tableY2b, length, "zp_b");
	X = dft(tableY2b, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, deltaT);
	max(length, Mprim, Fk);
	GenerateData(Fk, Mprim, length, "dtfb_zp");

	xtmp = 0;
	double* tableYc = new double[length];
	double* tableY2c = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYc[i] = za(ka_c, xtmp, fn);
		xtmp += deltaT;
	}
	GenerateData(tableX, tableYc, length, "za_c");
	X = dft(tableYc, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, deltaT);
	max(length, Mprim, Fk);
	GenerateData(Fk, Mprim, length, "dtfc_za");

	xtmp = 0;

	for (int i = 0; i < length; i++)
	{
		tableY2c[i] = zp(kp_c, xtmp, fn);
		xtmp += deltaT;
	}

	GenerateData(tableX, tableY2c, length, "zp_c");
	X = dft(tableY2c, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, deltaT);
	max(length, Mprim, Fk);
	GenerateData(Fk, Mprim, length, "dtfc_zp");
	/*
	Szerokosc pasma:1HZ
	Szerokosc pasma : 10HZ
	Szerokosc pasma : 7HZ
	Szerokosc pasma : 13HZ
	Szerokosc pasma : 7HZ
	Szerokosc pasma : 412HZ
	*/
}