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

double s(double t, double f, double fi) {
	return (sin(2 * M_PI * t * f + fi));
}
void squareFunction(int min, int max, double deltat) {

	int lenght = (max - min) / deltat;
	double* tableX = new double[lenght];
	double* tableY = new double[lenght];
	double xtmp = min;
	for (int i = 0; i <= lenght; i++)
	{
		tableX[i] = xtmp;
		tableY[i] = A * pow(xtmp, 2) + B * xtmp + C;
		xtmp = xtmp + deltat;
	}
	double delta = B ^ 2 - 4 * A * C;
	if (delta > 0) {
		double t1 = (B * (-1) + sqrt(delta)) / 2 * A;
		double t2 = (B * (-1) - sqrt(delta)) / 2 * A;
		cout << "Miejsca zerowe: t1=" << t1 << " t2=" << t2 << endl;
	}
	else if (delta == 0)
	{
		double t1 = (B / (2 * A)) * (-1);
		cout << "Miejsca zerowe: t1=" << t1 << endl;
	}
	else {
		cout << "Brak miejsc zerowych" << endl;
	}
}

double X(double x) {
	return(A * pow(x, 2) + B * x + C);
}

double Y(double x) {
	return(2 * pow(X(x), 2) + 12 * cos(x));
}

double Z(double x) {
	return(sin(2 * M_PI * 7 * x) * X(x) - 0.2 * log10(Y(x + M_PI)));
}

double U(double x) {
	return(sqrt(abs(Y(x) * Y(x) * Z(x)) - 1.8 * sin(0.4 * x * Z(x) * X(x))));
}

double V(double x) {
	if (x < 0.22) {
		return((1 - 7 * x) * sin((2 * M_PI * 10) / (x + 0.4)));
	}
	else if (x < 0.77)
	{
		return(0.63 * x * sin(125 * x));
	}
	else {
		return(pow(x, -0.662) + 0.77 * sin(8 * x));
	}
}

double P(double x, int N) {
	double suma = 0;
	for (int i = 1; i <= N; i++)
	{
		suma += (cos(12 * x * pow(i, 2)) + cos(16 * x * i)) / (pow(i, 2));
	}
	return suma;
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

int main() {
	
	int min = 0;
	int max = 765;
	double a = 7;
	double b = 6;
	double c = 5;
	double deltaT = 1 / 100.;
	int length = (max - min) / deltaT;
	double* tableX = new double[length];
	double* tableY = new double[length];
	double xtmp = min;

	for (int i = 0; i < length; i++)
	{
		tableX[i] = xtmp;
		tableY[i] = s(xtmp, b, (c * M_PI));
		xtmp += deltaT;
	}

	int size = pow(2, 12);
	complex<double>* X = new complex<double>[size];
	X = dft(tableY,size);
	double* M = new double[size];
	M = Mk(X, size);
	double* Fk = new double[size];

	for (int i = 0; i < size; i++)
	{
		Fk[i] = (double)i * (1 / deltaT) / size;
	}

	GenerateData(Fk, M, size, "dft.plt");

	return 0;
 }