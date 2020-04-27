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

string tobyte(string a, bool sw) {
	int length = a.length();
	string bit;
	for (int i = 0; i < length; i++)
	{
		bitset<8> b = (bitset<8>)(int(a[i]));
		if (sw == true)
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

double ASK(double x, double y, int A1, int A2, int N, double Tb,double fi) {
	int f = N / Tb;
	if (y == 0) {
		return (A1 * sin(2 * M_PI * f * x + fi));
	}
	else {
		return (A2 * sin(2 * M_PI * f * x+ fi));
	}
}

double FSK(double x, double y, int A, int N, double Tb, double fi) {

	if (y == 0) {
		return (A * sin(2 * M_PI * (N  / Tb) * x + fi));
	}
	else {
		return (A * sin(2 * M_PI * (N*10  / Tb) * x + fi));
	}
}

double PSK(double x, double y, int A, int N, double Tb) {
	int f = N / Tb;
	if (y == 0) {

		return (A * sin(2 * M_PI * f * x));
	}
	else {
		return (A * sin(2 * M_PI * f * x + M_PI));
	}
}

double sygnalNosny(double x, int A1, int N, double Tb, double fi) {
	int f = N / Tb;
	return (A1 * sin(2 * M_PI * f * x + fi));
}
double sygnalNosny2(double x, int A, int N, double Tb, double fi) {
	int f = N / Tb;
	return (A * sin(2 * M_PI * f * x + M_PI));
}
double sygnalNosny3(double x, int A, int N, double Tb, double fi,int y) {

	if (y == 0) {
		return (A * sin(2 * M_PI * (N / Tb) * x + fi));
	}
	else {
		return (A * sin(2 * M_PI * (N * 10 / Tb) * x + fi));
	}
}

double calka(double* y,int start,int end, double h) {
	double calka = 0;
	for (int i = start; i < end; i++)
	{
		calka += y[i];
	}
	return h * calka;
}
double* pt(double* y,int length) {
	for (int i = 0; i < length; i++)
	{
		if (y[i] > 0.7)
		{
			y[i] = 1;
		}
		else {
			y[i] == 0;
		}
	}
	return y;
}




double* pt2(double* y, int length) {
	for (int i = 0; i < length; i++)
	{
		if (y[i] > 0.8)
		{
			y[i] = 1;
		}
		else {
			y[i] = 0;
		}
	}
	return y;
}

int main()
{
	string byte = tobyte("haslo", true);
	double A = 4;
	double A1 = 0;
	double A2 = 4;
	double Tb = 0.1;
	double fi = 0;
	int N = 1;
	int czestotlowosc = 8000;
	int length = Tb * czestotlowosc * byte.length();
	double* tableX = new double[length];
	double* tableY = new double[length];
	double* tableYASK = new double[length];
	double* tableYFSK = new double[length];
	double* tableYPSK = new double[length];
	double xtmp = 0;
	for (int i = 0; i < byte.length(); i++)
	{
		int a = byte[i] - '0';
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableY[j] = a;
			tableX[j] = j * 1. / czestotlowosc;
		}
	}

	for (int i = 0; i < length; i++)
	{
		tableYASK[i] = ASK(tableX[i], tableY[i], A1, A2, N, Tb,fi);
	}


	for (int i = 0; i < length; i++)
	{
		tableYFSK[i] = FSK(tableX[i], tableY[i], A, N, Tb,fi);
	}


	for (int i = 0; i < length; i++)
	{
		tableYPSK[i] = PSK(tableX[i], tableY[i], A, N, Tb);
	}

	//LAB_07
	double* sygnal_nosny = new double[length];
	double* modulacja = new double[length];
	double* calkatab = new double[length];
	double* bit = new double[length];
	//for (int i = 0; i < length; i++)
	//{
	//	sygnal_nosny[i] = sygnalNosny(tableX[i], A, N, Tb, fi);
	//}

	//GenerateData(tableX, sygnal_nosny, length, "sygnalnosny");

	//for (int i = 0; i < length; i++)
	//{
	//	modulacja[i] = sygnal_nosny[i] * tableYASK[i];
	//}
	//GenerateData(tableX, modulacja, length, "mnozenie");

	//for (int i = 0; i < byte.length(); i++)
	//{
	//	/*double h1 = tableX[0];
	//	double h2 = tableX[1];
	//	double a = calka(modulacja, i * Tb * czestotlowosc, (i + 1) * Tb * czestotlowosc, h2-h1);*/
	//	double suma = 0;
	//	for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
	//	{
	//		suma += modulacja[j];
	//		calkatab[j] = suma;
	//		//calkatab[j] = a;
	//	}
	//}

	//GenerateData(tableX, calkatab, length, "calka");

	//bit = pt(calkatab, length);

	//GenerateData(tableX, bit, length, "m");

	//PSK
	for (int i = 0; i < length; i++)
	{
		sygnal_nosny[i] = sygnalNosny2(tableX[i], A, N, Tb, fi);
	}

	GenerateData(tableX, sygnal_nosny, length, "sygnalnosny");

	for (int i = 0; i < length; i++)
	{
		modulacja[i] = sygnal_nosny[i] * tableYPSK[i];
	}
	GenerateData(tableX, modulacja, length, "mnozenie");

	for (int i = 0; i < byte.length(); i++)
	{
		//double h1 = tableX[0];
		//double h2 = tableX[1];
		//double a = calka(modulacja, i * Tb * czestotlowosc, (i + 1) * Tb * czestotlowosc, h2 - h1);
		double suma = 0;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{

			suma += modulacja[j];
			calkatab[j] = suma;
			//calkatab[j] = a;
		}
	}

	GenerateData(tableX, calkatab, length, "calka");

	bit = pt2(calkatab, length);

	GenerateData(tableX, bit, length, "m");
	//FSK
	//double* modulacja2 = new double[length];
	//double* sygnalnosny2 = new double[length];
	//double* odejmowanie = new double[length];
	//for (int i = 0; i < length; i++)
	//{
	//	sygnal_nosny[i] = sygnalNosny3(tableX[i], A, N, Tb, fi,0);
	//}

	//GenerateData(tableX, sygnal_nosny, length, "sygnalnosny");

	//for (int i = 0; i < length; i++)
	//{
	//	sygnalnosny2[i] = sygnalNosny3(tableX[i], A, N, Tb, fi,1);
	//}

	//GenerateData(tableX, sygnalnosny2, length, "sygnalnosny2");

	//for (int i = 0; i < length; i++)
	//{
	//	modulacja[i] = sygnal_nosny[i] * tableYFSK[i];
	//}
	//GenerateData(tableX, modulacja, length, "mnozenie");

	//for (int i = 0; i < length; i++)
	//{
	//	modulacja2[i] = sygnalnosny2[i] * tableYFSK[i];
	//}
	//GenerateData(tableX, modulacja2, length, "mnozenie2");


	//for (int i = 0; i < length; i++)
	//{
	//	modulacja2[i] = sygnalnosny2[i] * tableYFSK[i];
	//}
	//GenerateData(tableX, modulacja2, length, "mnozenie2");

	//for (int i = 0; i < length; i++)
	//{
	//	odejmowanie[i] = modulacja2[i] - modulacja[i];
	//}


	//for (int i = 0; i < byte.length(); i++)
	//{
	//	//double h1 = tableX[0];
	//	//double h2 = tableX[1];
	//	//double a = calka(odejmowanie, i * Tb * czestotlowosc, (i + 1) * Tb * czestotlowosc, h2 - h1);
	//	double suma = 0;
	//	for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
	//	{
	//		suma += odejmowanie[j];
	//		calkatab[j] = suma;
	//		//calkatab[j] = a;
	//	}
	//}

	//GenerateData(tableX, calkatab, length, "calka");

	//bit = pt2(calkatab, length);

	//GenerateData(tableX, bit, length, "m");

}
