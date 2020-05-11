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

double ASK(double x, double y, int A1, int A2, int N, double Tb, double fi) {
	int f = N / Tb;
	if (y == 0) {
		return (A1 * sin(2 * M_PI * f * x + fi));
	}
	else {
		return (A2 * sin(2 * M_PI * f * x + fi));
	}
}

double FSK(double x, double y, int A, int N, double Tb, double fi) {

	if (y == 0) {
		return (A * sin(2 * M_PI * (N / Tb) * x + fi));
	}
	else {
		return (A * sin(2 * M_PI * (N * 10 / Tb) * x + fi));
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
double sygnalNosny3(double x, int A, int N, double Tb, double fi, int y) {

	if (y == 0) {
		return (A * sin(2 * M_PI * (N / Tb) * x + fi));
	}
	else {
		return (A * sin(2 * M_PI * (N * 10 / Tb) * x + fi));
	}
}

double calka(double* y, int start, int end, double h) {
	double calka = 0;
	for (int i = start; i < end; i++)
	{
		calka += y[i];
	}
	return h * calka;
}
double* pt(double* y, int length) {
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
	string byte = tobyte("ab", false);
	byte = "1010110010001100";
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
	double* tableYTTL = new double[length];
	double* tableYM = new double[length];
	double* tableYN = new double[length];
	double* tableYB = new double[length];
	double xtmp = 0;
	//ttl
	for (int i = 0; i < byte.length(); i++)
	{
		int a = byte[i] - '0';
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableYTTL[j] = a;
			tableX[j] = j * 1. / czestotlowosc;
		}
	}
	
	GenerateData(tableX, tableYTTL, length, "TTL");
	
	//clk
	for (int i = 0; i < byte.length(); i++)
	{
		int a = 1;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{

			if (j > (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2) {
				a = 0;
			}
			tableY[j] = a;
		}
	

	}
	GenerateData(tableX, tableY, length, "CLK");
	//manchester
	int a = 0;
	int b = byte[0] - '0';
	for (int i = 0; i < byte.length(); i++)
	{
		b = byte[i] - '0';
		if (i == 0) {
			a = 0;
		}
		else {
			if (b == 1) {
				a = 1;
			}
			else {
				a = -1;
			}
		}
		for (int j = i * Tb * czestotlowosc; j < (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j++)
		{
			tableYM[j] = a;
		}

		if (b == 1) {
			a = -1;
		}
		else {
			a = 1;
		}
	
		for (int j = (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableYM[j] = a;
		}
	}
	GenerateData(tableX, tableYM, length, "manchester");
	//nrzi
	a = 0;
	int flag = a;

	for (int i = 0; i < byte.length(); i++)
	{
		
		b = byte[i] - '0';
		if (i == 0) {
			a = 0;
		}
		for (int j = i * Tb * czestotlowosc; j < (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j++)
		{
			tableYN[j] = a;
		}
		if (i == 0 && b ==1) {
			a = 1;
		}
		if (b == 1) {
			a *= -1;
		}

		for (int j = (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableYN[j] = a;
		}
	}
	GenerateData(tableX, tableYN, length, "nrzi");

	//bami
	a = 0;
	flag = a;
	int flag1 = -1;
	for (int i = 0; i < byte.length(); i++)
	{

		b = byte[i] - '0';
		if (i == 0) {
			a = 0;
		}
		else {
			if (b == 1) {
				flag = flag1;
				flag1 *= -1;

			}
			else {
				flag = 0;
			}
			a = flag;
		}
		for (int j = i * Tb * czestotlowosc; j < (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j++)
		{
			tableYB[j] = a;
		}
		for (int j = (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableYB[j] = a;
		}
	}
	GenerateData(tableX, tableYB, length, "bami");

	//dekodowanie ttl
	int* ttl = new int[byte.length()];

	for (int i = 0; i < byte.length(); i++)
	{
		a = 0;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			a = tableYTTL[j-1];
		}
		ttl[i] = a;
	}
	cout << "TTL: ";
	for (int i = 0; i < byte.length(); i++)
	{
		cout << ttl[i];

	}
	cout << endl;

	//dekodowanie manchester
	int* manchester = new int[byte.length()];
	int o = Tb * czestotlowosc * 0.25;
	for (int i = 0; i < byte.length(); i++)
	{

		for (int j = i * Tb * czestotlowosc + o ; j < (i * Tb * czestotlowosc+o  + (i + 1) * Tb * czestotlowosc+o ) / 2; j++)
		{
			a = tableYM[j];
		}
		if (a == 1) {
			manchester[i] = 0;
		}
		else {
			manchester[i] = 1;
		}
	}
	cout << "Manchester: ";

	for (int i = 0; i < byte.length(); i++)
	{
		cout << manchester[i];
	}
	cout << endl;
	//dekodowanie nrzi
	int* nrzi = new int[byte.length()];

	for (int i = 0; i < byte.length(); i++)
	{

		for (int j = i * Tb * czestotlowosc; j < (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j++)
		{
			a = tableYN[j-1];
		}
		if (a >= 0) {
			nrzi[i] = 1;
		}
		else {
			nrzi[i] = 0;
		}
		for (int j = (i * Tb * czestotlowosc + (i + 1) * Tb * czestotlowosc) / 2; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			a = tableYN[j];
		}
	}
	cout << "NRZI: ";

	for (int i = 0; i < byte.length()-1; i++)
	{
		nrzi[i] = nrzi[i] ^ nrzi[i+1];
	}

	for (int i = 0; i < byte.length(); i++)
	{
		cout << nrzi[i];
	}
	cout << endl;

	//dekodowanie bami
	int* bami = new int[byte.length()];

	for (int i = 0; i < byte.length(); i++)
	{
		a = 0;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			a = tableYB[j - 1];
		}
		if (a != 0) {
			bami[i] = 1;
		}
		else {
			bami[i] = 0;
		}
	}
	cout << "Bami: ";
	for (int i = 0; i < byte.length(); i++)
	{
		cout << bami[i];

	}
}
