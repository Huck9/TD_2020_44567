#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <bitset>
#include <cstdlib>
#include <ctime>

using namespace std;
int A = 7;
int B = 6;
int C = 5;
double alfa = 0.01;

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

double ASK(double x, double y, int A1, int A2, int N, double Tb) {
	int f = N / Tb;
	if (y == 0) {
		return (A1 * sin(2 * M_PI * f * x));
	}
	else {
		return (A2 * sin(2 * M_PI * f * x));
	}
}

double FSK(double x, double y, int A, int N, double Tb) {

	if (y == 0) {
		return (A * sin(2 * M_PI * ((N) / Tb) * x));
	}
	else {
		return (A * sin(2 * M_PI * ((N * 10) / Tb) * x));
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

double* pt(double* y, int length, int prog) {
	for (int i = 0; i < length; i++)
	{
		if (y[i] > prog)
		{
			y[i] = 1;
		}
		else {
			y[i] = 0;
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

string negative(string a, int b) {
	if (a[b] == '1') {
		a[b] = '0';
	}
	else {
		a[b] = '1';
	}
	return a;
}

string hamming(string byte) {
	string out;
	int g[7][4] = { {1,1,0,1},{1,0,1,1},{1,0,0,0},{0,1,1,1 },{0,1,0,0},{0,0,1,0},{0,0,0,1} };
	int pakiety = byte.size() / 8;
	int size = 4;


	for (int j = 0; j < pakiety; j++)
	{
		int* d1 = new int[size];
		int* d2 = new int[size];
		int d1T[4][1];
		int d2T[4][1];
		int Ks[8][1];
		int Ks2[8][1];
		string k;
		int k2[8];
		int index = j * 8;
		for (int i = 0; i < size; i++)
		{
			d1[i] = byte[i + index] - '0';
			d2[i] = byte[i + index + size] - '0';
		}
		for (int i = 0; i < 4; i++)
		{
			d1T[i][0] = d1[i];
		}
		for (int i = 0; i < 4; i++)
		{
			d2T[i][0] = d2[i];
		}
		//zerowanie macierzy k
		for (int i = 0; i < 7; i++)
		{
			Ks[i][0] = 0;
		}
		//wyznaczanie macierzy k
		for (int i = 0; i < 7; i++)
		{
			for (int j = 0; j < 1; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					Ks[i][j] += g[i][k] * d1T[k][j];

				}
			}
		}
		for (int i = 0; i < 7; i++)
		{
			Ks[i][0] = Ks[i][0] % 2;
		}
		int sum = 0;
		for (int i = 0; i < 7; i++)
		{
			sum += Ks[i][0];
			k = to_string(Ks[i][0]);
			out.append(k);
		}
		k = to_string(sum % 2);
		out.append(k);

		for (int i = 0; i < 7; i++)
		{
			Ks2[i][0] = 0;
		}
		//wyznaczanie macierzy k
		for (int i = 0; i < 7; i++)
		{
			for (int j = 0; j < 1; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					Ks2[i][j] += g[i][k] * d2T[k][j];
				}
			}
		}
		for (int i = 0; i < 7; i++)
		{
			Ks2[i][0] = Ks2[i][0] % 2;
		}
		sum = 0;
		for (int i = 0; i < 7; i++)
		{
			sum += Ks2[i][0];
			k = to_string(Ks2[i][0]);
			out.append(k);
		}
		sum = sum % 2;
		k = to_string(sum);
		out.append(k);
	}

	return out;
}

string dekodowanie(string byte) {
	string out;
	int h[3][7] = { {1,0,1,0,1,0,1},{0,1,1,0,0,1,1,}, {0,0,0,1,1,1,1} };
	int Ks[8][1];
	int s[4][1];
	int length = byte.length() / 8;
	for (int i = 0; i < length; i++)
	{
		int index = i * 8;
		for (int i = 0; i < 8; i++)
		{
			Ks[i][0] = byte[i + index] - '0';
		}
		int sum = 0;
		for (int i = 0; i < 7; i++)
		{
			sum += Ks[i][0];
		}
		sum = sum % 2;
		if (sum != Ks[7][0])
		{
			for (int i = 0; i < 4; i++)
			{
				s[i][0] = 0;
			}

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 1; j++)
				{
					for (int k = 0; k < 7; k++)
					{
						s[i][j] += h[i][k] * Ks[k][j];

					}
				}
			}
			int blad = 0;

			for (int i = 0; i < 3; i++)
			{
				s[i][0] = s[i][0] % 2;
				blad += pow(2, i) * s[i][0];
			}

			if (blad != 0) {
				Ks[blad - 1][0] = (Ks[blad - 1][0] == 1) ? 0 : 1;
			}

			sum = 0;
			for (int i = 0; i < 7; i++)
			{
				sum += Ks[i][0];
			}
			sum = sum % 2;
			if (sum != Ks[7][0]) {
				out.append("x");
				out.append("x");
				out.append("x");
				out.append("x");
			}
			else {
				out.append(to_string(Ks[2][0]));
				out.append(to_string(Ks[4][0]));
				out.append(to_string(Ks[5][0]));
				out.append(to_string(Ks[6][0]));
			}
		}
		else {
			for (int i = 0; i < 4; i++)
			{
				s[i][0] = 0;
			}

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 1; j++)
				{
					for (int k = 0; k < 7; k++)
					{
						s[i][j] += h[i][k] * Ks[k][j];

					}
				}
			}
			int blad = 0;

			for (int i = 0; i < 3; i++)
			{
				s[i][0] = s[i][0] % 2;
				blad += pow(2, i) * s[i][0];
			}

			if (blad != 0) {
				Ks[blad - 1][0] = (Ks[blad - 1][0] == 1) ? 0 : 1;
			}

			out.append(to_string(Ks[2][0]));
			out.append(to_string(Ks[4][0]));
			out.append(to_string(Ks[5][0]));
			out.append(to_string(Ks[6][0]));
		}

	}

	return out;
}

string toasc(string byte) {
	string asc;
	int length = byte.length() / 8;
	cout << endl;
	for (int i = 0; i < length; i++)
	{
		int a = 0;
		int index = i * 8;
		int k = 0;
		for (int j = 7; j >= 0; j--)
		{
			int tmp = byte[j + index] - '0';
			a += tmp * pow(2, k);
			k++;
		}
		char b = char(a);
		asc += b;
	}
	return asc;
}

string modulacjaASK(string hammingOut) {
	string out;
	double A = 4;
	double A1 = 1;
	double A2 = 4;
	double Tb = 0.1;
	double fi = 0;
	int N = 1;
	int czestotlowosc = 1000;
	int length = Tb * czestotlowosc * hammingOut.length();
	double* tableX = new double[length];
	double* tableY = new double[length];
	double* tableYASK = new double[length];
	double* sygnal_nosny = new double[length];
	double* modulacja = new double[length];
	double* calkatab = new double[length];
	double* bit = new double[length];


	double xtmp = 0;
	srand(time(0));

	/* Setup constants */
	const static int q = 15;
	const static float c1 = (1 << q) - 1;
	const static float c2 = ((int)(c1 / 3)) + 1;
	const static float c3 = 1.f / c1;

	/* random number in range 0 - 1 not including 1 */
	float random = 0.f;

	/* the white noise */
	double noise = 0.f;
	double* szum = new double[length];
	for (int i = 0; i < length; i++)
	{
		random = ((float)rand() / (float)(RAND_MAX + 1));
		noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
		szum[i] = noise;
	}

	for (int i = 0; i < hammingOut.length(); i++)
	{
		int a = hammingOut[i] - '0';
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableY[j] = a;
			tableX[j] = j * 1. / czestotlowosc;
		}
	}

	GenerateData(tableX, tableY, length, "bit");

	//ask
	for (int i = 0; i < length; i++)
	{
		tableYASK[i] = ASK(tableX[i], tableY[i], A1, A2, N, Tb);
	}
	GenerateData(tableX, tableYASK, length, "ASK");
	//widmo ask

	complex<double>* X = new complex<double>[length];
	double* Mprim = new double[length];
	double* M = new double[length];
	double* Fk = new double[length];

	/*X = dft(tableYASK, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, 1. / czestotlowosc);
	GenerateData(Fk, Mprim, length, "dtf_za");*/

	for (int i = 0; i < length; i++)
	{
		tableYASK[i] = (tableYASK[i] * alfa) + (szum[i] * (1 - alfa));
	}
	/*
		X = dft(tableYASK, length);
		M = Mk(X, length);
		Mprim = Mp(M, length);
		Fk = FK(length, 1. / czestotlowosc);
		GenerateData(Fk, Mprim, length, "dtf_za_szum");
	*/
	for (int i = 0; i < length; i++)
	{
		sygnal_nosny[i] = sygnalNosny(tableX[i], A, N, Tb, fi);
	}

	for (int i = 0; i < length; i++)
	{
		modulacja[i] = sygnal_nosny[i] * tableYASK[i];
	}

	for (int i = 0; i < hammingOut.length(); i++)
	{
		double suma = 0;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			suma += modulacja[j];
			calkatab[j] = suma;
		}
	}
	
	bit = pt(calkatab, length, 200*alfa +5);

	for (int i = 0; i < hammingOut.length(); i++)
	{
		int index = i * Tb * czestotlowosc+28;
		
		out.append(to_string((int)bit[index]));
	}
	return out;
}

string modulacjaFKS(string hammingOut) {
	string out;
	double A = 4;
	double A1 = 1;
	double A2 = 4;
	double Tb = 0.1;
	double fi = 0;
	int N = 1;
	int czestotlowosc = 1000;
	int length = Tb * czestotlowosc * hammingOut.length();
	double* tableX = new double[length];
	double* tableY = new double[length];
	double* tableYASK = new double[length];
	double* tableYFSK = new double[length];
	double* tableYPSK = new double[length];
	double* sygnal_nosny = new double[length];
	double* modulacja = new double[length];
	double* calkatab = new double[length];
	double* bit = new double[length];
	double xtmp = 0;
	srand(time(0));

	/* Setup constants */
	const static int q = 15;
	const static float c1 = (1 << q) - 1;
	const static float c2 = ((int)(c1 / 3)) + 1;
	const static float c3 = 1.f / c1;

	/* random number in range 0 - 1 not including 1 */
	float random = 0.f;

	/* the white noise */
	double noise = 0.f;
	double* szum = new double[length];
	for (int i = 0; i < length; i++)
	{
		random = ((float)rand() / (float)(RAND_MAX + 1));
		noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
		szum[i] = noise;
	}

	for (int i = 0; i < hammingOut.length(); i++)
	{
		int a = hammingOut[i] - '0';
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableY[j] = a;
			tableX[j] = j * 1. / czestotlowosc;
		}
	}

	GenerateData(tableX, tableY, length, "bit");

	
	complex<double>* X = new complex<double>[length];
	double* Mprim = new double[length];
	double* M = new double[length];
	double* Fk = new double[length];

	//FSK
	for (int i = 0; i < length; i++)
	{
		tableYFSK[i] = FSK(tableX[i], tableY[i], A, N, Tb);
	}

	GenerateData(tableX, tableYFSK, length, "FSK");

	//widmo
	/*
	X = dft(tableYFSK, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, 1. / czestotlowosc);
	GenerateData(Fk, Mprim, length, "dtf_zf");
	*/
	for (int i = 0; i < length; i++)
	{
		tableYFSK[i] = (tableYFSK[i] * alfa) + (szum[i] * (1 - alfa));
	}
	/*
	X = dft(tableYFSK, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, 1. / czestotlowosc);
	GenerateData(Fk, Mprim, length, "dtf_zf_szum");
	*/
	double* modulacja2 = new double[length];
	double* sygnalnosny2 = new double[length];
	double* odejmowanie = new double[length];
	for (int i = 0; i < length; i++)
	{
		sygnal_nosny[i] = sygnalNosny3(tableX[i], A, N, Tb, fi, 0);
	}

	for (int i = 0; i < length; i++)
	{
		sygnalnosny2[i] = sygnalNosny3(tableX[i], A, N, Tb, fi, 1);
	}

	for (int i = 0; i < length; i++)
	{
		modulacja[i] = sygnal_nosny[i] * tableYFSK[i];
	}

	for (int i = 0; i < length; i++)
	{
		modulacja2[i] = sygnalnosny2[i] * tableYFSK[i];
	}

	for (int i = 0; i < length; i++)
	{
		modulacja2[i] = sygnalnosny2[i] * tableYFSK[i];
	}

	for (int i = 0; i < length; i++)
	{
		odejmowanie[i] = modulacja2[i] - modulacja[i];
	}


	for (int i = 0; i < hammingOut.length(); i++)
	{
		double suma = 0;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			suma += odejmowanie[j];
			calkatab[j] = suma;
		}
	}
	bit = pt2(calkatab, length);

	for (int i = 0; i < hammingOut.length(); i++)
	{
		int index = i * Tb * czestotlowosc + 28;
		out.append(to_string((int)bit[index]));
	}
	return out;
}

string modulacjaPSK(string hammingOut) {
	string out;
	double A = 4;
	double A1 = 1;
	double A2 = 4;
	double Tb = 0.1;
	double fi = 0;
	int N = 1;
	int czestotlowosc = 1000;
	int length = Tb * czestotlowosc * hammingOut.length();
	double* tableX = new double[length];
	double* tableY = new double[length];
	double* tableYASK = new double[length];
	double* tableYFSK = new double[length];
	double* tableYPSK = new double[length];
	double* sygnal_nosny = new double[length];
	double* modulacja = new double[length];
	double* calkatab = new double[length];
	double* bit = new double[length];
	double xtmp = 0;

	srand(time(0));
	
	/* Setup constants */
	const static int q = 15;
	const static float c1 = (1 << q) - 1;
	const static float c2 = ((int)(c1 / 3)) + 1;
	const static float c3 = 1.f / c1;

	/* random number in range 0 - 1 not including 1 */
	float random = 0.f;

	/* the white noise */
	double noise = 0.f;
	double* szum = new double[length];
	for (int i = 0; i < length; i++)
	{
		random = ((float)rand() / (float)(RAND_MAX + 1));
		noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
		szum[i] = noise;
	}

	for (int i = 0; i < hammingOut.length(); i++)
	{
		int a = hammingOut[i] - '0';
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			tableY[j] = a;
			tableX[j] = j * 1. / czestotlowosc;
		}
	}

	GenerateData(tableX, tableY, length, "bit");

	//psk
	for (int i = 0; i < length; i++)
	{
		tableYPSK[i] = PSK(tableX[i], tableY[i], A, N, Tb);
	}

	GenerateData(tableX, tableYPSK, length, "PSK");
	/*
	X = dft(tableYPSK, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, 1. / czestotlowosc);
	GenerateData(Fk, Mprim, length, "dtf_zp");
	*/
	for (int i = 0; i < length; i++)
	{
		tableYPSK[i] = (tableYPSK[i] * alfa) + (szum[i] * (1 - alfa));
	}
	/*
	X = dft(tableYPSK, length);
	M = Mk(X, length);
	Mprim = Mp(M, length);
	Fk = FK(length, 1. / czestotlowosc);
	GenerateData(Fk, Mprim, length, "dtf_zp_szum");
	*/
	for (int i = 0; i < length; i++)
	{
		sygnal_nosny[i] = sygnalNosny2(tableX[i], A, N, Tb, fi);
	}

	for (int i = 0; i < length; i++)
	{
		modulacja[i] = sygnal_nosny[i] * tableYPSK[i];
	}

	for (int i = 0; i < hammingOut.length(); i++)
	{
		double suma = 0;
		for (int j = i * Tb * czestotlowosc; j < (i + 1) * Tb * czestotlowosc; j++)
		{
			suma += modulacja[j];
			calkatab[j] = suma;
		}
	}
	bit = pt2(calkatab, length);

	for (int i = 0; i < hammingOut.length(); i++)
	{
		int index = i * Tb * czestotlowosc + 18;
		out.append(to_string((int)bit[index]));
	}

	return out;

}

void ber(string a, string b) {
	int ber = 0;
	for (int i = 0; i < a.length(); i++)
	{
		if (a[i] != b[i]) {
			ber += 1;
		}
	}
	cout << "ber: " << ber << "/" << a.length();
}

int main()
{
	string napis = "Ala ma kota";
	cout << "Napis:" << napis << endl;
	string byte = tobyte(napis, false);
	string hammingOut;
	string demodulacjaASK;
	string demodulacjaFSK;
	string demodulacjaPSK;
	string dekodowanieham;
	cout << "Wejsciowe bity:" << endl;
	cout << byte << endl;
	cout << "Zakodowany sygnal:" << endl;
	hammingOut = hamming(byte);
	cout << hammingOut << endl;

	cout << endl;
	cout << "Sygnal zdemodulowany ASK:" << endl;
	demodulacjaASK =  modulacjaASK(hammingOut);
	cout << demodulacjaASK << endl;
	ber(hammingOut, demodulacjaASK);
	cout << endl;
	cout << "Sygnal zdekodowany:" << endl;
	dekodowanieham = dekodowanie(demodulacjaASK);
	cout << dekodowanieham << endl;
	ber(byte, dekodowanieham);
	cout << "Napis zdekodowany: " << toasc(dekodowanieham) << endl;

	cout << endl;
	cout << "Sygnal zdemodulowany FSK:" << endl;
	demodulacjaFSK = modulacjaFKS(hammingOut);
	cout << demodulacjaFSK << endl;
	ber(hammingOut, demodulacjaFSK);
	cout << endl;
	cout << "Sygnal zdekodowany:" << endl;
	dekodowanieham = dekodowanie(demodulacjaFSK);
	cout << dekodowanieham<< endl;
	ber(byte, dekodowanieham);
	cout << "Napis zdekodowany: " << toasc(dekodowanieham) << endl;

	cout << endl;
	cout << "Sygnal zdemodulowany PSK:" << endl;
	demodulacjaPSK = modulacjaPSK(hammingOut);
	cout << demodulacjaPSK << endl;
	ber(hammingOut, demodulacjaPSK);
	cout << endl;
	cout << "Sygnal zdekodowany:" << endl;
	dekodowanieham = dekodowanie(demodulacjaPSK);
	cout << dekodowanieham << endl;
	ber(byte, dekodowanieham);
	cout << "Napis zdekodowany: " << toasc(dekodowanieham) << endl;

}

