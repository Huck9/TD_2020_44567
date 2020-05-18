#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <bitset>

using namespace std;

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
int main()
{
	
	int size = 3;
	int* d1 = new int[size];
	int* d2 = new int[size];
	int d1T[4][1];
	//int d2T[4][1];
	int K[7][1];
	int Ks[8][1];
	int s[4][1];
	int g[7][4] = { {1,1,0,1},{1,0,1,1},{1,0,0,0},{0,1,1,1 },{0,1,0,0},{0,0,1,0},{0,0,0,1} };
	int h[3][7] = { {1,0,1,0,1,0,1},{0,1,1,0,0,1,1,}, {0,0,0,1,1,1,1} };
	string byte = tobyte("C", false);
	byte = "1101000";
	cout << byte << endl;
	cout << "macierz G" << endl;
	for (int i = 0; i < 7; i++)
	{
		cout << "|";
		for (int j = 0; j < 4; j++)
		{
			cout << " " << g[i][j];
		}
		cout << "|" << endl;
	}
	cout << "macierz H" << endl;
	for (int i = 0; i < 3; i++)
	{
		cout << "|";
		for (int j = 0; j < 7; j++)
		{
			cout << " " <<  h[i][j];
		}
		cout << "|" << endl;
	}
	//print po zamianie
	cout << "Pakiet:";
	for (int i = 0; i < size+1; i++)
	{
		d1[i] = byte[i] - '0';
		d2[i] = byte[i + size + 1] - '0';
		cout << d1[i];
	}	
	
	//transponowanie
	cout << endl;
	cout << "Macierz transponowana:" << endl;
	for (int i = 0; i < 4; i++)
	{
		d1T[i][0] = d1[i];
		cout <<"|"<< d1T[i][0]<< "|" << endl;
	}
	/*cout << endl;
	for (int i = 0; i < 4; i++)
	{
		d2T[i][0] = d2[i];
		cout << d2T[i][0] << endl;
	}*/
	for (int i = 0; i < 7; i++)
	{
		K[i][0] = 0;
	}
	//mnożenie macierzy
	cout << endl;
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				K[i][j] += g[i][k] * d1T[k][j];
			
			}
		}
	}
	//modulo 2
	cout << "macierz K:" << endl;
	for (int i = 0; i < 7; i++)
	{
		K[i][0] = K[i][0] % 2;
		cout << "|" << K[i][0]<< "|" << endl;
	}
	//
	K[2][0] = (K[2][0] == 1) ? 0 : 1;
	cout << "macierz K przed naprawa:" << endl;
	for (int i = 0; i < 7; i++)
	{
		cout << "|" << K[i][0] << "|" << endl;
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			for (int k = 0; k < 7; k++)
			{
				s[i][j] += h[i][k] * K[k][j];

			}
		}
	}
	//wyznaczanie błędu
	int blad = 0;
	cout << "Macierz S:" << endl;
	for (int i = 0; i < 3; i++)
	{
		s[i][0] = s[i][0] % 2 * -1;
		blad += pow(2,i) * s[i][0];
		cout << "|" << s[i][0] << "|" << endl;
	}

	cout << "Miejsce blednego bitu: " << blad << endl;
	K[blad-1][0] = (K[blad-1][0] == 1) ? 0 : 1;
	cout << "macierz K po naprawie:" << endl;
	for (int i = 0; i < 7; i++)
	{
		cout << "|" << K[i][0] << "|" << endl;
	}
	cout <<"SECDEDS"<< endl;
	//SECDEDS
	int sum = 0;
	for (int i = 0; i < 7; i++)
	{
		sum += K[i][0];
		Ks[i][0] = K[i][0];
	}
	Ks[7][0] = sum%2;
	cout << "Macierz k z 8 bitem:" << endl;
	for (int i = 0; i < 8; i++)
	{
		cout << "|" << Ks[i][0] << "|" << endl;
	}
	Ks[2][0] = (Ks[2][0] == 1) ? 0 : 1;
	Ks[4][0] = (Ks[4][0] == 1) ? 0 : 1;
	cout << "macierz K przed naprawa:" << endl;
	sum = 0;
	for (int i = 0; i < 8; i++)
	{
		cout << "|" << Ks[i][0] << "|" << endl;
	}
	for (int i = 0; i < 7; i++)
	{
		sum += Ks[i][0];
	}

	if (sum != Ks[7][0])
	{
		cout << "Blad bitu" << endl;
	}
	else {
		cout << "Ostatni bit ok" << endl;
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
	//wyznaczanie błędu
	blad = 0;
	cout << "Macierz S:" << endl;
	for (int i = 0; i < 3; i++)
	{
		s[i][0] = s[i][0] % 2;
		blad += pow(2, i) * s[i][0];
		cout << "|" << s[i][0] << "|" << endl;
	}
	if (blad != 0) {
		cout << "Miejsce blednego bitu: " << blad << endl;
		Ks[blad - 1][0] = (Ks[blad - 1][0] == 1) ? 0 : 1;
		cout << "macierz K po naprawie:" << endl;
		for (int i = 0; i < 7; i++)
		{
			cout << "|" << Ks[i][0] << "|" << endl;
		}
	}
	sum = 0;
	for (int i = 0; i < 7; i++)
	{
		sum += Ks[i][0];
	}

	if (sum != Ks[7][0])
	{
		cout << "Blad bitu" << endl;
	}
	else {
		cout << "Ostatni bit ok" << endl;
	}
}
