#define _USE_MATH_DEFINES
#include<iostream>
#include "fstream"
#include <string>
#include <math.h>

using namespace std;
template<typename T>
void GenerateData(T* Xtable, T* Ytable, int length, string name)
{
	ofstream outfile;
	outfile.open(name + ".plt");
	outfile << "set xzeroaxis" << endl;
	outfile << "set yzeroaxis" << endl;
	outfile << "plot \""<<name <<".plt\" using 1:2 with dots" << endl;
	
	for (int i = 0; i < length; i++)
	{
		outfile << Xtable[i] << "\t" << Ytable[i] << endl;
	}
	outfile.close();
}
//44567 A = 7 B = 6 C = 5
int A = 7;
int B = 6;
int C = 5;

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
	double delta =B^2 -  4 * A * C;
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
	GenerateData(tableX, tableY, lenght + 1,"zadanie1");
}

double X(double x) {
	return( A * pow(x, 2) + B * x + C);
}

double Y(double x) {
	return(2 * pow(X(x), 2) + 12 * cos(x));
}

double Z(double x) {
	return(sin(2 * M_PI * 7 * x) * X(x) - 0.2 * log10(Y(x + M_PI)));
}

double U(double x) {
	return(sqrt(abs(Y(x)*Y(x)*Z(x))- 1.8*sin(0.4*x*Z(x)*X(x))));
}

double V(double x) {
	if (x < 0.22) {
		return((1 - 7 * x) * sin((2 * M_PI * 10) / (x + 0.4)));
	}
	else if(x < 0.77)
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

int main() {
	//Zadanie 1
	squareFunction(-10,10,(1/100.));

	//Zadanie 2
	//Y
	int min = 0;
	int max = 1;
	double deltaT = 1 / 22050.;
	int length = (max - min) / deltaT;
	double* tableX = new double[length];
	double* tableYy = new double[length];
	double xtmp = min;
	for (int i = 0; i < length; i++)
	{
		tableX[i] = xtmp;
		tableYy[i] = Y(xtmp);
		xtmp += deltaT;
	}

	GenerateData(tableX, tableYy,length, "zadanie2y");
	//Z
	double* tableYz = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYz[i] = Z(tableX[i]);
	}
	GenerateData(tableX, tableYz, length, "zadanie2z");
	//U
	double* tableYu = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYu[i] =U(tableX[i]);
	}
	GenerateData(tableX, tableYu, length, "zadanie2u");

	//V
	double* tableYv = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYv[i] = V(tableX[i]);
	}
	GenerateData(tableX, tableYv, length, "zadanie2v");

	//P2

	double* tableYp = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYp[i] = P(tableX[i],2);
	}
	GenerateData(tableX, tableYp, length, "zadanie2p");

	double* tableYp4 = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYp4[i] = P(tableX[i], 4);
	}
	GenerateData(tableX, tableYp4, length, "zadanie2p4");
	//A = 7 B= 6 AB = 42
	double* tableYpAB = new double[length];
	for (int i = 0; i < length; i++)
	{
		tableYpAB[i] = P(tableX[i], 42);
	}
	GenerateData(tableX, tableYpAB, length, "zadanie2pAB");
}


