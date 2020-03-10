#define _USE_MATH_DEFINES
#include<iostream>
#include "fstream"
#include <string>
#include <math.h>

using namespace std;
//A = 7 B = 6 C = 5 D = 4 E = 4
template<typename T>
void GenerateData(T* Xtable, T* Ytable, int length, string name)
{
	ofstream outfile;
	outfile.open(name + ".plt");
	outfile << "set xzeroaxis" << endl;
	outfile << "set yzeroaxis" << endl;
	outfile << "plot \"" << name << ".plt\" using 1:2 with dots" << endl;

	for (int i = 0; i < length; i++)
	{
		outfile << Xtable[i] << "\t" << Ytable[i] << endl;
	}
	outfile.close();
}


double s(double a, double t, double f, double fi ) {
	return (a * sin(2 * M_PI * t * f + fi));
}



int main() {
	int min = 0;
	int max = 7;
	double a = 7;
	double b = 6;
	double c = 5;
	double deltaT = 1 / 22050.;
	int length = (max - min) / deltaT;
	double* tableX = new double[length];
	double* tableY = new double[length];
	double xtmp = min;
	for (int i = 0; i < length; i++)
	{
		tableX[i] = xtmp;
		tableY[i] = s(a,xtmp,b,(c*M_PI));
		xtmp += deltaT;
	}
	GenerateData(tableX, tableY, length + 1, "probki");


}