	#define _USE_MATH_DEFINES
	#include<iostream>
	#include "fstream"
	#include <string>
	#include <math.h>

	using namespace std;
	//A = 7 B = 6 C = 5 D = 4 E = 4
	template<typename T>
	void GenerateData(T* Xtable, T* Ytable, int length, string name,double deltat)
	{
		ofstream outfile;
		outfile.open(name + ".plt");

		for (int i = 0; i < length; i++)
		{
			outfile << Xtable[i] << "\t" << Ytable[i] << endl;
		}
		outfile.close();
	}


	double s( double t, double f, double fi ) {
		return (sin(2 * M_PI * t * f + fi));
	}


	double quaquantize(double x, int q) {
		int sample_max = pow(2, q - 1) - 1;
		double y = floor(x * sample_max);
		y += pow(2, q - 1);
		return y;
	}


	int main() {
		int min = 0;
		int max = 7;
		double a = 7;
		double b = 6;
		double c = 5;
		double deltaT = 1 / 1000.;
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
		GenerateData(tableX, tableY, length + 1, "probkowanie", deltaT);

		for (int i = 0; i < length; i++)
		{
			tableY[i] = quaquantize(tableY[i], 16);
		}

		GenerateData(tableX, tableY, length + 1, "probkowanie2", deltaT);

		deltaT = 1 / 500.;
		length = (max - min) / deltaT;
		double* tableXq = new double[length];
		double* tableYq = new double[length];
		xtmp = min;
		for (int i = 0; i < length; i++)
		{
			tableXq[i] = xtmp;
			tableYq[i] = quaquantize(s(xtmp, b, (c * M_PI)),8);
			xtmp += deltaT;
		}
		GenerateData(tableXq, tableYq, length + 1, "probkowanie3", deltaT);

	}