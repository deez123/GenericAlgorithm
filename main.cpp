/*
 *  file: main.cpp
 *  autor: Alexandr Pavlovets
 */

#include <time.h>
#include "generic.h"

double f1(double* x) //max
{
	return 100*SQR(x[0]*x[0] - x[1]) + SQR(1-x[0]);
	// x1,x2 E [-2.048,2.048]
}

double f2(double* x) //max
{
	return (-2*pow(x[1],3) + 6*pow(x[1],2) + 6*x[1] + 10)*sin(log(x[0])*exp(x[1]));
 	//0.5 <= x1 <= 1.1
    //1.0 <= x2 <= 4.6
}

double f3(double* x) //min
{
	double sum = 0;
	int n = 2;
	for (int i=0; i<n; ++i)
	{
		sum += -x[i] * sin( sqrt (abs ( x[i] ) ) );
	}
	return sum;
	//-500<=x(i)<=500
}

double f4(double* x) //min
{
	double a = 1;
	double b = 5.1/(4*pow(M_PI,2));
	double c = 5/M_PI;
	double d = 6;
	double e = 10;
	double f0 = 1/(8*M_PI);
	return a*pow((x[1] - b*pow(x[0],2) + c*x[0] - d),2) + e*(1-f0)*cos(x[0]) + e;
    //-5<=x1<=10, 0<=x2<=15.
}

double f5(double* x) //min
{
	return (4-2.1*pow(x[0],2) + pow(x[0],4)/3)*pow(x[0],2) + x[0]*x[1] + (-4 + 4*pow(x[1],2))*pow(x[1],2);
	//-3<=x1<=3, -2<=x2<=2.
}


int main()
{	
	srand(time(0)%RAND_MAX);

	double (*f)(double*);
	f = f5;
	bool mmax = false;
	int maxrestart = 10;
	double mindeviation = 0.0000001;
	int maxiter = 500;
    int variable_number = 2;
	int population_size = 256;
	int bits_per_gen = 15;
	double Pcrossover = 0.95;
	double Pmutation = 0.01;
	double Pelitizm = 0.2;
	int mdeviter = 20;

	double* a = new double[variable_number];
	double* b = new double[variable_number];
	
	//--------f1-----------------
	//a[0] = -2.048; b[0] = 2.048;
	//a[1] = -2.048; b[1] = 2.048;

	//--------f2-----------------
	//a[0] = 0.5; b[0] = 1.1;
	//a[1] = 1.0; b[1] = 4.6;

	//--------f3-----------------
	//a[0] = -500; b[0] = 500;
	//a[1] = -500; b[1] = 500;

	//--------f4-----------------
	//a[0] = -5; b[0] = 10;
	//a[1] = 0; b[1] = 15;

	//--------f5-----------------
	a[0] = -3; b[0] = 3;
    a[1] = -2; b[1] = 2;

    std::cout<<"Init..."<<std::endl;
	MyGenetic MG(f,a,b,variable_number,population_size,bits_per_gen,Pcrossover,Pmutation,maxiter,mindeviation,mdeviter,Pelitizm,maxrestart,mmax);
	
	std::cout<<"Starting..."<<std::endl;
	double* solution = MG.Run();

	std::cout<<"Solution "<<std::endl;

	std::cout<<"X = { ";
	for (int j=0; j<variable_number; ++j)
	{
		std::cout<<solution[j]<<" ";
	}
	std::cout<<"}"<<std::endl;

	std::cout.precision(15);
	std::cout<<"F(X) = "<<f(solution)<<std::endl;

	//MG.ShowSolutionsBuffer();
	std::cout<<"Done."<<std::endl;
	getchar();
	return 0;
}