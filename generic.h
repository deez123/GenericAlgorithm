/*
 *  file: generic.h
 *  autor: Alexandr Pavlovets
 */

#pragma once

#pragma warning(disable:4996)
#pragma warning(disable:4244)

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>

#define SQR(x) pow((x),2.0)

class MyGenetic
{
private:
	bool find_max;//what we try to find max  - find_max=true, min - find_max=false

	int restartNo;//number of restarts
	int maxrestart;//max number of restarts algorithm

	double** solutions;//solutions per restart

	int iter;//current iteration
	int maxiter;//max iterations
	
	double deviation;//average deviation in population per restart
	double mindeviation;//min posible deviation in population
	double bestperson;//best person in previous generation
	int deviter;//current number of iteration to stop by deviation
	int mindeviter;//min number of iteration to stop by deviation

	double (*function)(double*);//function
	double* A;//left borders of variables
	double* B;//right borders of variables

	int N;//population size
	int sigma;//bits per gen
	int n;//number of gens(variables in function)
	int k;//number of segments per interval for variable
	int bslen;//length of chromosome

	char** S;//persons chromosoms
	double** X;//persons fenotips
	double* Fitness;//person fitness ("weight")
	double all_fitness_plus;//sum of fitness of all persons +
 	double all_fitness_minus;//sum of fitness of all persons -
	double max_minus;
	double max_plus;

	double Pc;//probability of crossover
	double Pm;//probability of mutation per bit
	double Pe;//percent of elitizm

	void dec2bin(long decimal, char *binary, int len);
	int bin2dec(const char *bin);

	double Ps(int i);//fitness function
	char* Random_binary_string();
	double* Encode_binary_string(char* bs);
	double Calculate_fitness(double* x);

	void InitializePopulation();
	void ShakePopulation();
	void Select_Roulette_Wheel();
	void Crossover_two_points();
	void Mutation();
	void Encode_all_strings();
	void Calculate_all_fitness();
	void CalculateDeviation();
	void SaveBestPerson();
	
public:	
	void ShowPopulation();
	void ShowSolutionsBuffer();
	double* GetBestPerson();
	double* Run();

	MyGenetic(double (*f)(double*), double* a, double* b, int variable_number, int population_size, int bits_per_gen, double Pcrossover, double Pmutation, int miter, double mdeviation, int mdeviter, double Pelitizm, int mrestart, bool mmax);
};

