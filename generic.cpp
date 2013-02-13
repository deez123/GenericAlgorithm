/*
 *  file: generic.cpp
 *  autor: Alexandr Pavlovets
 */

#include "generic.h"

void MyGenetic::dec2bin(long decimal, char *binary, int len)
{
  int  k1 = 0, n1 = 0;
  int  remain;
  char temp[80];

  do 
  {
    remain    = decimal % 2;
    decimal   = decimal / 2;
    temp[k1++] = remain + '0';
  } while (decimal > 0);

  // reverse the spelling
  while (n1<(len-k1))
	  binary[n1++] = '0';

  while (k1 >= 0)
    binary[n1++] = temp[--k1];
 
  binary[n1-1] = 0;// end with NULL
}

int MyGenetic::bin2dec(const char *bin)   
{
  int  b1, k1, m1, n1;
  int  len, sum = 0;
 
  len = strlen(bin) - 1;
  for(k1 = 0; k1 <= len; k1++) 
  {
    n1 = (bin[k1] - '0'); // char to numeric value
    if ((n1 > 1) || (n1 < 0)) 
    {
      return 0;
    }
    for(b1 = 1, m1 = len; m1 > k1; m1--) 
    {
      b1 *= 2;
    }
    sum = sum + n1 * b1;
  }
  return(sum);
}

MyGenetic::MyGenetic(double (*f)(double *), double* a, double* b, int variable_number, int population_size, int bits_per_gen, double Pcrossover, double Pmutation, int miter, double mdeviation,int mdeviter, double Pelitizm, int mrestart, bool mmax)
{
	if (bits_per_gen>15)//max_rand 32767 => 15 bits
	{
		std::cout<<"Max number of bits per gen 15!"<<std::endl;
		bits_per_gen = 15;
	}
	
	find_max = mmax;
	mindeviter = mdeviter;
	Pe = Pelitizm;
	restartNo = 0;
	maxrestart = mrestart;
	solutions = new double*[maxrestart];
	iter = 0;
	maxiter=miter;	
	mindeviation = mdeviation;
	deviation = mindeviation + 1;
	function = f;
	n = variable_number;
	N = population_size;
	sigma = bits_per_gen;
	Pc = Pcrossover;
	Pm = Pmutation;
	k = (int)pow(2.0,sigma) - 1;
	bslen = n*sigma;
	A = a;
	B = b;

	S = NULL;
	X = NULL;
	Fitness = NULL;

    InitializePopulation();
}

double MyGenetic::Ps(int i)
{
	if (find_max)
	{
		if (all_fitness_plus>0)
		{
			if (Fitness[i]<0) return 0;
			else
			{
				return Fitness[i]/all_fitness_plus;
			}
		}
		else
		{
			if (all_fitness_minus == 0)
			{
				deviation = -1;
				return 1;
			}
			return abs(Fitness[i])/all_fitness_minus;
		}
	}
	else
	{
		if (all_fitness_minus>0)
		{
			if (Fitness[i]>0) return 0;
			else
			{
				return abs(Fitness[i])/all_fitness_minus;
			}
		}
		else
		{
			if (all_fitness_plus == 0)
			{
				deviation = -1;
				return 1;
			}
			return Fitness[i]/all_fitness_plus;
		}
	}
}

void MyGenetic::Calculate_all_fitness()
{
	all_fitness_plus = 0;
	max_plus = 0;
	all_fitness_minus = 0;
	max_minus = 0;
	
	for (int i=0; i<N; ++i)
	{
		Fitness[i] = Calculate_fitness(X[i]);
		
		if (Fitness[i] > 0)
		{
		     all_fitness_plus += Fitness[i];
			 if (Fitness[i]>max_plus) max_plus = Fitness[i];
		}
		else
		{
			 all_fitness_minus += abs(Fitness[i]);
			 if (abs(Fitness[i])>max_minus) max_minus = abs(Fitness[i]);
		}
	}

	if ((all_fitness_plus==0)&&(find_max))
	{
		for (int i=0; i<N; ++i)
		{
			Fitness[i] = -max_minus - Fitness[i];
		}

		all_fitness_minus = 0;

		for (int i=0; i<N; ++i)
		{
			all_fitness_minus += abs(Fitness[i]);
		}
	}

	if ((all_fitness_minus==0)&&(!find_max))
	{
		for (int i=0; i<N; ++i)
		{
			Fitness[i] = max_plus - Fitness[i];
		}

		all_fitness_plus = 0;

		for (int i=0; i<N; ++i)
		{
			all_fitness_plus += Fitness[i];
		}
	}
}

void MyGenetic::InitializePopulation()
{
	if (S!=NULL) free(S);
	if (X!=NULL) free(X);
	if (Fitness!=NULL) free(Fitness);

	S = new char*[N];
	X = new double*[N];
	Fitness = new double[N];

	for (int i=0; i<N; ++i)
	{
		S[i] = new char[bslen];
		strcpy(S[i],Random_binary_string());
		
		X[i] = new double[n];		
	}

	this->Encode_all_strings();
	this->Calculate_all_fitness();
	this->CalculateDeviation();
}

char* MyGenetic::Random_binary_string()
{
	char* bstring = new char[bslen];
	strcpy (bstring,"");
	for (int i=0; i<n; ++i)
	{	
		char* bsubstring = new char[sigma];	
		int maxvalue = (int)pow(2.0,sigma) - 1;
		dec2bin(rand()%maxvalue,bsubstring,sigma);
		strcat(bstring,bsubstring);
	}
	return bstring;
}

double* MyGenetic::Encode_binary_string(char *bs)
{
	char* variable_bin = new char[sigma];
	int variable_dec;
	double *EncodeBS = new double[n];

	for (int i=0, i1=0; i<bslen; i+=sigma, ++i1)
	{
		for (int j=i, j1=0; j<(i+sigma); ++j,++j1)
		{
			variable_bin[j1] = bs[j];
		}
		
		variable_bin[sigma] = 0;

		variable_dec = bin2dec(variable_bin);
		double h = (B[i1] - A[i1])/k;
		EncodeBS[i1] = A[i1] + variable_dec*h;
	}

	return EncodeBS;
}

double MyGenetic::Calculate_fitness(double *x)
{
	return function(x);
}

void MyGenetic::ShowPopulation()
{
	if (S==NULL)
	{
		std::cout<<"Population not initialized!"<<std::endl<<std::endl;
	}
	else
	{
	   std::cout<<"RestartNo "<<restartNo<<" Iteration "<<iter<<std::endl;
	   for (int i=0; i<N; ++i)
	   {      
		  std::cout<<"S["<<i<<"] = "<<S[i]<<"  ";

		  std::cout<<"Fitness = "<<Fitness[i]<<"  ";

		  std::cout<<"Ps = "<<Ps(i)<<"  ";
		
	      /*
		  cout<<"X["<<i<<"] = { ";
		  for (int j=0; j<n; ++j)
		  {
		     cout<<X[i][j]<<" ";
		  }
		  cout<<"}";
          */

		  std::cout<<std::endl;
	   }

	  std::cout<<"Deviation "<<deviation<<std::endl;
	  std::cout<<std::endl;
	}
}

void MyGenetic::Select_Roulette_Wheel()
{
	char** newPopulation = new char*[N];

    //add elitizm
	int N1 = (int)N*Pe;
	
	//find best
	double psmax = Ps(0);
	int best = 0;
	for(int i=0; i<N; ++i)
	{
		if (Ps(i)>psmax)
		{
			psmax = Ps(i);
			best = i;
		}		
	}

	//add best to population
	for (int i=0; i<N1; ++i)
	{
	   newPopulation[i] = new char[bslen];
	   strcpy(newPopulation[i],S[best]);
	}


	int rand_val;
	//0<=Ps(i)<=1

	for (int i=N1; i<N; ++i)
	{
		newPopulation[i]=NULL;

		rand_val = rand();

		int check_val = 0;

		for (int j=0; j<N; ++j)
		{
			if ((rand_val>=check_val)&&(rand_val<(check_val+(int)RAND_MAX*Ps(j)))) 
			{
				newPopulation[i] = new char[bslen];
				strcpy(newPopulation[i],S[j]);
				break;
			}

			check_val += (int)RAND_MAX*Ps(j);
		}

		if (newPopulation[i]==NULL)
		{
			newPopulation[i] = new char[bslen];
			strcpy(newPopulation[i],S[N-1]);
		}
	}

	free(S);

	S = newPopulation;
}

void MyGenetic::ShakePopulation()
{
	for (int i=0; i<N; ++i)
	{
		int j1 = rand()%N;
		int j2 = rand()%N;

		if (j1!=j2)
		{
			char* temp = S[j1];
			S[j1] = S[j2];
			S[j2] = temp;
		}
	}
}

void MyGenetic::Crossover_two_points()
{
	int rand_val;

	for(int i=0; i<N; i+=2)
	{
		if ((i+1)<N)
		{		
			rand_val = rand();
		
			//do two points crossover		
			if ((rand_val>=0)&&(rand_val<(int)RAND_MAX*Pc))		
			{			
				int point1 = 0;			
				int point2 = -1;
			
				//while ((point1 == 0)||(point1 == bslen-1)){ point1 = rand()%bslen; }			
				//while ((point2 == 0)||(point2 == bslen-1)||(point2 == point1)){ point2 = rand()%bslen; }

				point1 = rand()%bslen;
				while ((point2 == -1)||(point2 == point1)){ point2 = rand()%bslen; }
			
				if (point1>point2)			
				{				
					int temp = point1;				
					point1 = point2;				
					point2 = temp;			
				}
			
				char* middle1 = new char[point2 - point1];
			
				for (int j=point1,j1=0; j<=point2; ++j,++j1)			
				{				
					middle1[j1] = S[i][j];			
				}
		
				for (int j=point1; j<=point2; ++j)			
				{				
					S[i][j] = S[i+1][j];			
				}

				for (int j=point1,j1=0; j<=point2; ++j,++j1)
				{
					S[i+1][j] = middle1[j1];
				}	
			}
		}
	}
}

void MyGenetic::Mutation()
{
	int rand_val;

	for (int i=0; i<N; ++i)
	{
		for(int j=0; j<bslen; ++j)
		{
		   rand_val = rand();
		   if ((rand_val>=0)&&(rand_val<(int)RAND_MAX*Pm))
		   {
			   if (S[i][j]=='0')
			   {
				   S[i][j] = '1';
			   }
			   else
			   {
				   S[i][j] = '0';
			   }
		   }
		}
	}
}

void MyGenetic::Encode_all_strings()
{
	for(int i=0; i<N; ++i)
	{		
		double* EX = Encode_binary_string(S[i]);
		for(int j=0; j<n; ++j)
		{
			X[i][j] = EX[j];		   
		}	
	}
}

void MyGenetic::CalculateDeviation()
{
	//find best
	double psmax = Ps(0);
	int best = 0;
	for(int i=0; i<N; ++i)
	{
		if (Ps(i)>psmax)
		{
			psmax = Ps(i);
			best = i;
		}		
	}

	//calculate deviation in near generations
	if (iter==0)
	{
		bestperson = Fitness[best];
		deviation = mindeviation + 1;
	}
	else
	{
		deviation = abs(abs(bestperson) - abs(Fitness[best]));
		bestperson = Fitness[best];
	}
}

void MyGenetic::SaveBestPerson()
{
	//find best
	double psmax = Ps(0);
	int best = 0;
	for(int i=0; i<N; ++i)
	{
		if (Ps(i)>psmax)
		{
			psmax = Ps(i);
			best = i;
		}		
	}

	//save best
	solutions[restartNo] = new double[n];

	for (int i=0; i<n; ++i)
	{
		solutions[restartNo][i] = X[best][i];
	}

}

double* MyGenetic::GetBestPerson()
{
	if (solutions[0])
	{
	   double f = function(solutions[0]);
	   int best = 0;

	for (int i=0; i<maxrestart; ++i)
	{
		if (solutions[i])
		{
		   double y = function(solutions[i]);

		   if (find_max)
		   {		  
			   if (y>f)
		   	   {		  
				   f = y;			  
				   best = i;		   
			   }
		   }
		   else
		   {
			   if (y<f)
		   	   {		  
				   f = y;			  
				   best = i;		   
			   }
		   }
		}
	}

	    return solutions[best];
	}
	else
	{
		return NULL;
	}
}

void MyGenetic::ShowSolutionsBuffer()
{
	std::cout<<"Solutions buffer"<<std::endl;
	if (restartNo>0)
	{
		for (int i=0; i<maxrestart; ++i)
		{
			if (solutions[i])
			{
			std::cout<<"Solution "<<i<<" ";

			std::cout<<"X = { ";
		    for (int j=0; j<n; ++j)
		    {
		      std::cout<<solutions[i][j]<<" ";
		    }
		    std::cout<<"} ";

			std::cout<<"F(X) = "<<function(solutions[i])<<" ";

			std::cout<<std::endl;
			}
		}

		std::cout<<std::endl;
	}
	else
	{
		std::cout<<"No solutions in buffer!"<<std::endl;
	}
}

double* MyGenetic::Run()
{
    for (restartNo=0; restartNo<maxrestart; ++restartNo)
	{
		this->InitializePopulation();
	    
		iter = 0;
		deviter = 0;
		//this->ShowPopulation();
	
		for(iter = 0; iter < maxiter; ++iter)
		{		
			this->Select_Roulette_Wheel();
		    this->ShakePopulation();
		    this->Crossover_two_points();
		    this->Mutation();
		    this->Encode_all_strings();
		    this->Calculate_all_fitness();

			if (deviation<0) break;

		    this->CalculateDeviation();

			if (deviation<mindeviation)
			{
				deviter++;
			}
			else
			{
				deviter = 0;
			}

			if (deviter>mindeviter) break;

			//this->ShowPopulation();

	    }

		this->SaveBestPerson();
	
		//this->ShowPopulation();

	}

	std::cout<<"Iteration number "<<iter<<std::endl;
	return this->GetBestPerson();
}
