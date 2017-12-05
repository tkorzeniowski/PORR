//============================================================================
// Name        : Firefly.cpp
// Authors     : Dr. Iztok Fister and Iztok Fister Jr.
// Version     : v1.0
// Created on  : Jan 23, 2012
//============================================================================

/* Classic Firefly algorithm coded using C/C++ programming language */

/* Reference Paper*/

/*I. Fister Jr.,  X.-S. Yang,  I. Fister, J. Brest, Memetic firefly algorithm for combinatorial optimization,
in Bioinspired Optimization Methods and their Applications (BIOMA 2012), B. Filipic and J.Silc, Eds.
Jozef Stefan Institute, Ljubljana, Slovenia, 2012 */

/*Contact:
Iztok Fister (iztok.fister@uni-mb.si)
*/
/*
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <memory.h>
#include <math.h> // round

#define DUMP	1
#define MAX_FFA	1000
#define MAX_D	1000

using namespace std;

int D = 6; //1000;			// dimension of the problem
int n = 20;			// number of fireflies
int MaxGeneration;		// number of iterations
int NumEval;			// number of evaluations
int Index[MAX_FFA];		// sort of fireflies according to fitness values

double ffa[MAX_FFA][MAX_D];	// firefly agents
double ffa_tmp[MAX_FFA][MAX_D]; // intermediate population
double f[MAX_FFA];		// fitness values
double I[MAX_FFA];		// light intensity
double nbest[MAX_FFA];          // the best solution found so far
double lb[MAX_D];		// upper bound
double ub[MAX_D];		// lower bound

double alpha = 0.9;	//0.5;		// alpha parameter
double betamin = 0.8; //0.2;           // beta parameter
double gama = 2.0; //1.0;		// gamma parameter

double fbest;			// the best objective function

double masa[] = {1.0, 2.0, 2.0, 3.0, 3.0, 5.0};
double zysk[] = {2.0, 3.0, 4.0, 8.0, 5.0, 50.0};
const double maxMasaPlecaka = 10.0;


typedef double (*FunctionCallback)(double sol[MAX_D]);

//benchmark functions
double cost(double sol[MAX_D]);
double sphere(double sol[MAX_D]);

//Write your own objective function
FunctionCallback function = &cost;

// optionally recalculate the new alpha value
double alpha_new(double alpha, int NGen)
{
	double delta;			// delta parameter
	delta = 1.0-pow((pow(10.0, -4.0)/0.9), 1.0/(double) NGen);
	return (1-delta)*alpha;
}

// initialize the firefly population
void init_ffa()
{
	int i, j;
	double r;

	//int x = 16 + rand() % (17); //random_number = _min + rand() % ((_max - _min ) + 1);
	//cout<<"moj x to: "<<x<<endl;

	// initialize upper and lower bounds
	for (i=0;i<D;i++)
	{
		lb[i] = 0.0;
		ub[i] = 1.0; // wektory binarne - zaokraglenie do najblizszego int-a
	}

	for (i=0;i<n;i++)
	{
		for (j=0;j<D;j++)
		{
			r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			ffa[i][j]=round(r*(ub[j]-lb[j])+lb[j]); // zmiana na wektory binarne - pozycja odpowiada przedmiotowi, ktory wlozymy do plecaka
			//if(j<5)
			cout<<ffa[i][j]<<" ";
		}
		cout<<endl;
		f[i] = 1.0;			// initialize attractiveness
		I[i] = f[i];
		//cout<<I[i]<<endl;
	}
}

void change_ffa()
{
	int i, j, k = round(0.2*n); // najgorsze 10% swietlikow
	double r;

	for (i=n-1;i<n-k-1;--i)
	{
		for (j=0;j<D;j++)
		{
			r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			ffa[i][j]=round(r*(ub[j]-lb[j])+lb[j]); // zmiana na wektory binarne - pozycja odpowiada przedmiotowi, ktory wlozymy do plecaka
			//if(j<5)
			cout<<ffa[i][j]<<" ";
		}
		cout<<endl;
		f[i] = 1.0;			// initialize attractiveness
		I[i] = f[i];
		//cout<<I[i]<<endl;
	}
}

// implementation of bubble sort
void sort_ffa()
{
	int i, j;

	// initialization of indexes
	for(i=0;i<n;i++)
		Index[i] = i;

	// Bubble sort
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
			if(I[i] < I[j]) // zmiana znaku nierównoœci na zgodnosc z artykulem
			{
				double z = I[i];	// exchange attractiveness
				I[i] = I[j];
				I[j] = z;
				z = f[i];			// exchange fitness
				f[i] = f[j];
				f[j] = z;
				int k = Index[i];	// exchange indexes
				Index[i] = Index[j];
				Index[j] = k;
			}
		}
	}
}

// replace the old population according the new Index values
void replace_ffa()
{
	int i, j;

	// copy original population to temporary area
	for(i=0;i<n;i++)
	{
		for(j=0;j<D;j++)
		{
			ffa_tmp[i][j] = ffa[i][j];
		}
	}

	// generational selection in sense of EA
	for(i=0;i<n;i++)
	{
		for(j=0;j<D;j++)
		{
			ffa[i][j] = ffa_tmp[Index[i]][j];
		}
	}
}

void findlimits(int k)
{
	int i;

	for(i=0;i<D;i++)
	{
		if(ffa[k][i] < lb[i])
			ffa[k][i] = lb[i];
		if(ffa[k][i] > ub[i])
			ffa[k][i] = ub[i];
	}
}

void move_ffa()
{
	int i, j, k;
	double scale;
	double r, beta;

	for(i=0;i<n;i++)
	{
		scale = abs(ub[i]-lb[i]);
		for(j=0;j<n;j++)
		{
			r = 0.0;
			for(k=0;k<D;k++)
			{
				r += (ffa[i][k]-ffa[j][k])*(ffa[i][k]-ffa[j][k]);
			}
			r = sqrt(r);
			if(I[i] < I[j])	// brighter and more attractive // zmiana znaku nierównoœci na zgodnosc z artykulem
			{
				double beta0 = 1.0;
				beta = (beta0-betamin)*exp(-gama*pow(r, 2.0))+betamin;
				for(k=0;k<D;k++)
				{
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					double tmpf = alpha*(r-0.5)*scale;
					ffa[i][k] = round(ffa[i][k]*(1.0-beta)+ffa_tmp[j][k]*beta+tmpf); // zamiast ulamka zmiana binarna
				}
			}
		}
		findlimits(i);
	}
}

void dump_ffa(int gen)
{
	cout << "Generacja= " << gen << ", f_celu= " << fbest << endl;
}

// display syntax messages
void help()
{
	cout << "Syntax:" << endl;
	cout << "  Firefly [-h|-?] [-l] [-p] [-c] [-k] [-s] [-t]" << endl;
	cout << "    Parameters: -h|-? = command syntax" << endl;
	cout << "				 -n = number of fireflies" << endl;
	cout << "				 -d = problem dimension" << endl;
	cout << "				 -g = number of generations" << endl;
	cout << "				 -a = alpha parameter" << endl;
	cout << "				 -b = beta0 parameter" << endl;
	cout << "				 -c = gamma parameter" << endl;
}

int main(int argc, char* argv[])
{
        int i;
        int t = 1;		// generation  counter

         // interactive parameters handling

        MaxGeneration = 50;
        // firefly algorithm optimization loop
        // determine the starting point of random generator
    //srand(1);
	srand(time(0));

	// generating the initial locations of n fireflies
	init_ffa();
#ifdef DUMP
	dump_ffa(t);
#endif

	while(t <= MaxGeneration)
	{
		// this line of reducing alpha is optional
		//alpha = alpha_new(alpha, MaxGeneration);

		if(t!=1){
			//change_ffa();
		}

		// evaluate new solutions
		for(i=0;i<n;i++)
		{
                        f[i] = function(ffa[i]);                        // obtain fitness of solution
			I[i] = f[i];					// initialize attractiveness
		}
		for(int d = 0; d<n; ++d){
					//cout<<I[d]<<endl;
				}
//cout<<"Koniec wyznaczania funkcji celu dla kazdego swietlika"<<endl;

		// ranking fireflies by their light intensity
		sort_ffa();

		for(int d = 0; d<n; ++d){
			//cout<<I[d]<<endl;
		}
//cout<<"Doszed³em tutaj"<<endl;
		// replace old population
		replace_ffa();
//cout<<"nbest: ";
		// find the current best
		for(i=0;i<D;i++){
			nbest[i] = ffa[0][i];
			cout<<nbest[i]<<" ";
		}
		cout<<endl;
		fbest = I[0];

		// move all fireflies to the better locations
		move_ffa();
#ifdef DUMP
		dump_ffa(t);
#endif
		t++;
	}

	cout << "End of optimization: fbest = " << fbest << endl;

	return 0;
}

// FF test function
double cost(double* sol) // obliczanie funkcji celu
{
	double sum = 0.0; // wartosc funkcji celu
	double masaPlecaka = 0.0; 	// sumaryczna masa plecaka dla swietlika (rozwiazania)

	for(int i = 0; i<D; i++){

		if(masa[i]*sol[i] + masaPlecaka <= maxMasaPlecaka){ // plecak po dodaniu przedmiotu bedzie co najwyzej pelny
			sum += sol[i]*zysk[i];
			//cout<<"masa "<<masa[i]<<" masaPlecaka "<<"masaPlecaka";
			masaPlecaka += sol[i]*masa[i]; // dodawaj je do plecaka
		}else{ // chcemy dodac za duzo
			//double x = masa[i]*sol[i];
			//cout<<"zle "<<i<< " "<< x <<" "<< masaPlecaka<<" ";
			sum= -1.0; // niedopuszczalne rozwiazanie
			//masaPlecaka = 0.0;
			break;
		}

	}
	//cout<<" masa "<<masaPlecaka<<", zysk "<<sum<<endl;
	return sum;
}

double sphere(double* sol) {
	int j;
	double top = 0;
	for (j = 0; j < D; j++) {
		top = top + sol[j] * sol[j];
	}
	return top;
}
*/
