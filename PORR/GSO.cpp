//-----------------GSO.cpp---------------------------------------------------
// Glowworm swarm optimization (GSO)
// Developed by K.N. Kaipa and D. Ghose in 2005
// This is the main front-end code
//-------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cmath>
#include <omp.h>
#include <random>
#include "GSO.h"
using namespace std;



/* Parameters
Constants
----------
n - Number of glowworms
r - Sensor range
rho - Luciferin decay constant
gama - Luciferin enhancement constant
beta - neighborhood ehancement constant
s - step size (agent speed)
nd - Desired number of neighbors
D - dimension of the search space
Variables
----------
Lc - Luciferin value
Rd - Neighborhood range
P - Probability matrix
Ld - Leader set
N - Neighborhood matrix
Na - Actual number of neighbors
X - Glowworm positions
W - Workspace size
*/



//#define n 100 //20 //1000
#define r 50.0//125.0
#define rho 0.4//0.4
#define gama 0.6//0.6
#define beta 0.8//0.08
#define s 5
//#define d 10//6//2
#define nd 50//5
#define PI 3.14159
#define W 3
//#define IterMax 100//500
#define epsilon 400

#define MAX_N	1000
#define MAX_D	1000

extern int D, n, MaxGeneration;
extern double *masa, *zysk, maxMasaPlecaka;

static int Ld[MAX_N], N[MAX_N][MAX_N], Na[MAX_N];//, randSeed = 1;
//static float s = 0.03;
static double X[MAX_N][MAX_D], Lc[MAX_N], Rd[MAX_N], P[MAX_N][MAX_N];//, Sol[] = {-PI/2, PI/2};



class RNG
{
public:
    typedef std::mt19937 Engine;
    typedef std::uniform_real_distribution<double> Distribution;

    RNG() : engines(), distribution(0.0, 1.0)
    {
        int threads = std::max(1, omp_get_max_threads());
        for(int seed = 0; seed < threads; ++seed)
        {
            engines.push_back(Engine(seed));
        }
    }

    double operator()()
    {
        int id = omp_get_thread_num();
        return distribution(engines[id]);
    }

    std::vector<Engine> engines;
    Distribution distribution;
};





double Distance(int i, int j) {
	double dis = 0;
	for(int k = 0; k < D; ++k)
			dis = dis + pow((X[i][k]-X[j][k]),2);
	//return sqrt(dis);
	return dis;
}
void DeployGlowworms(int liczbaWatkow) {
	int lb = 0, ub = 1;
	double rr;

	RNG rand232;
	omp_set_nested(1);
    #pragma omp parallel num_threads(liczbaWatkow) private(rr)
	{
	        #pragma omp for
	        for(int i = 0; i < n; i++) {
                for(int j = 0; j < D; ++j) {
                    rr = rand232();
                    X[i][j]=round(rr*(ub-lb)+lb); // zmiana na wektory binarne - pozycja odpowiada przedmiotowi, ktory wlozymy do plecaka
                   //  std::cout << X[i][j] << " " ;
                }
            // cout<<endl;
        }
	}
}

void UpdateLuciferin(int liczbaWatkow) {
	double sum = 0.0, masaPlecaka = 0.0;
	omp_set_nested(0);
    #pragma omp parallel num_threads(liczbaWatkow)
    {

        #pragma omp for
        for (int i = 0; i < n; ++i) {
            sum = 0.0; // wartosc funkcji celu
            masaPlecaka = 0.0; 	// sumaryczna masa plecaka dla swietlika (rozwiazania)

            for(int j = 0; j<D; ++j){

                if(masa[j]*X[i][j] + masaPlecaka <= maxMasaPlecaka){ // plecak po dodaniu przedmiotu bedzie co najwyzej pelny
                    sum += X[i][j]*zysk[j];
                    masaPlecaka += X[i][j]*masa[j]; // dodawaj je do plecaka
                }else{ // chcemy dodac za duzo
                    sum= 0.0; // niedopuszczalne rozwiazanie
                    break;
                }
            }

            //Lc[i] = (1-rho)*Lc[i] + gama*sum;//J;
            Lc[i] = sum;
        }
    }
}
void FindNeighbors(int liczbaWatkow) {
    omp_set_nested(0);
	#pragma omp parallel num_threads(liczbaWatkow)
	{
        #pragma omp for
        for(int i = 0; i < n; ++i) {
            N[i][i] = 0; Na[i] = 0;
            for(int j = 0; j < n; ++j){
                if (j!=i){
                    if ((Lc[i] < Lc[j]) && (Distance(i,j) < Rd[i])) N[i][j] = 1;
                    else N[i][j] = 0;

                    Na[i] = Na[i] + N[i][j];
                }
            }
        }
	}
}

void FindProbabilities(int liczbaWatkow) {
	omp_set_nested(0);
	#pragma omp parallel num_threads(liczbaWatkow)
	{
        #pragma omp for
        for(int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) sum = sum + N[i][j]*(Lc[j] - Lc[i]);

            for(int j = 0; j < n; ++j) {
                if (sum != 0) P[i][j] = N[i][j]*(Lc[j] - Lc[i])/sum;
                else P[i][j] = 0;
            }
        }
	}
}


void SelectLeader() {
	for (int i = 0; i < n; ++i) {
		double b_lower = 0;
		Ld[i] = i;
		double toss = rand()/(RAND_MAX + 1.0);
		for (int j = 0; j < n; ++j) {
			if (N[i][j] == 1) {
				double b_upper = b_lower + P[i][j];
				if ((toss >= b_lower) && (toss < b_upper)) {
					Ld[i] = j;
					break;
				} else b_lower = b_upper;
			}
		}
	}
}


void Move() {
	for (int i = 0; i < n; ++i) {
		if (Ld[i]!=i) { // nie jestem liderem
			int flag = 0;
			double temp[D];
			double dis = Distance(i,Ld[i]);
			//if (Na[i] > 15) s = 0.0001; else s = 0.03;
			for (int j = 0; j < D; ++j) {
				temp[j] = round(X[i][j] + s*(X[Ld[i]][j] - X[i][j])/dis);
				if (fabs(temp[j]) > W) {
					flag = 1;
					break;
				}
			}
			if (flag == 0) for (int j = 0; j < D; ++j) X[i][j] = temp[j];
		}
	}

}
void UpdateNeighborhood(int liczbaWatkow) {
    #pragma omp parallel num_threads(liczbaWatkow)
    {
        #pragma omp for
        for (int i = 0; i < n; ++i)
            Rd[i] = max(0.0, min(r, Rd[i] + beta*(nd - Na[i])));
    }
}




//int main (int argc, char * const argv[]) {
void GSO(int liczbaWatkow){

	srand(time(0));
	DeployGlowworms(liczbaWatkow);

	//double rozw = 0, foBest = 0;
	int licznik = 0;
	for (int i = 0; i < n; ++i) {
		Lc[i] = 5;//0;//5;
		Rd[i] = r;
	}

	int bestNa = -1;

	for (int t = 0; t < MaxGeneration; ++t) {
		UpdateLuciferin(liczbaWatkow);

		FindNeighbors(liczbaWatkow);

		FindProbabilities(liczbaWatkow);

		SelectLeader();

		Move();

		UpdateNeighborhood(liczbaWatkow);
		/*
		cout<<endl<<"Najlepszy swietlik"<<endl;
			for (int j=0; j<D; ++j){
				cout<<X[indexNa][j]<<" ";
			}
			*/
		int minNa = abs(Na[0]), indexNa = 0;
		for (int i=0; i<n; ++i){
			if(abs(Na[i])<=minNa && Lc[i]>Lc[indexNa]){
				minNa = Na[i];
				indexNa = i;
			}
			//cout<<Na[i]<<" ";
		}

		if(t>0 && indexNa == bestNa){
			//cout<<"Indeks siê nie zmieni³ od "<<licznik<<" iteracja "<<t<<endl;
			++licznik;
			if(licznik == epsilon){
				cout<<"Osi¹gnieto warunek stopu"<<endl;
				break;
			}
		}else{
			cout<<"Nast¹pi³a zmiana w iteracji "<<t<<endl;
			cout << Lc[indexNa] << endl;
			licznik = 0;
			bestNa = indexNa;
		}



	/*
		rozw=0;
		for(int i =0; i<n; ++i){
			if(Lc[i]>rozw){
				rozw=Lc[i];
			}
		}
	*/
		//cout<<abs(rozw-foBest);
	/*
		if(t!= 0 && abs(rozw-foBest) <= 10){
			++licznik;
			cout<<"Jestem w ifie, a moj licznik to "<<licznik<<endl;
			cout<<"rozw "<<rozw<<" i foBest "<<foBest<<endl;
			if(licznik == 100){
				break;
			}
		}else{
			cout<<"Zmieni³o siê rozwi¹zanie"<<endl;
			cout<<"stare rozw "<<rozw<<" i foBest "<<foBest<<endl;
			licznik=0;
		}
		foBest=rozw;
	*/



	}
/*
	cout<<endl<<"Macierz N"<<endl;
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			cout<<N[i][j]<<" ";
		}
		cout<<endl;
	}

	cout<<endl<<"Macierz Ld"<<endl;
	for (int i=0; i<n; ++i){
		cout<<Ld[i]<<" ";
	}
*/
	cout<<endl<<"Macierz Na"<<endl;
	//int minNa = abs(Na[0]), indexNa = 0;
	int minNa = abs(Na[0]), indexNa = 0;
	for (int i=0; i<n; ++i){
		if(abs(Na[i])<=minNa && Lc[i]>Lc[indexNa]){
			minNa = Na[i];
			indexNa = i;
		}
		cout<<Na[i]<<" ";
	}

	cout<<endl<<"Macierz x"<<endl;
	for (int i=0; i<n; ++i){
		for (int j=0; j<D; ++j){
			cout<<X[i][j]<<" ";
		}
		cout<<endl;
	}

	cout<<endl<<"Macierz Lc"<<endl;
	for (int i=0; i<n; ++i){
		cout<<Lc[i]<<" ";
	}


	// znalezione rozwiazanie
	cout<<endl<<"Najlepszy swietlik"<<endl;
	for (int j=0; j<D; ++j){
		cout<<X[indexNa][j]<<" ";
	}
	cout<<endl<<"Wartosc f_celu: "<<Lc[indexNa]<<" uzyskal swietlik "<<indexNa+1<<endl;

	cout<<endl<<"koniec GSO"<<endl;
}

