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
#define r 5.0//125.0
#define rho 0.4
#define gama 0.6
#define beta 0.08
#define s 1
//#define d 10//6//2
#define nd 2//5
#define PI 3.14159
#define W 3
//#define IterMax 100//500

#define MAX_N	1000
#define MAX_D	1000

extern int D, n, MaxGeneration;
extern double *masa, *zysk, maxMasaPlecaka;

static int Ld[MAX_N], N[MAX_N][MAX_N], Na[MAX_N];//, randSeed = 1;
//static float s = 0.03;
static double X[MAX_N][MAX_D], Lc[MAX_N], Rd[MAX_N], P[MAX_N][MAX_N];//, Sol[] = {-PI/2, PI/2};

double Distance(int i, int j) {
	double dis = 0;
	for(int k = 0; k < D; ++k)
			dis = dis + pow((X[i][k]-X[j][k]),2);
	return sqrt(dis);
}
void DeployGlowworms(float lim) {
	int lb = 0, ub = 1;
	double rr;

	for(int i = 0; i < n; i++) {
			for(int j = 0; j < D; ++j) {
				rr = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				X[i][j]=round(rr*(ub-lb)+lb); // zmiana na wektory binarne - pozycja odpowiada przedmiotowi, ktory wlozymy do plecaka
				std::cout << X[i][j] << " " ;
			}
		cout<<endl;
	}
}

void UpdateLuciferin() {
	double sum = 0.0, masaPlecaka = 0.0;

	for (int i = 0; i < n; ++i) {
		sum = 0.0; // wartosc funkcji celu
		masaPlecaka = 0.0; 	// sumaryczna masa plecaka dla swietlika (rozwiazania)

		for(int j = 0; j<D; ++j){

			if(masa[j]*X[i][j] + masaPlecaka <= maxMasaPlecaka){ // plecak po dodaniu przedmiotu bedzie co najwyzej pelny
				sum += X[i][j]*zysk[j];
				masaPlecaka += X[i][j]*masa[j]; // dodawaj je do plecaka
			}else{ // chcemy dodac za duzo
				sum= -1.0; // niedopuszczalne rozwiazanie
				break;
			}
		}

		//Lc[i] = (1-rho)*Lc[i] + gama*sum;//J;
		Lc[i] = sum;
	}
}
void FindNeighbors() {
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

void FindProbabilities() {
	for(int i = 0; i < n; ++i) {
		double sum = 0;
		for (int j = 0; j < n; ++j) sum = sum + N[i][j]*(Lc[j] - Lc[i]);

		for(int j = 0; j < n; ++j) {
			if (sum != 0) P[i][j] = N[i][j]*(Lc[j] - Lc[i])/sum;
			else P[i][j] = 0;
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
void UpdateNeighborhood() {
	for (int i = 0; i < n; ++i)
		Rd[i] = max(0.0, min(r, Rd[i] + beta*(nd - Na[i])));
}




//int main (int argc, char * const argv[]) {
void GSO(){

	srand(time(0));
	DeployGlowworms(W);

	for (int i = 0; i < n; ++i) {
		Lc[i] = 0;//5;
		Rd[i] = r;
	}

	for (int t = 0; t < MaxGeneration; ++t) {
		UpdateLuciferin();

		FindNeighbors();

		FindProbabilities();

		SelectLeader();

		Move();

		UpdateNeighborhood();
	}

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

	cout<<endl<<"Macierz Na"<<endl;
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

