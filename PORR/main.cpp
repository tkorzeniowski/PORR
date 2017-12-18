#include <iostream>
#include <string>
#include <fstream>
/*
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <memory.h>
#include <algorithm>
*/
#include "Firefly.h"
#include "GSO.h"
#include <Windows.h>

using namespace std;

int D; 				// dimension of the problem
int N, n; 			// number of fireflies
int MaxGeneration; 	// number of iterations



double *masa, *zysk;
double maxMasaPlecaka;

void readCostM(string fileName)
{
    ifstream file(fileName.c_str());
    if(file.is_open()){
    	file >> D;
    	file >> N;
    	n = N;
    	file >> MaxGeneration;
    	file >> maxMasaPlecaka;

    	masa = new double[D];
    	zysk = new double[D];

        for(int i = 0; i < D; ++i){
        	file >> masa[i];
			file >> zysk[i];



        	/*
        	for(int j=0; j<N; ++j)
        	{
        		file >> myArray[i][j];
        	}
        	*/
        }
    }
}


int main(){

    LARGE_INTEGER frequency;
    LARGE_INTEGER poczatek;
    LARGE_INTEGER koniec;
    double elapsedSeconds;
    int liczbaWatkow = 8;
    double wynikKoncowy = 0;
    for (int i = 0; i < 1; i++)
    {
        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&poczatek);
        /* Do stuff */
        readCostM("plik.txt");
     //   cout<<"Algorytm FA"<<endl;

        FA(liczbaWatkow);
     //   cout<<"Algorytm GSO"<<endl;
        // GSO(liczbaWatkow);

         QueryPerformanceCounter(&koniec);
         elapsedSeconds = (koniec.QuadPart - poczatek.QuadPart) / (double)frequency.QuadPart;
     //   cout << "\nWynik: " << elapsedSeconds << endl;
      //  cout<<"KONIEC"<<endl;
        wynikKoncowy += elapsedSeconds;
    }
    cout << "\n" << wynikKoncowy;
	return 0;
}

