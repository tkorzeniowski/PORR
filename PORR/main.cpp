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
	readCostM("plik.txt");
	cout<<"Algorytm FA"<<endl;
	FA();
	cout<<"Algorytm GSO"<<endl;
	GSO();

	cout<<"KONIEC"<<endl;
	return 0;
}
