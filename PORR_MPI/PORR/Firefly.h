#include <iostream>

double alpha_new(double alpha, int NGen);
void init_ffa(int liczbaWatkow);
void sort_ffa();
void replace_ffa(int liczbaWatkow);
void findlimits(int k, int liczbaWatkow);
void move_ffa(int liczbaWatkow);
//void FA(int liczbaWatkow);
double* FA(int liczbaWatkow);
double cost(double* sol);
