#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <cmath>
#include <algorithm>
#include<vector>

double getNorm(const std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &x);
double Norm(int n,int totalThreads,std::vector<std::vector<double>>& a,std::vector<std::vector<double>>& x, double* ans);
void* NormProcess(void* arg);
#endif

