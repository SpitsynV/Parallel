#include "functions.h"
double getNorm(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> x)
{
    int n = a.size();
    double el = 0.0, nm = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            el = 0;
            for (int k = 0; k < n; k++)
            {
                el += ((a[i][k]) * (x[k][j]));
            }
            if (i == j)
                el -= 1.0;
            nm += el * el;
        }
    }
    nm = sqrt(nm);
    return nm;
}
