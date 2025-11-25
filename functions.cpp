#include "functions.h"
#include <chrono>
struct threadDataN {
    int n;
    int myID;
    int p;//total threads
    std::vector<std::vector<double>> &a;
    std::vector<std::vector<double>> &x;
    double *ans;
    threadDataN(int n_, int myID_, int p_,std::vector<std::vector<double>>& a_,std::vector<std::vector<double>>& x_,double* ans_)
    :n(n_), myID(myID_), p(p_), a(a_), x(x_), ans(ans_) {}
};
static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
double Norm(int n,int totalThreads,std::vector<std::vector<double>>& a,std::vector<std::vector<double>>& x, double* ans){
    pthread_t threadPool[totalThreads];
    std::vector<threadDataN>Data;
    Data.reserve(totalThreads);
    for (int i = 0; i < totalThreads; i++)
    {
        
        Data.emplace_back(n,i,totalThreads,a,x,ans);
    }

    std::chrono::duration<double> elapsed;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < totalThreads; i++)
        if (pthread_create(&threadPool[i], 0, NormProcess, &Data[i]))
        {
            throw std::runtime_error("Error in creating thread");
            return -10;//error
        }
    for (int i = 0; i < totalThreads; i++)
        if (pthread_join(threadPool[i], 0))
        {
            throw std::runtime_error("Error in joining thread");
            return -20;//error
        }
    auto end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    auto time=elapsed.count();
    return time;
}
void* NormProcess(void* arg){
    threadDataN* data = static_cast<threadDataN*>(arg);

    int n=data->n;
    int myID=data->myID;
    int totalThreads=data->p;
    
    int rows_per_thread =n/totalThreads;//распределяем
    int remainder = n % totalThreads;
    int first_row = myID * rows_per_thread + std::min(myID, remainder);
    int last_row = first_row + rows_per_thread - 1 + (myID < remainder ? 1 : 0);
    double el = 0.0, nm = 0.0;
    for(int i=first_row; i<=last_row; i++){
        for(int j=0;j<n;j++){
            el=0;
            for(int k=0; k<n;k++){
                el += ((data->a[i][k]) * (data->x[k][j]));
            }
            if (i == j)
                el -= 1.0;
            nm += el * el;
        }
    }
    pthread_mutex_lock(&m);
    (*data->ans)+=nm;
    pthread_mutex_unlock(&m);
    return 0;
}


double getNorm(const std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &x)
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
