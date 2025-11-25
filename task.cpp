#include "task.h"
#include "matrix.h"
#include <stdexcept>
#include <cmath>
#include <pthread.h>
#include <algorithm>

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
struct threadData {
    int n;
    int myID;
    int p;//total threads
    std::vector<std::vector<double>> &A;
    std::vector<std::vector<double>> &inv;
    std::vector<int> &columnOrder;
    double *globalMaxVal; // Pointer to globalMaxVal for ALL threads
    int *globalPivotRow;
    int *globalPivotCol;
    int v;//version of programm: 0->row subdivision, 1->col subdivison
    threadData(int n_, int myID_, int p_,std::vector<std::vector<double>>& A_, std::vector<std::vector<double>>& inv_, std::vector<int>& columnOrder_, double* globalMaxVal_, int* globalPivotRow_, int* globalPivotCol_, int v_)
    :n(n_), myID(myID_), p(p_), A(A_), inv(inv_),columnOrder(columnOrder_),globalMaxVal(globalMaxVal_),globalPivotRow(globalPivotRow_), globalPivotCol(globalPivotCol_), v(v_) {}
};

void await(int totalThreads)
{
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t out = PTHREAD_COND_INITIALIZER;
    static int threadIn = 0;
    static int threadOut = 0;

    pthread_mutex_lock(&mutex);

    threadIn++;
    if (threadIn >= totalThreads)
    {
        threadOut = 0;
        pthread_cond_broadcast(&in);
    }
    else
    {
        while (threadIn < totalThreads)
            pthread_cond_wait(&in, &mutex);
    }

    threadOut++;
    if (threadOut >= totalThreads)
    {
        threadIn = 0;
        pthread_cond_broadcast(&out);
    }
    else
    {
        while (threadOut < totalThreads)
            pthread_cond_wait(&out, &mutex);
    }
    pthread_mutex_unlock(&mutex);
}


void* SolveProcess(void *threadArg)
{ //эту функцию собираюсь дать каждому потоку

    //распаковка
    threadData* data = static_cast<threadData*>(threadArg);
    int n=data->n;
    int myID=data->myID;
    int totalThreads=data->p;
    int v=data->v;
    static double tol=1e-19;
    if(v==0){

    int rows_per_thread =n/totalThreads;//распределяем
    int remainder = n % totalThreads;//как-то балансируем
    int first_row = myID * rows_per_thread + std::min(myID, remainder);
    int last_row = first_row + rows_per_thread - 1 + (myID < remainder ? 1 : 0);

    //--
    for (int step = 0; step < n; step++)
    {
        
           
            if((*data->globalPivotRow)<step &&(*data->globalPivotCol)<step){
            (*data->globalMaxVal) =((data->A)[step][step]);
            (*data->globalPivotRow)=step;
            (*data->globalPivotCol)=step;
            }
            
        
        
        int localPivotRow=-1;
        int localPivotCol=-1;//тк неизвестно лежит ли в подматрице [(step,step) ...]
        //pivot find loc//
        double localMaxVal = -1;
        for (int i = first_row; i <= last_row; i++)
        {
            for (int j = step; j < n; j++)
            {
                double current = std::fabs((data->A)[i][j]);
                if (current > localMaxVal &&i>=step &&j>=step)
                {
                    localMaxVal = current;
                    localPivotRow = i;
                    localPivotCol = j;

                }
            }
        }
        


        pthread_mutex_lock(&mutex);
        if (localMaxVal > (*data->globalMaxVal))
        {

            (*data->globalMaxVal) = localMaxVal;
            (*data->globalPivotRow) = localPivotRow;
            (*data->globalPivotCol) = localPivotCol;
        }

        pthread_mutex_unlock(&mutex);
        await(totalThreads);
        if(myID==0)
        {
            if ((*data->globalMaxVal)< tol){
                throw std::runtime_error("Matrix is close to singular.");
                return nullptr;
            }
            if ((*data->globalPivotRow) != step) {
                std::swap((data->A)[step], (data->A)[(*data->globalPivotRow)]);
                std::swap((data->inv)[step], (data->inv)[(*data->globalPivotRow)]);
                (*data->globalPivotRow) = step;
            }
            if ((*data->globalPivotCol) != step) {
                for (int i = 0; i < n; ++i) {
                    std::swap((data->A)[i][step], (data->A)[i][(*data->globalPivotCol)]);
                }
                std::swap((data->columnOrder)[step], (data->columnOrder)[(*data->globalPivotCol)]);
                (*data->globalPivotCol) = step;
            }
            //Normalize//
            double pivotVal = (data->A)[step][step]; //Доступ не в свою часть, но строго один
            for (int j = step; j < n; ++j) {
                (data->A)[step][j] /= pivotVal;
            }
            for (int j = 0; j < n; ++j) {
                (data->inv)[step][(data->columnOrder)[j]] /= pivotVal;
            }
        }
        await(totalThreads);
        //Eliminate part//
        for (int i = first_row; i <= last_row; ++i) {
            if (i == step)
                continue;
            double factor = (data->A)[i][step];
            for (int j = step; j < n; ++j) {
                (data->A)[i][j] -= factor * (data->A)[step][j];
            }
            for (int j = 0; j < n; ++j) {
                (data->inv)[i][(data->columnOrder)[j]] -= factor * (data->inv)[step][(data->columnOrder)[j]];
            }
        }

        await(totalThreads);

    }
    if(myID==0){
        std::vector<double> undo(n);
        for(int q=0;q<n;q++){
            undo[(data->columnOrder)[q]]=q;
        }
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++)
                (data->A)[i][j]=(data->inv)[undo[i]][j];
        }
    }
    return 0;
}
    if(v==1){//распределение по столбцам
        int col_per_thread =n/totalThreads;//распределяем
        int remainder = n % totalThreads;//как-то балансируем
        int first_col = myID * col_per_thread + std::min(myID, remainder);
        int last_col = first_col + col_per_thread - 1 + (myID < remainder ? 1 : 0);    
        for (int step = 0; step < n; step++)
    {
        
           
            if((*data->globalPivotRow)<step &&(*data->globalPivotCol)<step){
            (*data->globalMaxVal) =((data->A)[step][step]);
            (*data->globalPivotRow)=step;
            (*data->globalPivotCol)=step;
            }
            
        
        
        int localPivotRow=-1;
        int localPivotCol=-1;//тк неизвестно лежит ли в подматрице [(step,step) ...]
        //pivot find loc//
        double localMaxVal = -1;
        for (int i = step; i <= n-1; i++)
        {
            for (int j = first_col; j <= last_col; j++)
            {
                double current = std::fabs((data->A)[i][j]);
                if (current > localMaxVal &&i>=step &&j>=step)
                {
                    localMaxVal = current;
                    localPivotRow = i;
                    localPivotCol = j;

                }
            }
        }
        


        pthread_mutex_lock(&mutex);
        if (localMaxVal > (*data->globalMaxVal))
        {
            (*data->globalMaxVal) = localMaxVal;
            (*data->globalPivotRow) = localPivotRow;
            (*data->globalPivotCol) = localPivotCol;
        }

        pthread_mutex_unlock(&mutex);
        await(totalThreads);
        if(myID==0)
        {
            if ((*data->globalMaxVal)< tol){
                throw std::runtime_error("Matrix is close to singular.");
                return nullptr;
            }
            if ((*data->globalPivotRow) != step) {
                std::swap((data->A)[step], (data->A)[(*data->globalPivotRow)]);
                std::swap((data->inv)[step], (data->inv)[(*data->globalPivotRow)]);
                (*data->globalPivotRow) = step;
            }
            if ((*data->globalPivotCol) != step) {
                for (int i = 0; i < n; ++i) {
                    std::swap((data->A)[i][step], (data->A)[i][(*data->globalPivotCol)]);
                }
                std::swap((data->columnOrder)[step], (data->columnOrder)[(*data->globalPivotCol)]);
                (*data->globalPivotCol) = step;
            }
            //Normalize//
            double pivotVal = (data->A)[step][step]; //Доступ не в свою часть, но строго один
            for (int j = step; j < n; ++j) {
                (data->A)[j][step] /= pivotVal;
            }
            for (int j = 0; j < n; ++j) {
                (data->inv)[j][(data->columnOrder)[step]] /= pivotVal;
            }

        }
        await(totalThreads);
        //Eliminate part//
       for(int i= first_col; i<=last_col; ++i){// i- номер текущего СТОЛБЦА
        if(i == step)continue;
        double factor = (data->A)[step][i];
        //вычитаем столбцы//
        for(int j=0; j<n; ++j ){//j--номер СТРОКИ
            (data->A)[j][i]-=factor*(data->A)[j][step];
            (data->inv)[j][(data->columnOrder)[i]]-=factor * (data->inv)[j][(data->columnOrder)[step]];
        }
       }
        await(totalThreads);
    }
    if(myID==0){
        std::vector<double> undo(n);
        for(int q=0;q<n;q++){
            undo[(data->columnOrder)[q]]=q;
        }
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++)
                (data->A)[i][j]=(data->inv)[undo[i]][j];
        }
    }


    
    }
    return 0;
}

//Возвращаем время если удалось решить и число<0 если не удалось
double Solve(int n,int totalThreads,std::vector<std::vector<double>>& A,
    std::vector<std::vector<double>>& inv,std::vector<int> &columnOrder, std::vector<int> &undo, int v=0){

    pthread_t threadPool[totalThreads];//массив указателей на потоки
    std::vector<threadData>Data;//массив данных для потоков
    Data.reserve(totalThreads);
    inv.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        inv[i][i] = 1.0;

    // Initialize shared data
    
    for (int i = 0; i < n; ++i)
        columnOrder[i] = i;
    //Intialize shared memory
    double ptr=-1;
    int ptr1=-1;
    int ptr2=-1;
    //Create threads
    for (int i = 0; i < totalThreads; i++)
    {
        Data.emplace_back(n,i,totalThreads,A,inv,columnOrder,&ptr,&ptr1,&ptr2,v);
    }

    std::chrono::duration<double> elapsed;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < totalThreads; i++)
        if (pthread_create(&threadPool[i], 0, SolveProcess, &Data[i]))
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
