#include <iostream>
#include <chrono>
#include "matrix.h"
#include "functions.h"
#include "task.h"
#include <pthread.h>




int main(int argc, char* argv[]) {
    if (argc < 6 || argc > 7) {
        std::cerr << "Usage: " << argv[0] << " n m k p version filename" << std::endl; //проверь формат
        return 1;
    }

    int n = std::stoi(argv[1]);   	// matrix dimension
    int m = std::stoi(argv[2]);   	// amount of output values
    int k = std::stoi(argv[3]);   	// formula number
    int totalThreads = std::stoi(argv[4]);
    int v = std::stoi(argv[5]);
    if(v!=0 && v!=1){
        std::cerr << "Error:Invalid version(0 for row, 1 for  col)" << std::endl;
        return 1;
    }
    std::string filename="";
    if (m > n) {
        std::cerr << "Error:Number of values to be output is greater than the matrix dimension " << std::endl;
        return 1;
    }
    int err = 0;
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> inv;
    std::vector<int> columnOrder(n);
    std::vector<int> undo(n);
    A.resize(n, std::vector<double>(n));
    inv.resize(n, std::vector<double>(n));
    if (k == 0) {
        filename = argv[6];
        err = readMatrixFromFile(filename, A, n);
        if (!err) {
            std::cerr << "Error: Can't read from file " << filename << std::endl;
            return err;
        }
    }
    else {
        initializeMatrix(A, k, n);
    }

    std::cout << "Initial matrix A:" << std::endl;
    printMatrix(A, m);


    double time=-1;
    try {
        time=Solve(n,totalThreads,A,inv,columnOrder,undo,v);
        std::cout << "Inverse matrix:\n";
        printMatrix(A,m);
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 0;
    }
    // time
    std::cout << "Time to solve: " << time <<" seconds"<< std::endl;




    if(k==0){
        readMatrixFromFile(filename, inv, n);
    }else{
        initializeMatrix(inv, k, n);
    }
    
    /*
    double normError = 0.0;
    std::chrono::duration<double> elapsed;
    auto start = std::chrono::high_resolution_clock::now();

    normError=getNorm(A,inv);

    auto end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    auto t=elapsed.count();
    std::cout<<"Single thread norm time == "<<t<< std::endl;
    std::cout << "Норма невязки: " << std::scientific << normError << std::endl;
    */

    // results

   
   
    
    double ans=0;
    auto t2=Norm(n,totalThreads,A,inv,&ans);
    ans=sqrt(ans);
    std::cout << "Параллельная норма невязки: " <<ans<< std::endl;
    std::cout<< "Multithread norm time: "<<t2<<" seconds"<<std::endl;
    
    return 0;
}


