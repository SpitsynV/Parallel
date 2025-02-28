#include <iostream>
#include <chrono>
#include "matrix.h"
#include "functions.h"
#include "task.h"
#include <pthread.h>




int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Usage: " << argv[0] << " n m k p filename" << std::endl; //проверь формат
        return 1;
    }

    int n = std::stoi(argv[1]);   	// matrix dimension
    int m = std::stoi(argv[2]);   	// amount of output values
    int k = std::stoi(argv[3]);   	// formula number
    int totalThreads = std::stoi(argv[4]);
    std::string filename="";
    if (m > n) {
        std::cerr << "Error:Number of values to be output is greater than the matrix dimension " << std::endl;
        return 1;
    }
    int err = 0;
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> inv;
    A.resize(n, std::vector<double>(n));
    inv.resize(n, std::vector<double>(n));
    if (k == 0) {
        filename = argv[5];
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
        time=Solve(n,totalThreads,A,inv);
        std::cout << "Inverse matrix:\n";
        printMatrix(inv,m);
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 0;
    }
    // time
    std::cout << "Time to solve: " << time << std::endl;




    if(k==0){
        readMatrixFromFile(filename, A, n);
    }else{
        initializeMatrix(A, k, n);
    }
    double normError = 0.0;
    normError=getNorm(A,inv);
    // results

    std::cout << "Норма невязки: " << std::scientific << normError << std::endl;

    return 0;
}


