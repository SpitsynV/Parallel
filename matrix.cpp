#include "matrix.h"

int readMatrixFromFile(const std::string& filename, std::vector<std::vector<double>>& A, int n) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        return 0;
    }
    if (n == 0 || A[0].size() != n)
        return 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if(!(file >> A[i][j])){
                return 0;
            }
        }
    }
    file.close();
    return 1;
}

void printMatrix(const std::vector<std::vector<double>>& A, int m) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cout << std::setw(10) << std::setprecision(3) << std::scientific << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double>& vec, int m) {
    for (int i = 0; i < m; ++i) {
        std::cout << std::setw(10) << std::setprecision(3) << std::scientific << vec[i] << " ";
    }
    std::cout << std::endl;
}

double f(int k, int n, int i, int j) {
    switch (k) {
    case 1:
        return n - std::max(i, j) + 1;
    case 2:
        return std::max(i, j);
    case 3:
        return std::abs(i - j);
    case 4:
        return 1.0 / (i + j - 1);
    default:
        throw std::invalid_argument("Error: Wrong formula number");
    }
}
void initializeMatrix(std::vector<std::vector<double>>& A, int k, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            A[i][j]=f(k,n,i+1,j+1);//initially enumerate from 1 to n
        }
    }
}
