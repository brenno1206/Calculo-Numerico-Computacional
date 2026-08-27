#include "LSSolver.hpp"
#include <iostream>
#include <format>
#include <cmath>

LSSolver::LSSolver(const std::vector<std::vector<double>>& pAugMat) : augMat(pAugMat) {}

void LSSolver::printMatrix() {
    printMatrix(augMat);
}

void LSSolver::printMatrix(const std::vector<std::vector<double>>& matrix) {

    std::cout << " --     --\n";
    for (auto row : matrix) {
        std::cout << '|';
        for(auto n : row) {
            std::cout << std::format(" {} ", n);
        }
        std::cout << "|\n";
    }
    std::cout << " --     --";
}

double LSSolver::frobeniusForm() {
    return frobeniusForm(augMat);
}

double LSSolver::frobeniusForm(const std::vector<std::vector<double>> &matrix) {
    double total{0.0};
    for (auto row : matrix) {
        for(auto n : row) {
            total += (n * n);
        }
    }
    return(std::sqrt(total));
}

std::vector<std::vector<double>> LSSolver::zerosMatriz(int rows, int cols) {
    std::vector<std::vector<double>> matrix;
    for(int i{0}; i < rows; i++) {
        matrix.push_back({});
        for(int j{0}; j < cols; j++) {
        matrix.at(i).push_back(0);
        }
    }
    return matrix;
}

std::vector<std::vector<double>>  LSSolver::coefMatrix() {
    int rows = augMat.size();
    int cols = augMat.at(0).size();

    auto MC = LSSolver::zerosMatriz(rows, cols);

    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols - 1; j++) {
            MC[i][j] = this->augMat[i][j];
        }
    }

    return MC;
}

void LSSolver::pivotCheckAndSwap(int rowIndex) {
    auto rows = augMat.size();

    if(augMat[rowIndex][rowIndex] == 0) {
        for(int j{rowIndex+1}; j < rows; j++) {
            if(augMat[j][rowIndex] != 0) {
                auto aux = augMat[rowIndex];
                augMat[rowIndex] = augMat[j];
                augMat[j] = aux;
                return;
            }
        }
    }
    // LANCAR ERRO ValueError
    // Matriz singular! Nao eh possivel aplicar o metodo (pivo zero)
}

int gaussianElimination() {
    
}
