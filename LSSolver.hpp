#pragma once
#include <vector>
#include <optional>
#include <tuple>
#include <string>

class LSSolver {
private: 
    std::vector<std::vector<double>> augMat;

    void pivotCheckAndSwap(int row_index);
    
public:
    LSSolver(const std::vector<std::vector<double>>& pAugMat);

    void printMatrix();
    void printMatrix(const std::vector<std::vector<double>>& matrix);

    double frobeniusForm();
    double frobeniusForm(const std::vector<std::vector<double>> &matrix);

    static std::vector<std::vector<double>> zerosMatriz(int rows, int cols);

    std::vector<std::vector<double>> coefMatrix();

    int gaussianElimination();

    std::tuple<double, std::vector<std::vector<double>>> determinant();

    void gaussJordanElimination();

    std::vector<std::vector<double>> inverse();

    double conditionNumber();

    std::string validateSystem(double conditionThreshold=100);

};