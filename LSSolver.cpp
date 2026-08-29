#include "LSSolver.hpp"
#include <iostream>
#include <format>
#include <cmath>

LSSolver::LSSolver(const std::vector<std::vector<double>>& pAugMat) : augMat(pAugMat) {}

void LSSolver::printMatrix() {
    printMatrix(augMat);
}

// void LSSolver::printMatrix(const std::vector<std::vector<double>>& matrix) {

//     std::cout << " --     --\n";
//     for (auto row : matrix) {
//         std::cout << '|';
//         for(auto n : row) {
//             std::cout << std::format(" {} ", n);
//         }
//         std::cout << "|\n";
//     }
//     std::cout << " --     --";
// }

void LSSolver::printMatrix(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty() || matrix[0].empty()) return;

    int cols = matrix[0].size();
    // Calcula o espaçamento interno exato (9 char por número + 1 espaço = 10 por coluna)
    std::string cap(cols * 10 + 1, ' ');

    std::cout << " -" << cap << "-\n";
    for (const auto& row : matrix) {
        std::cout << "| ";
        for(auto n : row) {
            // {:>9.2f} garante 9 caracteres de largura, alinhado à direita, com 2 decimais
            std::cout << std::format("{:>9.2f} ", n);
        }
        std::cout << "|\n";
    }
    std::cout << " -" << cap << "-\n";
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
    auto rows{augMat.size()};
    auto cols{augMat.at(0).size()};

    auto MC{LSSolver::zerosMatriz(rows, cols -1)};

    for(int i{0}; i < rows; i++) {
        for(int j = 0; j < cols - 1; j++) {
            MC[i][j] = this->augMat[i][j];
        }
    }

    return MC;
}

bool LSSolver::pivotCheckAndSwap(int rowIndex) {
    auto rows{augMat.size()};

    if(augMat[rowIndex][rowIndex] == 0) {
        for(int j{rowIndex+1}; j < rows; j++) {
            if(augMat[j][rowIndex] != 0) {
                auto aux = augMat[rowIndex];
                augMat[rowIndex] = augMat[j];
                augMat[j] = aux;
                return true;
            }
        }
    }
    return false;
    // LANCAR ERRO ValueError
    // Matriz singular! Nao eh possivel aplicar o metodo (pivo zero)
}

int LSSolver::gaussianElimination() {
    auto n{augMat.size()};
    auto m{augMat[0].size()};
    int numSwaps{0};

    for(int i{0}; i < n; i++) {
        bool swapped{LSSolver::pivotCheckAndSwap(i)};
        if(swapped) {
            numSwaps++;
        }
        
        double pivot{augMat[i][i]};

        for(int j{i+1}; j < n; j++) {
            double coef{augMat[j][i]  / pivot};
            for(int k{i}; k < m; k++) {
                augMat[j][k] -= coef * augMat[i][k];
            }
        }
    }
    return numSwaps;
}

std::tuple<double, std::vector<std::vector<double>>> LSSolver::determinant() {
    auto coef{coefMatrix()};
    auto originalAugMat{augMat};
    augMat = coef;

    auto n{augMat.size()};
    for(auto row : augMat) {
        if(row.size() != n) {
            augMat = originalAugMat;
            // RAISE ValueError
            // Matriz de coeficientes nao eh quadrada determinante indefinido.
        }
    }

    auto numSwaps{gaussianElimination()};

    double det{1};
    for(int i{0}; i < n; i++) {
        det *=augMat[i][i];
    }
    if(numSwaps % 2 == 1) {
        det = -det;
    }

    auto triangular = augMat;

    augMat = originalAugMat;

    return std::make_tuple(det, triangular);
}

void LSSolver::gaussJordanElimination() {
    auto n{augMat.size()};
    auto m{augMat[0].size()};

    for(int i{0}; i < n; i++) {
        // fazer try except pq lanca excessao
        pivotCheckAndSwap(i);

        auto pivot{augMat[i][i]};

        for(int k{0}; k < m; k++) {
            augMat[i][k] /= pivot;
        }

        for(int j{0}; j < n; j++) {
            if(i != j) {
                auto coef{augMat[j][i]};
                for(int k{0}; k < m; k++) {
                    augMat[j][k] -= coef * augMat[i][k];
                }
            }
        }
    }
}

std::vector<std::vector<double>> LSSolver::inverse() {
    auto A{coefMatrix()};
    auto n{A.size()};

    for(auto row : A) {
        if(row.size() != n) {
            // RAISE ValueError
            // Matriz de coeficientes nao eh quadrada determinante indefinido.
        }
    }

    auto [det, _]{determinant()};
    if(det == 0) {
        // RAISE ERROR
        // Matriz nao eh inversivel (determinante zero)
    }

    auto I{LSSolver::zerosMatriz(n,n)};
    for(int i {0}; i < n; i++) {
        I[i][i] = 1.0;
    }
    std::vector<std::vector<double>> augmented = A;
    for (int i{0}; i < n; ++i) {
       augmented[i].insert(augmented[i].end(), I[i].begin(), I[i].end());
    }
    auto original_augMat{this->augMat};
    augMat=augmented;

    // fazer try except pq pode lancar excessao
    gaussJordanElimination();

    std::vector<std::vector<double>> inverseMat;

    for (const auto& row : this->augMat) {
        inverseMat.push_back(std::vector<double>(row.begin() + n, row.end()));
    }

    this->augMat = original_augMat;

    return inverseMat;

}

double LSSolver::conditionNumber() {
    auto A{this->coefMatrix()};
    auto normA{this->frobeniusForm(A)};
    // try except pq pode lancar excessao
    auto inverseA{this->inverse()};

    auto normInvA{this->frobeniusForm(inverseA)};

    return normA * normInvA;
}

std::string LSSolver::validateSystem(double conditionThreshold) {
    std::string report{"---System analysis---\n"};

    auto A{this->coefMatrix()};
    auto n{A.size()};
    for(const auto& row : A) {
        if(n != row.size()){
            // Raise error 
            // matriz nao eh quadrada, sistema mal posto
        }
    }

    auto [det, _]{this->determinant()};
    report += std::format("Determinant={:.2f}\n", det);
    if(det == 0) {
        // RAISE ERROR,
        // distema nao possui solucao unica
    } else {
        report += "Well posed system\n";
    }

    auto kappa{this->conditionNumber()};
    report += std::format("Condition number={:.2f}\n", kappa);

    if(kappa > conditionThreshold) {
        report += std::format("System sensitive to disturbances! high condition number (k = {:.2f} > {:.2f})\n", kappa, conditionThreshold);
    } else {
        report += std::format("numerically stable system! (k = {:.2f})\n", kappa);
    }

    return report;
}