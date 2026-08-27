#pragma once
#include <vector>
#include <optional>
#include <tuple>
#include <string>

/**
 * @brief Classe que resolve sistemas lineares por gauss jordan.
 * 
 * Alem do resultado do sistema, ainda gera o determinante, a inversa, norma de frobenius e numero de condicao.
 */
class LSSolver {
private: 
    std::vector<std::vector<double>> augMat;

    /**
     * @brief Verifica se o pivo na linha 'rowIndex' eh zero.
     * Lanca um erro se a matriz for singular.
     * @param rowIndex (int)
     */
    void pivotCheckAndSwap(int rowIndex);
    
public:
    /**
     * @brief Construtor da classe LSSolver.
     * @param pAugMart (vetor de vetores)
     */
    LSSolver(const std::vector<std::vector<double>>& pAugMat);

    /**
     * @brief Imprime a matriz da classe no terminal.
     */
     void printMatrix();
    
     /**
     * @brief Imprime a matriz no terminal.
     * @param matrix (vetor de vetores)
     */
    void printMatrix(const std::vector<std::vector<double>>& matrix);

    /**
     * @brief Calcula a norma de Frobenius da matriz da classe.
     * @return FrobeniusNumber (double)
     */
    double frobeniusForm();

    /**
     * @brief Calcula a norma de Frobenius da matriz fornecida.
     * @param matrix (vetor de vetores)
     * @return FrobeniusNumber (double)
     */
    double frobeniusForm(const std::vector<std::vector<double>> &matrix);

    /**
     * @brief Cria uma matriz preenchida por zeros.
     * @param rows numero de linhas (int)
     * @param cols numero de colunas (int)
     * @return matrix (vetor de vetores)
     */
    static std::vector<std::vector<double>> zerosMatriz(int rows, int cols);

    /**
     * @brief Recupera a matriz de coeficientes.
     * @return matrix (vetor de vetores)
     */
    std::vector<std::vector<double>> coefMatrix();

    /**
     * @brief Execura o metodo de eliminacao de gauss para triangular a matriz aumentada.
     * @return numSwaps (int)
     */
    int gaussianElimination();

    /**
     * @brief Calcula o determinante de augMat por meio do metodo de gauss.
     * @return tupla de determinante (double) e matriz triangular (vetor de vetores)
     */
    std::tuple<double, std::vector<std::vector<double>>> determinant();

    /**
     * @brief Aplica Gauss-Jordan para trnasformar a matriz aumentada em form reduzida por linhas.
     */
    void gaussJordanElimination();

    /**
     * @brief Calcula a inversa da matriz de coeficientes usando Gauss-Jordan.
     * @return matriz inversa (vetor de vetores)
     * @exception ValueError: se a matriz nao for inversivel
     */
    std::vector<std::vector<double>> inverse();

    /**
     * @brief Calcula o numero de condicao da matriz de coeficientes com base na norma de Frobenius e do valor da inversa.
     * @return numero de condicao (double)
     * @exception ValueError: se a matriz nao for inversivel
     */
    double conditionNumber();

    /**
     * @brief Verifica se o sistema linear eh bem posto e estavel numericamente.
     * 
     * Criterios:
     * 
     * - Matriz quadrada e determinante diferente de zero
     * 
     * - Numero de condicao nao excessivamente alto
     * @param conditionThreshold: limite aceitavel para o numero de condicao (default: 100)
     * @exception ValueError: se o sistema nao for bem posto
     * 
     */
    std::string validateSystem(double conditionThreshold=100);

};