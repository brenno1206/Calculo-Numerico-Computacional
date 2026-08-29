#include <iostream>
#include "LSSolver.hpp"

int main() {
    
    LSSolver lss{{{1, 1, 2, 9},
        {2, 4, -3, 1},
        {3, 6, -5, 0}}};
    lss.printMatrix();

    std::cout << "\nAnalisando sistema linear\n";

    std::cout << lss.validateSystem();

    std::cout << "Matriz aumentada original:\n";
    lss.printMatrix();

    lss.gaussJordanElimination();

    std::cout << "Matriz na forma reduzida:\n";
    lss.printMatrix();
    std::vector<double> solution;
    for(auto row : lss.augMat) {
        solution.push_back(row.back());
    }
    std::cout << "\n\n";

    for(auto answer : solution) {
        std::cout << std::format(" x^? ={:.4f} ", answer);
    }
    
    
    return 0;
}

// Rode com:
// g++ -std=c++20 code_13.cpp -o out.exe


/*
CONFIGURACAO DO VSCODE

Pelo arquivo de configurações (settings.json)  
1. Pressione Ctrl + Shift + P no VS Code.  
2. Digite e selecione: Preferences: Open User Settings (JSON) (ou abra as configurações do workspace).
3. Adicione ou altere a linha do padrão C++ para c++20 (ou c++23):  
    JSON"C_Cpp.default.cppStandard": "c++20"
4.Salve o arquivo e reinicie o VS Code (ou feche e abra a janela com Ctrl + Shift + P -> Developer: Reload Window).
*/