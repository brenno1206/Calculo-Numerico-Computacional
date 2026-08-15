import math
from typing import Optional, Any

class LSSolver:
    def __init__(self, augMat: list[list[float]]) -> None:
        """
        Inicializa o solucionador de sistemas lineares.
        
        Args:
            augMat (list[list[float]]): A matriz aumentada [A | b] representando o sistema.
        """
        self.augMat = augMat

    def print_matrix(self, matrix: Optional[list[list[float]]] = None, format_spec: str = "7.2") -> None:
        """
        Imprime a matriz no terminal no formato especificado.
        
        Args:
            matrix (list[list[float]], opcional): Matriz a ser impressa. Se None, imprime self.augMat.
            format_spec (str): Especificação de formatação das casas decimais. Padrão "7.2".
        """
        if matrix is None:
            matrix = self.augMat

        # Itera sobre cada linha e cada valor, formatando de acordo com format_spec
        for row in matrix:
            for value in row:
                print(f"{value:{format_spec}f}", end=" ")
            print()

    def frobenius_norm(self, matrix: Optional[list[list[float]]] = None) -> float:
        """
        Calcula a norma de Frobenius da matriz fornecida.
        A norma de Frobenius é a raiz quadrada da soma dos quadrados de todos os elementos.
        
        Args:
            matrix (list[list[float]], opcional): Matriz alvo. Se None, usa a matriz de coeficientes.
            
        Returns:
            float: Valor numérico da norma de Frobenius.
        """
        if matrix is None:
            matrix = self.coef_matrix()

        total: float = 0.0
        for row in matrix:
            for value in row:
                total += value ** 2  # Soma o quadrado de cada elemento

        return math.sqrt(total)  # Retorna a raiz quadrada da soma

    @staticmethod
    def zeros_matrix(rows: int, cols: int) -> list[list[float]]:
        """
        Cria uma matriz preenchida com zeros (útil para alocação inicial de memória).

        Args:
            rows (int): Número de linhas.
            cols (int): Número de colunas.
            
        Returns:
            list[list[float]]: Matriz preenchida com zeros.
        """
        # Cria uma lista de listas com zeros usando list comprehension
        return [[0.0 for _ in range(cols)] for _ in range(rows)]

    def coef_matrix(self) -> list[list[float]]:
        """
        Recupera a matriz de coeficientes A (exclui a coluna de termos independentes b).
        
        Returns:
            list[list[float]]: Matriz quadrada dos coeficientes do sistema.
        """
        # Obtém as dimensões da matriz aumentada atual
        rows = len(self.augMat)
        cols = len(self.augMat[0])

        # Cria uma nova matriz de zeros com uma coluna a menos
        MC = self.zeros_matrix(rows, cols - 1)

        # Copia apenas os valores da matriz de coeficientes (ignorando a última coluna)
        for i in range(rows):
            for j in range(cols - 1):
                MC[i][j] = self.augMat[i][j]

        return MC

    def _pivot_check_and_swap(self, row_index: int) -> bool:
        """
        Verifica se o pivô na linha atual é zero. Se for, procura uma linha abaixo 
        para trocar (Pivoteamento Parcial).
        
        Args:
            row_index (int): Índice da linha (e coluna) atual do pivô.
            
        Returns:
            bool: True se uma troca de linhas foi feita, False caso contrário.
            
        Raises:
            ValueError: Se a matriz for singular (todos os pivôs da coluna são zero).
        """
        n = len(self.augMat)

        # Verifica se o pivô principal é zero
        if self.augMat[row_index][row_index] == 0:
            # Procura por linhas abaixo onde o pivô não seja zero
            for j in range(row_index + 1, n):
                if self.augMat[j][row_index] != 0:
                    # Faz o swap (troca) das linhas no próprio array
                    self.augMat[row_index], self.augMat[j] = self.augMat[j], self.augMat[row_index]
                    return True  # Retorna True informando que ocorreu uma troca
                    
            # Se procurou em tudo e tudo é 0, o sistema não tem solução única
            raise ValueError("Matriz singular! Não é possível aplicar o método (pivô zero).")
            
        return False  # Retorna False se não precisou trocar

    def gaussian_elimination(self) -> int:
        """
        Executa a Eliminação de Gauss clássica para triangular a matriz.
        Modifica a propriedade self.augMat in-place (direto na memória).

        Returns:
            int: O número de vezes que as linhas foram trocadas (importante para o sinal do determinante).
        """
        n = len(self.augMat)
        m = len(self.augMat[0])
        num_swaps = 0

        for i in range(n):
            # Verifica o pivô e atualiza o contador se houve troca
            swapped = self._pivot_check_and_swap(i)
            if swapped:
                num_swaps += 1

            pivot = self.augMat[i][i]

            for j in range(i + 1, n):
                # Calcula o multiplicador (linha alvo dividida pelo pivô)
                coef = self.augMat[j][i] / pivot
                
                for k in range(i, m):
                    # CORREÇÃO CRÍTICA AQUI: Faltava multiplicar o 'coef' antes de subtrair
                    self.augMat[j][k] -= coef * self.augMat[i][k]

        return num_swaps

    def determinant(self) -> tuple[float, list[list[float]]]:
        """ 
        Calcula o determinante da matriz de coeficientes usando Eliminação de Gauss.

        Returns:
            tuple[float, list[list[float]]]: Uma tupla contendo o determinante (float) 
                                             e a matriz triangular resultante.
        """
        coef = self.coef_matrix()
        
        # Salva a matriz original para restaurar depois (já que o Gauss modifica in-place)
        original_augMat = self.augMat
        
        # Substitui a self.augMat por uma cópia profunda da matriz de coeficientes
        self.augMat = [row[:] for row in coef]
        n = len(self.augMat)
        
        # Garantir que é uma matriz quadrada antes de tentar calcular determinante
        if any(len(row) != n for row in self.augMat):
            self.augMat = original_augMat
            raise ValueError("Matriz de coeficientes não é quadrada. Determinante indefinido.")

        # Executa a triangularização e pega o número de trocas de linha
        num_swaps = self.gaussian_elimination()

        # O determinante de uma matriz triangular é o produto da sua diagonal principal
        det: float = 1.0
        for i in range(n):
            det *= self.augMat[i][i]

        # Se o número de trocas de linhas for ímpar, o determinante inverte de sinal
        if num_swaps % 2 == 1:
            det = -det

        # Salva a matriz triangular
        triangular = [row[:] for row in self.augMat]

        # Restaura a matriz original aumentada de volta ao objeto
        self.augMat = original_augMat

        # Retorna o determinante e a matriz triangular!
        return det, triangular

    def gauss_jordan_elimination(self) -> None:
        """
        Executa a eliminação de Gauss-Jordan para encontrar a forma RREF (Reduced Row Echelon Form).
        Transfoma os coeficientes na matriz Identidade para descobrir os valores das variáveis.
        Modifica self.augMat in-place.
        """
        n = len(self.augMat)
        m = len(self.augMat[0])
        
        for i in range(n):
            self._pivot_check_and_swap(i)
            pivot = self.augMat[i][i]

            # Divide toda a linha pelo pivô, transformando o pivô em 1
            for k in range(m):
                self.augMat[i][k] /= pivot

            # Zera todos os outros elementos da coluna, acima e abaixo do pivô
            for j in range(n):
                if i != j:
                    coef = self.augMat[j][i]
                    for k in range(m):
                        self.augMat[j][k] -= coef * self.augMat[i][k]

    def inverse(self) -> list[list[float]]:
        """
        Calcula a inversa da matriz de coeficientes usando Gauss-Jordan [A | I] -> [I | A^-1].
        
        Returns:
            list[list[float]]: A matriz inversa de A.
            
        Raises:
            ValueError: Se a matriz não for invertível ou não for quadrada.
        """
        A = self.coef_matrix()
        n = len(A)

        # Verifica se A é quadrada
        if any(len(row) != n for row in A):
            raise ValueError("A matriz de coeficientes não é quadrada. Inversa indefinida.")

        # Verifica se tem inversa (Determinante != 0)
        det, _ = self.determinant()
        if det == 0:
            raise ValueError("Matriz não é invertível (determinante = 0).")

        # Cria a matriz Identidade
        I = [[float(i == j) for j in range(n)] for i in range(n)]

        # Junta as duas matrizes criando [A | I]
        augmented = [A[i] + I[i] for i in range(n)]

        # Salva e substitui
        original_augMat = self.augMat
        self.augMat = [row[:] for row in augmented]

        try:
            self.gauss_jordan_elimination()
        except ValueError:
            self.augMat = original_augMat
            raise

        # Extrai somente o lado direito da matriz (onde agora reside a inversa)
        inverse_mat = [row[n:] for row in self.augMat]

        # Restaura original
        self.augMat = original_augMat

        return inverse_mat

    def condition_number(self) -> float:
        """
        Calcula o número de condição (Kappa) da matriz usando a norma de Frobenius.
        Indica a sensibilidade do sistema a pequenos erros numéricos.

        Returns:
            float: O número de condição.
        """
        norm_A = self.frobenius_norm() 

        try:
            # Calcula a matriz inversa
            inverse_A = self.inverse()
        except ValueError as e:
            raise ValueError(f"Não é possível calcular o número de condições: {e}")

        # Calcula a norma da matriz inversa
        norm_inv_A = self.frobenius_norm(inverse_A)

        # O Número de Condição é ||A|| * ||A^-1||
        return norm_A * norm_inv_A

    def validate_system(self, condition_threshold: float = 1e2) -> dict[str, Any]:
        """
        Gera um relatório de saúde do sistema linear. Verifica existência de 
        solução e estabilidade numérica.

        Args:
            condition_threshold (float): Limite máximo tolerado para o número de condição.

        Returns:
            dict: Relatório com determinante, número de condição, avisos e booleanos de status.
            
        Raises:
            ValueError: Se o sistema for mal posto ou o determinante for zero.
        """
        report: dict[str, Any] = {}

        # 1. Verifica consistência estrutural
        A = self.coef_matrix()
        n = len(A)
        if any(len(row) != n for row in A):
            raise ValueError("Matriz de coeficientes não é quadrada. Sistema mal posto.")

        # 2. Verifica existência de solução única
        det, _ = self.determinant()
        report["determinant"] = det
        if det == 0:
            raise ValueError("Determinante nulo. Sistema não possui solução única (mal posto).")
        else:
            report["well_posed"] = True

        # 3. Analisa estabilidade numérica
        kappa = self.condition_number()
        report["condition_number"] = kappa

        if kappa > condition_threshold:
            report["sensitivity_warning"] = True
            report["message"] = (
                f"Sistema é sensível a perturbações.\n"
                f"Número de condição alto (k = {kappa:.2f} > {condition_threshold:.0f})."
            )
        else:
            report["sensitivity_warning"] = False
            report["message"] = (
                f"Sistema considerado numericamente estável.\n"
                f"k = {kappa:.2f}."
            )

        return report