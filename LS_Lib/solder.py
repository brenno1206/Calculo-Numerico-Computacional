import math

class LSSolver:
    def __init__(self, augMat):
        self.augMat = augMat

    def print_matrix(self, matrix=None, format_spec="7.2"):
        """
        Imprime a matriz no terminal no formato especificado
        ou no default 7.2. Caso nenhuma matriz seja fornecida,
        a funcao imprime self.augMat
        """
        if matrix is None:
            matrix = self.augMat

        for row in matrix:
            for value in row:
                print(f"{value:{format_spec}f}", end=" ")
            print()

    def frobenius_norm(self, matrix=None):
        """
        Calcula a norma de Frobenius da matriz fornecida.
        Se nenhuma matriz for fornecida, calcula para self.augMat
        
        :param matrix: matriz (lista de lista) opcional
        :return : norma de Frobenius (float)
        """
        if matrix is None:
            matrix = self.coef_matrix()

        total = 0.0
        for row in matrix:
            for value in row:
                total += value ** 2

        return math.sqrt(total)

    @staticmethod
    def zeros_matrix(rows, cols):
        """
        Cria uma matriz preenchida com zeros.

        :param rows: numero de linhas
        :param cols: numero de colunas
        :return: matriz de zeros (list of lists)
        """
        return [[0.0 for _ in range(cols)] for _ in range(rows)]

    def coef_matrix(self):
        """
        Recupera a matriz de coeficientes (sem a coluan de termos independentes)
        
        :return: Matriz de coeficientes (list of lists)
        """
        # Dimensoes da matriz aumentada
        rows = len(self.augMat)
        cols = len(self.augMat[0])

        # Nova matriz de zeros para os coeficientes
        MC = self.zeros_matrix(rows, cols - 1)

        # Copiando APENAS os coeficientes
        for i in range(rows):
            for j in range(cols - 1):
                MC[i][j] = self.augMat[i][j]

        return MC

    def _pivot_check_and_swap(self, row_index):
        """
        Verifica se o pivo na linha 'row_index' eh zero.
        Se for, tenta trocar com uma linha abaixo que enha um pivo nao-nulo.
        Lanca um erro se a matriz for singular.
        """
        # numero de linhas da matriz
        n = len(self.augMat)

        if self.augMat[row_index][row_index] == 0:
            for j in range(row_index + 1, n):
                if self.augMat[j][row_index] != 0:
                    # troca as linahs quando encontra um pívo valido (!= 0)
                    self.augMat[row_index], self.augMat[j] = self.augMat[j], self.augMat[row_index]
                    return
            # lanca erro porque verificou todas as linhas da coluna row_index e todas elas sao 0.
            raise ValueError("Matriz singular! Nao eh possivel aplicar o metodo (pivo zero).")

    def gaussian_elimination(self):
        """
        Executa o metodo de Eliminacao de Gauss para triangular a matriz aumentada
        Modifica self.augMat, in-place.
        """
        n = len(self.augMat)
        m = len(self.augMat[0])
        num_swaps = 0

        for i in range(n):
            swapped = self._pivot_check_and_swap(i)
            if swapped:
                num_swaps +=1

            pivot = self.augMat[i][i]

            for j in range(i + 1, n):
                # divide pelo pivo pra calcular coeficiente correto
                coef = self.augMat[j][i] / pivot
                # O laco a seguir óde otimizar comecando de i, ja que antes eh zero
                for k in range(i, m):
                    self.augMat[j][k] -= self.augMat[i][k]

        return num_swaps