# pip install -r requirements.txt
import numpy as np

def build_permutation_matriz(ipvt):
    """
    Controi a matriz de permutacao P
    Args:
        ipvt(list[int]): linhas trocadas
    Returns:
        list[list[double]]: matriz P
    """
    n = len(ipvt);
    P = np.eye(n) # cria uma matriz identidade

    for k in range(n-1):
        i = k
        j = ipvt[k]

        if i != j:
            P[[i, j]] = P[[j, i]] # troca as linhas

    return P

def lu_factor(A):
    """"
    Realiza a decomposicao LU na matriz A
    Args:
        A(list[list[double]]): matriz
    Returns:
        list[list[double]]: matriz L
        list[list[double]]: matriz U
        list[int]: ipvt
        int: info
    """
    A = A.copy()
    n = A.shape[0]

    ipvt = np.arange(n);
    info = 0;

    for k in range(n - 1):
        # pivoteamento parcial
        pivot_row = k + np.argmax(np.abs(A[k:, k]))
        ipvt[k] = pivot_row

        # troca de linhas
        if pivot_row != k:
            A[[k, pivot_row]] = A[[pivot_row, k]]

        # verificacao do pivo
        if A[k,k] == 0:
            info = k
            continue

        # calcula os multiplicadores
        A[k+1:, k] /= A[k,k]

        # atualiza as colunas seguintes
        for j in range(k + 1, n):
            A[k+1:, j] -= A[k, j] * A[k+1:, k]

    # verificacao do ultimo pivo
    if A[n-1, n-1] == 0:
        info = n-1

    # separa L e U
    L = np.tril(A, k=-1) + np.eye(n, dtype=A.dtype)
    U = np.triu(A)

    return L, U, ipvt, info

if __name__ == "__main__":
    A=np.array([
        [2, 3, 1],
        [4, 7, 7],
        [6, 18, 22]
    ], dtype=float)

    L, U, ipvt, info = lu_factor(A)

    print("Matriz A original:")
    print(A)

    print("\nL (triangular inferior)")
    print(L)

    print("\nU (triangular superior)")
    print(U)

    print("\nPivot indices:")
    print(ipvt)

    print("\nInfo:")
    print(info)

    P = build_permutation_matriz(ipvt)
    print("\nP @ A:")
    print(P@A)

    print("\nL @ U:")
    print(L@U)

    print("\nErro (norma de frobenius): ", np.linalg.norm((P@A)-(L@U)))
