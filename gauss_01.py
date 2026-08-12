from LS_Lib.solder import LSSolver
if __name__ == "__main__":
    # Matriz aumentada: A | b
    matrix = [
        [1, 1, 2, 9],
        [2, 4, -3, 1],
        [3, 6, -5, 0]
    ]

    print("Analisando o sistema linear...\n")

    lss = LSSolver(matrix)