from LS_Lib.solder import LSSolver

if __name__ == "__main__":
    # Matriz aumentada: A | b
    matrix = [
        [1, 1, 2, 9],
        [2, 4, -3, 1],
        [3, 6, -5, 0]
    ]

    print("Analisando o sistema linear...\n")

    # Inicializa o nosso solver passando a matriz
    lss = LSSolver(matrix)

    try:
        # Etapa 1: Validar sistema 
        report = lss.validate_system()

        print("--- Avaliação do Sistema Linear Proposto ---")
        print(f"Determinante: {report['determinant']:.4f}")
        print(f"Número de condição (k): {report['condition_number']:.4f}")
        print(f"Aviso de sensibilidade: {report['sensitivity_warning']}")
        print(f"Mensagem: {report['message']}\n")

        # Etapa 2: Exibir a matriz aumentada original
        print("--- Matriz Aumentada Original ---")
        lss.print_matrix()
        print()

        # Etapa 3: Resolver o sistema com Gauss-Jordan
        lss.gauss_jordan_elimination()

        # Etapa 4: Exibir a matriz reduzida (RREF)
        print("--- Matriz na Forma Reduzida (RREF) ---")
        lss.print_matrix()
        print()

        # Etapa 5: Extrair e exibir a solução do sistema
        # A última coluna da matriz transformada armazena as respostas (b)
        solution = [row[-1] for row in lss.augMat]
        
        print("--- Solução do Sistema ---")
        for i, val in enumerate(solution):
            print(f"x^{i + 1} = {val:.4f}")

    except ValueError as e:
        print("\n[!] Erro durante a resolução do sistema:")
        print(str(e))