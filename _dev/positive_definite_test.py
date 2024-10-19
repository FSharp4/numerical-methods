import matrix

# TODO: Evaluate for submission 1B/C

if __name__ == "__main__":
    for i in range(2, 12):
        M = matrix.read(f"../matrices/Test_matrix_{i}.csv")
        print(f"{i}: Positive-Definite: {matrix.is_positive_definite(M)}")
        matrix.pprint_matrix(matrix.chol(M))
