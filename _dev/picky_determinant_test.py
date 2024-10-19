import matrix

if __name__ == "__main__":
    M = matrix.read("../matrices/Test_matrix_4.csv")
    print(f"4: {matrix.det(M)}")
