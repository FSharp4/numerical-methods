import matrix

if __name__ == "__main__":
    for size in range(2, 12):
        M = matrix.read(f"Test_matrix_{size}.csv")
        matrix.pprint_matrix(M)
