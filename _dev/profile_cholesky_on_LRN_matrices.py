import time

import matrix

if __name__ == "__main__":
    timing_matrix = matrix.new_matrix(14, 3)
    for i in range(2, 16):
        AYAT = matrix.read(f"AYAT{i}.csv")
        AJYE = matrix.mat_to_vec(matrix.read(f"AJYE{i}.csv"))
        T1 = time.time()
        matrix.chol_solve(AYAT, AJYE)
        T2 = time.time()
        matrix.chol_band_solve(AYAT, AJYE)
        T3 = time.time()
        timing_matrix[i-2] = [i, T2-T1, T3-T2]

    matrix.write(timing_matrix, "Isolated Timing Comparations.csv")
