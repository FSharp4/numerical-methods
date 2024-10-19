import time
import matrix

# TODO: Evaluate for submission 1B/C
import matrix_generator

if __name__ == "__main__":
    sizes = tuple(range(10, 150, 10))
    timing_matrix = matrix.new_matrix(len(sizes), 3)
    for i in range(len(sizes)):
        M = matrix_generator.generate_positive_definite_matrix(sizes[i])
        t1 = time.time()
        matrix.is_positive_definite(M)
        t2 = time.time()
        matrix.chol(M)
        t3 = time.time()
        timing_matrix[i] = (sizes[i], t3-t2, t2-t1)

    matrix.write(timing_matrix, "Cholesky_Timing.csv")
