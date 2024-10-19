import time
import matrix
from finite_difference_mesh import generate_fd_analysis_circuit
from matrix import new_matrix

# TODO: Evaluate for submission 2A/2B

if __name__ == "__main__":
    # # import timeit
    timing_matrix = new_matrix(14, 4)
    for i in range(2, 16):
        t1 = time.time()
        c, pos_node = generate_fd_analysis_circuit(i)
        t2 = time.time()
        v = c.voltages()
        R = (10000 * v[pos_node]) / (10 - v[pos_node])
        t3 = time.time()
        timing_matrix[i - 2] = [i, R, t3 - t2, t2 - t1]

    matrix.write(timing_matrix, "Timing_data.csv")
