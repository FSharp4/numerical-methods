import time

import matrix
from circuit import OptimizedCircuit
from finite_difference_mesh import generate_fd_analysis_circuit
from matrix import new_matrix

if __name__ == "__main__":
    # # import timeit
    timing_matrix = new_matrix(14, 5)
    for i in range(2, 16):
        t1 = time.time()
        c: OptimizedCircuit
        c, pos_node = generate_fd_analysis_circuit(i, optimized_circuit_class=True)
        # noinspection DuplicatedCode
        t2 = time.time()
        v = c.voltages()
        b: int = c.AYAT_b
        R = (10000 * v[pos_node]) / (10 - v[pos_node])
        t3 = time.time()
        timing_matrix[i - 2] = [i, R, t3 - t2, t2 - t1, b]

    matrix.write(timing_matrix, "Optimized_timing_data.csv")
