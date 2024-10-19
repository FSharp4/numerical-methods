import conductor_simulation
import matrix
from matrix import new_matrix

if __name__ == "__main__":
    # TODO: Evaluate for submission 3B
    info_matrix = new_matrix(10, 3)
    for i in range(10):
        w = 1 + 0.1 * i
        M, iterations = conductor_simulation.simulate_example_ECSE543_A1_Q3_monitored(w=w)
        info_matrix[i] = [w, iterations, conductor_simulation.locate_sentinel_node_potential(M)]

    matrix.write(info_matrix, "Basic SOR Analysis Info.csv")
