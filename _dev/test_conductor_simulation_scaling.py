import conductor_simulation
import matrix
from matrix import new_matrix

if __name__ == "__main__":
    # TODO: Evaluate for submission 3C
    w = 1.3
    n = 5
    test_matrix = new_matrix(n, 3)
    for h_factor in range(1, 1 + n):
        print(f"Running H Factor {h_factor}")
        M, iterations = conductor_simulation.simulate_example_ECSE543_A1_Q3_monitored(w=w, h_factor=h_factor)
        test_matrix[h_factor - 1] = [conductor_simulation.recover_node_spacing(M),
                                     iterations,
                                     conductor_simulation.locate_sentinel_node_potential(M)]

    matrix.write(test_matrix, "Simulation Scaling.csv")
"""
H Factor 1: iterations: 58, sentinel potential (0.06, 0.04) = 41.081536300442416
H Factor 2: iterations: 212, sentinel potential (0.06, 0.04) = 40.68872313203535
H Factor 3: iterations: 740, sentinel potential (0.06, 0.04) = 40.54257828991054
H Factor 4: iterations: 2511, sentinel potential (0.06, 0.04) = 40.46793718724143
"""