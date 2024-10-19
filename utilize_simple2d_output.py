import math

import conductor_simulation
import construct_simple2d_input
import matrix

S_CONJ_SQUARE = [[1, -0.5, 0, -0.5], [-0.5, 1, -0.5, 0], [0, -0.5, 1, -0.5], [-0.5, 0, -0.5, 1]]

EPS_0 = 8.854e-12

if __name__ == "__main__":
    increment = 0.02
    conductor_matrix = conductor_simulation.simulate_example_ECSE543_A1_Q3()
    m = len(conductor_matrix)
    n = len(conductor_matrix[0])
    conductor_matrix, high_dirichlet_matrix = construct_simple2d_input.process_conductor_nodes(conductor_matrix, 110)
    _, low_dirichlet_matrix = construct_simple2d_input.process_conductor_nodes(conductor_matrix, 0)
    labelled_dirichlet_matrix = construct_simple2d_input.construct_labelled_dirichlet_matrix(
        [low_dirichlet_matrix, high_dirichlet_matrix])
    position_matrix, triangle_matrix, dirichlet_position_matrix, position_mapping = construct_simple2d_input.construct_input_matrices(
        conductor_matrix, labelled_dirichlet_matrix, increment_size=0.02)
    simple2d_voltages = matrix.mat_to_vec(matrix.read("Voltages.csv"))

    voltage = math.nan
    for idx in range(len(simple2d_voltages)):
        if position_matrix[idx][1] == 0.06 and position_matrix[idx][2] == 0.04:
            print(idx)
            print(conductor_matrix[3][2])
            print(simple2d_voltages[idx])
            break

    num_nodes = len(simple2d_voltages)
    # voltage_matrix = matrix.new_matrix(num_nodes, 4)
    for i in range(num_nodes):
        # voltage_matrix[i] = [position_matrix[0], position_matrix[1], position_matrix[2], simple2d_voltages[i]]
        x_index = int(position_matrix[i][1] // increment)
        y_index = int(position_matrix[i][2] // increment)
        conductor_matrix[x_index][y_index] = simple2d_voltages[i]

    system_energy = 0
    for i in range(m - 1):
        for j in range(n - 1):
            idx = [position_mapping[i + 1][j], position_mapping[i][j], position_mapping[i][j + 1],
                   position_mapping[i + 1][j + 1]]
            if -1 in idx:
                continue

            v = matrix.transpose([[conductor_matrix[i + 1][j], conductor_matrix[i][j], conductor_matrix[i][j + 1],
                                   conductor_matrix[i + 1][j + 1]]])

            u = matrix.mat_to_vec(v)
            uSu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2] + u[3] * u[3] - u[1] * u[2] - u[2] * u[3] - u[3] * u[0] - u[
                0] * u[1]

            system_energy += 0.5 * EPS_0 * uSu

            # system_energy += 0.5 * EPS_0 * matrix.multiply(matrix.transpose(v), matrix.multiply(S_CONJ_SQUARE, v))

    system_energy *= 4  # Get all quadrants
    print(f"W = {system_energy} J")
    print(f"C/L = {2 * system_energy / (110 * 110)} F/m")
