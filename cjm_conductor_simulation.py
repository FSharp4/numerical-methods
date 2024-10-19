from typing import Tuple

import conductor_simulation
import conjugate_gradient
import matrix


def within_inner_conductor(h, i, j):
    x = h * i
    y = h * j
    return 0.06 <= x <= 0.14 and 0.08 <= y <= 0.12


def on_outer_conductor(h, i, j):
    x = h * i
    y = h * j
    return x == 0 or x == 0.2 or y == 0 or y == 0.2


def next_to_conductor(h, i, j):
    x = h * i
    y = h * j
    if within_inner_conductor(h, i, j):
        return False

    if 0.06 <= x <= 0.14:
        return y == 0.06 or y == 0.14
    elif 0.08 <= y <= 0.12:
        return x == 0.04 or x == 0.16
    else:
        return False


CTR = 0

def construct_conductor_system():
    global CTR
    h = 0.02
    n = int(0.1 // h) + 1
    coordinate_map = matrix.new_matrix(n * n, 2, init_val=-1)  # Preallocate
    ctr = 0
    for i in range(1, n):
        for j in range(1, n):  # Starting the matrices at 1 removes the need to check for outside conductor
            if within_inner_conductor(h, i, j):  #'On' == 'In' here
                continue

            if i * h == 0.06 and j * h == 0.04:
                print(f"Sentinel CTR: {ctr}")
                CTR = ctr

            coordinate_map[ctr] = [i,j]
            ctr += 1

    coordinate_map = coordinate_map[:ctr]

    num_nodes = len(coordinate_map)
    system_matrix = matrix.new_square_matrix(num_nodes)
    b = matrix.new_vector(num_nodes)

    def find(i: int, j: int) -> int:
        for a in range(num_nodes):
            if coordinate_map[a] == [i, j]:
                return a

        return -1

    for index in range(num_nodes):
        entry = coordinate_map[index]
        e_ptr = 0
        if entry[0] != n - 1:
            # Not on rightside neumann
            e_ptr += 2
            if entry[0] - 1 > 0:  # Check if node should exist
                hit_index = find(entry[0] - 1, entry[1])
                if hit_index != -1:
                    system_matrix[hit_index][index] += 1

            if entry[0] + 1 <= n:  # Check if node should exist
                hit_index = find(entry[0] + 1, entry[1])
                if hit_index != -1:
                    system_matrix[hit_index][index] += 1

        if entry[1] != n - 1:
            # Not on topside neumann
            e_ptr += 2
            if entry[1] - 1 > 0:  # Check if node should exist
                hit_index = find(entry[0], entry[1] - 1)
                if hit_index != -1:
                    system_matrix[hit_index][index] += 1

            if entry[1] + 1 <= n:  # Check if node should exist
                hit_index = find(entry[0], entry[1] + 1)
                if hit_index != -1:
                    system_matrix[hit_index][index] += 1

        system_matrix[index][index] -= 4
        if next_to_conductor(h, entry[0], entry[1]):
            b[index] = -110

    matrix.pprint_matrix(coordinate_map)
    return system_matrix, b


# noinspection DuplicatedCode
def adjusted_conductor_matrix() -> Tuple[matrix.Matrix, matrix.Vector]:
    h = 0.02
    n = int(0.2 // h) + 1
    coordinate_map = matrix.new_matrix(n * n, 2, init_val=-1)  # Preallocate
    ctr = 0
    for i in range(1, n - 1):
        for j in range(1, n - 1):  # Starting the matrices at 1 removes the need to check for outside conductor
            if within_inner_conductor(h, i, j):  # 'On' == 'In' here
                continue

            if i*h == 0.06 and j*h == 0.04:
                print(f"Sentinel Index: {ctr}")

            coordinate_map[ctr] = [i, j]
            ctr += 1

    coordinate_map = coordinate_map[:ctr]

    b = matrix.new_vector(len(coordinate_map))

    num_nodes = len(coordinate_map)
    system_matrix = matrix.new_square_matrix(num_nodes)

    def find(i: int, j: int) -> int:
        for a in range(num_nodes):
            if coordinate_map[a] == [i, j]:
                return a

        return -1

    for index in range(num_nodes):
        entry = coordinate_map[index]
        e_ptr = 4
        if entry[0] - 1 > 0:  # Check if node should exist
            hit_index = find(entry[0] - 1, entry[1])
            if hit_index != -1:
                system_matrix[hit_index][index] += 1

        if entry[0] + 1 <= n:  # Check if node should exist
            hit_index = find(entry[0] + 1, entry[1])
            if hit_index != -1:
                system_matrix[hit_index][index] += 1

        if entry[1] - 1 > 0:  # Check if node should exist
            hit_index = find(entry[0], entry[1] - 1)
            if hit_index != -1:
                system_matrix[hit_index][index] += 1

        if entry[1] + 1 <= n:  # Check if node should exist
            hit_index = find(entry[0], entry[1] + 1)
            if hit_index != -1:
                system_matrix[hit_index][index] += 1

        system_matrix[index][index] -= e_ptr
        if next_to_conductor(h, entry[0], entry[1]):
            b[index] = 110  # / (h * h)

    # matrix.pprint_matrix(coordinate_map)
    # system_matrix = matrix.scalar_multiply(system_matrix, 1 / (h * h))
    return system_matrix, b




if __name__ == "__main__":
    mat, b = construct_conductor_system()
    matrix.pprint_matrix(mat)
    matrix.illustrate_sparseness(mat)
    print(matrix.is_positive_definite(mat, use_cholesky=True))

    matT = matrix.transpose(mat)
    matTmat = matrix.multiply(matT, mat)
    adj_b = matrix.mat_to_vec(matrix.multiply(matT, matrix.vec_to_mat(b)))

    V = matrix.chol_solve(matTmat, adj_b)
    V2 = conjugate_gradient.solve(matTmat, adj_b)
    results = [V, matrix.mat_to_vec(V2)]

    matrix.pprint_matrix(results)

    print(f"Value of Conjugate Gradient at (0.06, 0.04): {V2[11]}")
    print(f"Value of Cholesky Solution at (0.06, 0.04): {V[11]}")

    full_matrix, fullb = adjusted_conductor_matrix()
    full_matrix = matrix.scalar_multiply(full_matrix, -1)
    print(matrix.is_positive_definite(full_matrix))
    fV = matrix.chol_solve(full_matrix, fullb)
    fV2 = conjugate_gradient.solve(full_matrix, fullb)

    matrix.pprint_matrix([fV, matrix.mat_to_vec(fV2)])
    print(f"Value of Conjugate Gradient at (0.06, 0.04): {fV2[19]}")
    print(f"Value of Cholesky Solution at (0.06, 0.04): {fV[19]}")


    sor_results = conductor_simulation.simulate_example_ECSE543_A1_Q3(w=1.3)
    sentinel_sor = conductor_simulation.locate_sentinel_node_potential(sor_results)

    print(f"Value of SOR solution at (0.06, 0.04): {sentinel_sor}")
