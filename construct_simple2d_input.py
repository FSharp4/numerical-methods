import copy
import math
from typing import Tuple

import conductor_simulation
import matrix
from matrix import Matrix


def process_conductor_nodes(mat: Matrix, volt: int) -> Tuple[Matrix, Matrix]:
    """
    Does two things:
    1. Obscures nodes in a conductor simulation matrix that correspond to being inside a conductor at a constant voltage
    level.
    2. Constructs a boolean matrix where only node indices that are part of the conductor dirichlet boundary in the
    conductor_matrix are set as true.

    This method only obscures such nodes for one conductor voltage level at a time. Nodes along the boundary of the
    conductor are not obscured.

    To "obscure" a node in this instance refers to setting its voltage to NaN.
    :return:
    :param volt: Constant voltage level inside desired conductor
    :param mat: Conductor simulation matrix
    :return: Obscured conductor matrix and dirichlet boundary matrix
    """

    mask = matrix.mask_int_equals(mat, volt)
    m = len(mask)
    n = len(mask[0])
    obscure_matrix = copy.deepcopy(mat)
    dirichlet_matrix = matrix.new_matrix(m, n, init_val=math.nan)

    for i in range(m):
        for j in range(n):
            if not mask[i][j]:
                continue

            boundary_count = 0

            if i > 0:
                if not mask[i - 1][j]:
                    boundary_count += 1

            if j > 0:
                if not mask[i][j - 1]:
                    boundary_count += 1

            if i < m-1:
                if not mask[i + 1][j]:
                    boundary_count += 1

            if j < n-1:
                if not mask[i][j + 1]:
                    boundary_count += 1

            if boundary_count == 0:
                obscure_matrix[i][j] = math.nan
            else:
                dirichlet_matrix[i][j] = volt

    for i in range(m):
        for j in range(n):
            if not math.isnan(obscure_matrix[i][j]):
                continue

            if i > 0 and j > 0:
                if math.isnan(dirichlet_matrix[i - 1][j - 1]):
                    obscure_matrix[i][j] = mat[i][j]
                    continue

            if i > 0 and j < n - 1:
                if math.isnan(dirichlet_matrix[i - 1][j + 1]):
                    obscure_matrix[i][j] = mat[i][j]
                    continue

            if i < m - 1 and j > 0:
                if math.isnan(dirichlet_matrix[i + 1][j - 1]):
                    obscure_matrix[i][j] = mat[i][j]
                    continue

            if i < m - 1 and j < n - 1:
                if math.isnan(dirichlet_matrix[i + 1][j + 1]):
                    obscure_matrix[i][j] = mat[i][j]
                    dirichlet_matrix[i][j] = volt
                    continue

    return obscure_matrix, dirichlet_matrix


def construct_labelled_dirichlet_matrix(dirichlet_matrices: list[Matrix]) -> Matrix:
    """
    Constructs a matrix where all dirichlet boundaries are given unique integer identifier labels.
    :param dirichlet_matrices: Individual boundary boolean matrices
    :return: Labelled Dirichlet Matrix
    """
    m = len(dirichlet_matrices[0])
    n = len(dirichlet_matrices[0][0])
    labelled_dirichlet_matrix = matrix.new_matrix(m, n, init_val=math.nan)
    for boundary in range(len(dirichlet_matrices)):
        # dirichlet_matrix = labelled_dirichlet_matrix[boundary]
        for i in range(m):
            for j in range(n):
                if not math.isnan(dirichlet_matrices[boundary][i][j]):
                    labelled_dirichlet_matrix[i][j] = dirichlet_matrices[boundary][i][j]

    return labelled_dirichlet_matrix


def construct_input_matrices(obscured_matrix: Matrix, labelled_dirichlet_matrix: Matrix, increment_size, write=False):
    m = len(obscured_matrix)
    n = len(obscured_matrix[0])
    position_matrix = matrix.new_matrix(m * n, 3, init_val=-1)
    position_mapping = matrix.new_matrix(m, n, init_val=-1)
    dirichlet_position_matrix = []
    triangle_matrix = matrix.new_matrix((m - 1) * (n - 1) * 2, 4, init_val=-1)
    ctr = 1
    for i in range(m):
        for j in range(n):
            if math.isnan(obscured_matrix[i][j]):
                continue

            position_matrix[ctr - 1] = [ctr, i * increment_size, j * increment_size]
            position_mapping[i][j] = ctr
            if not math.isnan(labelled_dirichlet_matrix[i][j]):
                dirichlet_position_matrix.append([ctr, labelled_dirichlet_matrix[i][j]])

            ctr += 1

    dirichlet_position_matrix = matrix.col_sort(dirichlet_position_matrix, 1)
    dpm_ctr = ctr
    ctr = 0
    for i in range(m - 1):
        for j in range(n - 1):
            n1 = position_mapping[i][j + 1]
            n2 = position_mapping[i][j]
            n3 = position_mapping[i + 1][j]
            n4 = position_mapping[i + 1][j + 1]
            if n1 == -1 or n2 == -1 or n3 == -1 or n4 == -1:
                continue
            triangle_matrix[ctr] = [n1, n2, n3, 0]
            ctr += 1
            triangle_matrix[ctr] = [n4, n1, n3, 0]
            ctr += 1

    if write:
        with open("file1.dat", "w") as file:
            ptr = 0
            for row in position_matrix:
                index = row[0]
                if index == -1:
                    ptr -= 1
                    break

                x = '{0:.3f}'.format(row[1])
                y = '{0:.3f}'.format(row[2])

                s = f"{index}\t{x}\t{y}\n"
                file.write(s)
                ptr += 1

            integer_length = len(str(ptr))
            file.write("\n")
            template = "{" + f":<{integer_length}d" + "}"
            for row in triangle_matrix:
                n1 = row[0]
                if n1 == -1:
                    break

                n2 = row[1]
                n3 = row[2]
                file.write(f"{template.format(n1)} {template.format(n2)} {template.format(n3)}\t0.000\n")

            file.write("\n")
            for row in dirichlet_position_matrix:
                file.write(f"{row[0]}\t{'{0:.3f}'.format(row[1])}\n")

    return position_matrix, triangle_matrix, dirichlet_position_matrix, position_mapping


if __name__ == "__main__":
    conductor_matrix = conductor_simulation.simulate_example_ECSE543_A1_Q3()
    conductor_matrix, high_dirichlet_matrix = process_conductor_nodes(conductor_matrix, 110)
    _, low_dirichlet_matrix = process_conductor_nodes(conductor_matrix, 0)
    construct_input_matrices(conductor_matrix, labelled_dirichlet_matrix=construct_labelled_dirichlet_matrix(
        [low_dirichlet_matrix, high_dirichlet_matrix]), increment_size=0.02, write=True)
