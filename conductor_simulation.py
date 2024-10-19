import math
import warnings
from copy import deepcopy

from matrix import Matrix, new_square_matrix

# TODO: Evaluate for submission 3A


def simulate_example_ECSE543_A1_Q3(h_factor: int = 0, w: float = 1.1, use_gauss_seidel=False) -> Matrix:
    return simulate_example_ECSE543_A1_Q3_monitored(h_factor, w)[0]


def locate_sentinel_node_potential(m: Matrix) -> float:
    n = len(m)
    h_coefficient = int((n - 1) / 5)
    h = 0.02 / h_coefficient
    corresponding_i = int(0.06 / h)
    corresponding_j = int(0.04 / h)
    return m[corresponding_i][corresponding_j]


def recover_node_spacing(m: Matrix) -> float:
    n = len(m)
    h_coefficient = int((n - 1) / 5)
    return 0.02 / h_coefficient


def simulate_example_ECSE543_A1_Q3_monitored(h_factor: int = 0, w: float = 1.1, use_gauss_seidel=False) \
        -> tuple[Matrix, int]:
    """
    Simulates the potential over open air of a system cross-section consisting of two conductors:
    - An outer square conductor of length/height 0.2 m, grounded
    - A centred rectangular conductor of height 0.04 m, width 0.08 m, held at 110 V
    :param h_factor:
    :param w:
    :return:
    """
    if w == 1:
        use_gauss_seidel = True
    h_coefficient: int = int(math.pow(2, h_factor))
    n = 5 * h_coefficient + 1
    h = 0.02 / h_coefficient
    M = new_square_matrix(n)
    for i in range(n):
        for j in range(n):
            if within_inner_conductor(h, i, j):
                M[i][j] = 110

    max_residual = 0.0
    index_queue: Matrix = construct_diagonal_visitation_queue(n)
    ctr = 0
    while ctr < len(index_queue):
        i, j = index_queue[ctr]
        if M[i][j] != 0:
            del [index_queue[ctr]]
        else:
            ctr += 1

    def calculate_potential(i: int, j: int) -> float:
        # Assume NOT WITHIN CONDUCTOR REGION
        a, b, d = 0, 0, 0
        if i < n - 1:
            a += M[i - 1][j] + M[i + 1][j]
            d += 2
        if j < n - 1:
            b += M[i][j - 1] + M[i][j + 1]
            d += 2

        if use_gauss_seidel:
            return 1 / d * (a + b)
        else:
            return (1 - w) * M[i][j] + w / d * (a + b)

    continue_condition: bool = True
    continue_threshold = 10e-5
    error_threshold = 10e3
    iteration = 0
    while continue_condition:
        max_residual = 0.0
        for index_pair in index_queue:
            i, j = index_pair
            updated_potential = calculate_potential(i, j)
            residual = abs(updated_potential - M[i][j])
            M[i][j] = updated_potential
            if max_residual < residual:
                max_residual = residual

        if max_residual < continue_threshold:
            continue_condition = False
        elif max_residual > error_threshold:
            raise ArithmeticError(f"SOR Diverges with w={w}")

        iteration += 1

    return M, iteration


def construct_diagonal_visitation_queue(n) -> Matrix:
    index_queue = []
    for norm in range(2, 2 * n - 1):
        i = min(n, norm) - 1
        j = norm - i
        while i >= 1 and j < n:
            index_queue += [[i, j]]
            i -= 1
            j += 1

    return index_queue

def within_inner_conductor(h, i, j):
    return i * h >= 0.06 and j * h >= 0.08