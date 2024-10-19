"""
Matrix Generation

The goal is to generate real, symmetric, positive-definite matrices for usage with Choleski Decomposition
"""

import random

import matrix
from matrix import new_square_matrix
import sys

__HIGH__ = 5
__LOW__ = -5
__THRESHOLD__ = sys.float_info.epsilon * 256


def set_bounds(high, low):
    """
    Sets bounds on values generated from random generator.
    Values should be such that low < 0 < high; it will work for low < high, but will produce strange matrices.
    :param high:
    :param low:
    :return:
    """
    if low < high:
        raise ArithmeticError("high < low not allowed")
    global __HIGH__
    global __LOW__
    __HIGH__ = high
    __LOW__ = low


def set_seed(seed=None):
    """
    Sets rng seed (for testing/debug)
    :param seed:
    :return:
    """
    random.seed(seed)


def _next_random_nonzero(flag_int: bool = True):
    """
    Retrieves next random number
    :param flag_int:
    :return:
    """
    factor = 0
    if flag_int:
        while factor == 0:
            factor = random.randint(int(__LOW__), int(__HIGH__))
    else:
        while abs(factor) < __THRESHOLD__:
            factor = (random.random() - 0.5) * (__HIGH__ - __LOW__)

    return factor

def _next_random(flag_int: bool):
    if flag_int:
        return random.randint(int(__LOW__), int(__HIGH__))
    else:
        return (random.random() - 0.5) * (__HIGH__ - __LOW__)


def generate_positive_definite_matrix(n: int, flag_int: bool = True):
    """
    Generates Positive-Definite matrices using lower-triangular factors.
    Matrices generated using this method tend to have nonuniform distribution (with magnitudes increasing with row
    and column) and additionally have their Cholesky factor computed as an intermediary step.
    :param n:
    :param flag_int:
    :return:
    """
    L = new_square_matrix(n)
    for i in range(n):
        for j in range(i):
            L[i][j] = _next_random_nonzero(flag_int)

        L[i][i] = abs(_next_random_nonzero(flag_int))

    return matrix.chol_compose(L)


def _generate_small_positive_definite_matrix(n: int):
    """
    Creates small positive-definite matrices.
    These matrices have small magnitudes but still have nonuniform distribution.
    :param n:
    :return:
    """
    L = new_square_matrix(n)
    for i in range(n):
        for j in range(i):
            L[i][j] = random.random()

        while abs(L[i][i]) < __THRESHOLD__:
            L[i][i] = random.random()

    return matrix.chol_compose(L)


def safely_generate_positive_definite_matrix(n: int):
    """
    Creates positive-definite matrices without calculating the lower triangular choleski factor of the output.
    These matrices still tend to have nonuniform distribution.
    :param n: size
    :return:
    """
    N = generate_positive_definite_matrix(n, True)
    M = _generate_small_positive_definite_matrix(n)
    return matrix.multiply(N, matrix.multiply(M, N))


def generate_random_matrix(n: int, m: int = 0, flag_int: bool = True) -> matrix.Matrix:
    """
    Generate random matrix. There is no limit as to what this matrix is, other than dimension
    :param n: # Rows
    :param m: Optional # columns (if omitted, matrix is square)
    :param flag_int: Optional flag determining if elements should be integers
    :return: Random matrix
    """
    if m == 0:
        m = n

    M = matrix.new_matrix(n, m)
    for i in range(n):
        for j in range(m):
            M[i][j] = _next_random(flag_int)

    return M


def generate_real_symmetric_matrix(n: int, flag_int: bool = True) -> matrix.Matrix:
    """
    Creates real-symmetric matrices.
    :param n: size
    :param flag_int: Flag for if matrix elements are floating-point (false) or integer (true).
    :return:
    """
    global __HIGH__, __LOW__
    M = new_square_matrix(n)
    __HIGH__ *= n
    __LOW__ *= n
    for i in range(n):
        for j in range(i):
            M[i][j] = _next_random_nonzero(flag_int)
            M[j][i] = M[i][j]

        M[i][i] = abs(_next_random_nonzero(flag_int))

    __HIGH__ /= n
    __LOW__ /= n
    return M