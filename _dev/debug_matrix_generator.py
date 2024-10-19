import matrix
from matrix import new_square_matrix
from matrix_generator import _next_random_nonzero


def _debug_positive_definite(n: int, flag_int: bool = True):
    """
    Debug version of main positive definite routine
    Returns the lower triangular choleski factor alongside the positive-definite matrix
    :param n:
    :param flag_int:
    :return: (L, M) where M = L * transpose(L)
    """
    L = new_square_matrix(n)
    for i in range(n):
        for j in range(i):
            L[i][j] = _next_random_nonzero(flag_int)

        L[i][i] = abs(_next_random_nonzero(flag_int))

    return L, matrix.chol_compose(L)
