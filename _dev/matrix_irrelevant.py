from copy import deepcopy

from matrix import Matrix, _swap


def row_echelon_form(M: Matrix):
    """
    Gaussian Elimination algorithm for transforming arbitrary matrix into row echelon form
    :param M: Matrix
    :return: ref(M)
    """
    """
    Variable naming follows common Gaussian Elimination notation, with the caveat that pythonic arrays start at 0
    """
    A = deepcopy(M)
    m = len(A)  # Rows
    n = len(A[0])  # columns
    h = 0  # pivot row
    k = 0  # pivot column
    polarity_ctr = 0

    while h < m and k < n:
        # Find k-th pivot
        i_max = h
        i_max_value = abs(A[i_max][k])
        for i in range(h, m):
            value = abs(A[i][k])
            if i_max_value < value:
                i_max_value = value
                i_max = i

        if A[i_max][k] == 0:
            # Column has no pivot
            k += 1
        else:
            _swap(A, h, i_max)  # polarity is reversed
            polarity_ctr += 1
            for i in range(h + 1, m):
                f = A[i][k] / A[h][k]
                # Fill with zero at lower part of pivot column
                A[i][k] = 0
                for j in range(k + 1, n):
                    A[i][j] -= A[h][j] * f
                    # This is a multiply-accumulate row operation; determinant polarity preserved.

            h += 1
            k += 1

    if polarity_ctr % 2 == 1:
        for j in range(n):
            A[0][k] *= -1
    return A
