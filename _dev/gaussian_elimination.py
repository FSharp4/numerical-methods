from copy import deepcopy

from matrix import Matrix, Vector


def _swap(M: Matrix, r1: int, r2: int):
    M[r1], M[r2] = M[r2], M[r1]


def _multiply(M: Matrix, row: int, coefficient):
    # noinspection PyTypeChecker
    for j in len(M[0]):
        M[row][j] *= coefficient


def _add(M: Matrix, target: int, source: int, coefficient):
    # noinspection PyTypeChecker
    for j in len(M[0]):
        M[target][j] += M[source][j] * coefficient


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


def pivots(M: Matrix) -> Vector:
    """
    Retrieves the pivots of a matrix in row-echelon form
    If the matrix is not in row-echelon form, use gaussian_elimination.row_echelon first.
    :param M: ref matrix
    :return: pivots of ref matrix (array of the first nonzero element in each row, starting at pythonic index 0)
    """
    n = len(M[0])
    pivot_buffer = [0.0] * n
    count = 0

    for j in range(n):
        for i in range(len(M)):
            if M[i][j] != 0:
                pivot_buffer[j] = M[i][j]
                count += 1
                break

    P = [0.0] * count
    pivot_ptr = 0
    for element in pivot_buffer:
        if element != 0:
            P[pivot_ptr] = element
            pivot_ptr += 1

    return P


def pivot_test(M: Matrix) -> bool:
    n = len(M[0])
    if n != len(M):
        return False

    A = deepcopy(M)
    row_ptr = 0
    for column in range(n):
        # find nonzero column
        col_zero_flag = True
        for row in range(row_ptr, n):
            if A[row][column] != 0:
                col_zero_flag = False
                # Found nonzero column
                f = A[row][column]
                for _column in range(n):
                    A[row][_column] /= f  # makes row have leading one; determinant is preserved

                for _row in range(row + 1, n):
                    # make all other column entries below this a zero
                    f = A[_row][column]
                    if f != 0:
                        for _column in range(n):
                            A[_row][_column] -= A[row][_column] * f  # ri -= rj * k

                if row != row_ptr:
                    # make the next available top row the pivot for this column, and eliminate the original pivot row
                    f = A[row_ptr][column] - 1
                    if f != 0:
                        for _column in range(n):
                            A[row_ptr][_column] -= A[row][_column] * f
                            A[row][_column] -= A[row_ptr][_column]

                row_ptr += 1
                break

        if col_zero_flag:
            return False

    for k in range(n):
        if A[k][k] < 0:
            return False

    return True
