import copy
import csv
from copy import deepcopy
from math import sqrt
from typing import List

import matrix

# TODO: Separate cholesky material into file and evaluate for submission 1A

Matrix = List[List[int | float]]
Vector = List[int | float]


def write(matrix: Matrix, filename: str) -> None:
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',', quotechar="|", quoting=csv.QUOTE_MINIMAL)
        for row in matrix:
            writer.writerow(row)


def read(filename: str) -> Matrix:
    mat: Matrix = []
    with open(filename, 'r', newline='') as file:
        reader = csv.reader(file, delimiter=',', quotechar="|", quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            n = len(row)
            vector = [0.0] * n
            for i in range(n):
                vector[i] = float(row[i].replace('\ufeff', ''))

            mat.append(vector)

    return mat


def new_matrix(rows: int, cols: int, init_val: int | float = 0) -> Matrix:
    return [[init_val] * cols for _ in range(rows)]


def new_square_matrix(size: int, init_val: int | float = 0) -> Matrix:
    return new_matrix(size, size, init_val)


def new_vector(size: int, init_val: float | int = 0) -> Vector:
    # At this point I don't think I'll move from using lists to arrays, but just in case, here's a convenience method
    return [init_val] * size


def diag(entries: Vector) -> Matrix:
    mat = new_square_matrix(len(entries))
    for index in range(len(entries)):
        mat[index][index] = entries[index]

    return mat


def pprint_matrix(matrix: Matrix) -> None:
    print()
    for row in matrix:
        s = ""
        for element in row:
            s += "{:10.4f}\t".format(element)
        print(s)

    print()


def pretty_vector(vector: Vector) -> str:
    s = "["
    for index in range(len(vector)):
        element = round_if_int(vector[index])
        if type(element) is int:
            s += f"{element}"
        else:
            s += "{:10.4f}".format(element).strip()
        if index != len(vector) - 1:
            s += ", "

    s += "]"
    return s


def matrix_to_string(matrix) -> List[str]:
    n = len(matrix)
    out = [""] * n
    for index in range(n):
        for element in matrix[index]:
            out[index] += "{:10.4f}\t".format(element)

    return out


def transpose(matrix: Matrix) -> Matrix:
    rows = len(matrix)
    cols = len(matrix[0])
    t = new_matrix(cols, rows)
    for i in range(rows):
        for j in range(cols):
            t[j][i] = matrix[i][j]

    return t


def is_symmetric(matrix: Matrix) -> bool:
    return fp_matrix_equals(matrix, transpose(matrix))


def dot(v1: Vector, v2: Vector) -> float:
    """
    Computes dot product of two vectors
    :param v1: Vector, size n
    :param v2: Vector, size n
    :return: Scalar, v1 • v2
    """
    n = len(v1)
    if n != len(v2):
        raise Exception(f"Vectors must be of equal length: Expected {n} in v2 but got {len(v2)}")

    product = 0
    for i in range(n):
        product += v1[i] * v2[i]

    return product


def multiply(m1: Matrix, m2: Vector | Matrix) -> Matrix | float | int:
    """
    Multiplies two matrices together
    :param m1: First matrix, mxk
    :param m2: Second matrix, kxn
    :return: Product matrix, mxn
    """
    if type(m2[0]) == int | float:
        m2: Matrix = matrix.vec_to_mat(m2)
    m1_rows = len(m1)
    m1_cols = len(m1[0])
    m2_rows = len(m2)
    m2_cols = len(m2[0])

    if m1_cols != m2_rows:
        raise Exception(f"Incompatible matrix sizes: expected {m1_cols} rows in m2 but got {m2_rows}")

    product = new_matrix(m1_rows, m2_cols)
    m2_transpose = transpose(m2)
    for i in range(m1_rows):
        for j in range(m2_cols):
            product[i][j] = dot(m1[i], m2_transpose[j])

    if len(product) == 1:
        if len(product[0]) == 1:
            return product[0][0]

    return product


def is_positive_definite(matrix: Matrix, use_cholesky: bool = False) -> bool:
    """
    TODO: COMPLETE

    Tests if real-valued matrix is positive definite
    :param use_cholesky: Flag for determining if cholesky algorithm can be used to determine P-D (can't set to true for
                         assignment).
    :param matrix: Real-valued matrix
    :return: Positive definiteness
    """

    if not is_symmetric(matrix):  # Tests if matrix is Hermitian
        return False

    if use_cholesky:
        # Cholesky algorithm is very efficient at determining if matrices are p-d.
        try:
            chol(matrix)  # If it works, this is p-d
            return True
        except ValueError:
            return False

    else:
        """
        Cholesky decomposition is commonly used for determining if matrices are positive-definite (P-D), but it is not 
        the only way to do so.
        This algorithm is an adaptation of Gaussian Elimination and relies on upper left sub-matrix determinant values; 
        the number of sign changes in that sequence are equal to the number of negative eigenvalues. 
        Positive-definite matrices, by definition, have strictly positive eigenvalues.
        Thus, the idea here is to transform the matrix into a triangular matrix using row operations that do not
        change the value of the determinant (namely, by adding multiples of rows to other rows).
        This is equivalent to finding the row-echelon form of the entire matrix by only using the "add" row 
        operation.
        This is done iteratively as the iteration for transforming the matrix matches nicely with checking the
        upper left sub-matrices for positive determinants.
        We start by considering the top left element; if it is not positive, then the matrix cannot be P-D.
        At each iteration, we add another "row" to our consideration, considering the next sub-matrix 
        simultaneously.
        Considering iteration i:
        - We first subtract multiples of prior rows from row i to make column elements 0...(i-1) in row i equal to 0
        - Then, we check to make sure that the bottom right element in the sub-matrix is positive, breaking if it is 
          not. This is equivalent to checking the determinant positivity because the determinant of a triangular matrix 
          is equal to its trace, and the prior elements of the trace were confirmed to be strictly positive in prior 
          iterations.
        
        Gaussian Elimination is O(n^3) and has known numerical instability issues, but is simple and sufficient enough 
        for our matrices.
        """
        if matrix[0][0] == 0:
            return False

        mat: Matrix = deepcopy(matrix)
        for i in range(1, len(mat)):
            for j in range(i):
                factor = mat[i][j] / mat[j][j]
                if factor != 0:
                    _add(mat, i, j, -factor)

            if mat[i][i] <= 0:
                return False

        return True


def chol(mat: Matrix) -> Matrix:
    """
    Performs cholesky decomposition on matrix A
    Implementation is naive and not in-place (A is preserved).
    No safety checking is performed, but program has type hints; algorithm assumes matrix is positive-definite
    :param mat: Positive-definite matrix M = l_factor * l_factor^T
    :return: l_factor: Lower-triangular matrix decomposed from A
    """
    n = len(mat)
    l_factor = new_matrix(n, n)

    for j in range(n):
        l_sum = 0
        for i in range(j):
            l_sum += l_factor[j][i] * l_factor[j][i]

        l_factor[j][j] = sqrt(mat[j][j] - l_sum)
        for i in range(j + 1, n):
            l_sum = 0.0
            for k in range(j):
                l_sum += l_factor[i][k] * l_factor[j][k]

            l_factor[i][j] = 1 / l_factor[j][j] * (mat[i][j] - l_sum)

    return l_factor


def chol_compose(l_factor: Matrix):
    return multiply(l_factor, transpose(l_factor))


# noinspection DuplicatedCode
def chol_solve(mat: Matrix, c: Vector) -> Vector:
    """
    Performs cholesky decomposition on system Mx = c
    No safety checking is performed, but program has type hints; algorithm assumes system is solvable.
    :param mat: Real, symmetric, positive-definite matrix of size [n,n]
    :param c: Real-valued vector of size n
    :return: x: Solution to System of Linear Equations, real-valued vector of size n
    """
    """
    cholesky Decomposition:
    1. Factor sys_mat = L L^T where L is lower-triangular

        sys_mat[1,1] = L[1,1]^2                       -> L[1,1] = sqrt(sys_mat[1,1])
        sys_mat[i,1] = L[i,1]L[1,1]                   -> L[i,1] = sys_mat[i,1]/L[1,1]

        sys_mat[j,j] = sum_{i = 1}^j [L[j, i]^2]      -> L[j,j] = sqrt(sys_mat[j,j] - sum_{i = 1}^{j - 1} [L[j, i]^2]
        sys_mat[i,j] = sum_{k = 1}^j [L[i,k]L[j,k]]   -> L[i,j] = \frac{sys_mat[i,j] - sum_{k=1}^{j-1} [L[i,k]L[j,k]]

    2. Solve Ly = b for y
    3. Solve L^T x = y for x
    """
    sys_mat = deepcopy(mat)
    b = deepcopy(c)
    n = sys_mat.__len__()
    for j in range(n):
        sys_mat[j][j] = sqrt(sys_mat[j][j])  # Throws math domain error if less than 0
        b[j] /= sys_mat[j][j]
        for i in range(j + 1, n):
            sys_mat[i][j] /= sys_mat[j][j]
            b[i] -= sys_mat[i][j] * b[j]
            for k in range(j + 1, i + 1):
                sys_mat[i][k] -= sys_mat[i][j] * sys_mat[k][j]

    # Now, b is y and sys_mat is L
    x = [0] * n
    for i in range(n - 1, -1, -1):
        temp = 0
        for j in range(i + 1, n):
            temp += sys_mat[j][i] * x[j]
        # noinspection PyTypeChecker
        x[i] = (b[i] - temp) / sys_mat[i][i]

    return x


# noinspection DuplicatedCode
def chol_band_solve(mat: Matrix, c: Vector, profile_bandwidth=False) -> tuple[Vector, int] | Vector:
    """
    Performs cholesky decomposition on system Mx = c where M is "smoothly" banded.
    This requires that, if the leftmost nonzero element of row i occurs in column j, then there can be no row i+k where
    column j-l (k, l > 0) is nonzero. Violating this constraint will produce incorrect results.
    This is an aggressive optimization of chol_solve for the ECSE 543 A1 assignment.
    :param mat: Real, symmetric, positive-definite matrix of size [n,n]
    :param c: Real-valued vector of size n
    :param profile_bandwidth: Optional flag to return bandwidth of M
    :return: x: Solution to System of Linear Equations, real-valued vector of size n, or a tuple of that plus the
                bandwidth of the matrix.
    """
    sys_mat = deepcopy(mat)
    b = deepcopy(c)
    n = sys_mat.__len__()

    column_bounds_exclusive = new_vector(len(sys_mat[0]), init_val=n)
    row_ptr = 0
    for column_ptr in range(n):
        if row_ptr == n:
            break

        for row_ptr in range(row_ptr, n):
            if sys_mat[row_ptr][column_ptr] == 0:
                column_bounds_exclusive[column_ptr] = row_ptr
                break

    for j in range(n):
        # Modify diagonal
        sys_mat[j][j] = sqrt(sys_mat[j][j])  # Throws math domain error if less than 0
        # Adjust constant vector towards b
        b[j] /= sys_mat[j][j]
        for i in range(j + 1, column_bounds_exclusive[j]):
            # Propagate lookahead increment below diagonal entry
            sys_mat[i][j] /= sys_mat[j][j]
            b[i] -= sys_mat[i][j] * b[j]
            for k in range(j + 1, i + 1):
                # Propagate lookahead increment rightwards of element below diagonal back to it.
                sys_mat[i][k] -= sys_mat[i][j] * sys_mat[k][j]

    # Now, b is y and sys_mat is L
    x = new_vector(n, init_val=0)
    for i in range(n - 1, -1, -1):
        temp = 0
        for j in range(i + 1, column_bounds_exclusive[i]):
            temp += sys_mat[j][i] * x[j]
        x[i] = (b[i] - temp) / sys_mat[i][i]

    if profile_bandwidth:
        band: int = max(column_bounds_exclusive[j] - j for j in range(n))
        return x, band

    return x


def _swap(mat: Matrix, r1: int, r2: int) -> None:
    mat[r1], mat[r2] = mat[r2], mat[r1]


def _multiply(mat: Matrix, row: int, coefficient) -> None:
    # noinspection PyTypeChecker
    for j in range(len(mat[0])):
        mat[row][j] *= coefficient


def _add(mat: Matrix, target: int, source: int, coefficient) -> None:
    # noinspection PyTypeChecker
    for j in range(len(mat[0])):
        mat[target][j] += mat[source][j] * coefficient


# noinspection GrazieInspection
def round_if_int(f: float) -> int | float:
    """
    Casts floating-point quantities to integer if they are deemed to be equal to integers.
    Only really works for small numbers and is meant to assist pretty print routines rather than serious applications.
    :param f: FP quantity to check
    :return: round(FP) if FP 'is an integer', else the FP quantity
    """
    atol = 1e-8
    i = round(f)
    difference = abs(f - i)
    return i if difference < atol else f


# noinspection GrazieInspection
def fp_equals(f: float, i: int | float, rtol=1e-5, atol=1e-8) -> bool:
    """
    More forgiving equality check for looking between two floats, or a float and an int
    :param f: Float number
    :param i: Number to compare against (float or int)
    :param rtol: Optional Relative tolerance (wrt i)
    :param atol: Optional Absolute tolerance
    :return: If f == i is considered true
    """
    difference = abs(f - i)
    return difference < (atol + rtol * i)


def fp_matrix_equals(mat1: Matrix, mat2: Matrix, rtol=1e-5, atol=1e-8) -> bool:
    """
    Returns if M1 is sufficiently close to M2 to be considered equal under the constraints of finite-precision
    arithmetic.
    Idea and tolerance parameters taken from numpy.allclose
    :param mat1: Matrix 1
    :param mat2: Matrix 2, of equal size
    :param rtol: Relative tolerance (calculated as proportion of abs(M2)
    :param atol: Absolute tolerance margin
    :return: |M1 - M2| < atol + rtol * |M2|
    """
    n = len(mat1)
    m = len(mat1[0])
    square_sum_diff: float = 0
    square_sum_mat2: float = 0
    for i in range(n):
        for j in range(m):
            square_sum_diff += (mat1[i][j] - mat2[i][j]) * (mat1[i][j] - mat2[i][j])
            square_sum_mat2 += mat2[i][j] * mat2[i][j]

    absolute_value_difference = sqrt(square_sum_diff)
    absolute_value_mat2 = sqrt(square_sum_mat2)
    threshold = atol + rtol * absolute_value_mat2
    return absolute_value_difference < threshold


def vec_to_mat(v: Vector) -> Matrix:
    mat = new_matrix(len(v), 1)
    for i in range(len(v)):
        mat[i][0] = v[i]

    return mat


def mat_to_vec(m: Matrix) -> Vector:
    vec = [0.0] * len(m) * len(m[0])
    ptr = 0
    for i in range(len(m)):
        for j in range(len(m[0])):
            vec[ptr] = m[i][j]
            ptr += 1

    return vec


def subtract(mat1: Matrix, mat2: Matrix) -> Matrix:
    mat: Matrix = new_matrix(len(mat1), len(mat1[0]))
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            mat[i][j] = mat1[i][j] - mat2[i][j]

    return mat


def add(mat1: Matrix | Vector, mat2: Matrix) -> Matrix:
    if type(mat1[0]) == int or type(mat1[0]) == float:
        mat1 = matrix.transpose([mat1])

    mat: Matrix = new_matrix(len(mat1), len(mat1[0]))
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            mat[i][j] = mat1[i][j] + mat2[i][j]

    return mat

def scalar_multiply(mat: Matrix, sca: float | int) -> Matrix:
    m: Matrix = new_matrix(len(mat), len(mat[0]))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            m[i][j] = sca * mat[i][j]

    return m


def norm1_int(v: Vector) -> int:
    """
    Returns the first norm of a vector
    :param v: Vector
    :return: First norm of vector (sum of the absolute values of all elements)
    """
    return sum(abs(element) for element in v)


def norm2(v: Vector) -> float:
    return sum(element * element for element in v)


def illustrate_sparseness(mat: Matrix) -> None:
    s = ""
    for vector in mat:
        for element in vector:
            if fp_equals(element, 0):
                s += "   "
            else:
                s += " * "
        s += '\n'
    print(s)


def mask_equals(mat: Matrix, level: float | int) -> Matrix:
    bool_mat = new_matrix(len(mat), len(mat[0]))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            bool_mat[i][j] = int(fp_equals(mat[i][j], level))

    return bool_mat


def mask_int_equals(mat: Matrix, level: int) -> Matrix:
    bool_mat = new_matrix(len(mat), len(mat[0]))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            bool_mat[i][j] = (type(mat[i][j]) == int) and mat[i][j] == level

    return bool_mat

def rot180(mat: Matrix) -> Matrix:
    m = len(mat)
    n = len(mat[0])
    rot_mat = new_matrix(m, n)
    for i in range(m):
        for j in range(n):
            rot_mat[m - i - 1][n - j - 1] = mat[i][j]

    return rot_mat

def col_sort(mat: Matrix, col: int) -> Matrix:
    def col_key(entry: Vector):
        return entry[col]
    sort_mat = copy.deepcopy(mat)
    sort_mat.sort(key=col_key)
    return sort_mat

def norm_inf(vec: Matrix | Vector) -> int | float:
    if type(vec[0]) == list:
        vec = mat_to_vec(vec)
    return max(vec)