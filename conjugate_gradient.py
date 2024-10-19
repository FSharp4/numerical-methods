import matrix
from matrix import Matrix, Vector


def solve(mat: Matrix, vec: Matrix | Vector, threshold: int = 100, monitor: bool = False):
    if type(vec) == Matrix:
        vec: Vector = matrix.mat_to_vec(vec)
    m = len(mat)
    n = len(mat[0])
    x = matrix.new_vector(n)

    r = matrix.subtract(matrix.vec_to_mat(vec), matrix.multiply(mat, matrix.vec_to_mat(x)))
    p = r

    k = 0

    norm_2 = [matrix.norm2(matrix.mat_to_vec(r))]
    norm_inf = [matrix.norm_inf(matrix.mat_to_vec(r))]

    while k < threshold and matrix.norm2(matrix.mat_to_vec(r)) > 10e-8:  # TODO: Verify update condition
        alpha = matrix.multiply(matrix.transpose(p), r) / matrix.multiply(matrix.transpose(p), matrix.multiply(mat, p))
        x = matrix.add(x, matrix.scalar_multiply(p, alpha))
        r = matrix.subtract(matrix.vec_to_mat(vec), matrix.multiply(mat, x))
        beta = -matrix.multiply(matrix.transpose(p), matrix.multiply(mat, r)) \
               / matrix.multiply(matrix.transpose(p), matrix.multiply(mat, p))
        p = matrix.add(r, matrix.scalar_multiply(p, beta))
        k += 1
        norm_2.append(matrix.norm2(matrix.mat_to_vec(r)))
        norm_inf.append(matrix.norm_inf(matrix.mat_to_vec(r)))

    if monitor:
        return x, norm_2, norm_inf
    return x

# if __name__ == "__main__":
