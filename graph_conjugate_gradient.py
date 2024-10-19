"""
THIS IS THE SOLE EXCEPTION TO THE "NO LIBRARIES CLAUSE". Matplotlib is used ONLY within this file, ONLY for generating
graphs.
"""
from matplotlib import pyplot as plt

import conjugate_gradient
import matrix
from cjm_conductor_simulation import adjusted_conductor_matrix

if __name__ == "__main__":
    full_matrix, b = adjusted_conductor_matrix()
    full_matrix = matrix.scalar_multiply(full_matrix, -1)
    print(matrix.is_positive_definite(full_matrix))
    V = matrix.chol_solve(full_matrix, b)

    V2, norm_2, norm_inf = conjugate_gradient.solve(full_matrix, b, monitor=True)
    iterations = list(range(len(norm_2)))
    plt.figure()
    plt.semilogy(iterations, norm_2)
    # ax[1].semilogy(iterations, norm_inf)
    plt.title("Logplot of Norm2 of residual over iteration count")
    plt.show()

    plt.figure()
    plt.semilogy(iterations, norm_inf)
    plt.title("Logplot of Norm_inf of residual over iteration count")
    plt.show()

    plt.figure()
