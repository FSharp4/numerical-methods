import _dev.debug_matrix_generator
import matrix
import matrix_generator

SEED = 0

if __name__ == "__main__":
    matrix_generator.set_seed(SEED)
    for iteration in range(100):
        # harness = _dev.debug_matrix_generator._debug_positive_definite(5, flag_int=True)
        M = matrix_generator.generate_positive_definite_matrix(5, flag_int=True)
        # preL = harness[0]
        # print("Generating Factor: ")
        # matrix.pprint_matrix(preL)
        print("Composed Matrix: ")
        matrix.pprint_matrix(M)
        print("Recovered Choleski Factor: ")
        matrix.pprint_matrix(matrix.chol(M))
