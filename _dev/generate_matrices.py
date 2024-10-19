import matrix
import matrix_generator

# TODO: Evaluate for submission 1B

if __name__ == "__main__":
    """
    This script should generate the set of test matrices we use for the assignment.
    At small size 2, the matrices are floating-point not inherently spoiled by knowing pre-emptively the L factor.
    At mid-small sizes 3-5, the matrices are integer for human readability
    At mid-large sizes 6-9, the matrices are floating point but bounded reasonably
    At large sizes 10-11, the matrices are floating point and not inherently spoiled as previously described.
    """
    matrix_generator.set_seed(5)
    for size in range(2, 12):
        M: matrix.Matrix
        if size == 2:
            M = matrix_generator.safely_generate_positive_definite_matrix(size)
        elif size < 6:
            M = matrix_generator.generate_positive_definite_matrix(size, flag_int=True)
        elif size < 10:
            M = matrix_generator.generate_positive_definite_matrix(size, flag_int=False)
        else:
            M = matrix_generator.safely_generate_positive_definite_matrix(size)

        if (matrix.is_positive_definite(M)):
            matrix.write(M, f"../matrices/Test_matrix_{size}.csv")
        else:
            raise ArithmeticError(f"Matrix {size} is not positive definite!")
