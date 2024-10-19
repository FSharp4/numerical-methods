import matrix
import matrix_generator
from gaussian_elimination import pivot_test

matrix_generator.set_seed(0)

for size in range(2, 11):
    M = matrix_generator.generate_positive_definite_matrix(size, flag_int=True)
    pd = True
    if not pivot_test(M):
        print(f"Not PD: size {size} @ check 2")

