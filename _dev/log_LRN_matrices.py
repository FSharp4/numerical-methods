from circuit import Circuit
from matrix import diag, multiply, vec_to_mat, subtract, write, transpose

if __name__ == "__main__":
    for i in range(2, 16):
        circuit = Circuit.read_from_csv(f"incidence{i}.csv")
        Y = diag(circuit.R)
        for j in range(len(circuit.R)):
            Y[j][j] = 1 / Y[j][j]

        AYAT = multiply(circuit.A, multiply(Y, transpose(circuit.A)))
        AJYE = multiply(circuit.A, subtract(vec_to_mat(circuit.J), multiply(Y, vec_to_mat(circuit.E))))
        write(AYAT, f"AYAT{i}.csv")
        write(AJYE, f"AJYE{i}.csv")
