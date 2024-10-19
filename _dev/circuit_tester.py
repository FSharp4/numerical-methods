import circuit
from matrix import pretty_vector

# TODO: Evaluate for submission 1D

if __name__ == "__main__":
    for i in range(1, 6):
        c = circuit.read_from_csv(f"../circuits/circuit{i}.csv")
        print(f"Circuit {i} Voltages: {pretty_vector(c.voltages())}")
