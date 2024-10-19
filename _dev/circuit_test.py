import circuit
from circuit import Circuit

CIRCUIT = 1

if __name__ == "__main__":
    circuit = circuit.read_from_csv(f"../circuits/circuit{CIRCUIT}.csv")
    print(str(circuit))
