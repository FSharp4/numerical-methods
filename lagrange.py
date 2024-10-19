import copy

import matrix
from matrix import Vector
from typing import Self, List  # Very convenient


class LagrangianPolynomial:
    def __init__(self):
        self.main_node_position = -1
        self.nodes: Vector = []

    def fit(self, main_node_position: int, nodes: Vector) -> Self:
        self.main_node_position = main_node_position
        self.nodes = copy.deepcopy(nodes)
        return self

    def evaluate(self, x: int):
        y = 1
        for j in range(len(self.nodes)):
            if j != self.main_node_position:
                y *= (x - self.nodes[j]) / \
                    (self.nodes[self.main_node_position] - self.nodes[j])

        return y


class FullDomainLagrangeInterpolator:
    def __init__(self):
        self.polynomials: List[LagrangianPolynomial] = []
        self.values: Vector = []

    def fit(self, x: Vector, y: Vector) -> Self:
        self.values = copy.deepcopy(y)
        self.polynomials = []
        for i in range(len(x)):
            self.polynomials.append(LagrangianPolynomial().fit(i, x))

        return self

    def predict(self, x: float | int) -> float:
        y = 0
        for i in range(len(self.values)):
            y += self.values[i] * self.polynomials[i].evaluate(x)

        return y
