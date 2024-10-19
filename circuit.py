import csv

import matrix
from matrix import Matrix, diag, Vector


class Circuit:
    def __init__(self, J: Vector, R: Vector, E: Vector, A: Matrix):
        self.J = J
        self.R = R
        self.E = E
        self.A = A
        self._v: Vector | None = None

    def shape(self) -> (int, int):
        """
        :return: #branches by #nodes
        """
        return len(self.J), len(self.A)

    def __str__(self):
        s = "Branches:\n\n"
        for index in range(len(self.J)):
            s += f"{self.J[index]}, {self.R[index]}, {self.E[index]}\n"

        s += "\nIncidence:\n\n"
        for entry in self.A:
            v = ""
            for number in entry:
                v += f"{str(number)}, "

            s += f"{v[:-2]}\n"

        return s

    # noinspection DuplicatedCode
    def voltages(self):
        if self._v is not None:
            return self._v

        Y = diag(self.R)
        for i in range(len(self.R)):
            Y[i][i] = 1 / Y[i][i]

        AYAT = matrix.multiply(self.A, matrix.multiply(Y, matrix.transpose(self.A)))
        AJYE = matrix.multiply(self.A,
                               matrix.subtract(matrix.vec_to_mat(self.J),
                                               matrix.multiply(Y, matrix.vec_to_mat(self.E))))
        self._v = matrix.chol_solve(AYAT, matrix.mat_to_vec(AJYE))
        return self._v

    @staticmethod
    def read_from_csv(file: str):
        J: Vector = []
        R: Vector = []
        E: Vector = []
        A: Matrix = []
        with open(file, newline='', encoding='utf-8-sig') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            jre_mode = True
            for row in csv_reader:
                if jre_mode:
                    test = row[0]
                    if test != "":
                        J.append(float(row[0]))
                        R.append(float(row[1]))
                        E.append(float(row[2]))
                    else:
                        jre_mode = False
                else:
                    A_row: Vector = []
                    for entry in row:
                        if entry != "":
                            A_row.append(int(entry))
                        else:
                            break

                    A.append(A_row)

        return Circuit(J, R, E, A)

    def load_to_csv(self, file: str) -> None:
        with open(file, 'w', newline='', encoding='utf-8-sig') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            for i in range(len(self.J)):
                row = [self.J[i], self.R[i], self.E[i]]
                writer.writerow(row)

            writer.writerow(['', '', ''])

            for node in self.A:
                writer.writerow(node)


class OptimizedCircuit(Circuit):
    def __init__(self, J: Vector, R: Vector, E: Vector, A: Matrix):
        super().__init__(J, R, E, A)
        self.AYAT_b: int | None = None

    # noinspection DuplicatedCode
    def voltages(self):
        if self._v is not None:
            return self._v

        Y = diag(self.R)
        for i in range(len(self.R)):
            Y[i][i] = 1 / Y[i][i]

        AYAT = matrix.multiply(self.A, matrix.multiply(Y, matrix.transpose(self.A)))
        AJYE = matrix.multiply(self.A,
                               matrix.subtract(matrix.vec_to_mat(self.J),
                                               matrix.multiply(Y, matrix.vec_to_mat(self.E))))
        self._v, self.AYAT_b = matrix.chol_band_solve(AYAT, matrix.mat_to_vec(AJYE), profile_bandwidth=True)
        return self._v
