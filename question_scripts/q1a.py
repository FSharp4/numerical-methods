import matrix
from matrix import Matrix, Vector
import settings
from lagrange import FullDomainLagrangeInterpolator

RESOLUTION = 100


def run():
    global RESOLUTION
    data: Matrix = matrix.transpose(matrix.read(f"{settings.PROJECT_ROOT}/m19_bh.csv")[:6])
    x: Vector = data[0]
    y: Vector = data[1]

    interpolator: FullDomainLagrangeInterpolator = FullDomainLagrangeInterpolator().fit(x, y)
    step_size = (max(x) - min(x)) / RESOLUTION
    x_test: Vector = matrix.arange(min(x), max(x), step_size, end_inclusive=True)
    y_test: Vector = list(map(interpolator.predict, x_test))
    test_data: Matrix = [x_test, y_test]
    matrix.write(matrix.transpose(test_data), "outputs/q1a.csv")


if __name__ == "__main__":
    run()
