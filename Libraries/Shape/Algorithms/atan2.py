from numpy import pi
from math import fma

def atan(z: float) -> float:
    raise NotImplementedError
    # implemented dynamically in main()


def atan_(y: float, x: float) -> float:
    if abs(y) <= abs(x):
        return atan(y / x)

    z = x / y
    if z >= 0:
        return pi / 2 - atan(z)
    else:
        return -pi / 2 - atan(z)


def atan2(y: float, x: float) -> float:
    if x == 0:
        if y > 0:
            return pi / 2
        if y < 0:
            return -pi / 2
        return 0

    if x > 0:
        return atan_(y, x)

    if y >= 0:
        return atan_(y, x) + pi
    else:
        return atan_(y, x) - pi



def random_test():
    import numpy as np

    # axes and quadrants
    X, Y = np.meshgrid([-1, 0, 1], [-1, 0, 1])

    actual = np.vectorize(atan2)(X, Y)
    expected = np.atan2(X, Y)

    assert np.allclose(actual, expected, rtol=1e-8)

    # random data
    X, Y = np.random.normal(size=(2, 1_000_000))

    actual = np.vectorize(atan2)(X, Y)
    expected = np.atan2(X, Y)

    assert np.allclose(actual, expected, rtol=1e-14)
    print('maximum absolute error:', np.max(np.abs(actual - expected)))
    print('maximum relative error:', np.max(np.abs(1 - actual / expected)))

def main():
    import numpy as np
    import io
    import contextlib

    N = 19
    assert N & 1, 'odd'

    K = np.arange(N)
    nodes = np.cos((2 * K + 1) / (2 * N) * np.pi)
    M = np.vander(nodes)
    Y = np.atan(nodes)
    X = np.linalg.solve(M, Y)

    print('Coefficients:')
    X = X[1::2]
    print(*X, sep='\n', end='\n\n')

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        print('def atan(z: float) -> float:')
        print('    zz = z*z')
        print('    return', '(' * len(X), end='')
        for c in X[:-1]:
            print(f'{c})*zz + ', end='')
        c = X[-1]
        print(f'{c})*z')
    defn = out.getvalue()

    print(defn)
    exec(defn, globals=globals())

if __name__ == '__main__':
    main()
    random_test()
