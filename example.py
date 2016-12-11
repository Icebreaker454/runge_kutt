"""
This is an example usage of the Runge-Kutt differential equation solving
algorithm
"""

import math

from runge_kutt import RungeKuttMethod
from runge_kutt.functions import DoubleVectorFunction, VectorFunction


def main():
    # Defining single-element vector function for one differential equation.
    # u' = 2u
    F = DoubleVectorFunction([
        lambda t, u: (
            2. * u[0]), ])
    # Defining corresponding exact solution for mismatch checking
    U = VectorFunction([
        lambda t: 10 * math.exp(2 * t),
    ])

    method = RungeKuttMethod(F, 0., 1., [10.], 1E-4, 0.1, 1E-6, U)
    method.solve()
    method.output()

    # Defining a vector-function for a differential equation system.
    FS = DoubleVectorFunction([
        lambda t, U: U[0] - U[1],
        lambda t, U: U[0] * U[1]
    ])

    method = RungeKuttMethod(FS, 0., 3., [0., 0.5], 1E-4, 0.3, 1E-6)
    method.solve()
    method.output()


if __name__ == '__main__':
    main()
