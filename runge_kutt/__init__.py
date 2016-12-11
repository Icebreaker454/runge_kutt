"""
This module contains the Runge-Kutt algorithm for solving differential
equations and differential equation systems.
"""
import math
from operator import sub
from copy import copy

from tabulate import tabulate


class RungeKuttMethod(object):

    COEF = {
        'c': [0, 0.5, 1.],
        'a': [
            [0, 0, 0],
            [0.5, 0, 0],
            [-1., 2., 0]],
        'b': [1. / 6., 4. / 6., 1. / 6.],
        'p': 3.
    }
    RGC = math.pow(2., COEF['p']) - 1.

    def __init__(self, F, t_0, T, y_0, eps, tau_0, eps_M, U=None, **kwargs):
        """
        :param F: - the vector function describing the equation to be solved
        :param t_0: - initial time
        :param T: - end of the time interval
        :param y_0: - initial function vector values
        :param eps: - desired precision
        :param tau_0: - initial integration step
        :param eps_M: - "machine noise" value. Lowest posible value delta
        :param U: - exact solution (if a comparison is needed)
        """
        self.F = F
        self.t_0 = t_0
        self.T = T
        self.y_0 = y_0
        self.eps = eps
        self.tau_0 = tau_0
        self.eps_M = eps_M
        self.U = U

    @classmethod
    def norm(cls, a):
        return max(map(abs, a))

    def solve(self):

        calls = 0

        steps = []

        t = self.t_0
        y = self.y_0
        tau = self.tau_0

        step = [t, y]
        if self.U:
            step.append(self.U(t))
            step.append(self.norm(map(sub, y, self.U(t))))
            self.e_max = 0.
        steps.append(step)

        while abs(self.T - t) >= self.eps_M:
            if t + tau > self.T:
                tau = self.T - t
            v = copy(y)
            t_1 = t
            while True:
                kf = 0
                while True:
                    if kf != 1:
                        k_1 = self.F(t, y)
                        calls += 1

                    k_2 = self.F(
                        t + self.COEF['c'][1] * tau,
                        map(lambda x, y: x + tau * self.COEF['a'][1][0] * y,
                            y, k_1)
                    )
                    k_3 = self.F(
                        t + self.COEF['c'][2] * tau,
                        map(
                            lambda x, y, z:
                                x + y * tau * self.COEF['a'][2][0] +
                                z * tau * self.COEF['a'][2][1],
                            y, k_1, k_2))
                    y = map(
                        lambda x, y, z, t:
                            x + tau * (
                                y * self.COEF['b'][0] +
                                z * self.COEF['b'][1] +
                                t * self.COEF['b'][2]),
                        y, k_1, k_2, k_3)

                    calls += 3

                    if kf == 2:
                        break

                    if kf == 1:
                        t = t + tau
                        kf = 2

                    if kf == 0:
                        w = copy(y)
                        y = copy(v)
                        tau /= 2.
                        kf = 1

                E = self.norm(
                    map(sub, y, w)) / self.RGC / max(1., self.norm(y))

                tau_H = 2. * tau * min(
                    5.,
                    max(0.1, 0.9 * math.pow(
                        (self.eps / E), 1. / (1. + self.COEF['p']))))
                if E <= self.eps:
                    t = t + tau
                    y = map(lambda x, y: x + (x - y) / self.RGC,
                            y, w)
                    tau = tau_H

                    step = [t, y]
                    if self.U:
                        step.append(self.U(t))
                        step.append(self.norm(map(sub, y, self.U(t))))
                        if self.e_max < step[-1]:
                            self.e_max = step[-1]
                    steps.append(step)
                    break

                y = copy(v)
                t = t_1
                tau = tau_H
        self.steps = steps
        self.calls = calls

    def output(self):
        print tabulate(self.steps, headers=['t', 'y', 'U(t)', '|y - U(t)|'])
        if self.U:
            print 'e_max = %s' % self.e_max
        print 'calls : %s' % self.calls
