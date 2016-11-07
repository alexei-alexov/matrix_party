#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import factorial

import numpy as np
from numpy import array
from numpy.linalg import inv


class MatrixParty(object):

    # PREDEFINED
    # v, m, l, c, N
    c = 5
    l = 4
    m = 2
    v = 3
    N = 51

    def generate_array(self, size):
        """
        Generate square array filled with 0
        """
        print "MATRIX GENERATION. SIZE: ", size
        result = []
        for i in xrange(size):
            result.append([0] * size)
        return result

    def get_one(self):
        """
        Return ¯1¯

         - OK
        """
        return [1] * (self.c+1)

    def get_e(self, i):
        """
        Return ei = (di0, di1, di2, ... dic)T

         - OK
        """
        result = np.array([1 if i == j else 0 for j in xrange(self.c+1)])
        print "GENERATION E: ", result
        return result.T

    def get_fii(self, j, i):
        return self.l + i * self.m + j * self.v

    def get_fiimo(self, j, i):
        """
        Return f[j][i;i-1]
        """
        return -self.l

    def get_fipo(self, j, i):
        """
        Return f[j][i;i+1]
        """
        return -((i+1)*self.m)

    def get_fcc(self, j):
        return self.l + self.c * self.m

    def get_fccmo(self, j):
        return -(self.l + j * self.v)

    def get_fck(self, j, k):
        return -j*self.v

    def get_f(self, j):
        """
        This function return F(j)

         - OK
        """
        print "GET F CALLED. J:", j
        # creating array filled with 0
        f = self.generate_array(self.c+1)

        f[0][0] = self.l + j * self.v
        f[self.c][self.c] = self.l + self.c * self.m
        f[0][1] = -self.m
        f[self.c][self.c-1] = -(self.l + j * self.v)
        f[self.c-1][self.c-1] = self.l + (self.c-1) * self.m + j * self.v
        f[self.c-1][self.c-2] = -self.l
        f[self.c-1][self.c] = -self.c*self.m

        for i in xrange(1, self.c - 1):
            f[i][i] = self.l + i * self.m + j * self.v
            f[i][i-1] = -self.l
            f[i][i+1] = -(i+1) * self.m
            f[self.c][i] = -j * self.v

        print "RESULT: ", f
        return f

    def get_b(self):
        """
        This function return B
        Example:
         0 0 0
         1 0 0
         0 1 1

          - OK
        """
        b = self.generate_array(self.c+1)
        b[self.c][self.c] = 1
        for i in xrange(1, self.c+1):
            b[i][i-1] = 1
        print "GET B RESULT: ", b
        return b

    def get_d(self):
        """
        This function return D array
        """
        d = self.generate_array(self.c+1)

        d[0][0] = 1

        f = self.get_f(self.N)

        for i in xrange(1, self.c+1):
            for k in xrange(0, self.c+1):
                d[i][k] = f[i-1][k]

        print "GET D CALLED. RESULT: ", d
        return d

    def get_w(self, j):
        """
        This function return W(j)

        """
        result = inv(np.dot(self.get_f(j), self.get_b()))
        print "GET W CALLED. RESULT: ", result
        return result

    def get_bigw(self, j, r):
        """
        This return multiplyed matrices
        """
        print "BIG W ", j, " ", r
        result = self.get_f(j)
        for k in xrange(j+1, r+1):
            print "K: ", k
            print "BEFORE:\n", result
            result = np.dot(result, self.get_f(k))
            print "AFTER:\n", result

        return result

    def delta_j(self, j):
        print ">>>>>> DELTA J >>>>>>>>>>>"
        first_part = 1 / (factorial(j) * self.v**j)
        # D inversed
        dmo = inv(self.get_d())
        print dmo, type(dmo)
        # e with zero
        ewz = self.get_e(0)
        second_part_up = self.get_bigw(j, self.N-1)
        second_part_up = np.dot(second_part_up, dmo)
        second_part_up = np.dot(second_part_up, array(ewz))
        print second_part_up, type(second_part_up)
        second_part_down = np.dot(
            np.dot(np.dot(ewz.T, self.get_bigw(0, self.N-1)), dmo), ewz)

        return first_part * (second_part_up / second_part_down)

    # def delta_n(self):
    #     first_part = 1 / (factorial(self.N) * self.v**self.N)
    #     # D inversed
    #     dmo = inv(self.get_d())
    #     # e with zero
    #     ewz = self.get_e(0)
    #     second_part_up = dmo * ewz
    #     second_part_down = ewz.T * self.get_bigw(0, self.N-1) * dmo * ewz
    #     return first_part * (second_part_up / second_part_down)

    def get_p(self):
        result = self.get_one() * self.delta_j(1)
        for j in xrange(1, self.N+1):
            result += (self.get_one() * self.delta_j(j))

        return inv(result)


if __name__ == "__main__":
    m = MatrixParty()
    # print m.get_d()
    print m.get_p()

    # m1 = array([[0, -2147483648, 0, -2147483648, 0, 0],
    #             [0, -2147483648, 0, -2147483648, 0, 0],
    #             [0, -2147483648, 0, 0, 0, 0],
    #             [0, -2147483648, 0,   0, 0, 0],
    #             [-2147483648, -1073741824, 0, 1073741824, 0, 0],
    #             [-2147483648, -1073741824, 0, 1073741824, 0, 0]])
    # print m1.T
    #
    # m2 = array([[130, -2, 0, 0, 0, 0],
    #             [-4, 132, -4, 0,    0,    0],
    #             [0, -4,  134,   -6,    0,    0],
    #             [0, 0,   -4,  136,   -8,    0],
    #             [0, 0,    0,   -4,  138,  -10],
    #             [0, -126, -126, -126, -130,   14]])
    # print m2.T
    #
    # print np.dot(m1.T.tolist(), m2.T.tolist())

    # print get_d(3, get_f(N))
    # print get_e(1, 5)
