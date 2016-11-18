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

    def getElement(self, m1, m2, row, column):
        result = 0
        for i in xrange(len(m1)):
            result += m1[row][i]*m2[i][column]

        return result

    def matrix_mult(self, m1, m2):
        """
        This method multiply two same sized square matrices
        """
        length = len(m1)
        return [[self.getElement(m1, m2, row, column)
                 for column in xrange(length)] for row in xrange(length)]

    def vector_mult(self, m, v):
        """
        This method multiply matrix on vector
        """
        result = []
        length = len(m)
        for row in xrange(length):
            row_result = 0
            for column in xrange(length):
                row_result += m[row][column]*v[column][0]
            result.append([row_result])
        return result

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
        result = [1 if i == j else 0 for j in xrange(self.c+1)]
        print "GENERATION E: ", result
        return result

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
        # creating array filled with 0
        f = self.generate_array(self.c+1)

        f[0][0] = self.l + j * self.v
        f[self.c][self.c] = self.l + self.c * self.m
        f[self.c-1][self.c-2] = -1 * self.l
        f[self.c][self.c-1] = -(self.l + j * self.v)
        f[0][1] = -1 * self.m

        f[self.c-1][self.c-1] = self.l + (self.c-1) * self.m + j * self.v
        f[self.c-1][self.c] = -self.c*self.m

        f[self.c][0] = -j * self.v

        for i in xrange(1, self.c - 1):
            f[i][i] = self.l + i * self.m + j * self.v
            f[i][i-1] = -1 * self.l
            f[i][i+1] = -1 * (i+1) * self.m
            f[self.c][i] = -j * self.v

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
        result = inv(self.matrix_mult(self.get_f(j), self.get_b())).tolist()
        print "GET W CALLED. RESULT: ", result
        return result

    def get_bigw(self, j, r):
        """
        This return multiplyed matrices
        """
        result = self.get_f(j)
        for k in xrange(j+1, r+1):
            f = self.get_f(k)
            result = self.matrix_mult(result, f)
        return result

    def delta_j(self, j):
        print ">>>>>> DELTA J >>>>>>>>>>>"
        first_part = 1 / (factorial(j) * self.v**j)
        # D inversed
        dmo = inv(self.get_d()).tolist()
        print dmo, type(dmo)
        # e with zero
        ewz = self.get_e(0)
        second_part_up = self.get_bigw(j, self.N-1)
        second_part_up = self.matrix_mult(second_part_up, dmo)
        second_part_up = self.vector_mult(second_part_up, array(ewz))
        print second_part_up, type(second_part_up)
        second_part_down = self.vector_mult(
            self.matrix_mult(self.matrix_mult(ewz, self.get_bigw(0, self.N-1)), dmo), ewz)

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

        return inv(result).tolist()


if __name__ == "__main__":
    m = MatrixParty()

    # print m.get_d()
    # m1 = [[1, 4, 44, 1], [3, 5, 14, 51], [15, 31, 8, 4], [0, 2, 33, 6]]
    # m2 = [[1, 4, 44, 1], [3, 5, 14, 51], [15, 31, 8, 4], [0, 2, 33, 6]]
    # v = [5, 5, 5, 5]
    # print m.vector_mult(m1, v)
    # print m.get_bigw(0, 51)
    for j in xrange(51):
        print j
        print m.get_f(j)

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
    # print self.matrix_mult(m1.T.tolist(), m2.T.tolist())

    # print get_d(3, get_f(N))
    # print get_e(1, 5)
