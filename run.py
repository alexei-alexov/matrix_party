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

    def matrix_mult_matrix(self, m1, m2):
        """
        This method multiply two same sized square matrices
        """
        length = len(m1)
        return [[self.getElement(m1, m2, row, column)
                 for column in xrange(length)] for row in xrange(length)]

    def matrix_mult_vector(self, m, v):
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

    def vector_mult_matrix(self, v, m):
        result = []
        length = len(v)
        for column in xrange(length):
            column_result = 0
            for row in xrange(length):
                print v[row], " * ", m[row][column]
                column_result += v[row] * m[row][column]
            result.append(column_result)
        return result

    def generate_array(self, size):
        """
        Generate square array filled with 0
        """
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
        result = [[1] if i == j else [0] for j in xrange(self.c+1)]
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

        return d

    def get_w(self, j):
        """
        This function return W(j)

        """
        result = self.matrix_mult_matrix(inv(self.get_f(j)).tolist(), self.get_b())
        return result

    def get_bigw(self, j, r):
        """
        This return multiplyed matrices
        """
        result = self.get_w(j)
        for k in xrange(j+1, r+1):
            f = self.get_w(k)
            result = self.matrix_mult_matrix(result, f)
        return result

    def delta_j(self, j):
        print ">>>>>> DELTA J >>>>>>>>>>>"
        first_part = (factorial(j) * self.v**j)
        # D inversed
        dmo = inv(self.get_d()).tolist()
        print "D:\n", dmo
        # e with zero
        ewz = self.get_e(0)
        e_t = (array(self.get_e(0)).T).tolist()[0]
        print "Et: ", e_t
        second_part_up = self.get_bigw(j, self.N-1)
        second_part_up = self.matrix_mult_matrix(second_part_up, dmo)
        second_part_up = self.matrix_mult_vector(second_part_up, ewz)
        print "SECOND_PART_UP:\n", second_part_up
        second_part_down = self.vector_mult_matrix(e_t, self.get_bigw(0, self.N-1))
        second_part_down = self.vector_mult_matrix(second_part_down, dmo)
        print "VECTOR_MULT_MATRIX(2): ", second_part_down
        second_part_down_result = 0
        for i in xrange(len(ewz)):
            second_part_down_result += second_part_down[i] * ewz[i][0]
        print "SECOND_PART_DOWN:\n", second_part_down_result
        down = (first_part * second_part_down_result)
        result = []
        for i in xrange(len(second_part_up)):
            result.append([second_part_up[i][0] / down])
        return result

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
    print m.delta_j(1)
    # print m.get_d()
    # v = [1, 2, 3]
    # m1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    # print m.vector_mult_matrix(v, m1)
    # print inv(m1).tolist()
    # m2 = [[1, 4, 44, 1], [3, 5, 14, 51], [15, 31, 8, 4], [0, 2, 33, 6]]

    # print m.matrix_mult_vector(m1, v)
    # print m.get_bigw(0, 51)
    # for j in xrange(51):
    #     print j
    #     print m.get_f(j)

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
    # print self.matrix_mult_matrix(m1.T.tolist(), m2.T.tolist())

    # print get_d(3, get_f(N))
    # print get_e(1, 5)
