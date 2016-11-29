#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from math import factorial

import numpy as np
from numpy import array
from numpy.linalg import inv


class MatrixParty(object):

    # PREDEFINED
    # v, m, l, c, N
    c = 2
    l = 3
    m = 1
    v = 2
    N = 2

    def __init__(self, c, l, m, v, N):
        self.c = int(c)
        self.l = int(l)
        self.m = int(m)
        self.v = int(v)
        self.N = int(N)

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
                column_result += v[row] * m[row][column]
            result.append(column_result)
        return result

    def vector_mult_number(self, m, n):
        length = len(m)
        for row in xrange(length):
            m[row][0] *= n
        return m

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
        # print "DELTA J[", j, "]"
        first_part = 1.0 / (factorial(j) * self.v**j)
        # print "FIRST PART: ", first_part
        # D inversed
        dmo = inv(self.get_d()).tolist()
        # print "DMO: ", dmo
        # e with zero
        ewz = self.get_e(0)
        # print "E0: ", ewz
        e_t = (array(self.get_e(0)).T).tolist()[0]
        # print "E0t: ", e_t
        if j == self.N:
            second_part_up = self.matrix_mult_vector(dmo, ewz)
            # print "J == N, UP PART: ", second_part_up
        else:
            second_part_up = self.get_bigw(j, self.N-1)
            second_part_up = self.matrix_mult_matrix(second_part_up, dmo)
            second_part_up = self.matrix_mult_vector(second_part_up, ewz)
            # print "J != N, UP PART: ", second_part_up

        second_part_down = self.vector_mult_matrix(e_t, self.get_bigw(0, self.N-1))
        second_part_down = self.vector_mult_matrix(second_part_down, dmo)
        # print "SECOND PART DOWN: ", second_part_down
        second_part_down_result = 0
        for i in xrange(len(ewz)):
            second_part_down_result += second_part_down[i] * ewz[i][0]
        # print "SECOND PART DOWN RESULT: ", second_part_down_result
        down = (first_part * (1 / second_part_down_result))
        # print "DOWN: ", down
        return self.vector_mult_number(second_part_up, down)

    def pi_zero(self):
        result = 0
        for j in xrange(0, self.N+1):
            delta_j = self.delta_j(j)
            # print "J: ", j, " = ", delta_j
            for i in xrange(len(delta_j)):
                result += delta_j[i][0]
        return 1 / result

    def pi(self, j):
        return self.vector_mult_number(self.delta_j(j), self.pi_zero())

    def get_p(self):
        result = self.delta_j(0)
        for j in xrange(1, self.N+1):
            result += self.delta_j(j)

        return 1 / result

    def get_big_pi(self):
        result = []
        for j in xrange(self.N+1):
            result.append(str(self.pi(j)))
        return "\n".join(result)

if __name__ == "__main__":
    if sys.argv[1] == "--help":
        print "This is matrix party program."
        print "ARGUMENTS: c l m v N"
    else:
        if len(sys.argv) != 6:
            print "No enought arguments!"
        else:
            sys.argv.pop(0)
            m = MatrixParty(*sys.argv)
            print m.get_big_pi()
