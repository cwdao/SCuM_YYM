import numpy as np
import matplotlib.pyplot as plt
import sympy
import sys
from sklearn.linear_model import LinearRegression
from matplotlib.ticker import MultipleLocator
from sympy import Symbol, solve

SetPerGroup = 32
SetPerGroup2 = SetPerGroup*SetPerGroup
freq_Test = 2502


x1 = 2
y1 = (2298+2287)/2
x2 = 16
y2 = (2474+2461)/2
x3 = 31
y3 = (2716+2699)/2

y_matrix = np.matrix([[y2-y1],[y3-y1]])
x_matrix = np.matrix([[pow(x2,2) - pow(x1,2), x2 - x1], [pow(x3,2) - pow(x1,2), x3 -x1]])

print(x_matrix)
print(x_matrix.I)
print(y_matrix)

para_matrix = x_matrix.I*y_matrix
Inv_equ_a = para_matrix[0][0]
Inv_equ_b = para_matrix[1][0]
Inv_equ_c = y2-Inv_equ_a*pow(x2,2)-Inv_equ_b*x2
Inv_equ = [Inv_equ_a.tolist()[0][0], Inv_equ_b.tolist()[0][0], Inv_equ_c.tolist()[0][0]]
print(Inv_equ)

x = Symbol('x')
fn = Inv_equ[0]*x**2 + Inv_equ[1]*x + Inv_equ[2] - freq_Test
index_Inv = sympy.solve(fn, x)[1]
print(index_Inv)

p11 = 2468
p12 = 2483
p21 = 2481
p22 = 2496
p11x = 16
p12x = 32
p21x = 0
p22x = 16
kp1 = (p12-p11)/(p12x-p11x)
print(kp1)
b = p11 - kp1*p11x
print(b)
setf1 = (freq_Test-b)/kp1
print(setf1) 

kp2 = (p22-p21)/(p22x-p21x)
print(kp2)
b = p21 - kp2*p21x
print(b)
setf2 = (freq_Test-b)/kp2
print(setf2)

