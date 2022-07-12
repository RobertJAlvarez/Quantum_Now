from PythonFunctions import *
from sys import exit
from copy import deepcopy

"""
Author: Robert Alvarez
Date:   July 11th, 2022
"""

Matrix = list[list[float]]

def print_matrix(X: Matrix) -> None:
  """ Print a matrix row by row in scientific notation with 8 digits precision """
  for row in X:
    print(str(["{0:14.7e}".format(num) for num in row]).replace("'",""))
  pass

def INVERSE(AA: Matrix) -> Matrix:
  n = len(AA)
  A = [[0.0]*(2*n) for _ in range(n)]
  CC = [[0.0]*n for _ in range(n)]

  # Copy matrix AA in first n columns and an identity matrix after it
  for i in range(n):
    A[i][i+n] = 1.0
    A[i][:n] = deepcopy(AA[i][:n])
  print('\nMatrix A:')
  print_matrix(A)

  # Find invert matrix
  for i in range(n):
    # Find largest value of column i
    MAX = ABS(A[i][i])
    k = i
    for j in range(i+1,n):
      if ABS(A[j][i]) > MAX:
        k = j
        MAX = ABS(A[j][i])

    MAX = A[k][i]
    if MAX == 0.0:
      exit('A row is linearly dependent of one or more other rows')

    # Swap row with highest value in column i
    if k != i:
      for j in range(i,2*n):
        A[i][j], A[k][j] = A[k][j], A[i][j]

    # Normalize matrix
    for j in range(i,2*n):
      A[i][j] = DIV(A[i][j],MAX)

    # Subtract value A(j,i) to every column in row i
    for j in range(n):
      temp = A[j][i]
      if j != i:
        for k in range(i,2*n):
          A[j][k] -= temp*A[i][k]

  # Copy inverse matrix
  BB = [[0]*n for _ in range(n)]
  for i in range(n):
    BB[i][:n] = deepcopy(A[i][n:2*n])
  print('\nInvert matrix:')
  print_matrix(BB)

  # Dot product between AA and BB is the identity matrix
  for i in range(n):
    for j in range(n):
      for k in range(n):
        CC[i][j] += AA[i][k]*BB[k][j]
  print('\nIdentity matrix:')
  print_matrix(CC)
  return BB

# if A = [[d,e], [f,g]] 
# det(A) = d*g-f*e
# An eigenvalue satisfies:
# (d-E)*(g-E) - e*f = 0 <=> E^2 - E*(d+g) + (d*g-f*e) = 0
# From quadratic formula: b^2-4*a*c ->
# (D+g)^2 - 4*1*(d*g-f*e) <=> (d-g)^2 + 4*e*f
# E = answers of quadratic formula
def J2x2(H: Matrix, E: Matrix, O: Matrix) -> None:
  """ Diagonalize a 2x2 matrix using Jacobean 2 by 2 analytic diagonalization """
  rad = SQR((H[1][1] - H[2][2])*(H[1][1] - H[2][2]) + 4.0*(H[1][2]*H[2][1]))
  trc = H[1][1] + H[2][2]
  E[1] = 0.5*(trc+rad)
  E[2] = 0.5*(trc-rad)

  #Explain ortagonal matrix:
  O[1][1] = -H[1][2]
  O[2][1] =  H[1][1] - E[1]
  O[1][2] =  H[2][2] - E[2]
  O[2][2] = -H[1][2]

  for i in range(2):
    dot = SQR(O[1][i]*O[1][i] + O[2][i]*O[2][i])
    O[1][i] = DIV(O[1][i],dot)
    O[2][i] = DIV(O[2][i],dot)
  pass

def JAC2BY2GEN():
  pass

def DIAGDVR():
  pass

def sort():
  pass

def MergeIdx():
  pass

def DIAGNxN():
  pass

def LEASTSQUARE():
  pass

