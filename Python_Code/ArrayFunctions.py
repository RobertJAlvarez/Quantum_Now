from PythonFunctions import *
from sys import exit
from copy import deepcopy

def print_matrix(X):
  for row in X:
    print(str(["{0:9.6f}".format(i) for i in row]).replace("'",""))
  pass

def INVERSE(AA):
  n = len(AA)
  A = [[0.0]*(2*n) for _ in range(n)]
  CC = [[0.0]*n for _ in range(n)]

  # Copy matrix AA in first n columns and an identity matrix after it
  for i in range(n):
    A[i][i+n] = 1.0
    A[i][:n] = deepcopy(AA[i][:n])
  print('\nMatrix A')
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
  print('\nInvert matrix')
  print_matrix(BB)

  # Dot product between AA and BB is the identity matrix
  for i in range(n):
    for j in range(n):
      for k in range(n):
        CC[i][j] += AA[i][k]*BB[k][j]
  print('\nIdentity matrix')
  print_matrix(CC)
  return BB

def J2x2():
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

