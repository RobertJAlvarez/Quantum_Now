from PythonFunctions import ABS, DIV, SQR
from sys import exit

"""
Author: Robert Alvarez
Date:   July 11th, 2022
"""

Vector = list[float]
Matrix = list[list[float]]

def print_mtx(X: Matrix) -> None:
  """ Print a matrix row by row in scientific notation with 8 digits precision """
  for row in X:
    print(str(["{0:14.7e}".format(num) for num in row]).replace("'",""))
  pass

def INVERSE(AA: Matrix) -> Matrix:
  n = len(AA)

  # Copy matrix AA in first n columns and an identity matrix after it
  A = [[num for num in row]+[1 if i == j else 0 for j in range(n)] for i,row in enumerate(AA)]
  print('\nMatrix A:')
  print_mtx(A)

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
    if MAX == 0.:
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
  BB = [[num for num in row[n:]] for row in A]
  print('\nInvert matrix:')
  print_mtx(BB)

  # Multiplication of A and A inverse = identity matrix
  CC = [[0.]*n for _ in range(n)]
  for i in range(n):
    for j in range(n):
      for k in range(n):
        CC[i][j] += AA[i][k]*BB[k][j]

  print('\nIdentity matrix:')
  print_mtx(CC)
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
  rad = SQR((H[1][1] - H[2][2])*(H[1][1] - H[2][2]) + 4.*(H[1][2]*H[2][1]))
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

def JAC2BY2GEN(H: Matrix, O: Matrix, V: Matrix, E: Vector) -> None:
  A = O(1,1)*O(2,2) - O(1,2)*(2,1)
  if A < 1.E-15:
    exit('Non positive overlap matrix')

  B = -( H(1,1)*O(2,2) - O(1,2)*H(2,1) + O(1,1)*H(2,2) - H(1,2)*O(2,1) )
  C = H(1,1)*H(2,2) - H(1,2)*H(2,1)

  trc = DIV(-B,A)
  rad = DIV(SQR(B*B - 4.*A*C),A)

  E[1] = 0.5*(trc + rad)
  E[2] = 0.5*(trc - rad)

  print("Eigenvalues: {0:14.7e} {0:14.7e}".format(E[1],E[2]))

  #Calculate eigenvectors
  T = [[0.]*2 for _ in range(2)]  #Eigenvectors
  for k in range(2):
    for i in range(2):
      for j in range(2):
        T[i][j] = H[i][j] - E[k]*O[i][j]
    V[1][k] = -T[k][2]
    V[2][k] =  T[k][1]

  # <V_1 | V_2> = 0 ?
  D = [[0.]*2 for _ in range(2)]
  for iTry in range(3):
    for i in range(2):
      for k in range(2):
        D[i][j] = 0.
        for k in range(2):
          for l in range(2):
            D[i][j] = D[i][j] + V[k][i]*V[l][j]*(O[k][l] if iTry <= 1 else H[k][l])

    if iTry == 0:
      for i in range(2):
        for k in range(2):
          V[k][i] = DIV(V[k][i],SQR(D[i][i]))
  pass

def SortIdex(mtx: Matrix, PRD: Matrix) -> None:
  if len(mtx[0]) > 1:
    mid = len(mtx[0])//2  #Find middle of columns

    L = [[num for num in row[:mid]] for row in mtx]
    R = [[num for num in row[mid:]] for row in mtx]

    SortIdex(L, PRD)  #Sort first half
    SortIdex(R, PRD)  #Sort second half

    i = j = k = 0

    while i < len(L[0]) and j < len(R[0]):
      if ABS(PRD[L[0][i],L[1][i]]) > ABS(PRD[R[0][j],R[1][j]]):
        mtx[0][k], mtx[1][k] = L[0][i], L[1][i]
        i += 1
      else:
        mtx[0][k], mtx[1][k] = R[0][j], R[1][j]
        j += 1
      k += 1

    while i < len(L[0]):
      mtx[0][k], mtx[1][k] = L[0][i], L[1][i]
      i += 1
      k += 1

    while j < len(L[0]):
      mtx[0][k], mtx[1][k] = R[0][j], R[1][j]
      j += 1
      k += 1
  pass

def DIAGNxN():
  pass

def LEASTSQUARE():
  pass

