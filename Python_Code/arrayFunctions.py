from pythonFunctions import ABS, DIV, SQR
import numpy as np
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
    print(str(["{0:10.7f}".format(num) for num in row]).replace("'",""))
  pass

def print_EV(E: Vector, O: Matrix) -> None:
  """ Print Eigenvalues (E) and eigenvectors (V) with format E_i: V_i """
  print("Eigenvalues and eigenvectors")
  for i in range(len(E)):
    print(" {0:10.7f}: ".format(E[i]) + str(["{0:10.7f}".format(num) for num in O[:,i]]).replace("'","")[1:-1:1])
  pass

def INVERSE(AA: Matrix) -> Matrix:
  """ Receives matrix A and returns A inverse """
  n = len(AA)

  # Copy matrix AA in first n columns and an identity matrix after it
  A = np.append(AA,np.identity(n),axis=1)
  print('\nMatrix A:')
  print_mtx(A)

  # Find invert matrix
  for i in range(n):
    # Find largest value of column i
    k = i
    MAX = ABS(A[k,i])
    for j in range(i+1,n):
      if ABS(A[j,i]) > MAX:
        k = j
        MAX = ABS(A[k,i])

    MAX = A[k,i]
    if MAX < 1E-15:
      exit('A row is linearly dependent of one or more other rows')

    if k != i:  # Swap row with highest value in column i
      A[[i,k],i:] = A[[k,i],i:]

    # Normalize matrix
    A[i,i:] = [DIV(A[i,j],MAX) for j in range(i,2*n)]

    # Subtract value A[j,i] to every column in row >= i
    for j in range(n):
      temp = A[j,i]
      if j != i:
        A[j,i:] -= temp*A[i,i:]

  # Copy inverse matrix
  BB = np.array(A[:,n:])
  print('\nInvert matrix:')
  print_mtx(BB)

  # Multiplication of A and A inverse = identity matrix
  CC = np.zeros(shape=(n,n))
  for i in range(n):
    for j in range(n):
      for k in range(n):
        CC[i,j] += AA[i,k]*BB[k,j]

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
def J2X2(H: Matrix, E: Matrix, O: Matrix) -> None:
  """ Diagonalize a 2x2 matrix using Jacobean 2 by 2 analytic diagonalization """
  rad = SQR((H[0,0] - H[1,1])*(H[0,0] - H[1,1]) + 4.*(H[0,1]*H[1,0]))
  trc = H[0,0] + H[1,1]
  E[0], E[1] = 0.5*(trc+rad), 0.5*(trc-rad)

  #Explain ortagonal matrix:
  O[0,0] = -H[0,1]
  O[1,0] =  H[0,0] - E[0]
  O[0,1] =  H[1,1] - E[1]
  O[1,1] = -H[0,1]

  for i in range(2):
    dot = SQR(O[0,i]*O[0,i] + O[1,i]*O[1,i])
    O[0,i], O[1,i] = DIV(O[0,i],dot), DIV(O[1,i],dot)
  pass

def JAC2BY2GEN(H: Matrix, O: Matrix, V: Matrix, E: Vector) -> None:
  A = O[0,0]*O[1,1] - O[0,1]*O[1,0]
  if A < 1.E-15:
    exit('Non positive overlap matrix')

  B = -( H[0,0]*O[1,1] - O[0,1]*H[1,0] + O[0,0]*H[1,1] - H[0,1]*O[1,0] )
  C = H[0,0]*H[1,1] - H[0,1]*H[1,0]

  trc = DIV(-B,A)
  rad = DIV(SQR(B*B - 4.*A*C),A)

  E[0], E[1] = 0.5*(trc+rad), 0.5*(trc-rad)

  print("Eigenvalues: {0:10.7f} {1:10.7f}".format(E[1],E[2]))

  #Calculate eigenvectors
  T = np.zeros(shape=(2,2)) #Eigenvectors
  for k in range(2):
    for i in range(2):
      for j in range(2):
        T[i,j] = H[i,j] - E[k]*O[i,j]
    V[0,k], V[1,k] = -T[k,1], T[k,0]

  # <V_1 | V_2> = 0 ?
  D = np.zeros(shape=(2,2))
  for iTry in range(3):
    for i in range(2):
      for k in range(2):
        D[i,j] = 0.
        for k in range(2):
          for l in range(2):
            D[i,j] = D[i,j] + V[k,i]*V[l,j]*(O[k,l] if iTry <= 1 else H[k,l])

    if iTry == 0:
      for i in range(2):
        for k in range(2):
          V[k,i] = DIV(V[k,i],SQR(D[i,i]))
  pass

def SortIdx(mtx: Matrix, PRD: Matrix) -> None:
  """ Sort pair of indices in mtx[:,i] from largest to smallest by using the value at PRD[mtx[0,i], mtc[1,i]]."""
  if len(mtx[0]) > 1:
    mid = len(mtx[0])//2  #Find middle of columns

    L = np.array(mtx[:,:mid])
    R = np.array(mtx[:,mid:])

    SortIdx(L, PRD)  #Sort first half
    SortIdx(R, PRD)  #Sort second half

    i = j = k = 0

    while i < len(L[0]) and j < len(R[0]):
      if ABS(PRD[L[0,i],L[1,i]]) >= ABS(PRD[R[0,j],R[1,j]]):
        mtx[0,k], mtx[1,k] = L[0,i], L[1,i]
        i += 1
      else:
        mtx[0,k], mtx[1,k] = R[0,j], R[1,j]
        j += 1
      k += 1

    while i < len(L[0]):
      mtx[0,k], mtx[1,k] = L[0,i], L[1,i]
      i += 1
      k += 1

    while j < len(L[0]):
      mtx[0,k], mtx[1,k] = R[0,j], R[1,j]
      j += 1
      k += 1
  pass

def DIAGNxN(HAM: Matrix, UMT: Matrix) -> None:
  """ Receive matrix HAM (Hamiltonian matrix) and UMT (Unitary matrix). HAM gets diagonalize and UMT keep the eigenvectors"""
  NBS = len(HAM[0])
  np.fill_diagonal(UMT,1.)  # Make UMT a unitary matrix
  PRD = np.array(HAM)       # Copy HAM into PRD

  MXIT = NBS*NBS*2
  for iTry in range(MXIT):
    ERROLD = 0.
    idxAll = [[] for _ in range(2)]

    for i in range(NBS):
      ERROLD += sum(PRD[i+1:,i]*PRD[i+1:,i])
      for j in range(i+1,NBS):
        if ABS(PRD[i,j]) > 1.E-10:
          idxAll[0].append(i)
          idxAll[1].append(j)

#    print('All indexes before sorting:')
#    for i in range(len(idxAll[0])):
#      print("{0:2d} {1:2d} {2:14.7e}".format(idxAll[0,i], idxAll[1,i], PRD[idxAll[0,i],idxAll[1,i]]))

    idxAll = np.array(idxAll)
    SortIdx(idxAll, PRD)

#    print('All indexes:')
#    for i in range(len(idxAll[0])):
#      print("{0:2d} {1:2d} {2:14.7e}".format(idxAll[0,i], idxAll[1,i], PRD[idxAll[0,i],idxAll[1,i]]))

    useIdx = set()
    idx = [[] for _ in range(2)]
    for j in range(len(idxAll[0])):
      if idxAll[0,j] not in useIdx and idxAll[1,j] not in useIdx:
        useIdx.add(idxAll[0,j])
        useIdx.add(idxAll[1,j])
        idx[0].append(idxAll[0,j])
        idx[1].append(idxAll[1,j])
        if len(idx[0]) >= len(useIdx):
          break

#    print('Non repetitive indexes with highest values:')
#    for i in range(len(idx[0])):
#      print("{0:2d} {1:2d} {2:14.7e}".format(idx[0][i], idx[1][i], PRD[idx[0][i],idx[1][i]]))

    k, l = idx[0][0], idx[1][0]

    H = np.array([[PRD[k,k], PRD[k,l]], [PRD[l,k], PRD[l,l]]])
    E = np.zeros(shape=(2))
    O = np.zeros(shape=(2,2))
    J2X2(H, E, O)

    #print_EV(E, O)

    PRD = np.array(UMT)
    PRD[:,k] = UMT[:,k]*O[0,0] + UMT[:,l]*O[1,0]
    PRD[:,l] = UMT[:,k]*O[0,1] + UMT[:,l]*O[1,1]

    for i in range(NBS):
      for j in range(NBS):
        UMT[i,j] = PRD[i,j]

    SPC = np.identity(NBS)  # Start with identity matrix
    SPC[k,k], SPC[l,k], SPC[l,l], SPC[k,l] = O[0,0], O[1,0], O[1,1], O[0,1]
    for i in range(NBS):
      for j in range(NBS):
        SPC[j,i] = 0.
        for k in range(NBS):
          SPC[j,i] += UMT[k,i]*HAM[k,j]

    PRD = np.zeros(shape=(NBS,NBS))
    for i in range(NBS):
      for j in range(NBS):
        for k in range(NBS):
          PRD[j,i] += UMT[k,j]*SPC[k,i]

    #print_mtx(PRD)

    ERRNW = 0.
    for i in range(NBS):
      ERRNW += sum(PRD[i+1:,i]*PRD[i+1:,i])
    print("{0:2d} {1:7.5e} {2:7.5e}".format(iTry, ERRNW, ERROLD))

    if ERRNW < 1.E-15:
      break

  if iTry == MXIT:
    print("Warning: No convergence")

  print("{0:2d} {1:2d} {2:7.5f} Diag Eff".format(iTry, NBS, DIV(iTry, NBS*NBS)))

  for i in range(NBS):
    for j in range(NBS):
      HAM[i,j] = PRD[i,j]
  pass

def LEASTSQUARE() -> None:
  pass

