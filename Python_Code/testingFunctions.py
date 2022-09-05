import numpy as np
from math import sin, sqrt
from time import time
from random import random, uniform

from pythonFunctions import PI, MOD, DIV, SQR, SINE
from retireFunctions import Class_DIV
from arrayFunctions import print_mtx, INVERSE, J2X2, JAC2BY2GEN, DIAGNxN, LEASTSQUARE
from applications import STARKDVR, RINGDVR, HMODVR, BOXDVR

"""
  file:   PythonFunctions.py
  brief:  Tester driver

  This contains the basics to test every function used by the Quantum Now software

  Author: Robert Alvarez
  Date:   August 5th, 2022
  bug:    No known bugs
"""

"""
  Testing for pythonFunctions.py
"""
def modulus_test() -> None:
  print('modulus_test')

  n_nums = set()
  D = PI
  for i in range(100000):
    #Python modulus operator can't handle negative numerators so we stay with positives
    N = random()*720.
    n_nums.add(N)
    if abs(N%D-MOD(N,D)) > 1.E-12:
      print('\nNN = {}, DD = {}'.format(N,D))
      print('N%D = {}, MOD(N,D) = {}'.format(N%D, MOD(N,D)))
  print("{} different N's were used for testing".format(len(n_nums)))
  pass

def auto_div(f: callable, low_b: float, up_b: float) -> float:
  D = uniform(low_b,up_b)
  N = random()
  return abs(N/D - f(N,D))

def auto_sqr(f: callable, low_b: float, up_b: float) -> float:
  num = uniform(low_b,up_b)
  return abs(sqrt(num) - SQR(num))

def func_comp(f: list[callable], auto_f: callable) -> None:
  bounds = [(0.0,0.1), (0.1,1.0), (1.0,1000)]
  nTimes = 100000

  for i in range(len(bounds)):
    for j in range(len(f)):
      print('\nFunction: {} - range {} < D < {}'.format(f[j].__name__, bounds[i][0], bounds[i][1]))
      error = 0.
      t = time()
      for _ in range(nTimes):
        error += auto_f(f[j], bounds[i][0], bounds[i][1])
      print('Time:  {}'.format(time()-t))
      print('Error: {}'.format(error))
  pass

def trigTest() -> None:
  print('SINE/COSINE Testing:')
  add = PI/36.
  angle = -PI
  print('       angle        SINE         SIN       error')
  for _ in range(72):
    print('{:12.8f}{:12.8f}{:12.8f}{:12.8f}'.format(angle,SINE(angle),sin(angle),abs(sin(angle)-SINE(angle))))
    angle += add
  pass

"""
  Testing for arrayFuncitons.py
"""
def inverseTest() -> None:
  n = int(input('What is the size of the matrix?'))
  A = np.zeros(shape=(n,n))
  for i in range(len(A)):
    print('Enter {} numbers, separated by a white space and in a single line, for row {} as a single input'.format(n,i+1))
    A[i] = np.array([float(j) for j in input().strip().split(" ")])
  B = INVERSE(A)
  pass

def solve_2by2() -> None:
  pass

def solve_gen() -> None:
  pass

def DIAGDVR() -> None:
  NBS = 7
  G = -1.0  # Ground
  X = -0.5  # Exited
  P =  0.1  # Perturbation
  TXR = 0.2 # Transfer
  UMT = np.zeros(shape=(NBS,NBS))
  HAM = np.zeros(shape=(NBS,NBS))

  HAM[0,0] = HAM[3,3] = HAM[6,6] = G
  HAM[1,1] = HAM[2,2] = HAM[4,4] = HAM[5,5] = X

  HAM[0,1] = HAM[2,3] = HAM[3,4] = HAM[5,6] = P
  HAM[1,2] = HAM[4,5] = TXR

  HAM[1,0] = HAM[3,2] = HAM[4,3] = HAM[6,5] = P
  HAM[2,1] = HAM[5,4] = TXR

  print('Original Hamiltonian:')
  print_mtx(HAM)

  DIAGNxN(HAM, UMT)

  print('Updated Hamiltonian:')
  print_mtx(HAM)
  pass

def LSA_test() -> None:
  pass

def menu() -> int:
  print('\nYou can test for:')
  print('PythonFunctions: 1. Modulus, 2. Division, 3. Square root, 4. Sin')
  print('ArrayFunctions: 5. Inverse, 6. J2X2, 7. JAC2BY2GEN, 8. DIAGNxN, 9. LEASTSQUARE')
  print('Applications: 10. STARKDVR, 11. RINGDVR, 12. BOXDVR, 13. HMODVR')
  print('Anything else to exit')
  return int(input())

if __name__ == "__main__":
  choose = menu()

  while (1 <= choose <= 12):
    if choose == 1:
      modulus_test()
    elif choose == 2:
      print('Division functions comparison:')
      func_comp([Class_DIV, DIV], auto_div)
    elif choose == 3:
      func_comp([SQR], auto_sqr)
    elif choose == 4:
      trigTest()
    elif choose == 5:
      inverseTest()
    elif choose == 6:
      solve_2by2()
    elif choose == 7:
      solve_gen()
    elif choose == 8:
      DIAGDVR()
    elif choose == 9:
      LSA_test()
    elif choose == 10:
      STARKDVR()
    elif choose == 11:
      RINGDVR()
    elif choose == 12:
      BOXDVR()
    elif choose == 13:
      HMODVR()
    else:
      print('Have a nice day:)')
      break

    if 10 <= input <= 13:
      draw_plot()
      open_plot()

    choose = menu()

