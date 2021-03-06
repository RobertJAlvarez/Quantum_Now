from PythonFunctions import PI, MOD, DIV, SQR, SINE
from retireFunctions import Class_DIV
from ArrayFunctions import INVERSE, J2x2, JAC2BY2GEN, DIAGNxN, LEASTSQUARE, print_matrix
from Applications import STARKDVR, RINGDVR, HMODVR, BOXDVR
from math import sin, sqrt
from time import time
from random import random, uniform

"""
  file:   PythonFunctions.py
  breif:  Tester driver

  This contains the basics to test every function used by the Quantum Now software

  Author: Robert Alvarez
  bug:    No known bugs
"""

"""
  Testing for PythonFunctions.py
"""

def modulus_test() -> None:
  print('modulus_test')

  n_nums = set()
  D = PI
  for i in range(100000):
    #Python modulus operator can't handle negative numerators so we stay with positives
    N = random()*720.
    n_nums.add(N)
    if abs(N%D-MOD(N,D)) > 1.E-9:
      print()
      print('NN = {}, DD = {}'.format(N,D))
      print(N%D)
      print(MOD(N,D))
  print("{} different N's were used for testing".format(len(n_nums)))
  pass

def auto_div(f: callable, low_b: float, up_b: float) -> float:
  D = uniform(low_b,up_b)
  N = random()
  return abs(N/D - f(N,D))

def func_comp(f: list[callable]) -> None:
  bounds = [(0.0,0.1), (0.1,1.0), (1.0,1000)]
  nTimes = 100000

  for i in range(len(bounds)):
    for j in range(len(f)):
      print('\nFunction: {} - range {} < D < {}'.format(f[j].__name__, bounds[i][0], bounds[i][1]))
      error = 0.
      t = time()
      for _ in range(nTimes):
        error += auto_div(f[j], bounds[i][0], bounds[i][1])
      print('Time:  {}'.format(time()-t))
      print('Error: {}'.format(error))
  pass

def trigTest() -> None:
  print('SINE/COSINE Testing:')
  add = PI/36.
  angle = -PI
  print('       angle        SINE         SIN       error')
  for i in range(72):
    print('{:12.8f}{:12.8f}{:12.8f}{:12.8f}'.format(angle,SINE(angle),sin(angle),abs(sin(angle)-SINE(angle))))
    angle += add
  pass

"""
  Testing for ArrayFuncitons.py
"""

def inverseTest() -> None:
  n = int(input('What is the size of the matrix?'))
  A = [[0.]*n for _ in range(n)]
  for i in range(len(A)):
    print('Enter {} numbers for row {} as a single input'.format(n,i+1))
    A[i] = [float(j) for j in input().strip().split(" ")]
  B = INVERSE(A)
  pass

def J2x2_test() -> None:
  pass

def gen_J2x2_test() -> None:
  pass

def DIAGDVR() -> None:
  pass

def LSA_test() -> None:
  pass

def menu() -> int:
  print('\nYou can test for:')
  print('PythonFunctions: 1. Modulus, 2. Division, 3. Sin')
  print('ArrayFunctions: 4. Inverse, 5. J2x2, 6. JAC2BY2GEN, 7. DIAGNxN, 8. LEASTSQUARE')
  print('Applications: 9. STARKDVR, 10. RINGDVR, 11. BOXDVR, 12. HMODVR')
  print('Anything else to exit')
  return int(input())

if __name__ == "__main__":
  choose = menu()

  while (1 <= choose <= 12):
    if choose == 1:
      modulus_test()
    elif choose == 2:
      print('Division functions comparison:')
      func_comp([Class_DIV, DIV])
    elif choose == 3:
      trigTest()
    elif choose == 4:
      inverseTest()
    elif choose == 5:
      J2x2_test()
    elif choose == 6:
      gen_J2x2_test()
    elif choose == 7:
      DIAGDVR()
    elif choose == 8:
      LSA_test()
    elif choose == 9:
      STARKDVR()
    elif choose == 10:
      RINGDVR()
    elif choose == 11:
      BOXDVR()
    elif choose == 12:
      HMODVR()
    else:
      print('Have a nice day:)')
      break
    choose = menu()

