from PythonFunctions import PI, MOD, DIV, SQR, SINE
from retireFunctions import Class_DIV
from ArrayFunctions import INVERSE
from math import sin, sqrt
from time import time
from random import random, uniform

def modulus_test():
  print('modulus_test')

  n_nums = set()
  D = PI
  N = -10.0
  for i in range(100000):
    N = newAngle(N)
    #Python modulus function can't handle a negative numerator
    if N<0.0:
      continue
    n_nums.add(N)
    if abs(N%D-MOD(N,D)) > 1.0E-9:
      print()
      print('NN = {}, DD = {}'.format(N,D))
      print(N%D)
      print(MOD(N,D))
  print("{} different N's were used for testing".format(len(n_nums)))
  pass

def auto_div(f, low_b, up_b):
  D = uniform(low_b,up_b)
  N = random()
  return abs(N/D - f(N,D))

def func_comp(f):
  bounds = [(0.0,0.1), (0.1,1.0), (1.0,1000)]

  for i in range(len(bounds)):
    for j in range(len(f)):
      print('\nFunction: {} - range {} < D < {}'.format(f[j].__name__, bounds[i][0], bounds[i][1]))
      error = 0.0
      nTimes = 100000
      t = time()
      for _ in range(nTimes):
        error += auto_div(f[j], bounds[i][0], bounds[i][1])
      print('Time:  {}'.format(time()-t))
      print('Error: {}'.format(error))
  pass

def newAngle(N):
  x = random()*720.0
  return -x if N > 0.0 else x

def trigTest():
  print('SINE/COSINE Testing:')
  add = PI/36.0
  angle = -PI
  print('       angle        SINE         SIN       error')
  for i in range(72):
    print('{:12.8f}{:12.8f}{:12.8f}{:12.8f}'.format(angle,SINE(angle),sin(angle),abs(sin(angle)-SINE(angle))))
    angle += add
  pass

def inverseTest():
  n = int(input('What is the size of the matrix?'))
  A = [[0]*n for _ in range(n)]
  for i in range(len(A)):
    print('Enter {} numbers for row {} as a single input'.format(n,i+1))
    A[i] = [float(j) for j in input().strip().split(" ")]
  B = INVERSE(A)
  pass

def menu():
  print('\nYou can test for:')
  print('PythonFunctions: 1. Modulus, 2. Division, 3. Sin')
  print('ArrayFunctions: 4. Inverse')
  print('Anything else to exit')
  return int(input())

if __name__ == "__main__":
  choose = menu()

  while (1 <= choose <= 4):
    if choose == 1:
      modulus_test()
    elif choose == 2:
      print('Division functions comparison:')
      func_comp([Class_DIV, DIV])
    elif choose == 3:
      trigTest()
    elif choose == 4:
      inverseTest()
    else:
      print('Have a nice day:)')
      break
    choose = menu()

