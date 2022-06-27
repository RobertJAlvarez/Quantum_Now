from PythonFunctions import *
from math import *
from time import time
import random

def modulus_test():
  print('modulus_test')

  n_nums = set()
  D = PI
  N = -10.0
  for i in range(5):
    N = newAngle(N)
    n_nums.add(N)
    if abs(N%D-MOD(N,D)) > 1.0E-12:
      print()
      print('NN = {}, DD = {}'.format(N,D))
      print(N%D)
      print(MOD(N,D))
  print("{} different N's were used for testing".format(len(n_nums)))
  pass

def auto_div(f, low_b, up_b):
  j = 1
  error = 0.0
  nTimes = 100000
  t = time()
  for i in range(nTimes):
    D = random.uniform(low_b,up_b)
    N = random.random()
    error += abs(N/D - f(N,D))
  print('Time:  {}'.format(time()-t))
  print('Error: {}'.format(error))
  pass

def div_comp():
  print('Division functions comparition:')
  
  f_names = ['Class_DIV', 'DIV']
  f = [Class_DIV, DIV]
  bounds = [(0.0,0.1), (0.1,1.0), (1.0,1000)]

  for i in range(3):
    for j in range(2):
      print('\nFunction: {} - range {} < D < {}'.format(f_names[j], bounds[i][0], bounds[i][1]))
      auto_div(f[j], bounds[i][0], bounds[i][1])

  pass

def newAngle(N):
  x = random.random()*720.0
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
  #
  pass

def menu():
  print('You can test for:')
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
      div_comp()
    elif choose == 3:
      trigTest()
    elif choose == 4:
      inverseTest()
    else:
      print('Have a nice day:)')
      break
    choose = menu()

