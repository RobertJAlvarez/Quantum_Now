from sys import exit

"""
Author: Robert Alvarez
Date:   Jun 26th, 2022
"""

PI = 3.14159265358979

def MOD(NN: float, DD: float) -> float:
  """ NN = numerator, DD = denominator. Calculate and return |N|%|D|. """

  #Make N and D positive
  N, fact = (-NN, -1.) if NN<0. else (NN, 1.)
  D, fact = (-DD, -fact) if DD<0. else (DD, fact)

  #Calculate residual of N/D
  while N > D:
    N -= D
  return fact * N

def ABS(num: float) -> float:
  """ Calculate and return |number| """
  return -num if num < 0. else num

def DIV(NN: float, DD: float) -> float:
  """ Use Goldschmidt division to calculate and return NN/DD,
  where NN = numerator, DD = denominator """
  if abs(DD) < 1.E-15:
    exit("Can't divide by 0")

  #If DD < 0 multiply D and N by -1 so D > 0
  N, D = (-NN, -DD) if DD < 0. else (NN,DD)

  while D > 1.:  #Scale N and D so 0<D<1
    N *= 0.1
    D *= 0.1

  while D + 1.E-15 <= 1.: #Until D = 1
    F = 2. - D  #Set a factor F so D approaches 1.0
    N *= F
    D *= F
  return N

def SQR(num: float) -> float:
  """ Approximation and return the square root of a number """
  xn = 1.
  for _ in range(25):
    xn += DIV(num - xn*xn, 2.*xn)
  return xn

def SINE(num: float) -> float:
  """ Use weights in a first and second degree equation in p where
  p=x(180-x)/8100 to approximate sine(x) """
  x = MOD(ABS(num),PI)

  temp = DIV(x,PI)
  temp = temp*(temp-1.)
  ans = temp*(2.21652*(temp - DIV(31.,36.)) - DIV(1.5372,1.25+temp))

  ans += 2.563E-5  #Shift the graph down by 2.6E-5 to minimize error
  #Adjust for negative angles
  if num > 0.:
    return -ans if MOD(num,2.*PI) > PI else ans
  return -ans if MOD(num,2.*PI) > -PI else ans

def COSINE(num: float) -> float:
  """ Return cosine(x) by using sine(x) function """
  return SINE(0.5*PI - num)

