from sys import exit

PI = 3.14159265358979

def MOD(NN,DD):
  """ NN = numerator, DD = denominator. Calculate and return |N|%|D|. """
  #Make N and D positive
  if NN<0.0 or DD<0.0:
    if NN<0.0 and DD<0.0: #If both are negative
      N = -NN
      D = -DD
      fact = 1.0
    elif NN<0.0:   #If N is negative and D positive
      N = -NN
      D = DD
      fact = -1.0
    else:
      N = NN
      D = -DD
      fact = -1.0
  else:           #If both are positive
    N = NN
    D = DD
    fact = 1.0

  #Calculate residual of N/D
  while N > D:
    N -= D

  return fact * N

def ABS(num):
  """ Calculate and return |number| """
  return -num if num < 0.0 else num

def DIV(NN,DD):
  """ NN = numerator, DD = denominator. Calculate and return N/D """
  if DD == 0.0:
    try:
      ans = NN/DD
    except ZeroDivisionError:
      exit("Can't divide by 0")

  #If DD < 0 multiply D and N by -1 so D > 0
  if DD < 0.0:
    N, D = -NN, -DD
  else:
    N, D = NN, DD

  while D > 1.0:  #Scale N and D so 0<D<1
    N *= 0.1
    D *= 0.1

  while D + 1.0E-15 <= 1.0: #Until D = 1
    #Set a factor F so D approaches 1.0
    F = 2.0 - D if D > 0.1 else 10.0
    N *= F
    D *= F
  return N

# Bhaskara approx:  temp = (PI-x)*x
#                   sin(x) = (16*temp) / (5*PI*PI - 4*temp)
# second approx:    temp = (x/PI)*(36*temp - 31)
#                   sin(x) = (temp/10)*(36*temp - 31)
# Weight average: Bhaskara -> 0.385 Second -> 0.615
#
# sin(x) approx with weight average: temp = (x/PI)*(x/PI - 1)
#   sin(x) = temp*(2.21652(temp - 31/36) - 1.5372/(1.25 + temp))
def SINE(num):
  """ Use weights in a first and second degree equation in p where p=x(180-x)/8100 to approximate sine(x) """
  x = MOD(ABS(num),PI)

  temp = DIV(x,PI)
  temp = temp*(temp-1.0)
  ans = temp*(2.21652*(temp - DIV(31.0,36.0)) - DIV(1.5372,1.25+temp))

  ans += 2.6E-5  #Shift the graph down by 2.6E-5 to minimize error
  #Adjust for negative angles
  if num > 0.0:
    return -ans if MOD(num,2.0*PI) > PI else ans
  return -ans if MOD(num,2.0*PI) > -PI else ans

def COSINE(num):
  """ Return cosine(x) by using sine(x) function """
  return SINE(DIV(PI,2.0) - num)

def SQR(num):
  """ Return the approximation  """
  xn = 1.0
  xn_1 = -0.1
  xn_2 = -0.2

  for i in range(100):
    EPS = DIV(num - xn*xn, 2.0*xn)
    xn += EPS
    if xn == xn_1:
      break
    elif xn == xn_2:
      xn = DIV(xn_1 + xn_2, 2.0)
      break
    xn_2 = xn_1
    xn_1 = xn
  return xn

"""
Extra functions for reference but never used
"""

def SINB(num):
  """ Bhaskara approximation: temp = (PI-x)x
  sin(x) = (16*temp) / (5*PI*PI - 4*temp)"""
  x = MOD(ABS(num),PI)
  temp = (PI-x)*x
  ans = (16.0*temp) / (5.0*PI*PI - 4.0*temp)

  if num > 0.0:
    if MOD(num,2.0*PI) > PI:
      ans = -ans
  else:
    if MOD(num,2.0*PI) < -PI:
      ans = -ans
  return ans

#DIVIDE (1/(1+x)) 0 < x < 0.5
def DIVIDE(num):
  n = 50
  recid = 0.0
  p = 1.0

  for i in range(n):
    recid += p
    p *= (-num)
  return recid

#DIVIDER (1/(1+x)) 0 < x < 1
def DIVIDER(num):
  if num <= 0.5:
    return DIVIDE(num)

  TWOTHIRDS = 0.666666666666667
  y = TWOTHIRDS*(num - 0.5)
  return DIVIDE(y)*TWOTHIRDS

# 0.0 < num <= 0.1
def TNYDIVIDE(num):
  P = 1.0
  R = 1.0
  recidual = 1.0
  for i in range(64):
    if num > R and num < R*2.0:
      x = P*(num-R)
      recidual = P*DIVIDER(num)
    R *= 0.5
    P *= 2.0
  return recidual

# 0.1 < num < 1.0
def LT1DIVIDE(num):
  n = 1000
  recid = 0.0
  p = 1.0

  for i in range(n):
    recid += p
    p *= num

  return recid

# Only DIV call this function and num exist between (0.0,1.0)
def RECIPROCAL(num):
  return TNYDIVIDE(num) if num <= 0.1 else LT1DIVIDE(num)

def Class_DIV(NN,DD):
  """ Approximate and return NN/DD and return """
  if DD == 1.0: #If DD = 1.0, NN/DD = NN
    return NN

  if DD == 0.0: #Check that DD /= 0.0
    try:
      ans = NN/DD
    except ZeroDivisionError:
      exit("Can't divide by 0")

  #If DD < 0 multiply D and N by -1 so D > 0
  if DD < 0.0:
    N = -NN
    D = -DD
  else:
    N = NN
    D = DD

  while D > 1.0:
    N *= 0.1
    D *= 0.1

  return N * RECIPROCAL(D)
