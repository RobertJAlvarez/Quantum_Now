from pythonFunctions import PI, MOD, ABS, DIV

"""
Extra functions for reference but never used
"""
def SINB(num: float) -> float:
  """ Bhaskara approximation: temp = (PI-x)x; sin(x) = (16*temp) / (5*PI*PI - 4*temp)"""
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

def Class_SQR(num: float) -> float:
  """ Approximation the square root of a number and return it """
  xn = 1.0
  xn_1 = -0.1
  xn_2 = -0.2

  for _ in range(100):
    xn += DIV(num - xn*xn, 2.0*xn)
    #Exit loop if values are cycling
    if xn == xn_1:
      break
    elif xn == xn_2:
      xn = DIV(xn_1 + xn_2, 2.0)
      break
    #Update x_ns
    xn_2 = xn_1
    xn_1 = xn
  return xn

#DIVIDE (1/(1+x)) 0 < x < 0.5
def DIVIDE(num: float) -> float:
  n = 50
  recid = 0.0
  p = 1.0

  for i in range(n):
    recid += p
    p *= (-num)
  return recid

#DIVIDER (1/(1+x)) 0 < x < 1
def DIVIDER(num: float) -> float:
  if num <= 0.5:
    return DIVIDE(num)

  TWOTHIRDS = 0.666666666666667
  y = TWOTHIRDS*(num - 0.5)
  return DIVIDE(y)*TWOTHIRDS

# 0.0 < num <= 0.1
def TNYDIVIDE(num: float) -> float:
  P = 1.0
  R = 1.0
  recid = 0.0
  for i in range(64):
    if num > R and num < R*2.0:
      x = (num-R)*P
      recid = P*DIVIDER(x)
    R *= 0.5
    P *= 2.0
  return recid

# 0.1 < num < 1.0
def LT1DIVIDE(num: float) -> float:
  n = 1000
  recid = 0.0
  p = 1.0

  for i in range(n):
    recid += p
    p *= num

  return recid

# Only Class_DIV call this function and num exist between (0.0,1.0)
def RECIPROCAL(num: float) -> float:
  return TNYDIVIDE(num) if num <= 0.1 else LT1DIVIDE(1.0-num)

def Class_DIV(NN: float, DD: float) -> float:
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
