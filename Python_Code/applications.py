import numpy as np
from pythonFunctions import PI, ABS, COSINE, SINE, DIV, SQR
from arrayFunctions import Matrix, print_mtx, print_EV, DIAGNxN

def STARKDVR() -> None:
  EFIELD = float(input("Welcome to start driver, EFIELD = "))

  NBS = 3
  DIP = np.zeros(shape=(NBS,NBS))
  DIP[0,2] = 0.01
  DIP[1,2] = 0.03
  for i in range(NBS):
    for j in range(i+1,NBS):
      DIP[j,i] = DIP[i,j]

  HAM = np.zeros(shape=(NBS,NBS))
  HAM[0,0] = -0.5
  HAM[1,1] = -0.125
  HAM[2,2] = -0.125

  HAM[0,2] = DIP[0,2]*EFIELD
  HAM[1,2] = DIP[1,2]*EFIELD
  for i in range(NBS):
    for j in range(i+1,NBS):
      HAM[j,i] = HAM[i,j]

  UMT = np.zeros(shape=(NBS,NBS))
  DIAGNxN(HAM, UMT)
  print_diag_mtx_info(HAM, UMT)

#At time t=0, occupy the 2s function:
#|phi_k(t)> = sum_j exp(ie(j) t)*OVR(k,j)|phi_j> = sum_jl exp(iEjt)OVR(j,k)UMT(k,l)|phi_l where E_j = lamda_j
  tau = DIV(8.*PI,ABS(HAM[0,0]))
  f = open("PLOT", "w")
  for j in range(1000):
    t = float(j)*DIV(tau,1.E3)
    calc_something(HAM, UMT, DIP, t, f)
  f.close()
  pass

def RINGDVR() -> None:
  NBS = 11

  RINGSZ = float(input("How large is your ring?"))
  iAP = int(input("Electric(1) or Magnetic(2)"))
  E = float(input("How big the field?"))

  aM = -5.
  HAM = np.zeros(shape=(NBS,NBS))
  for i in range(NBS):
    HAM[i,i] = DIV(aM*aM,RINGSZ*RINGSZ)
    aM += 1.
  
  DIP = np.zeros(shape=(NBS,NBS))
  if iAP != 2:
    for i in range(1,NBS):
      DIP[i-1,i] = DIV(E*RINGSZ,2.)
      DIP[i,i-1] = DIP[i-1,i]
  else:
    i = -1
    for j in range(-5,5):
      i += 1
      DIP[i,i] = float(j)*E

  print("Initial HAM:")
  print_mtx(HAM)

  DIP += HAM

  UMT = np.zeros(shape=(NBS,NBS))
  DIAGNxN(HAM, UMT)
  print_diag_mtx_info(HAM, UMT)

#At time t=0, occupy the 2s function:
#|phi_k(t)> = sum_j exp(ie(j) t)*OVR(k,j)|phi_j> = sum_jl exp(iEjt)OVR(j,k)UMT(k,l)|phi_l where E_j = lamda_j
  tau = DIV(8.*PI,ABS(HAM[0,0]))
  f = open("PLOT", "w")
  for j in range(-100,10000):
    if j == 0:
      HAM = np.copy(DIP)
      DIAGNxN(HAM, UMT)
      print(str(["{0:10.7f}".format(HAM[i,i]) for i in range(NBS)]).replace("'",""))
    t = float(j)*DIV(tau,5.E1)
    calc_something(HAM, UMT, DIP, t, f)
  f.close()
  pass

def BOXDVR() -> None:
  NBS = 9
  BOXSZ = float(input("Welcome to box driver, how large is your box?"))

  i = -1
  HAM = np.zeros(shape=(NBS,NBS))
  for j in range(-4,4):
    i += 1
    TWOM = 2.*float(j)*PI
    HAM[i,i] = DIV(TWOM,BOXSZ)*DIV(TWOM,BOXSZ)

  print("Initial HAM:")
  print_mtx(HAM)

  UMT = np.zeros(shape=(NBS,NBS))
  DIAGNxN(HAM, UMT)
  print_diag_mtx_info(HAM, UMT)

#At time t=0, occupy the 2s function:
#|phi_k(t)> = sum_j exp(ie(j) t)*OVR(k,j)|phi_j> = sum_jl exp(iEjt)OVR(j,k)UMT(k,l)|phi_l where E_j = lamda_j
  tau = DIV(1.,ABS(DIV(HAM[0,0],4.)))
  f = open("PLOT", "w")
  for j in range(1000):
    t = float(j)*DIV(tau,1.E3)
    calc_something(HAM, UMT, DIP, t, f)
  f.close()
  pass
  pass

def HMODVR() -> None:
  pass

def print_diag_mtx_info(HAM: Matrix, UMT: Matrix) -> None:
  """ Print Hamiltonian, Eigenvalues and Eigenvectors and 2s wavefunction in terms of eigenstates """
  NBS = len(HAM)

  print("Update HAM")
  print_mtx(HAM)

  #Prove that the inverse of UMT is the transpose of UMT
  print_EV([HAM[i,i] for i in range(NBS)], UMT)
  OVR = UMT.transpose()

  print("2s Wavefunction in terms of new eigenstates")
  print(str(["{0:10.7f}".format(OVR[1,j]) for j in range(len(OVR[0]))]).replace("'",""))
  pass

def calc_something(HAM: Matrix, UMT: Matrix, DIP: Matrix, t: float, f) -> None:
  NBS = len(HAM)
  PRD = np.zeros(shape=(NBS,NBS))
  for k in range(NBS):
    for j in range(NBS):
      ar = sum([COSINE(HAM[i,i]*t)*UMT[k,i]*UMT[j,i] for i in range(NBS)])
      ai = sum([SINE(HAM[i,i]*t)*UMT[k,i]*UMT[j,i] for i in range(NBS)])
      PRD[k,j] = ar*ar + ai*ai

# |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
# <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
  dr = di = 0.
  for i in range(NBS):
    for j in range(NBS):
      phs = (HAM[i,i] - HAM[j,j])*t
      dr += COSINE(phs)*UMT[1,i]*UMT[1,j]*DIP[i,j]
      di += SINE(phs)*UMT[1,i]*UMT[1,j]*DIP[i,j]

  temp = [t] + [num for num in PRD[:,1]] + [dr, di]
  print(str(["{0:10.7f}".format(num) for num in temp]).replace("'",""))
  temp.pop()
  temp.pop()
  temp.append(SQR(dr*dr +di*di)*100)  #dipole
  f.write(str(["{0:10.7f}".format(num) for num in temp]).replace("'",""))
  pass

def draw_plot():
  pass

def open_plot():
  pass

