import numpy as np
from pythonFunctions import PI, ABS, COSINE, SINE, DIV, SQR
from arrayFunctions import Matrix, print_mtx, DIAGNxN

def STARKDVR() -> None:
  pass

def RINGDVR() -> None:
  pass

def BOXDVR() -> None:
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

def calc_something(HAM: Matrix, UMT: Matrix, DIP: Matrix, t: float) -> None:
  NBS = len(HAM)
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
      dr += COSINE(phs)*UMT(1,i)*UMT(1,j)*DIP[i,j]
      di += SINE(phs)*UMT(1,i)*UMT(1,j)*DIP[i,j]

  print('Finish printing information')
  pass

def draw_plot():
  pass

def open_plot():
  pass

