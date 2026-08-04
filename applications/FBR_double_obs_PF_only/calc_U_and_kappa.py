import numpy as np

def Ufunc(gamma, l):
  return 2*gamma/(np.pi*l)

def kappafunc(gamma, l):
  return 4*gamma*l/np.pi

def main():
  #Interfacial energy in J/m^2 (also MPa * microns)
  gamma = 0.258
  #Interfacial width (microns)
  l = 5
  U = Ufunc (gamma, l)
  kappa = kappafunc(gamma, l)
  print("U (MPa): ", U)
  print("kappa (MPa * micron^2): ", kappa)

if __name__ == "__main__":
    main()