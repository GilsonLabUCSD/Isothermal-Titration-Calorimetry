import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import sys

if len(sys.argv) < 7:
  print "ARGUMENT SYNTAX:\n [1]SYRINGE CONCENTRATION(mM)\n [2]CELL CONCENTRATION(mM)\n [3]INJECTION VOLUME(uL)\n [4]NUMBER OF INJECTIONS\n [5]BINDING ENTHALPY(kcal/mol)\n [6]EQUILIBRIUM CONSTANT\n [7]CONCENTRATION PRINT; IF 1 = PRINT CONCENTRATION AT EACH TITRATION STEP AND NORMALIZED INJECTIONS HEATS, IF 0 = PRINT RAW INJECTION HEATS ONLY"

V0 = 1.4614/1000  # Cell Volume (L)
X0 = float(sys.argv[1])/1000  # Injectant Concentration
M0 = float(sys.argv[2])/1000  # Cell Concentration
dV = float(sys.argv[3])/1000000  # Injection Size
Ninj = int(sys.argv[4])  # Number of injections

dH = float(sys.argv[5])*1000 # Binding Enthalpy (kcal/mol)
K = float(sys.argv[6])  # Equilibrium Constant (M^-1)
N = 1  # Binding Stoich

CONCPRNT = float(sys.argv[7]) # If 0, print raw heats only; if 1, print raw heats with concentrations at each titration step

QQ = np.zeros([28], np.float64)
Q = np.zeros([Ninj+1], np.float64) # Total Heat Array Initialized
Q[0]=0.0 # Total Heat before injections is zero
for i in range(1,Ninj+1): # Loop over injections
  X = (i*dV*X0/V0)*( (1 - i*dV/(2*V0) ) ) #New Injectant Concentration 
  M = M0*( (1 - i*dV/(2*V0) ) / (1 + i*dV/(2*V0) ) ) # New Cell Molecule Concentration
  Q[i] = (N*M*dH*V0/2)*( 1 + (X/(N*M)) + 1/(N*K*M) - np.sqrt( ((1 + X/(N*M) + 1/(N*K*M))**2) - 4*X/(N*M) ) ) # Total Heat
  dQ = (Q[i] + (dV/V0)*( (Q[i] + Q[i-1])/2 ) - Q[i-1])/(dV*X0) # Change in heat normalized by amount of injectant
  if CONCPRNT == 0:
    print dQ*(dV*X0)
  if CONCPRNT == 1:
    print X/M, dQ
