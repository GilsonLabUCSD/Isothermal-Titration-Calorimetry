import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import sys
#print(sys.version)

if len(sys.argv) < 10:
  print "ARGUMENT SYNTAX:\n [1]INJECTION VOLUME(uL)\n [2]NUMBER OF FITTING CYCLES\n [3]SYRINGE CONCENTRATION(mM)\n [4]CELL CONCENTRATION(mM)\n [5]SYRINGE CONCENTRATION PERCENT UNCERTAINTY\n [6]CELL CONCENTRATION PERCENT UNCERTAINTY\n [7]HEAT PERCENT UNCERTAINTY\n [8]BASELINE UNCERTAINTY (ucal) [9]N FIXED (TRUE = 1, FALSE = 0)\n [10]FILE CONTAINING RAW INTEGRATED INJECTION HEATS IN ONE COLUMN\n"
  sys.exit()

### ITC Settings
V0 = 1.4614/1000  # Cell Volume (L)
dV = float(sys.argv[1])/1000000  # Injection Size (L)
Rinj = 0  # Number of injections removed from the beginning of the ITC data
CYC = int(sys.argv[2])

### Build empty arrays for concentrations drawn off distribution (based on inputted % uncertainty)
Xx = np.zeros([CYC+1], np.float64)
Mm = np.zeros([CYC+1], np.float64)

### Define Input Parameters for experimental data
Xo = float(sys.argv[3])/1000 # Syringe concentration used for experimental data
Mo = float(sys.argv[4])/1000 # Cell concentration used for experimental data 
Syr_Err = float(sys.argv[5])/100 # Percent uncertainty in syringe concentration
Cel_Err = float(sys.argv[6])/100 # Percent uncertainty in cell concentration
Heat_Err = float(sys.argv[7])/100 # Percent uncertainty in heat error (peak to peak error)
Base_Err = float(sys.argv[8])/1000000 #Uncertainty in baseline (ucal)
FIX = float(sys.argv[9]) ### Checks if fixing N is true or not

### Generate Distribution of Concentrations Based on Inputted Uncertainties
for i in range(CYC):
  Xx[i]=np.random.normal(Xo,abs(Syr_Err*Xo))
  Mm[i]=np.random.normal(Mo,abs(Cel_Err*Mo))

### Read In ITC data
ITC = np.loadtxt(sys.argv[10], dtype=np.float64)#, delimiter=' ', unpack=True)
Ninj = len(ITC) + Rinj # Get number of injections from data file.  Add number of injections that were removed

### Fitting Function
if FIX == 0: ### IF N IS ALLOWED TO FLOAT
   def func(XMa, dH, K, N):
     Q = np.zeros([Ninj+1], np.float64) # Total Heat Array Initialized
     Q[0]=0.0 # Total Heat before injections is zero
     dQ = np.zeros([Ninj-Rinj], np.float64) # Heat Change Array (but not for injections without data!)
     for i in range(1,Ninj+1): # Loop over injections
         X = XMa[0,i]  # New Injectant Concentration
         M = XMa[1,i]  # New Cell Molecule Concentration
         Q[i] = (N*M*dH*V0/2)*( 1 + (X/(N*M)) + 1/(N*K*M) - np.sqrt( ((1 + X/(N*M) + 1/(N*K*M))**2) - 4*X/(N*M) ) ) # Total Heat
         if i > Rinj: # If initial injections were omitted from the ITC data, omit them here too
           dQ[i-1-Rinj] = (Q[i] + (dV/V0)*( (Q[i] + Q[i-1])/2 ) - Q[i-1])/(dV*X0) # Change in heat normalized by amount of injectant
     return dQ

if FIX == 1: ### IF N IS FIXED
  def func(XMa, dH, K):
    Q = np.zeros([Ninj+1], np.float64) # Total Heat Array Initialized
    Q[0]=0.0 # Total Heat before injections is zero
    dQ = np.zeros([Ninj-Rinj], np.float64) # Heat Change Array (but not for injections without data!)
    for i in range(1,Ninj+1): # Loop over injections
       X = XMa[0,i]  # New Injectant Concentration
       M = XMa[1,i]  # New Cell Molecule Concentration
       Q[i] = (M*dH*V0/2)*( 1 + (X/(M)) + 1/(K*M) - np.sqrt( ((1 + X/(M) + 1/(K*M))**2) - 4*X/(M) ) ) # Total Heat
       if i > Rinj: # If initial injections were omitted from the ITC data, omit them here too
         dQ[i-1-Rinj] = (Q[i] + (dV/V0)*( (Q[i] + Q[i-1])/2 ) - Q[i-1])/(dV*X0) # Change in heat normalized by amount of injectant
    return dQ

### Define parameters for Monte Carlo cycles
Ncyc=CYC #Number of times Wiseman plots are resampled
varITC=np.zeros([len(ITC)], np.float64) #Empty array for a resampled Wiseman plot
parITC=np.zeros([len(ITC)], np.float64) #Empty array for a resampled Wiseman plot
varDQ=np.zeros([len(ITC),Ncyc], np.float64) #Empty array for heat
vardH=np.zeros([Ncyc], np.float64) #Empty array for fitted dH
varK=np.zeros([Ncyc], np.float64) #Empty array for fitted K
varN=np.zeros([Ncyc], np.float64) #Empty array for fitted N
varSS=np.zeros([Ncyc], np.float64) #Empty array for sum of squares
covs=np.zeros([3], np.float64) #Empty array for diagonals of covariance matrices
covmat=np.zeros([Ncyc], np.float64) #Empty array for LS error from each fit (irrelevant in default mode)

### Begin resampling
for n in range(Ncyc):
  ### Initial Guesses, Drawing Parameters for Either Floating or Fixed N
  if FIX == 0:
    p0=np.zeros([3], np.float64) # Guess Array
    p0[0]=-6000.0  # Guess dH
    p0[1]=1000  # Guess K
    p0[2]=1  # Guess N
  if FIX == 1:
    p0=np.zeros([2], np.float64) # Guess Array
    p0[0]=-6000.0  # Guess dH
    p0[1]=1000  # Guess K

  ###	Take concentration from	previously bootstrapped	distribution
  X0 = Xx[n]
  M0 = Mm[n]
  
  ### Find concentrations after each injection
  XM = np.zeros([2,Ninj+1], np.float64)
  for i in range(1,Ninj+1): # Loop over injections
    XM[0,i] = (i*dV*X0/V0)*( (1 - i*dV/(2*V0) ) ) # New Injectant Concentration
    XM[1,i] = M0*( (1 - i*dV/(2*V0) ) / (1 + i*dV/(2*V0) ) ) # New Cell Molecule Concentration
  
  ### Generate new Wiseman plot based on heat uncertainty, re-scale heats based on each re-sampled syringe concentration
  for i in range(len(ITC)):
    parITC[i]=np.random.normal(ITC[i],abs(np.sqrt((ITC[i]*Heat_Err)**2)+((Base_Err)**2))) # Add heat error
    varITC[i]=(parITC[i]/(dV*X0)) # Scale nominal Wiseman plot by new bootstrapped syringe concentration
  ### Fit the data allowing N to float
  if FIX == 0:
    ### Fit the data
    popt, pcov = curve_fit(func, XM, varITC, p0, maxfev = 1000000000 )
    dH = popt[0]
    K = popt[1]
    N = popt[2]
  
    ### (Print) Original Data, Fit, and find SumSqr
    fitdQ = func(XM, dH, K, N)
    SumSqr = 0.0
    for i in range(len(ITC)):
      SumSqr += (varITC[i] - fitdQ[i])**2
      varDQ[i,n]=fitdQ[i]
    covs = np.sqrt(np.diag(pcov))
    vardH[n]=dH
    varK[n]=K
    varN[n]=N
    varSS[n]=SumSqr
    covmat[n] = covs[0]
    print dH, covs[0]
  ### Fit the data while fixing N      
  if FIX == 1:
    ### Fit the data
    popt, pcov = curve_fit(func, XM, varITC, p0, maxfev = 1000000000 )
    dH = popt[0]
    K = popt[1]
  
    ### (Print) Original Data, Fit, and find SumSqr
    fitdQ = func(XM, dH, K)
    SumSqr = 0.0
    for i in range(len(ITC)):
      SumSqr += (varITC[i] - fitdQ[i])**2
      varDQ[i,n]=fitdQ[i]
    covs = np.sqrt(np.diag(pcov))
    vardH[n]=dH
    varK[n]=K
    varSS[n]=SumSqr
    covmat[n] = covs[0]

#if FIX == 0:
  #print "# dH= ",np.mean(vardH),"+/-",np.std(vardH)/65,np.mean(covmat)/65,"K= ",np.mean(varK),"+/-",np.std(varK),"N= ",np.mean(varN),"+/-",np.std(varN)*100

#if FIX == 1:
  #print "# dH= ",np.mean(vardH),"+/-",np.std(vardH)/65,np.mean(covmat)/65,"K= ",np.mean(varK),"+/-",np.std(varK)
