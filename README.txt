README


There are two primary files included here: PYTC.py, the fitting software, and genit.py, a Python script that generates artificial raw heats based
on user inputted thermodynamic parameters.


#################
genit.py
#################

To use genit.py, you must supply seven (7) arguments. 

[1] Syringe concentration, in mM
[2] Cell concentration, in mM
[3] Injection volume, in microliters
[4] Number of injections
[5] Binding Enthalpy, in kcal/mol
[6] Equilibrium constant (M^-1)
[7] 0 if you want raw heats only, 1 if you want normalized heats with X/M (Wiseman plot)

#################
PYTC.py
#################

To use PYTC.py, you must supply ten (10) arguments. 

[1] Injection volume used in experiment, in microliters
[2] Number of Monte Carlo cycles to use
[3] Syringe concentration used in experiment, in mM
[4] Cell concentration used in experiment, in mM
[5] Percent uncertainty in syringe concentration
[6] Percent uncertainty in cell concentration
[7] Percent uncertainty in heat measurment (peak-to-peak)
[8] Baseline uncertainty (sigma_irreducible)
[9] Check if N is fixed at one (TRUE = 1, FALSE = 0)
[10] File containing raw heats, in units of calories


An example file (Example.dat) with nominal dH = -9 kcal/mol, K = 36000, and N = 1 is included. Below is the process to run the file:

python PYTC.py 10 2000 15 1.5 0.6 0.6 1 0.13 0 Example.dat

# dH=  -9002.41702359 +/- 63.2095073764 K=  36017.0903671 +/- 877.313523474 N=  0.999560306112 +/- 0.00848798071207

This yields an output with a 0.6% error in both syringe and cell concentrations, 1% error in heat error, and 0.13 irreducible heat error, with N
allowed to flow.

Output is dH in cal/mol, +/- the resampling S.D.(cal/mol), K (M^-1 units) and resampling S.D., and N (unitless) and resampling S.D..
