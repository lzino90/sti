# sti

Code used for the simulations in "A time-varying network model for sexually transmitted infections accounting for behavior and control actions" by K. Frieswijk, L. Zino, and M. Cao. International Journal of Robust and Nonlinear Control, 2021, doi: 10.1002/rnc.5930

INSTRUCTIONS AND SYSTEM REQUIREMENTS:

The code for the simulation is written in MATLAB, using version 2021a and requires no add-ons. All the functions needed are included in the corresponding folders. Hence, to run the code it is sufficient to open the corresponding folder with MATLAB.

The code has been tested on a PC with 16GB RAM and CPU 1.9 GHz with OS Windows 10 (64-bit). Maximum run time is approximately 10 seconds, for a population of n=10,000 individuals.

ORGANIZATION:

The project contains two main files, and some auxiliary functions. The main files are STI1s and STI2s, which can be used to simulate the epidemics with a 1-stage and a 2-stage negotiation process, respectively. The details on the parameters are reported in the paper and in the MATLAB code. The function parnters_net(n,p) can be used to generate a partners network with n individuals and a fraction p of them involved in dyadic relations.
