import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import itertools


namFiles = False
classFiles = True


# Nam's Files
if namFiles:
    ClsCDM = np.loadtxt('data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.txt')
    ClsFDM = np.loadtxt('data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.txt')

# CLASS Files
if classFiles:
    ClsCDM = np.loadtxt('/Users/sayanmandal/Documents/ProjectN/CLASS/FromCosmoGPU/hub67CDM2cl.dat')
    ClsFDM = np.loadtxt('/Users/sayanmandal/Documents/ProjectN/CLASS/FromCosmoGPU/hubWDM_AmandaParam1cl.dat')


# Loading the Bin Edges
ellBinEdges = np.load('data/190404_experiment_0.25arc_0.5uk_2000_reionKSZ_2.912603855229585sqdeg_lbin_edges_dl300.npy')
ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2


#Counter To Run Through the ellMids
c = 0


# These Are To Store The Shortened Version Of The Spectra
ClsCDMTrim = []
ClsFDMTrim = []


# Trimming the Spectra
for i in range(len(ClsFDM)):
    if ClsCDM[i,0] == ellMids[c]:
        c += 1

        # Temporary Arrays To Store Each Row
        tempCDM = []
        tempFDM = []

        # ell
        tempCDM.append(ClsCDM[i,0])
        tempFDM.append(ClsFDM[i,0])

        if namFiles:
            tempCDM.append(ClsCDM[i,1])
            tempFDM.append(ClsFDM[i,1])

        if classFiles:
            t1 = 0
            t1 = ClsCDM[i,0] * (ClsCDM[i,0] + 1) * ClsCDM[i,5]
            tempCDM.append(t1)

            t2 = 0
            t2 = ClsFDM[i,0] * (ClsFDM[i,0] + 1) * ClsFDM[i,5]
            tempFDM.append(t2)

        # Appending Rows
        ClsCDMTrim.append(tempCDM)
        ClsFDMTrim.append(tempFDM)
        
    if c == len(ellMids):
        break


# Convert List to Array
ClsCDMTrimArr = np.array(ClsCDMTrim)
ClsFDMTrimArr = np.array(ClsFDMTrim)


# Values
beam = 0.25
fskys = [0.5,0.25,0.1]
noise = 0.5


# Printing The Header
print 'beam,fsky,noise,S/N'


for fsky in fskys:
    # Loading The CovMat
    cov = np.load('data/190404_experiment_0.25arc_0.5uk_2000_reionKSZ_2.912603855229585sqdeg_covmat_dl300.npy')

    # Scaling and Inverting
    cov = cov*2.912603855229585/41253./fsky # scale with the right fsky
    covinv = inv(cov)
    
    
    # Calculate SNR squared
    sn2 = 0
    listmax = min(len(ClsCDMTrim),len(ellMids))

    for i in range (listmax):
        for j in range (listmax):
            sn2 += (ClsCDMTrimArr[i,1]-ClsFDMTrimArr[i,1]) * covinv[i,j] * (ClsCDMTrimArr[j,1]-ClsFDMTrimArr[j,1])
    
    
    print beam,fsky,noise,np.round(np.sqrt(np.sum(sn2)),0)
