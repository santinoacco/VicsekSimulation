# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:39:08 2018

@author: Santiago
"""
#import sys to conect path to dir
import sys
sys.path.insert(0, 'C:/Users/Santiago/Desktop/Facultad/Fisica/Simulaciones')
import ModSimPy as MSP
import matplotlib.pyplot as plt
import numpy as np
import math





nroFilesA=0
nroFilesV=0

InFile=[]
AngInFiles=[]
VecInFiles=[]
nroFilesA=0
for j in range(10,15):
    namearch  = 'VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'.dat'
    acorrFile='VAvsV_variandoN\Acorr_VA'+str(j)+'.dat'
    PhiFile='VAvsV_variandoN\CrudosVicsekA'+str(j)+'.dat'
    PhiSqrFile='VAvsV_variandoN\CrudosVicsekA'+str(j)+'_Sqr.dat'
    AngInFiles.append([acorrFile,PhiFile,PhiSqrFile])
    nroFilesA+=1
    InFile.append(namearch)
 
G=0
decG = math.modf(G)[0]
for j in range(10,15,2):
    namearch  = 'VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'.dat'
    acorrFile='VAvsV_variandoN\Acorr_VV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
#    print(acorrFile)
    PhiFile='VAvsV_variandoN\CrudosVicsekV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    PhiSqrFile='VAvsV_variandoN\CrudosVicsekvV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    VecInFiles.append([acorrFile,PhiFile,PhiSqrFile])
    nroFilesV+=1

    

N=[]
v0=[]
rho=[]
MidoCada=[]



VA_OBJ=[]
parameters=[]
for j in range(nroFilesA):
    MSP.ReadParameters(InFile[j],N,v0,rho,MidoCada,True)    
    parameters.append([N[j],v0[j],rho[j],MidoCada[j]])

tACTANG=[]
ACTANG=[]
for j in range(nroFilesA):
    
    VA_OBJ.append(MSP.AnalysisVicsek(j,AngInFiles[j],'A',parameters[j]))
    
    tACTANG.append(VA_OBJ[j].acorrTime)
    ACTANG.append(VA_OBJ[j].acorr)
#    AngInFiles[j]


eta=VA_OBJ[0].vicsekNoise


labels = [r'$\eta = %s$' % eta[j] for j in range(len(eta))]
graficoACT=[]
formatACT=['k.-','r.-','m.-','c.-','b.-','g.-','y.-','gv-','kv-','rv-']
for k in range(nroFilesA):
    tituloAcorr='$N = %s$'% N[k]
    graficoACT.append(
            [[tACTANG[k],ACTANG[k][j],None,formatACT[j],labels[j]]
            for j in range(len(eta))]
            )
    MSP.Graficador(graficoACT[k],'$t [MCS]$','$Acorr$',tituloAcorr)
    
# =============================================================================
# for k in range(nroFilesA): 
#     tituloAcorr='$N = %s$'% N[k]
#     [plt.errorbar(tACTANG[k],
#                   ACTANG[k][j],
#                   yerr=None,
#                   fmt='.-',
#                   label=labels[j]) for j in range(len(eta))]
#     plt.xlabel('$t [MCS]$',fontsize = 20)
#     plt.ylabel('$Acorr$',fontsize = 20)
#     plt.title(tituloAcorr)
#     plt.legend()
#     plt.show()
# =============================================================================

'''para ajustar'''    
def f(x,tau,b,y0): 
    return b*np.exp(-x/tau) + y0

# armar array con las semillas para cada archivo.
initGuessA=[[[MidoCada[k],ACTANG[k][j][0],0] for j in range(len(eta))]for k in range(nroFilesA)]#semilla para el ajuste  
# obtener el tiempo de eq. para cada curva de ACT
tauEq_A=[VA_OBJ[k].get_tauEq(f,initGuessA[k]) for k in range(nroFilesA)]

initGuessV=[[[MidoCada[k],ACTVEC[k][j][0],0] for j in range(len(eta))]for k in range(nroFilesV)]#semilla para el ajuste  
# obtener el tiempo de eq. para cada curva de ACT
tauEq_V=[VA_OBJ[k].get_tauEq(f,initGuessV[k]) for k in range(nroFilesV)]

# utilizar tauEq_A para calcular los valores medios

PhiMean_A=[VA_OBJ[k].get_OP_Mean(tauEq_A[k]) for k in range(nroFilesA)]
PhiVar_A=[VA_OBJ[k].get_OP_Var(tauEq_A[k]) for k in range(nroFilesA)]
PhiCB_A=[VA_OBJ[k]._calc_X_BinderCumulant(tauEq_A[k]) for k in range(nroFilesA)]

StatsLabels=[r'$N=%s$'%N[k] for k in range(nroFilesA)]
formatos=['yo--','go--','bo--','ro--','mo--']

graficoOP=[[eta,PhiMean_A[k],None,formatos[k],StatsLabels[k]] for k in range(nroFilesA)]
graficoSuc=[[eta,PhiVar_A[k],None,formatos[k],StatsLabels[k]]for k in range(nroFilesA)]
graficosCB=[[eta[:8],PhiCB_A[k][:8],None,formatos[k],StatsLabels[k]]for k in range(nroFilesA)]

MSP.Graficador(graficoOP,'$\eta$','<|$\phi$|>',None)
MSP.Graficador(graficoSuc,'$\eta$','Suceptibilidad',None)
MSP.Graficador(graficosCB,'$\eta$','CumulanteBinder',None)


""" extras"""

