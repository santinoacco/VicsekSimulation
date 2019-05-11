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
#nroFilesA=0
tEqA=[5000,12000,30000,5000,10000]
tEqV=[5000,5000,5000,5000,5000,]
for j in range(10,15):
    namearch  = 'VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'.dat'
    acorrFile='VAvsV_variandoN\Acorr_VA'+str(j)+'_tEq'+str(tEqA[j-10])+'.dat'
    PhiFile='VAvsV_variandoN\CrudosVicsekA'+str(j)+'.dat'
    PhiSqrFile='VAvsV_variandoN\CrudosVicsekA'+str(j)+'_Sqr.dat'
    AngInFiles.append([acorrFile,PhiFile,PhiSqrFile])
    nroFilesA+=1
    InFile.append(namearch)
 
G=0
decG = math.modf(G)[0]
for j in range(10,15,2):
    namearch  = 'VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'.dat'
    acorrFile='VAvsV_variandoN\Acorr_VV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'_tEq'+str(tEqV[j-10])+'.dat'
#    print(acorrFile)
    PhiFile='VAvsV_variandoN\CrudosVicsekV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    PhiSqrFile='VAvsV_variandoN\CrudosVicsekvV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    VecInFiles.append([acorrFile,PhiFile,PhiSqrFile])
    nroFilesV+=1

    

N=[]
v0=[]
rho=[]
MidoCada=[]




parameters=[]
for j in range(nroFilesA):
    MSP.ReadParameters(InFile[j],N,v0,rho,MidoCada,True)    
    parameters.append([N[j],v0[j],rho[j],MidoCada[j]])

VA_OBJ=[]
tACTANG=[]
ACTANG=[]
ACTVEC=[]
for j in range(nroFilesA):
    
    VA_OBJ.append(MSP.AnalysisVicsek(j,AngInFiles[j],'A',parameters[j]))
    
    tACTANG.append(VA_OBJ[j].acorrTime)
    ACTANG.append(VA_OBJ[j].acorr)
    
#    AngInFiles[j]


eta=VA_OBJ[0].vicsekNoise


labels = [r'$\eta = %s$' % eta[j] for j in range(len(eta))]
formatACT=['k.-','r.-','m.-','c.-','b.-','g.-','y.-','gv-','kv-','rv-']

graficoACT_A=[]
for k in range(nroFilesA):
    tituloAcorr=r'Anuglar - $N = %s$, $\rho = %s$' %(N[k],rho[k])
    graficoACT_A.append(
            [[tACTANG[k],ACTANG[k][j],None,formatACT[j],labels[j]]
            for j in range(0,len(eta),3)]
            )
    plt.figure(k)
    MSP.Graficador(graficoACT_A[k],'$t [MCS]$','$Acorr$',tituloAcorr)
    
'''para ajustar'''    
def f(x,tau,b,y0): 
    return b*np.exp(-x/tau) + y0
# =============================================================================
# def f(x,tau,b): 
#     return b*np.exp(-x/tau)
# =============================================================================

# armar array con las semillas para cada archivo.
initGuessA=[[[MidoCada[k],ACTANG[k][j][0],0] for j in range(len(eta))]for k in range(nroFilesA)]#semilla para el ajuste  
#initGuessA=[[[200,1] for j in range(len(eta))]for k in range(nroFilesA)]
# obtener el tiempo de eq. para cada curva de ACT
tauEq_A=[VA_OBJ[k].get_tauEq(f,initGuessA[k]) for k in range(nroFilesA)]

# =============================================================================
# initGuessV=[[[MidoCada[k],ACTVEC[k][j][0],0] for j in range(len(eta))]for k in range(nroFilesV)]#semilla para el ajuste  
# # obtener el tiempo de eq. para cada curva de ACT
# tauEq_V=[VA_OBJ[k].get_tauEq(f,initGuessV[k]) for k in range(nroFilesV)]
# =============================================================================

# utilizar tauEq_A para calcular los valores medios

PhiMean_A=[VA_OBJ[k].get_OP_Mean(tauEq_A[k]) for k in range(nroFilesA)]
PhiErr_A=[VA_OBJ[k].ErrOP(tauEq_A[k]) for k in range(nroFilesA)]
PhiVar_A=[VA_OBJ[k].get_OP_Var(tauEq_A[k]) for k in range(nroFilesA)]
PhiCB_A=[VA_OBJ[k]._calc_X_BinderCumulant(tauEq_A[k]) for k in range(nroFilesA)]

StatsLabels=[r'$N=%s$'%N[k] for k in range(nroFilesA)]
#formatos=['yo--','go--','bo--','ro--','mo--']
formatos=['y.--','g.--','b.--','r.--','m.--']

graficoOP=[[eta,PhiMean_A[k],PhiErr_A[k],formatos[k],StatsLabels[k]] for k in range(nroFilesA)]
graficoSuc=[[eta,PhiVar_A[k],None,formatos[k],StatsLabels[k]]for k in range(nroFilesA)]
graficosCB=[[eta[:],PhiCB_A[k][:],None,formatos[k],StatsLabels[k]]for k in range(nroFilesA)]

MSP.Graficador(graficoOP,'$\eta$','<|$\phi$|>',None)
MSP.Graficador(graficoSuc,'$\eta$','Suceptibilidad',None)
MSP.Graficador(graficosCB,'$\eta$','CumulanteBinder',None)


""" extras"""
# =============================================================================
# #probamos calcualr ACT de otra manera, usamos archivos de ACT_Barkema
# 
# #Phi=[]
# Phi=VA_OBJ[4].get_OrderParam
# time = VA_OBJ[4].get_OPtime
# 
# plt.plot(time,Phi[9],'')
# 
# #Phi[9]/=Phi[9][0]
# Phi[9]-=np.mean(Phi[9])
# 
# #print(Phi[9])
# plt.plot(time,Phi[9],'r')
# 
# plt.show()
# 
# #calculo ACTc con np.
# def autocorr(x):
#     result = np.correlate(x, x, mode='full')
#     return result[int(result.size/2):]
# ACTc=autocorr(Phi[9])
# 
# ent0=int(1000/MidoCada[0])
# ACTc2=autocorr(Phi[9][ent0::])
# ent1=int(5000/MidoCada[0])
# ACTc3=autocorr(Phi[9][ent1::])
# ent2=int(10000/MidoCada[0])
# ACTc4=autocorr(Phi[9][ent2::])
# ent3=int(20000/MidoCada[0])
# ACTc5=autocorr(Phi[9][ent3::])
# #normalizo ACTc
# ACTc/=ACTc[0]
# ACTc2/=ACTc2[0]
# ACTc3/=ACTc3[0]
# ACTc4/=ACTc4[0]
# ACTc5/=ACTc5[0]
# 
# ACTANG[4][9]/=ACTANG[4][9][0]
# 
# 
# tfin=int(len(time)/2)
# plt.plot(time[:tfin],ACTc[:tfin],'k-',
#          label = r'ACTc calcualda en Scipy para $\Delta t$, graf hasta $\Delta t / 2$')
# plt.plot(time[:tfin],ACTc2[:tfin],'g-',
#          label = r'ACTc calcualda en Scipy para $\Delta t$, graf hasta $\Delta t / 2$')
# plt.plot(time[:tfin],ACTc3[:tfin],'-',
#          label = r'ACTc calcualda en Scipy para $\Delta t$, graf hasta $\Delta t / 2$')
# plt.plot(time[:tfin],ACTc4[:tfin],'-',
#          label = r'ACTc calcualda en Scipy para $\Delta t$, graf hasta $\Delta t / 2$')
# plt.plot(time[:tfin],ACTc5[:tfin],'-',
#          label = r'ACTc calcualda en Scipy para $\Delta t$, graf hasta $\Delta t / 2$')
# plt.plot(time[:tfin],ACTANG[4][9],'r--',
#          label = r'ACTc calcualda con rutina Barkema')
# plt.plot(time[:tfin],[0]*tfin,'b')
# plt.legend()
# plt.title(r'Comparación ACTc Angular para N:%s, $\eta$:%s'%(4,eta[9]))
# plt.show()
# #
# #plt.plot(time[:tfin],[0]*tfin,'b')
# #plt.show()
# plt.plot(time,ACTc,'k-',label = r'graf hasta $\Delta t$')
# plt.plot(time,0*time,'b')
# plt.title(r'ACTc Angular calculada Scipy para N:%s, $\eta$:%s'%(4,eta[9]))
# plt.show()
# 
# """ lo que veo es que cuando calculo ACTc para todo el periódo de medida,
#     la curva efectivamente oscila en torno al 0.
#     si calculo para la mitad del periodo entonces hay un desfasaje del 0.
# """
# =============================================================================
