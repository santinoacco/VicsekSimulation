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
tEqV=[5000,5000,5000,5000,5000]
for j in range(10,15):
    namearch  = '..\VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'.dat'
    acorrFile='..\VAvsV_variandoN\Acorr_VA'+str(j)+'_tEq'+str(tEqA[j-10])+'.dat'
    PhiFile='..\VAvsV_variandoN\CrudosVicsekA'+str(j)+'.dat'
    PhiSqrFile='..\VAvsV_variandoN\CrudosVicsekA'+str(j)+'_Sqr.dat'
    AngInFiles.append([acorrFile,PhiFile,''])
    nroFilesA+=1
    InFile.append(namearch)
    
for j in range(10,15):
    namearch  = '..\VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'b.dat'
    acorrFile='..\VAvsV_variandoN\Acorr_VA2'+str(j)+'_tEq'+str(tEqA[j-10])+'.dat'
    PhiFile='..\VAvsV_variandoN\CrudosVicsek2A'+str(j)+'.dat'
    AngInFiles.append([acorrFile,PhiFile,''])

G=0
decG = math.modf(G)[0]
for j in range(10,15,2):
    namearch  = '..\VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'.dat'
    acorrFile='..\VAvsV_variandoN\Acorr_VV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'_tEq'+str(tEqV[j-10])+'.dat'
#    print(acorrFile)
    PhiFile='..\VAvsV_variandoN\CrudosVicsekV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    PhiSqrFile='..\VAvsV_variandoN\CrudosVicsekV'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    VecInFiles.append([acorrFile,PhiFile,PhiSqrFile])
    nroFilesV+=1

for j in range(10,15,2):
    namearch  = '..\VAvsV_variandoN\Entrada_Parametros_Vicsek'+str(j)+'b.dat'
    acorrFile='..\VAvsV_variandoN\Acorr2_VV'+str(j)+'_tEq5000.dat'
    PhiFile='..\VAvsV_variandoN\CrudosVicsek2V'+str(j)+'_G'+str(int(G))+'p'+str(int(decG*10))+'.dat'
    PhiSqrFile=''
    VecInFiles.append([acorrFile,PhiFile,PhiSqrFile])
    
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
for j in range(nroFilesA):
    
    VA_OBJ.append(MSP.AnalysisVicsek(j,AngInFiles[j],'A',parameters[j]))
    
    tACTANG.append(VA_OBJ[j].acorrTime)
    ACTANG.append(VA_OBJ[j].acorr)
    
VAextra=[]
for j in range(nroFilesA):
    VAextra.append(MSP.AnalysisVicsek(j,AngInFiles[5+j],'A',parameters[j]))


VV_OBJ=[]
tACTVEC=[]
ACTVEC=[]
parametersV=[parameters[k] for k in range(0,len(parameters),2)]

for j in range(nroFilesV):
    VV_OBJ.append(MSP.AnalysisVicsek(j,
                                     VecInFiles[j],
                                     'V',
                                     parametersV[j]))
for j in range(nroFilesV):
    tACTVEC.append(VV_OBJ[j].acorrTime)
    ACTVEC.append(VV_OBJ[j].acorr)


VVextra=[]
for j in range(nroFilesV):
    VVextra.append(MSP.AnalysisVicsek(j,
                                     VecInFiles[3+j],
                                     'V',
                                     parametersV[j]))
    

eta=VA_OBJ[0].vicsekNoise


labels = [r'$\eta = %s$' % eta[j] for j in range(len(eta))]
formatACT=['k.-','r.-','m.-','c.-','b.-','g.-','y.-','gv-','kv-','rv-']

'''para ajustar'''    
def f(x,tau,b,y0): 
    return b*np.exp(-x/tau) + y0

# armar array con las semillas para cada archivo.
initGuessA=[[[MidoCada[k],ACTANG[k][j][0],0] for j in range(len(eta))]for k in range(nroFilesA)]#semilla para el ajuste  
#initGuessA=[[[200,1] for j in range(len(eta))]for k in range(nroFilesA)]
# obtener el tiempo de eq. para cada curva de ACT
tauEq_A=[VA_OBJ[k].get_tauEq(f,initGuessA[k]) for k in range(nroFilesA)]

initGuessV=[[[MidoCada[k],ACTVEC[k][j][0],0] for j in range(len(eta))]for k in range(nroFilesV)]#semilla para el ajuste  
# obtener el tiempo de eq. para cada curva de ACT
tauEq_V=[VV_OBJ[k].get_tauEq(f,initGuessV[k]) for k in range(nroFilesV)]

graficoACT_A=[]
for k in range(nroFilesA):
    tituloAcorr=r'Anuglar - $N = %s$, $\rho = %s$' %(N[k],rho[k])
    graficoACT_A.append(
            [[tACTANG[k],ACTANG[k][j],None,formatACT[j],labels[j]]
            for j in range(0,len(eta),3)]
            )
    plt.figure(k)
    MSP.Graficador(graficoACT_A[k],'$t [MCS]$',r'$C(t)_c$',tituloAcorr)
#    plt.plot(tACTANG[3],f(tACTANG[3],tauEq_A[3][9],1,0),'b-')

graficoACT_V=[]
for k in range(nroFilesV):
    tituloAcorr=r'Vectorial - $N = %s$, $\rho = %s$' %(N[2*k],rho[k])
    graficoACT_V.append(
            [[tACTVEC[k],ACTVEC[k][j],None,formatACT[j],labels[j]]
            for j in range(0,len(eta),3)]
            )
    plt.figure(5+k)
    MSP.Graficador(graficoACT_V[k],'$t [MCS]$',r'$C(t)_c$',tituloAcorr)

# =============================================================================
# def f(x,tau,b): 
#     return b*np.exp(-x/tau)
# =============================================================================
#%%
plt.plot(eta,tauEq_A[3],'o:')
plt.xlabel(r'$\eta$',size=20)
plt.ylabel(r'$\tau_{eq}$',size=20)
plt.title(r'N=%s'%N[3])
plt.show()
#%%
# =============================================================================
# [plt.plot(eta,tauEq_A[k],'o') for k in range(nroFilesA)]
# plt.show()
# =============================================================================

# =============================================================================
# kTsorted = sorted(kT)
# print(kTsorted)
# 
# Esorted = [x for _,x in sorted(zip(kT,Emean))]
# =============================================================================

# utilizar tauEq_A para calcular los valores medios
PhiMean_A=[VA_OBJ[k].get_OP_Mean(tauEq_A[k]) for k in range(nroFilesA)]
PhiErr_A=[VA_OBJ[k].ErrOP(tauEq_A[k]) for k in range(nroFilesA)]
PhiVar_A=[VA_OBJ[k].get_OP_Var(tauEq_A[k]) for k in range(nroFilesA)]
PhiCB_A=[VA_OBJ[k]._calc_X_BinderCumulant(tauEq_A[k]) for k in range(nroFilesA)]

PhiMean_A_extra=[VAextra[k].get_OP_Mean(tauEq_A[k]) for k in range(nroFilesA)]
PhiErr_A_extra=[VAextra[k].ErrOP(tauEq_A[k]) for k in range(nroFilesA)]
PhiVar_A_extra=[VAextra[k].get_OP_Var(tauEq_A[k]) for k in range(nroFilesA)]
PhiCB_A_extra=[VAextra[k]._calc_X_BinderCumulant(tauEq_A[k]) for k in range(nroFilesA)]

etaAextra=VAextra[0].vicsekNoise
etaA2=np.concatenate((eta,etaAextra),axis=0)
etaA2S=sorted(etaA2)

PhiMean_V=[VV_OBJ[k].get_OP_Mean(tauEq_V[k]) for k in range(nroFilesV)]
PhiErr_V=[VV_OBJ[k].ErrOP(tauEq_V[k]) for k in range(nroFilesV)]
PhiVar_V=[VV_OBJ[k].get_OP_Var(tauEq_V[k]) for k in range(nroFilesV)]
PhiCB_V=[VV_OBJ[k]._calc_X_BinderCumulant(tauEq_V[k]) for k in range(nroFilesV)]

PhiMean_V_extra=[VVextra[k].get_OP_Mean(tauEq_V[k]) for k in range(nroFilesV)]
PhiErr_V_extra=[VVextra[k].ErrOP(tauEq_V[k]) for k in range(nroFilesV)]
PhiVar_V_extra=[VVextra[k].get_OP_Var(tauEq_V[k]) for k in range(nroFilesV)]
PhiCB_V_extra=[VVextra[k]._calc_X_BinderCumulant(tauEq_V[k]) for k in range(nroFilesV)]

etaextra=VVextra[0].vicsekNoise
eta2=np.concatenate((eta,etaextra),axis=0)
eta2S=sorted(eta2)


#PhiMean_V[2]=VV_OBJ[2].get_OP_Mean(tauEq_V[2])
#PhiErr_V[2]=VV_OBJ[2].ErrOP(tauEq_V[2])
#PhiVar_V[2]=VV_OBJ[2].get_OP_Var(tauEq_V[2])

#unir dos arrays que pertenecen a mismo parámetros
for k in range(nroFilesV):
    PhiMean_V[k]=np.concatenate((PhiMean_V[k],PhiMean_V_extra[k]),axis=0)
    PhiErr_V[k]=np.concatenate((PhiErr_V[k],PhiErr_V_extra[k]),axis=0)
    PhiVar_V[k]=np.concatenate((PhiVar_V[k],PhiVar_V_extra[k]),axis=0)
    PhiCB_V[k]=np.concatenate((PhiCB_V[k],PhiCB_V_extra[k]),axis=0)
    
for k in range(nroFilesA):
    PhiMean_A[k]=np.concatenate((PhiMean_A[k],PhiMean_A_extra[k]),axis=0)
    PhiErr_A[k]=np.concatenate((PhiErr_A[k],PhiErr_A_extra[k]),axis=0)
    PhiVar_A[k]=np.concatenate((PhiVar_A[k],PhiVar_A_extra[k]),axis=0)
    PhiCB_A[k]=np.concatenate((PhiCB_A[k],PhiCB_A_extra[k]),axis=0)

#PhiMean_V[1]=np.concatenate((PhiMean_V[1],PhiMean_V_extra[1]),axis=0)

#ordenar los arrays segun el ruido para graficarlos
PhiMean_V=[[x for _,x in sorted(zip(eta2,PhiMean_V[k]))] for k in range(nroFilesV)]
PhiErr_V=[[x for _,x in sorted(zip(eta2,PhiErr_V[k]))] for k in range(nroFilesV)]
PhiVar_V=[[x for _,x in sorted(zip(eta2,PhiVar_V[k]))] for k in range(nroFilesV)]
PhiCB_V=[[x for _,x in sorted(zip(eta2,PhiCB_V[k]))] for k in range(nroFilesV)]

PhiMean_A=[[x for _,x in sorted(zip(etaA2,PhiMean_A[k]))] for k in range(nroFilesA)]
PhiErr_A=[[x for _,x in sorted(zip(etaA2,PhiErr_A[k]))] for k in range(nroFilesA)]
PhiVar_A=[[x for _,x in sorted(zip(etaA2,PhiVar_A[k]))] for k in range(nroFilesA)]
PhiCB_A=[[x for _,x in sorted(zip(etaA2,PhiCB_A[k]))] for k in range(nroFilesA)]

plt.plot(eta2S,PhiVar_V[0],'ro-')
plt.plot(eta2S,PhiVar_V[1],'bo-')
plt.plot(eta2S,PhiVar_V[2],'go-')
plt.show()

plt.plot(etaA2S,PhiVar_A[0],'ro-')
plt.plot(etaA2S,PhiVar_A[1],'bo-')
plt.plot(etaA2S,PhiVar_A[2],'go-')
plt.show()

#%%
StatsLabels=[r'$N=%s$'%N[k] for k in range(nroFilesA)]
#formatos=['yo--','go--','bo--','ro--','mo--']
formatos=['y.--','g.--','b.--','r.--','m.--']

StatsLabels_V=[r'$N=%s$'%N[2*k] for k in range(nroFilesV)]
formatosV=['y*--','g*--','r*--','m*--','b*--',]

graficoOP=[[etaA2S,PhiMean_A[k],PhiErr_A[k],formatos[k],StatsLabels[k]] for k in range(nroFilesA)]
graficoSuc=[[etaA2S,PhiVar_A[k],None,formatos[k],StatsLabels[k]]for k in range(nroFilesA)]
graficosCB=[[etaA2S[0:10],PhiCB_A[k][0:10],None,formatos[k],StatsLabels[k]]for k in range(nroFilesA)]

graficoOP_V=[[eta2S,PhiMean_V[k],PhiErr_V[k],formatosV[k],StatsLabels_V[k]] for k in range(nroFilesV)]
graficoSuc_V=[[eta2S,PhiVar_V[k],None,formatosV[k],StatsLabels_V[k]]for k in range(nroFilesV)]
graficosCB_V=[[eta2S[6:13],PhiCB_V[k][6:13],None,formatosV[k],StatsLabels_V[k]]for k in range(nroFilesV)]


plt.figure(9)
MSP.Graficador(graficoOP,'$\eta$',r'$<|\phi|>$',r'Limite termodinámico, Angular, $\rho = %s$'%rho[0])
plt.figure(10)
MSP.Graficador(graficoSuc,'$\eta$',r'<$\sigma$>',r'Limite termodinámico, Angular, $\rho = %s$'%rho[0])
plt.figure(11)
MSP.Graficador(graficosCB,'$\eta$',r'$U_4$',r'Limite termodinámico, Angular, $\rho = %s$'%rho[0])

plt.figure(9)
MSP.Graficador(graficoOP_V,'$\eta$','<|$\phi$|>',r'Limite termodinámico, Vectorial, $\rho = %s$'%rho[0])
plt.figure(12)
MSP.Graficador(graficoSuc_V,'$\eta$',r'<$\sigma$>',r'Limite termodinámico, Vectorial, $\rho = %s$'%rho[0])
plt.figure(13)
MSP.Graficador(graficosCB_V,'$\eta$',r'$U_4$',r'Limite termodinámico, Vectorial, $\rho = %s$'%rho[0])

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
