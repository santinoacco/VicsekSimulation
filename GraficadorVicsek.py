# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:39:08 2018

@author: Santiago
"""
import matplotlib.pyplot as plt
#import csv
import numpy as np
#import sys to conect path to dir
import sys
sys.path.insert(0, 'C:/Users/Santiago/Desktop/Facultad/Fisica/Simulaciones')
import ModSimPy as MSP



nroFiles=5


etaA=[]
tA=[[]]*nroFiles
tACT_A=[[]]*nroFiles
tau_Eq_A=[[]]*nroFiles

ACT_A=[[]]*nroFiles

PhiA=[[]]*nroFiles
PhiSqrA=[[]]*nroFiles
Phi_Mean_A=[[]]*nroFiles
PhiSqr_Mean_A=[[]]*nroFiles

varPhiA=[[]]*nroFiles
CB_A=[[]]*nroFiles


ErrPhiMean=[[]]*nroFiles
#VarPhi_Mean_A=[[],[],[]]
#VarPhiSqr_Mean_A=[[],[],[]]
#ErrVarPhiMean=[[],[],[]] # no la puedo estimar directamente.. necesito usar blocking method ponele



#variables Vicsek Vectorial
'CrudosVicsekV08_G0p0.dat'

"""archivos analisis VA_variandoN"""
archIn_crudosAng=['VAvsV_variandoN\CrudosVicsekA10.dat',
                  'VAvsV_variandoN\CrudosVicsekA11.dat',
                  'VAvsV_variandoN\CrudosVicsekA12.dat',
                  'VAvsV_variandoN\CrudosVicsekA13.dat',
                  'VAvsV_variandoN\CrudosVicsekA14.dat',
                  ]
archIn_crudosvarAng=['VAvsV_variandoN\CrudosVicsekA10_Sqr.dat',
                     'VAvsV_variandoN\CrudosVicsekA11_Sqr.dat',
                     'VAvsV_variandoN\CrudosVicsekA12_Sqr.dat',
                     'VAvsV_variandoN\CrudosVicsekA13_Sqr.dat',
                     'VAvsV_variandoN\CrudosVicsekA14_Sqr.dat'
                     ]

archIn_ACTcAng=['VAvsV_variandoN\Acorr_VA10.dat',
                'VAvsV_variandoN\Acorr_VA11.dat',
                'VAvsV_variandoN\Acorr_VA12.dat',
                'VAvsV_variandoN\Acorr_VA13.dat',
                'VAvsV_variandoN\Acorr_VA14.dat'
                     ]
printing=[True,False,False,False,False]

for j in range(nroFiles):
    etaA,tA[j],PhiA[j]=MSP.LecturaTablasSim(archIn_crudosAng[j],printing[j])
    etaA=[]
    tA[j]=[]
    etaA,tA[j],PhiSqrA[j] = MSP.LecturaTablasSim(archIn_crudosvarAng[j],False)
    etaA=[]
    etaA,tACT_A[j],ACT_A[j] = MSP.LecturaTablasSim(archIn_ACTcAng[j],False)



#InFile='Entrada_Parametros_Vicsek.dat'
#InFile=['Entrada_Parametros_Vicsek.dat','Entrada_Parametros_Vicsek.dat']
InFile=['VAvsV_variandoN\Entrada_Parametros_Vicsek10.dat',
        'VAvsV_variandoN\Entrada_Parametros_Vicsek11.dat',
        'VAvsV_variandoN\Entrada_Parametros_Vicsek12.dat',
        'VAvsV_variandoN\Entrada_Parametros_Vicsek13.dat',
        'VAvsV_variandoN\Entrada_Parametros_Vicsek14.dat'
        ]
#def ReadParameters(InFile,printing):

#    return N,rho,v0,MidoCada
N=[]
v0=[]
rho=[]
MidoCada=[]


for j in range(nroFiles):
#    ReadParameters(InFile[j],N[j],v0[j],rho[j],MidoCada[j],True)
    MSP.ReadParameters(InFile[j],N,v0,rho,MidoCada,True)



""" 
    A PARTIR DEL GRAF. phi vs t estimo el tiempo de equilibrio,
    luego tomo el promedio a partir de dicho dato.
"""
#Introducir tiempo de equilibración
#tEqA = 5000
"""1er Paso: estudiar ACTc y calcular tau_Eq"""
#tau_Eq_A=[]

#tEqA=[10000,12000,5000]

'''para ajustar'''    
def f(x,tau,b,y0): 
    return b*np.exp(-x/tau) + y0
def fnorm(x,tau): 
    return 1.0*np.exp(-x/tau)


tau_Eq_A=[[]]*nroFiles
#print(tA[0])




#print(tACT_A[j])
initGuess=[[]]*nroFiles
for j in range(nroFiles):
    #NOTA Importante: por como guardo los tiempos en ACT debo multiplicar el
    # vector tiempo por MidoCada para que pueda hacer la comparación tau vs MidoCada
    for l in range(len(tACT_A[j])):
        elem = tACT_A[j][l]*MidoCada[j]
        tACT_A[j][l]=elem

    
#    #También tengo que normalizar a 1 ACTc_A
#    for noise in range(len(etaA)):
#        elem0=ACT_A[j][noise][0]
#        for l in range(len(ACT_A[j][noise])):
#            elem = ACT_A[j][noise][l]/elem0
#            ACT_A[j][noise][l]=elem
#        initGuess[j]=[MidoCada,ACT_A[j][noise][0],0] #semilla para el ajuste
#    tau_Eq_A[j]=ExtractTauEq(etaA,f,tACT_A[j],ACT_A[j][noise],initGuess[j],MidoCada[j])
    

    initGuess[j]=[MidoCada[j],ACT_A[j][0],0]#semilla para el ajuste
    tau_Eq_A[j]=MSP.ExtractTauEq(etaA,f,tACT_A[j],ACT_A[j],initGuess[j],MidoCada[j])
    
print(tau_Eq_A)



plt.plot(etaA,tau_Eq_A[0], 'o')
plt.plot(etaA,tau_Eq_A[1], 'ok')
plt.xlabel(r'$\eta$',fontsize = 20)
plt.ylabel(r'$\tau_{acorr}$',fontsize = 20)
plt.legend()
plt.show()

dim_fin=len(etaA)
#V0=[0.5,1,1.5]
#MidoCada=[200,200,1000]



#graficoACT=[[]]*nroFiles

graficoACT=[]
formatos=['k.-','b.-','g.-','c.-','r.-','k*-','b*-','g*-','c*-','r*-']
for j in range(nroFiles):
    graficoACTaux=[]
    for temp in range(0,len(etaA),1):
        labelACT= r'$\Delta t_{mido}$=%s, $\eta$ %s' % (MidoCada[j],etaA[temp])
        graficoACTaux.append([tACT_A[j],ACT_A[j][temp::dim_fin],None,formatos[temp],labelACT])
    graficoACT.append(graficoACTaux)

for j in range(nroFiles):
    tituloACT=r'Angular v0=%s, variación v0, N = %s, $\rho$ = %s' % (v0[j],N[j],rho[j])
    MSP.Graficador(graficoACT[j][1::2],'t [MCS]','Acorr(|$\phi$|)',tituloACT)
""""##########################################""" 

"""" VALORES MEDIOS """ 
def SquareList(lista):
    ListaSqr=[]
    for item in lista:
        ListaSqr.append(item**2)
    return ListaSqr

#grafico=[[graf1][graf2][graf3]]
#graf1=[x,y,yEr,'g-',label]

""" ErrPhi = sqrt((<Phi**2>-<Phi>**2)/n) = sqrt(Sucep/n)"""
""" Sucep = varPhi = <Phi**2>-<Phi>**2 """
""" CumulanteBinder = CB_A = 1- <Phi^4>/3*<Phi^2>^2"""

#Calcular valores Medios
#pregunta?? corresponde que tome el valor medio sujeto a cada tau de eta distinto?
#rta: en principio sí, pero en aquellos en los que no necesitas calcular
n=[[]]*nroFiles
Phi4th=[[]]*nroFiles
for j in range(nroFiles):
    PhiSqrA[j]=SquareList(PhiA[j])
    Phi4th[j]=SquareList(PhiSqrA[j])
# =============================================================================
#     CB_A[j]= kurtosis(PhiA[j],fisher=True)
#     print(CB_A[j])
# =============================================================================
    
Phi_Mean_ASqr=[[]]*nroFiles
PhiSqrSqr_A=[[]]*nroFiles

varPhiA=[]
n=[]
ErrPhiMean=[]
CB_A=[]

#print(varPhiA)
for j in range(nroFiles):
    
    
    Phi_Mean_A[j] = MSP.ValorMedio(PhiA[j],etaA,tau_Eq_A[j],MidoCada[j]) #-> <phi>
    PhiSqr_Mean_A[j] = MSP.ValorMedio(PhiSqrA[j],etaA,tau_Eq_A[j],MidoCada[j]) #-> <phi**2>
    
        
    Phi_Mean_ASqr[j]=SquareList(Phi_Mean_A[j])  #-> <phi>**2
    
    PhiSqrSqr_A[j]=SquareList(PhiSqr_Mean_A[j]) #-> <phi**2>**2
    
#    print(Phi_Mean_A[j])
    
    

    varPhiAaux=[]
    naux=[]
    ErrPhiMeanaux=[]
    CB_Aaux=[]
    for l in range(len(etaA)):
        elemvar=PhiSqr_Mean_A[j][l] - Phi_Mean_ASqr[j][l]
        elemn=tau_Eq_A[j][l]/MidoCada[j]
        varPhiAaux.append(elemvar)
        naux.append(elemn)
    
        CB_Aaux.append(1 - Phi4th[j][l]/(3.0*PhiSqrSqr_A[j][l]))
        
        #falta distingir entre las distintas formas de computar el error
        if tau_Eq_A[j][l]<=MidoCada[j]:
            ErrPhiMeanaux.append(np.sqrt(elemvar/(elemn)))
        else:
            tmax=max(tACT_A[j])-tau_Eq_A[j][l]
            ErrPhiMeanaux.append(np.sqrt(2*tau_Eq_A[j][l]*elemvar/tmax))
        
    varPhiA.append(varPhiAaux)
    n.append(naux)
    ErrPhiMean.append(ErrPhiMeanaux)
    CB_A.append(CB_Aaux)
# =============================================================================
#     for l in range(len(etaA)):
#         varPhiA[j].append(PhiSqr_Mean_A[j][l] - Phi_Mean_ASqr[j][l]) #notar que no es valor medio, sino que estimo 1 p/c ruido
#         n[j].append(tau_Eq_A[j][l]/MidoCada[j])
#         
# #       calculo el ERROR cuando medidas son independientes como:
# #       if tauEq > tmido -> corrección Barkema
# #        ErrPhiMean[j].append(np.sqrt(varPhiA[j][l]/(n[j]-1)))
# #        ErrPhiMean[j].append(np.sqrt(varPhiA[j][l]/(n[j][l]-1)))
#         ErrPhiMean[j].append(np.sqrt(varPhiA[j][l]/(n[j][l])))
# 
# 
# 
#         CB_A[j].append(1 - Phi4th[j][l]/(3.0*PhiSqrSqr_A[j][l]))
# =============================================================================
#    print(varPhiA[:])
#    print(ErrPhiMean[j])






# =============================================================================
# """ esto sería tomando promedios con cada tau particular"""
# for j in range(3):
#     for eta in range(len(etaA)):
#         Phi_Mean_A[j] = MSP.ValorMedio(PhiA[j],etaA,tauEq[j][eta],MidoCada[j])
#         VarPhi_Mean_A[j] = MSP.ValorMedio(varPhiA[j],etaA,tauEq[j][eta],MidoCada[j])
# =============================================================================


tituloG1 = r'Comparación Angular distintos $\tau$, $\rho = %s$' % (rho)
#tituloSuc = 'Comparación Angular distintos V0, N = %s, L = %s' % (N,L)

formatos=['g.--','b.--','r.--','m.--','k.--']

grafico=[]
grafSuc=[]
grafCumBind=[]

""" NOTA:
    a la hora de calcular los errores no puedo usar simplemente np.sqrt(VarPhi_Mean_A[j])
    hay que hacerlo como dice Barkema, teniendo en cuenta medidas independientes.
    para ello primero hay que calcular el tEq de cada Ruido.
    luego usar la expresión (3.39) y (3.40) del Barkema
    
    ALARGAR la SERIE.. la corrección solo para el caso azul, para el resto que veo que es menor a 500
    puedo decir que son independientes
    """


for j in range(nroFiles):
    labelAng= r'$N=%s$' % N[j] #, $\tau_{eq}$ = %s'% (N[j],tau_Eq_A[j])
#    labelAng=''
    
    grafico.append([etaA,Phi_Mean_A[j],None,
                    formatos[j],labelAng
                    ])
# =============================================================================
#     grafico.append([etaA,Phi_Mean_A[j],ErrPhiMean[j],
#                     formatos[j],labelAng
#                     ])
# =============================================================================
    grafSuc.append([etaA,varPhiA[j],None,
                    formatos[j],labelAng
                    ])
    grafCumBind.append([etaA,CB_A[j],None,
                    formatos[j],labelAng
                    ])


MSP.Graficador(grafico,'$\eta$','$<|\phi|>$',tituloG1)
MSP.Graficador(grafSuc,'$\eta$','Suceptibilidad',tituloG1)
MSP.Graficador(grafCumBind,'$\eta$','Cumulante de Binder',tituloG1)




""" extras"""

