# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:31:26 2019

@author: Santiago
"""

import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import kurtosis

# MODULO PARA SIMULACIONES
# Nomenclatura:
#   PO: Parámetro de Orden (Order Parameter)
#   t: tiempo (time)
#   X: variable que evoluciona con el tiempo y tiene como parámetro PO. (variable X[t,PO])

def ReadParameters(InFile,N,v0,rho,MidoCada,printing):
    with open(InFile) as f_in_names:
    
        plots=csv.reader(f_in_names,delimiter = ' ')
        Naux=int(f_in_names.readline().split()[1])
        N.append(Naux)
        L=int(f_in_names.readline().split()[1])
        v0.append(float(f_in_names.readline().split()[1]))
        for i in range(5):
            next(plots)
        MidoCada.append(int(f_in_names.readline().split()[1]))
        rho.append(Naux/L**2)
        if printing == True:
            print('N = %s' % N)
            print('rho = %s' % rho)
            print( r'$\Delta$ t = %s' % MidoCada)
    
    return 
# función para leer las tablas .dat


def Autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[int(result.size/2):]


#creamos una rutina que calcula el tiempo de equilibración
#distingue si es mayor a un dado límite o no
def ExtractTauEq(noise,f_fit,time,ACT,initGuessParam,lowerBound):
    TauEq=[]
    for j in range(len(noise)):
#        tau_temp=curve_fit(f_fit,time,ACT[j],p0=initGuessParam)[0][0] #el elem [0] contiene el ajuste de parametros
        tau_temp=curve_fit(f_fit,time,ACT[j],p0=initGuessParam[j],maxfev = 2000)[0][0]
        if tau_temp <= lowerBound:
            TauEq.append(lowerBound)
        else:
            TauEq.append(tau_temp)
        
            
    return TauEq


""" La función Graficador """
# -------------- ACLARACIONES --------------
# el parámetro 'grafico' es una lista cuyos elementos
# contienen la siguiente información, en el siguiente orden:
# grafico[0]=[x,y,yerr,fmt,label]
# x: variable indepentiente
# y: variable dependiente
# yerr: error de variable dependiente
# fmt: formato, es de tipo string
# label: etiquieta, tipo string
#-------------- EJEMPLO -------------
#grafico=[[graf1][graf2][graf3]]
#graf1=[x,y,yEr,'g-',label]
#----------------------------
def Graficador(grafico,xlabel,ylabel,title):
# =============================================================================
#     print("grafico=[[graf1][graf2][graf3]]")
#     print("graf1=[x,y,yEr,'g-',label]")
# =============================================================================
    assert type(grafico[0])==list, 'missing brakets in "grafico" '
    for item in grafico:
        plt.errorbar(
                item[0],
                item[1],
                yerr = item[2],
                fmt=item[3],
                ecolor='k',
                label = item[4] )
    
    plt.xlabel(xlabel,fontsize = 20)
    plt.ylabel(ylabel,fontsize = 20)
    plt.title(title)
    plt.legend()
    plt.show()
    return
""" END función Graficador """


class AnalysisVicsek(object):
    
    
    def __init__(self,runNbr,inFiles,typeAnalysis,parameters):
        self.runNbr = runNbr
        self.inFiles = inFiles
        assert len(inFiles)==3, 'Must privide a 3 length list'
        self.typeAnalysis = typeAnalysis
        self.parameters = parameters
        self.M_Period = parameters[3] #period of messurements
#        self.tauEq = 0

        return


    #from inFiles read Parameters and Variables  
    def getVariables(self,fileToRead):
        with open(fileToRead, 'r') as archivo:  #notar que se abre con 'r'
            #read noise & convert to np.array/float
            noise=np.asarray(archivo.readline().split()[1::],float)
            datos = pd.read_csv(
                    archivo,
                    delimiter = ' ',
                    index_col = 1,
                    header=0,
                    )
            
            #set time as index & save in variable
            time=datos.index
    
            # create a list to store each column from data
            X=np.array([datos[j] for j in datos.keys()[1:]],float)
            
        return noise, time, X
    
    # Set noise
    @property
    def vicsekNoise(self):
        return self.getVariables(self.inFiles[0])[0]
    
    def _calc_acorr(self,tEq):
        print('tEq is a presumed time of equilibration different than tauEq')
        first=int(tEq/self.M_Period)
        Acorr=[]
        OP=self.get_OrderParam
        for k in range(len(self.vicsekNoise)):
            # estimate Acorr from an stabilized array
            Acorr.append(Autocorr(OP[k][first::]))
    
        return np.array(Acorr,float)
    
    @property
    def acorrTime(self,t_f=250):
        # we save only up to len(OP)/2 because after that we have little
        # amount of points to estimate Acorr with few error.
#        t_f=250
        time=self.getVariables(self.inFiles[0])[1][:t_f]

        return np.asarray(time,float)
    
    @property
    def acorr(self):
        Acorr = self.getVariables(self.inFiles[0])[2][:,:250]
        return np.asarray(Acorr,float)
    
    
#    @property
    #Now we need to extract the equilibrium time 'tauEq'.
    def get_tauEq(self,FittingFunction,FittingSeed):
        #check dimensions
        assert len(FittingSeed)==len(self.vicsekNoise),"Must pass an array same size 'noise'"
        #check the way the function is defined
        print('CAREFUL: FittingFunction must be like F(x,tau,A,y0,**p)')
        #falta proteger la variable tauEq
        tauEq=ExtractTauEq(
                self.vicsekNoise,
                FittingFunction,
                self.acorrTime,
                self.acorr,
                FittingSeed,
                self.M_Period
                )
#        self.tauEq = np.array(tauEq,float)
        print(tauEq)
        return np.array(tauEq,float)
    
#    _tauEq = get_tauEq
    
    #get order parameter time array
    @property
    def get_OPtime(self):    
        return self.getVariables(self.inFiles[1])[1]
    #get order parameter array
    @property
    def get_OrderParam(self):
        return self.getVariables(self.inFiles[1])[2]
    #get order parameter squared array
    @property
    def get_OPSqr(self):
        return self.getVariables(self.inFiles[2])[2]
    
    
    # calculate the mean,var and kurtosis of the variable X from an initial time tauEq
    def _calc_X_Stats(self,X,tauEq):   
        X_mean=[]
        X_var=[]
        X_kurt=[]
        for k in range(len(self.vicsekNoise)):
            tauEq[k]/=self.M_Period
            X_mean.append(np.mean(X[k][int(tauEq[k]):]))
            X_var.append(np.var(X[k][int(tauEq[k]):]))
            X_kurt.append(kurtosis(X[k][int(tauEq[k]):],fisher=True))
            tauEq[k]*=self.M_Period #idk how to avoid de method to change the original value.. so re calc
        X_mean=np.asarray(X_mean,float)
        X_var=np.asarray(X_var,float)
        X_kurt=np.asarray(X_kurt,float)
        return X_mean,X_var,X_kurt
       
    def get_OP_Mean(self,tauEq):
        OPMean = self._calc_X_Stats(self.get_OrderParam,tauEq)[0]
        return OPMean
    
    def get_OP_Var(self,tauEq):
        OPVar = self._calc_X_Stats(self.get_OrderParam,tauEq)[1]
        return OPVar
    
    """ CumulanteBinder = CB_A = 1- <Phi^4>/3*<Phi^2>^2"""
    def _calc_X_BinderCumulant(self,tauEq):
        BindCum = 1 -self._calc_X_Stats(self.get_OrderParam,tauEq)[2]/3
        return BindCum
    
# =============================================================================
#     """Este método utiliza dos archivos y calcula 'a mano' var """
#     def compute_OP_var_TOW(self,tauEq):
#         #compute <OP>^2
#         OPMeanSqr=self.get_OP_Mean(tauEq)**2
#         #compute <OP^2>
#         OPSqrMean=self._calc_X_Mean(self.get_OPSqr(),tauEq)
#         #
#         var = OPSqrMean - OPMeanSqr
#         return var
# =============================================================================
# =============================================================================
#     def get_CB_TOW(self,tauEq):
#         X4th= self.get_OrderParam**4
#         X4th_mean=[]
#         for k in range(len(self.vicsekNoise)):
#             tauEq[k]/=self.M_Period
#             X4th_mean.append(np.mean(X4th[k][int(tauEq[k]):]))
#             tauEq[k]*=self.M_Period
#         CB = 1 - X4th_mean/(3*self.get_OP_Var(tauEq))
#         return CB
# =============================================================================
    def ErrOP(self,tauEq):
        errOP=[]
        for k in range(len(self.vicsekNoise)):
            if tauEq[k]<=self.M_Period:
                n=len(self.get_OPtime)
                errOP.append(np.sqrt(self.get_OP_Var(tauEq)[k]/(n-1)))
            else:
                tmax=max(self.get_OPtime) - tauEq[k]
                n=2*tauEq[k]/tmax
                errOP.append(np.sqrt(n*self.get_OP_Var(tauEq)[k]))
    #            errOP = np.sqrt(self.get_OP_Var/(n-1))
        
        return np.array(errOP,float)
"""END OF CLASS AnalysisVicsek"""
# =============================================================================
# class AnalysisAcorr(AnalysisVicsek):
#     
#     def __init__(self,runNbr,inFiles,typeAnalysis,parameters,tEq):
#         super().__init__(runNbr,inFiles,typeAnalysis,parameters)
#         self.tEq=tEq
#     
#     def _calc_acorr(self):
#         print('tEq is a presumed time of equilibration different than tauEq')
#         first=int(self.tEq/self.M_Period)
#         Acorr=[]
#         OP=self.get_OrderParam
#         for k in range(len(self.vicsekNoise)):
#             # estimate Acorr from an stabilized array
#             Acorr.append(Autocorr(OP[k][first::]))
#     
#         return np.array(Acorr,float)
#     
#     @property
#     def acorrTime(self):
#         # we save only up to len(OP)/2 because after that we have little
#         # amount of points to estimate Acorr with few error.
#         t_f=int(len(self.get_OrderParam)/2)
#         time=self.getVariables(self.inFiles)[1][:t_f]
# 
#         return np.asarray(time,float)
#     
#     @property
#     def acorr(self):
#         return self.getVariables(self.inFiles)[2]
#     
#     
# #    @property
#     #Now we need to extract the equilibrium time 'tauEq'.
#     def get_tauEq(self,FittingFunction,FittingSeed):
#         #check dimensions
#         assert len(FittingSeed)==len(self.vicsekNoise),"Must pass an array same size 'noise'"
#         #check the way the function is defined
#         print('CAREFUL: FittingFunction must be like F(x,tau,A,y0,**p)')
#         #falta proteger la variable tauEq
#         tauEq=ExtractTauEq(
#                 self.vicsekNoise,
#                 FittingFunction,
#                 self.acorrTime,
#                 self.acorr,
#                 FittingSeed,
#                 self.M_Period
#                 )
# #        self.tauEq = np.array(tauEq,float)
#         print(tauEq)
#         return np.array(tauEq,float)
#     
# class AnalysisOP(AnalysisVicsek):
#     pass    
# =============================================================================




""" old routines """ 
# =============================================================================
# def ValorMedio(X,ruido,tEq,MidoCada):
# # =============================================================================
# #     tEq/=MidoCada
# #     XMean=[]
# #     for entero in range(len(ruido)):
# #         XMean.append(np.mean(X[entero::len(ruido)][int(tEq):]))
# #     return XMean
# # =============================================================================
#     
#     XMean=[]
#     for entero in range(len(ruido)):
#         tEq[entero]/=MidoCada
#         XMean.append(np.mean(X[entero::len(ruido)][int(tEq[entero]):]))
#     return XMean
# =============================================================================

# =============================================================================
# #def LecturaTablasSim(tableName,noise,time,X,printing):
# def LecturaTablasSim(tableName,printing):
# 
#     time=[]
#     Xaux=[]
#     with open(tableName,'r') as archivo:
#         plots=csv.reader(archivo,delimiter = ' ')  #abro el archivo: me guarda la información en una lista(array).    
#         
#         noise=archivo.readline().split()[1::]
#         for i in range(len(noise)):
#             noise[i] = float(noise[i])
#         
#         if printing == True:
#             print('El Ruido es:')
#             print(noise)
# #        time=[]
#         
# 
# #        Xaux=[]
#         next(plots)
#         for line in plots: #leo las tres columnas\n
#             time.append(float(line[1]))
#             for entero in range(len(noise)):
#                 Xaux.append(float(line[2::][entero])) # el segundo [] te permite elegir la columna 
# #                Xaux.append(float(line[2+entero]))
#                 
#         #Tranformamos a X una lista de len(noise) componentes, c/u guarda
#         # la columna X(noise).
#         X=[[]]*len(noise)
#         for entero in range(len(noise)):
#             X[entero].append(Xaux[entero::len(noise)])
# 
# 
# #    return noise,time,X
#     return noise,time,Xaux
# #    return
# =============================================================================

