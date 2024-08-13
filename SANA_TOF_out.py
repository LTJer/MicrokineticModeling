import xlrd
import xlsxwriter
from RateConstantFormation import RateConstant
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy import errstate
import pandas as pd
from colorama import Fore, Back, Style
from functools import wraps
import time
from ase import units
from EnergyTool import  EaTransformation
from Entropy_partition_function import NNHgasEntropy, Entropy_HNNH_Z, \
    Entropy_N2H4, Entropy_N2, Entropy_H2, Entropy_NH3,Entropy_Natom, Entropy_Hatom,\
    Entropy_NH




def boltzmann_distN2(nlevels, Tv):
      

        ws = np.array([2372.45, 18.1017])
        cm1eV = 0.00012398426
        EN_vibs = cm1eV * ws

        levels = np.arange(nlevels)
        Ei_harmonic = levels * EN_vibs[0]  # energy of ith vibrational state
        Ei_anharmonic = EN_vibs[0] * levels - EN_vibs[1] * levels ** 2
        P_B = np.exp(-Ei_harmonic / units.kB / Tv)
        P_B = P_B / sum(P_B)
        return P_B
    
    
    
def boltzmann_distH2(nlevels, Tv):
      

        ws = np.array([4342])
        cm1eV = 0.00012398426
        EN_vibs = cm1eV * ws

        levels = np.arange(nlevels)
        Ei_harmonic = levels * EN_vibs[0]  # energy of ith vibrational state
        # Ei_anharmonic = EN_vibs[0] * levels - EN_vibs[1] * levels ** 2
        P_B = np.exp(-Ei_harmonic / units.kB / Tv)
        P_B = P_B / sum(P_B)
        return P_B
    
    
    
def k_Ey(T, TS_S, TS_H): #Eyring equation
    """
    Calculate reaction rate constant for a surface reaction in s-1

    T       - Temperature in K
    nu      - Pre-exponential factor in s^-1
    TS_S     - Entropies difference between TS and reactant in J/(K*mol)
    TS_H    - Activation energy in J/mol
    """
    R = 8.3144598 # gas constant J/(K*mol)
    kb = 1.38064852e-23 # boltzmann constant
    h = 6.62607004e-34  # planck constant

    # if TandSprint :  
    #     print(f'Apply entropy:{TS_S} and Temperature:{self.T} K')
        
    TS_G = TS_H - T*TS_S
    
    if TS_G < 0:
        TS_G = 0
        
    k =  kb * T / h  * np.exp(-TS_G / (R * T)) 
    ksci  = np.float64("{:.1e}".format(k))
    # k2 = np.format_float_scientific(k, exp_digits=1)
    
    return ksci

def treanor_distN2(nlevels, Tv, Tg):

        
        ws = np.array([2372.45, 18.1017])
        cm1eV = 0.00012398426
        EN_vibs = cm1eV * ws

        levels = np.arange(nlevels)
        Ei_harmonic = levels * EN_vibs[0]  # energy of ith vibrational state
        Ei_anharmonic = EN_vibs[0] * levels - EN_vibs[1] * levels ** 2

        P_T = np.exp(-Ei_harmonic / units.kB / Tv
                     + (Ei_harmonic - Ei_anharmonic) / units.kB / Tg)

        P_T = P_T / sum(P_T)
        return  P_T



def treanor_distH2(nlevels, Tv, Tg):

        
        ws = np.array([4342, 18.1017])
        cm1eV = 0.00012398426
        EN_vibs = cm1eV * ws

        levels = np.arange(nlevels)
        Ei_harmonic = levels * EN_vibs[0]  # energy of ith vibrational state
        Ei_anharmonic = EN_vibs[0] * levels - EN_vibs[1] * levels ** 2

        P_T = np.exp(-Ei_harmonic / units.kB / Tv
                     + (Ei_harmonic - Ei_anharmonic) / units.kB / Tg)

        P_T = P_T / sum(P_T)
        return  P_T


def nonnegative(f):
    @wraps(f)
    def wrapper(t, C, *args, **kwargs):
        negatives = C < 0
        C = np.maximum(C, np.zeros(np.shape(C)))

        recalculate = f(t, C, *args, **kwargs)
            
        recalculate[negatives] = np.maximum(recalculate[negatives], 0)

        return recalculate

    return wrapper

 
 
class Reaction_ODE() :
    
    def __init__(self,T = 400, metal = 'Fe', Exlfile='MKMInputSheet.xlsx',Exlsheet ='EnergyInput',\
                 N2vib = 'off', AAds = 'off', ER1 = 'off', SubDso = 'off', InitialAds = 'simple'\
                  ,LHW = 'off', ER2 = 'off', Dissociation = 'off', Desorption = 'off', \
                  StickingProbability=1, ConcRatio=1, DRC_label = None,TSR=1) :
          
        self.metal = metal
        self.Exlfile=Exlfile
        self.Exlsheet =Exlsheet
        #def temperature and pressures
        
        self.ConcRatio = ConcRatio
        self.T = T
        self.p0 = 1 # bar
        self.pa = 101325 #pressure in pa from bar
        self.pN2 = 0.683405 * 0.63656464 * self.p0  # pressure in bar
        self.pH2 = 0.300839 * 0.875 * self.p0 
        self.pNH3 = 0 * self.p0 
        self.pN = 0.00019 * self.p0 * self.ConcRatio 
        self.pH = 0.015148 * self.p0 * self.ConcRatio 
        self.pNH = 5.14e-7 * self.p0 * self.ConcRatio
 #############################
        #Excited species' pressures
        vibratio = 1
        self.pN2_v1 = 0.683405 * 0.21481658 * self.p0*vibratio #pressure in bar
        self.pN2_v2 = 0.683405 * 0.08035864 * self.p0 *vibratio
        self.pN2_v3 = 0.683405 * 0.03332245 * self.p0*vibratio
        self.pN2_v4 = 0.683405 * 0.01531724 * self.p0*vibratio
        self.pN2_v5 = 0.683405 * 0.00780484 * self.p0*vibratio
        self.pN2_v6 = 0.683405 * 0.00440846 * self.p0 *vibratio
        self.pN2_v7 = 0.683405 * 0.00276025 * self.p0 *vibratio
        self.pN2_v8 = 0.683405 * 0.0019158 * self.p0*vibratio
        self.pN2_v9 = 0.683405 * 0.00147398 * self.p0*vibratio
        self.pN2_v10 =0.683405 * 0.00125711 * self.p0*vibratio 
        
        
        self.pH2_v1 = 0.300839 * 0.109 * self.p0*vibratio
        self.pH2_v2 = 0.300839 * 0.0136 * self.p0 *vibratio
        self.pH2_v3 = 0.300839 * 0.00169 * self.p0*vibratio
        self.pH2_v4 = 0.300839 * 0.00021 * self.p0*vibratio
        self.pH2_v5 = 0.300839 * 2.63229847e-05 * self.p0*vibratio
        self.pH2_v6 = 0.300839 * 3.28068101e-06 * self.p0*vibratio
        self.pH2_v7 = 0.300839 * 4.08877185e-07 * self.p0 *vibratio
        self.pH2_v8 = 0.300839 * 5.09591002e-08 * self.p0*vibratio
        self.pH2_v9 = 0.300839 * 6.35112447e-09 * self.p0*vibratio
        self.pH2_v10 = 0.300839 * 7.91552086e-10 * self.p0*vibratio
        #Entropy manipulations
        AdsRatio = 1/3
        DesRatio = 2/3
        self.SN2 = round(Entropy_N2(self.T)*AdsRatio,0)
        self.SH2 = round(Entropy_H2(self.T)*AdsRatio,0)
        self.SNH3 = round(Entropy_NH3(self.T)*DesRatio,0)
        self.SHNNH = round(Entropy_HNNH_Z(self.T)*DesRatio,0)
        self.SH2NNH2 = round(Entropy_N2H4(self.T)*DesRatio,0)
        ##
    
        # print(Entropy_Natom(self.T)*AdsRatio,0)
        self.SN = round(Entropy_Natom(self.T)*AdsRatio,0)
        self.SH = round(Entropy_Hatom(self.T)*AdsRatio,0)
        self.SNH = round(Entropy_NH(self.T)*AdsRatio,0)
        self.SNNH = round(NNHgasEntropy(self.T,self.pa)*DesRatio,0)
        #Reading file
        self.Exlfile=Exlfile
        self.Exlsheet =Exlsheet
        ##Module Switch 
        self.InitialAds = InitialAds
        self.N2vib = N2vib
        self.AAds = AAds
        self.ER1 = ER1
        self.SubDso = SubDso
        self.LHW = LHW
        self.ER2 = ER2
        self.Dissociation = Dissociation
        self.Desorption = Desorption
        #Other Parameters
        self.SP = StickingProbability
        self.TSR = TSR
        self.Mratio  = 0
        
        print(f'fetch ConcRatio = {self.ConcRatio}')
        print(f'fetch StickingProb = {self.SP}')
        print(f'pN con.c = {self.pN}')
        print(f'pH con.c = {self.pH}')
        
        print(f'Temperature for rate constant at {self.T} K')
        RC_model = RateConstant(Exlfile=self.Exlfile,Exlsheet=self.Exlsheet,\
                                metal=self.metal,T=self.T,SN2=self.SN2,SH2=self.SH2\
                                    ,SN = self.SN,SH=self.SH,SNH3=self.SNH3,\
                                    SHNNH=self.SHNNH,SN2H4=self.SH2NNH2,SNH=self.SNH,SNNH=self.SNNH)
        
        
        if 'sim' in self.InitialAds.lower():
            print("\033[1;35m" + 'use simple initial adsorption'+'\033[39m')
            self.kS1f,self.kS1b,self.kS2f,self.kS2b = RC_model.k_ThermalSimple()
       
        if 'com' in self.InitialAds.lower():
            print("\033[1;35m"+'use complx initial adsorption'+'\033[39m')
            self.kC1f,self.kC1b,self.kC2f,self.kC2b,self.kC3f,self.kC3b,self.kC4f,self.kC4b\
            = RC_model.k_ThermalComplex()
            
        if 'mix' in self.InitialAds.lower(): 
            print("\033[1;35m" + 'use mix initial adsorption'+'\033[39m')
            self.kS1f,self.kS1b,self.kS2f,self.kS2b = RC_model.k_ThermalSimple()
            self.kC1f,self.kC1b,self.kC2f,self.kC2b,self.kC3f,self.kC3b,self.kC4f,self.kC4b\
            = RC_model.k_ThermalComplex()
            
        self.kH1f,self.kH1b,self.kH2f,self.kH2b,self.kH3f,self.kH3b=\
            RC_model.k_ThermalHydrogenation()
        
        if self.N2vib.lower() == 'on' :
            self.kNv1f,self.kNv1b,self.kNv2f,self.kNv2b,self.kNv3f,self.kNv3b,self.kNv4f,self.kNv4b,self.kNv5f,self.kNv5b, \
            self.kNv6f,self.kNv6b,self.kNv7f,self.kNv7b,self.kNv8f,self.kNv8b,self.kNv9f,self.kNv9b,self.kNv10f,self.kNv10b,\
            self.kHv1f,self.kHv1b,self.kHv2f,self.kHv2b,self.kHv3f,self.kHv3b,self.kHv4f,self.kHv4b,self.kHv5f,self.kHv5b,\
            self.kHv6f,self.kHv6b,self.kHv7f,self.kHv7b,self.kHv8f,self.kHv8b,self.kHv9f,self.kHv9b,self.kHv10f,self.kHv10b,\
                =RC_model.k_DissociativeAds_vib()
            
        if self.AAds.lower() == 'on' :
            self.kAd1f,self.kAd1b,self.kAd2f,self.kAd2b,self.kAd3f,self.kAd3b = RC_model.k_ActiveAdsorptions()
            
        if self.ER1.lower() == 'on' :
            self.kER11f,self.kER11b,self.kER12f,self.kER12b,self.kER13f,self.kER13b,self.kER14f,self.kER14b,\
                self.kER15f,self.kER15b,self.kER16f,self.kER16b = RC_model.k_ER_Reactions1()
                
        if self.SubDso.lower() == 'on':
            self.kSD1f,self.kSD1b,self.kSD2f,self.kSD2b = RC_model.k_SubsurfaceDissolution()

        if self.LHW.lower() == 'on':
            self.kLHW1f,self.kLHW1b,self.kLHW2f,self.kLHW2b,self.kLHW3f,self.kLHW3b,self.kLHW4f,self.kLHW4b,self.kLHW5f,self.kLHW5b,self.kLHW6f\
            ,self.kLHW6b,self.kLHW7f,self.kLHW7b,self.kLHW8f,self.kLHW8b,self.kLHW9f,self.kLHW9b,self.kLHW10f,self.kLHW10b,self.kLHW11f,self.kLHW11b\
            ,self.kLHW12f,self.kLHW12b = RC_model.k_LangmuirHinshelwood()

        if self.ER2.lower() == 'on' :
            self.kER21f,self.kER21b,self.kER22f,self.kER22b,self.kER23f,self.kER23b,self.kER24f,self.kER24b,self.kER25f,self.kER25b,\
                    self.kER26f,self.kER26b,self.kER27f,self.kER27b,self.kER28f,self.kER28b,self.kER29f,self.kER29b,self.kER30f,self.kER30b,\
                        self.kER31f,self.kER31b,self.kER32f,self.kER32b,self.kER33f,self.kER33b,self.kER34f,self.kER34b,self.kER35f,self.kER35b\
                            = RC_model.k_ER_Reactions2()
                            
        if self.Dissociation.lower() == 'on' :
            self.kDISO1f,self.kDISO1b,self.kDISO2f,self.kDISO2b,self.kDISO3f,self.kDISO3b\
                ,self.kDISO4f,self.kDISO4b,self.kDISO5f,self.kDISO5b,self.kDISO6f,self.kDISO6b\
                    ,self.kDISO7f,self.kDISO7b,self.kDISO8f,self.kDISO8b = RC_model.k_SurfaceDissociation()

        if self.Desorption.lower() == 'on' :
            self.kDes1f,self.kDes1b,self.kDes2f,self.kDes2b,self.kDes3f,self.kDes3b,\
                self.kDes4f,self.kDes4b = RC_model.k_Desorption()

        
        self.negate = 0
        ##################### General sticking probability  for ER 
        # self.SP_H = 8E-3
        # self.SP_N = 1E-2
        
        self.SP_H = self.SP
        self.SP_N = self.SP
        
        self.SP_Ads = 1#self.SP
        
        ##########
        self.kER11f = self.kER11f * self.SP_N
        self.kER12f = self.kER12f * self.SP_H
        self.kER13f = self.kER13f * self.SP_H
        self.kER14f = self.kER14f * self.SP_H
        self.kER15f = self.kER15f * self.SP_N
        self.kER16f = self.kER16f * self.SP_H
        
        self.kER21f = self.kER21f * self.SP_H
        self.kER22f = self.kER22f * self.SP_H
        self.kER23f = self.kER23f * self.SP_H
        self.kER24f = self.kER24f * self.SP_H
        self.kER25f = self.kER25f * self.SP_H
        self.kER26f = self.kER26f * self.SP_H
        self.kER27f = self.kER27f * self.SP_H
        self.kER28f = self.kER28f * self.SP_H
        self.kER29f = self.kER29f * self.SP_H
        self.kER30f = self.kER30f * self.SP_H
        self.kER31f = self.kER31f * self.SP_H
        self.kER32f = self.kER32f * self.SP_H
        self.kER33f = self.kER33f * self.SP_N
        self.kER34f = self.kER34f * self.SP_N
        ##########
        self.kAd1f  = self.kAd1f *  self.SP_Ads
        self.kAd2f  = self.kAd2f *  self.SP_Ads
        ##############
        #self.kER35f = self.kER35f * self.SP_N 
        ####################
        
        if StickingProbability != 1:
            #zero test
            self.kER35f = StickingProbability
            
        ###DRC section
        if DRC_label != None:
            
            delHpercent = 1 #%
            k_name     = f'k{DRC_label}'
            H          = EaTransformation.Ea_FetchFromExcel(self.Exlfile,\
                         self.Exlsheet, [DRC_label], self.metal)[:,1].astype(float)[0]
                
            if H <= 0 :
                H_add_delH = 1142.427981*delHpercent
            else:
                H_add_delH = H*(1+(0.01*delHpercent))
                
            TS_S = 0
            k_value = k_Ey(self.T, TS_S, H_add_delH)
            setattr(self, k_name, k_value) #equivalent to: self.k_name= k_value
            
            self.del_H = H_add_delH - H
            # print(k_name,H)
            
        else:
            self.del_H = 0
            
    def InitializeModel(self, t, C):
        
        '''
        The order of the term of the differential equation can change the 
        final result. Ex: the NH3 production is different when the adsorptions 
        are in the head or tail of the ODE.
        '''
        # surface coverage
        CN = C[0]
        CH = C[1]
        CNH = C[2]
        CNH2 = C[3]
        CNH3 = C[4]
        CNsub = C[5]
        CHsub = C[6]
        CN2 = C[7]
        CH2 = C[8]
        CNNH = C[9]
        CNNH2 = C[10]
        CHNNH = C[11]
        CNNH3 = C[12]
        CHNNH2 = C[13]
        CHNNH3 = C[14]
        CH2NNH2 = C[15]
        CH2NNH3 = C[16]
        
        
        
        CNsub = CNsub/self.TSR 
        CHsub = CHsub/self.TSR 

        Cacsub = 1.0 - CNsub - CHsub
        Cac = 1 -CN-CH-CNH-CNH2-CNNH-CNNH2-CHNNH-CNNH3-CHNNH2-CHNNH3-CH2NNH2-CH2NNH3-CN2-CH2-CNH3
        if Cac < 0 :
           self.negate += 1
           # print("\033[1;31m"+'Negative Cac duing integration'+'\033[39m')
        dCNdt = dCHdt = dCNHdt = dCNH2dt = dCNH3dt = 0
        dCNsubdt = dCHsubdt = 0
        dCN2dt = dCH2dt = 0
        dCNNHdt = dCNNH2dt = dCHNNHdt = dCNNH3dt = 0
        dCHNNH2dt = dCHNNH3dt = dCH2NNH2dt = dCH2NNH3dt = 0 
      
        '''define rate equations in s-1
        '''
        ##ThermalAdsSimple
        if 'sim' in self.InitialAds.lower():
            
            rS1f = self.kS1f * self.pN2 * (Cac**2) 
            rS1b = self.kS1b * (CN**2)  
            rS2f = self.kS2f * self.pH2 * (Cac**2)  
            rS2b = self.kS2b * (CH**2) 
            
            dCNdt += (2*(rS1f-rS1b))
            dCHdt += (2*(rS2f-rS2b))
            
        ##ThermalAdsComplex
        if 'com' in self.InitialAds.lower():
            
            rC1f = self.kC1f * self.pN2 * Cac 
            rC1b = self.kC1b * CN2
            rC2f = self.kC2f * self.pH2 * Cac
            rC2b = self.kC2b * CH2
            rC3f = self.kC3f * CN2 * Cac
            rC3b = self.kC3b * (CN**2)
            rC4f = self.kC4f * CH2 * Cac
            rC4b = self.kC4b * (CH**2)
            
            dCN2dt += rC1f - rC1b - rC3f + rC3b
            dCH2dt += rC2f - rC2b - rC4f + rC4b
            dCNdt += (2*rC3f) - (2*rC3b) 
            dCHdt += (2*rC4f) - (2*rC4b) 
        
        if 'mix' in self.InitialAds.lower():
            
            rC1f = self.kC1f * self.pN2 * Cac 
            rC1b = self.kC1b * CN2
            rC3f = self.kC3f * CN2 * Cac
            rC3b = self.kC3b * (CN**2)
            
            rS2f = self.kS2f * self.pH2 * (Cac**2)  
            rS2b = self.kS2b * (CH**2)

            
            dCN2dt += rC1f - rC1b - rC3f + rC3b
            dCNdt += (2*rC3f) - (2*rC3b) 
            dCHdt += (2*(rS2f-rS2b))
            
        if self.SubDso.lower() == 'on':
                
            rSD1f = self.kSD1f * Cacsub * CN
            rSD1b = self.kSD1b * CNsub * Cac
            rSD2f = self.kSD2f * Cacsub * CH
            rSD2b = self.kSD2b * CHsub * Cac
            
            dCNdt += - rSD1f + rSD1b
            dCHdt += - rSD2f + rSD2b
            dCNsubdt += rSD1f - rSD1b
            dCHsubdt += rSD2f - rSD2b
            
                
        if self.Desorption.lower() == 'on' :
        ##Desorption
  
            rDes1f = self.kDes1f * CNH3 
            # rDes1b = self.kDes1b * self.pNH3 * Cac
            rDes2f = self.kDes2f * CHNNH  
            # rDes2b = self.kDes2b * self.pHNNH * Cac
            rDes3f = self.kDes3f * CH2NNH2 
            # rDes3b = self.kDes3b * self.pH2NNH2 * Cac
            rDes4f = self.kDes4f * CNNH 
            # rDes4b = self.kDes4b * self.pNNH * Cac
            
            dCHNNHdt += -rDes2f
            dCH2NNH2dt += -rDes3f
            dCNNHdt += -rDes4f
            dCNH3dt += -rDes1f
        else:
            self.kDes2f = self.kDes3f = self.rDes4f = self.rDes1f = 0
            
        ##ThermalHydrogenation
           
        rH1f = self.kH1f * CN * CH  
        rH1b = self.kH1b * CNH * Cac  
        rH2f = self.kH2f * CNH * CH  
        rH2b = self.kH2b * CNH2 * Cac  
        #NH2*+ H* ⟶ NH3* + *, NH3 desorbe directly
        rH3f = self.kH3f * CNH2 * CH 
        rH3b = self.kH3b* CNH3 * Cac

        dCNdt += -rH1f + rH1b #+ 2*(rS1f-rS1b)
        dCHdt += -rH1f + rH1b - rH2f + rH2b - rH3f + rH3b #+ 2*(rS2f-rS2b)
        dCNHdt += rH1f - rH1b - rH2f + rH2b
        dCNH2dt += rH2f - rH2b - rH3f + rH3b
        dCNH3dt += rH3f - rH3b

        if self.N2vib.lower() == 'on' :
        ##DissociativeAdsoprion(Vib)  
        #N2

            rNv1f = self.kNv1f * self.pN2_v1 * (Cac**2)  
            rNv1b = self.kNv1b * (CN**2)  
            rNv2f = self.kNv2f * self.pN2_v2 * (Cac**2)  
            rNv2b = self.kNv2b * (CN**2)  
            rNv3f = self.kNv3f * self.pN2_v3 * (Cac**2)  
            rNv3b = self.kNv3b * (CN**2)
            rNv4f = self.kNv4f * self.pN2_v4 * (Cac**2)  
            rNv4b = self.kNv4b * (CN**2) 
            rNv5f = self.kNv5f * self.pN2_v5 * (Cac**2)  
            rNv5b = self.kNv5b * (CN**2) 
            rNv6f = self.kNv6f * self.pN2_v6 * (Cac**2)  
            rNv6b = self.kNv6b * (CN**2)  
            rNv7f = self.kNv7f * self.pN2_v7 * (Cac**2)  
            rNv7b = self.kNv7b * (CN**2)  
            rNv8f = self.kNv8f * self.pN2_v8 * (Cac**2)  
            rNv8b = self.kNv8b * (CN**2)
            rNv9f = self.kNv9f * self.pN2_v9 * (Cac**2)  
            rNv9b = self.kNv9b * (CN**2) 
            rNv10f = self.kNv10f * self.pN2_v10 * (Cac**2)  
            rNv10b = self.kNv10b * (CN**2)
            
            #H2
            rHv1f = self.kHv1f * self.pH2_v1 * (Cac**2)  
            rHv1b = self.kHv1b * (CH**2)  
            rHv2f = self.kHv2f * self.pH2_v2 * (Cac**2)  
            rHv2b = self.kHv2b * (CH**2)   
            rHv3f = self.kHv3f * self.pH2_v3 * (Cac**2)  
            rHv3b = self.kHv3b * (CH**2) 
            rHv4f = self.kHv4f * self.pH2_v4 * (Cac**2)  
            rHv4b = self.kHv4b * (CH**2)
            rHv5f = self.kHv5f * self.pH2_v5 * (Cac**2)  
            rHv5b = self.kHv5b * (CH**2)
            rHv6f = self.kHv6f * self.pH2_v6 * (Cac**2)  
            rHv6b = self.kHv6b * (CH**2)  
            rHv7f = self.kHv7f * self.pH2_v7 * (Cac**2)  
            rHv7b = self.kHv7b * (CH**2)   
            rHv8f = self.kHv8f * self.pH2_v8 * (Cac**2)  
            rHv8b = self.kHv8b * (CH**2) 
            rHv9f = self.kHv9f * self.pH2_v9 * (Cac**2)  
            rHv9b = self.kHv9b * (CH**2)
            rHv10f = self.kHv10f * self.pH2_v10 * (Cac**2)  
            rHv10b = self.kHv10b * (CH**2)
                
            dCNdt += (2*rNv1f) - (2*rNv1b) + (2*rNv2f) - (2*rNv2b) + (2*rNv3f) - (2*rNv3b) + (2*rNv4f) - (2*rNv4b) + (2*rNv5f) - (2*rNv5b)\
                   + (2*rNv6f) - (2*rNv6b) + (2*rNv7f) - (2*rNv7b) + (2*rNv8f) - (2*rNv8b) + (2*rNv9f) - (2*rNv9b) + (2*rNv10f) - (2*rNv10b)
                   
            dCHdt += (2*rHv1f) - (2*rHv1b) + (2*rHv2f) - (2*rHv2b) + (2*rHv3f) - (2*rHv3b) + (2*rHv4f) - (2*rHv4b) + (2*rHv5f) - (2*rHv5b)\
                   + (2*rHv6f) - (2*rHv6b) + (2*rHv7f) - (2*rHv7b) + (2*rHv8f) - (2*rHv8b) + (2*rHv9f) - (2*rHv9b) + (2*rHv10f) - (2*rHv10b)
        
        if self.AAds.lower() == 'on' :
        ##ActiveAdsorptions
    
            rAd1f = self.kAd1f * self.pN * Cac 
            rAd1b = self.kAd1b * CN  
            rAd2f = self.kAd2f * self.pH * Cac 
            rAd2b = self.kAd2b * CH 
            # rAd3f =0#self.kAd3f * self.pNH * Cac 
            # rAd3b =0#self.kAd3b * CNH
            
            dCNdt += rAd1f - rAd1b
            dCHdt += rAd2f - rAd2b
            # dCNHdt += rAd3f - rAd3b
        
        if self.ER1.lower() == 'on' :
            #Basic ER reactions
            rER11f = self.kER11f * self.pN * CH 
            rER11b = self.kER11b * CNH 
            rER12f = self.kER12f * self.pH * CN
            rER12b = self.kER12b * CNH 
            rER13f = self.kER13f * self.pH * CNH 
            rER13b = self.kER13b * CNH2 
            rER14f = self.kER14f * self.pH * CNH2 
            rER14b = self.kER14b * CNH3
            rER15f = self.kER15f * self.pN * CN 
            rER15b = self.kER15b * CN2  # can change to gas species
            rER16f = self.kER16f * self.pH * CH  
            rER16b = self.kER16b * CH2  # can change to gas species
        
            dCNdt += -rER12f + rER12b - rER15f + rER15b
            dCHdt += -rER11f + rER11b - rER16f + rER16b
            dCNHdt += +rER11f - rER11b + rER12f - rER12b - rER13f + rER13b
            dCNH2dt += + rER13f - rER13b - rER14f + rER14b
            dCNH3dt += rER14f - rER14b
            dCN2dt +=  +rER15f - rER15b
            dCH2dt +=  +rER16f - rER16b
        else:
            self.kER14f = 0
        
        

        
        if self.LHW.lower() == 'on':
    
            rLHW1f = self.kLHW1f * CN2 * CH
            rLHW1b = self.kLHW1b * CNNH * Cac
            rLHW2f = self.kLHW2f * CNNH * CH
            rLHW2b = self.kLHW2b * CNNH2 * Cac
            rLHW3f = self.kLHW3f * CNNH * CH
            rLHW3b = self.kLHW3b * CHNNH * Cac
            rLHW4f = self.kLHW4f * CNNH2 * CH
            rLHW4b = self.kLHW4b * CNNH3 * Cac
            rLHW5f = self.kLHW5f * CNNH2 * CH
            rLHW5b = self.kLHW5b * CHNNH2 * Cac
            rLHW6f = self.kLHW6f * CHNNH * CH
            rLHW6b = self.kLHW6b * CHNNH2 * Cac 
            rLHW7f = self.kLHW7f * CNNH3 * CH
            rLHW7b = self.kLHW7b * CHNNH3 * Cac
            rLHW8f = self.kLHW8f * CHNNH2 * CH
            rLHW8b = self.kLHW8b * CHNNH3 * Cac
            rLHW9f = self.kLHW9f * CHNNH2 * CH
            rLHW9b = self.kLHW9b * CH2NNH2 * Cac
            rLHW10f = self.kLHW10f * CHNNH3 * CH
            rLHW10b = self.kLHW10b * CH2NNH3 * Cac
            rLHW11f = self.kLHW11f * CH2NNH2 * CH
            rLHW11b = self.kLHW11b * CH2NNH3 * Cac
            rLHW12f = self.kLHW12f * CH2NNH3 * CH
            rLHW12b = self.kLHW12b * (CNH3 ** 2)  
            
            dCHdt += - rLHW1f + rLHW1b - rLHW2f + rLHW2b - rLHW3f + rLHW3b 
            -rLHW4f + rLHW4b - rLHW5f + rLHW5b -rLHW6f + rLHW6b -rLHW7f + rLHW7b 
            -rLHW8f + rLHW8b -rLHW9f + rLHW9b -rLHW10f + rLHW10b -rLHW11f + rLHW11b -rLHW12f + rLHW12b
            
            dCN2dt += -rLHW1f + rLHW1b
            dCNNHdt += rLHW1f - rLHW1b - rLHW2f + rLHW2b - rLHW3f + rLHW3b
            dCNNH2dt += rLHW2f - rLHW2b -rLHW4f + rLHW4b -rLHW5f + rLHW5b
            dCHNNHdt += rLHW3f - rLHW3b -rLHW6f + rLHW6b
            dCNNH3dt += rLHW4f - rLHW4b - rLHW7f + rLHW7b 
            dCHNNH2dt += rLHW5f - rLHW5b + rLHW6f - rLHW6b - rLHW8f + rLHW8b -rLHW9f + rLHW9b
            dCHNNH3dt += rLHW7f - rLHW7b + rLHW8f - rLHW8b -rLHW10f + rLHW10b 
            dCH2NNH2dt += rLHW9f - rLHW9b -rLHW11f + rLHW11b
            dCH2NNH3dt += rLHW10f - rLHW10b + rLHW11f - rLHW11b -rLHW12f + rLHW12b
            dCNH3dt += 2*rLHW12f-2*rLHW12b
            
        else:
            self.kLHW12f = 0
        
        if self.ER2.lower() == 'on' :

            rER21 = self.kER21f*self.pH*CN2  - self.kER21b*CNNH
            rER22 = self.kER22f*self.pH*CNNH   - self.kER22b*CNNH2
            rER23 = self.kER23f*self.pH*CNNH   - self.kER23b*CHNNH
            rER24 = self.kER24f*self.pH*CNNH2 - self.kER24b*CNNH3
            rER25 = self.kER25f*self.pH*CNNH2 - self.kER25b*CHNNH2
            rER26 = self.kER26f*self.pH*CNNH3 - self.kER26b*CHNNH3
            rER27 = self.kER27f*self.pH*CHNNH -self.kER27b*CHNNH2
            rER28 = self.kER28f*self.pH*CHNNH2 - self.kER28b*CHNNH3
            rER29 = self.kER29f*self.pH*CHNNH2 -self.kER29b*CH2NNH2
            rER30 = self.kER30f*self.pH*CH2NNH2  -self.kER30b*CH2NNH3
            rER31 = self.kER31f*self.pH*CHNNH3  -self.kER31b*CH2NNH3
            rER32 = self.kER32f*self.pH*CH2NNH3 -self.kER32b*(CNH3 ** 2)
            rER33 = self.kER33f*CNH*self.pN - self.kER33b*CNNH
            rER34 = self.kER34f*CNH2*self.pN - self.kER34b*CNNH2
            rER35 =0# self.kER35f*CH*self.pN2 - self.kER35b*CNNH
            
            #dCHdt += -rER35
            dCNHdt += -rER33
            dCNH2dt += -rER34         
            dCN2dt += -rER21
            dCNNHdt += rER21 - rER22 - rER23 + rER33 #+ rER35
            dCNNH2dt += rER22 - rER24 -rER25 +rER34
            dCNNH3dt += rER24 - rER26
            dCHNNH2dt += rER25 - rER28 -rER29 + rER27
            dCHNNHdt += rER23-rER27
            dCHNNH3dt += rER26 + rER28 - rER31
            dCH2NNH2dt += rER29 - rER30
            dCH2NNH3dt += rER30 + rER31 -rER32
            dCNH3dt += 2*rER32
            # print(rER32)
        else:
            self.kER32f = 0
            
            
        if self.Dissociation.lower() == 'on' :   
            
            rDISO1f = self.kDISO1f * CNNH *Cac
            rDISO1b = self.kDISO1b * CN * CNH
            rDISO2f = 0#self.kDISO2f * CNNH * Cac
            rDISO2b = 0#self.kDISO2b * CN2 * CH
            rDISO3f = self.kDISO3f * CNNH2 * Cac
            rDISO3b = self.kDISO3b * CN * CNH2
            rDISO4f = self.kDISO4f * CNNH3 * Cac
            rDISO4b = self.kDISO4b * CN * CNH3
            rDISO5f = self.kDISO5f * CHNNH * Cac
            rDISO5b = self.kDISO5b * (CNH ** 2)
            rDISO6f = self.kDISO6f * CH2NNH2 * Cac
            rDISO6b = self.kDISO6b * (CNH2 **2 )
            rDISO7f = self.kDISO7f * CHNNH3 * Cac
            rDISO7b = self.kDISO7b * CNH * CNH3
            rDISO8f = self.kDISO8f * CH2NNH3 * Cac
            rDISO8b = self.kDISO8b * CNH2 * CNH3
        
            dCNdt += rDISO1f - rDISO1b + rDISO3f - rDISO3b + rDISO4f- rDISO4b
            dCN2dt += rDISO2f - rDISO2b
            dCHdt += rDISO2f - rDISO2b
            dCNHdt += rDISO1f - rDISO1b + 2*(rDISO5f-rDISO5b) + rDISO7f -rDISO7b
            dCNH2dt += rDISO3f -rDISO3b + 2*(rDISO6f-rDISO6b) + rDISO8f -rDISO8b
            dCNH3dt += rDISO4f - rDISO4b + rDISO7f - rDISO7b + rDISO8f - rDISO8b
            dCNNHdt += -rDISO1f + rDISO1b -rDISO2f + rDISO2b 
            dCNNH2dt += -rDISO3f + rDISO3b
            dCNNH3dt += -rDISO4f + rDISO4b
            dCHNNHdt += -rDISO5f + rDISO5b
            dCH2NNH2dt += -rDISO6f + rDISO6b
            dCHNNH3dt += -rDISO7f + rDISO7b
            dCH2NNH3dt += -rDISO8f + rDISO8b
    
        else:
            
            self.kDISO4f = self.kDISO7f = self.kDISO8f = 0
            
            
            
        return [ dCNdt, dCHdt, dCNHdt, dCNH2dt, dCNH3dt, dCNsubdt, dCHsubdt,dCN2dt, dCH2dt,\
                dCNNHdt, dCNNH2dt,dCHNNHdt, dCNNH3dt , dCHNNH2dt , dCHNNH3dt , dCH2NNH2dt ,dCH2NNH3dt]
            
        
    def AnalyzeData (self, Conc) :
            
            CN = Conc[:,0]
            CH = Conc[:,1] 
            CNH = Conc[:,2]
            CNH2 = Conc[:,3]
            CNH3 = Conc[:,4]
            CNsub = Conc[:,5]
            CHsub = Conc[:,6]
            CN2 = Conc[:,7]
            CH2 = Conc[:,8]
            CNNH = Conc[:,9]
            CNNH2 = Conc[:,10]
            CHNNH = Conc[:,11]
            CNNH3 = Conc[:,12]
            CHNNH2 = Conc[:,13]
            CHNNH3 = Conc[:,14]
            CH2NNH2 = Conc[:,15]
            CH2NNH3 = Conc[:,16]
            Cac = 1 -CN-CH-CNH-CNH2-CNNH-CNNH2-CHNNH-CNNH3-CHNNH2-CHNNH3-CH2NNH2-CH2NNH3-CN2-CH2-CNH3
            Cacsub = 1 - CNsub - CHsub
     
            print(f'{self.negate} negative Cac duing integration')
            print("\033[41m"+f'{len(Conc[Conc<0])} negative coverages after integration'+'\033[0m')
            print("\033[41m"+f'{len(Cac[Cac<0])} negative Cac after analyzer'+'\033[0m')
            
            
            # else:
            #     print('Cac is positive')
            # rH3f = self.kH3f * CNH2 * CH 
            # rH3b = self.kH3b* CNH3 * Cac
            # rER14f = self.kER14f * self.pH * CNH2
            # rER14b = self.kER14b * CNH3
            # rLHW12f = self.kLHW12f * CH2NNH3 * CH
            # rLHW12b = self.kLHW12b * (CNH3 ** 2)  
            # rER32f = self.kER32f*self.pH*CH2NNH3
            # rER32b = self.kER32b*(CNH3 ** 2)  
            # rDISO4f = self.kDISO4f * CNNH3 * Cac
            # rDISO4b = self.kDISO4b * CN * CNH3
            # rDISO7f = self.kDISO7f * CHNNH3 * Cac
            # rDISO7b = self.kDISO7b * CNH * CNH3
            # rDISO8f = self.kDISO8f * CH2NNH3 * Cac
            # rDISO8b = self.kDISO8b * CNH2 * CNH3  
            ######
            # rH3 = rH3f[-1]-rH3b[-1]
            # rER14 = rER14f[-1]-rER14b[-1]
            # rER32 = rER32f[-1]-rER32b[-1]
            # rLHW12 = rLHW12f[-1]-rLHW12b[-1]
            # rDISO4 = rDISO4f[-1]-rDISO4b[-1]
            # rDISO7 = rDISO7f[-1]-rDISO7b[-1]
            # rDISO8 = rDISO8f[-1]-rDISO8b[-1]
            # print(f'ressible checking rH3 : {rH3}')
            # print(f'ressible checking rER14 : {rER14}')
            # print(f'ressible checking rER32 : {rER32}')
            # print(f'ressible checking rLHW12 : {rLHW12}')
            # print(f'ressible checking rDISO4 : {rDISO4}')
            # print(f'ressible checking rDISO7 : {rDISO7}')
            # print(f'ressible checking rDISO8 : {rDISO8}')
            ##Desorptions
            rDes1f = self.kDes1f * CNH3 
            rDes2f_HNNH = self.kDes2f * CHNNH  
            rDes3f_H2NNH2 = self.kDes3f * CH2NNH2 
            rDes4f_NNH = self.kDes4f * CNNH 
            #Radical Ads
            # rAd1f = self.kAd1f * self.pN * Cac 
            # rAd1b = self.kAd1b * CN  
            # print('Con.c',CN[-1], Cac[-1])
            # print('rate constant',self.kAd1f,self.kAd1b)
            # print('rate',rAd1f[-1],rAd1b[-1])
            ###
            # rNv1f = self.kNv1f * self.pN2_v1 * (Cac**2)  
            # rNv2f = self.kNv2f * self.pN2_v2 * (Cac**2)  
            # rNv3f = self.kNv3f * self.pN2_v3 * (Cac**2)  
            # rC1f = self.kC1f * self.pN2 * Cac
            # print('RateConstant',self.kC1f,self.kNv1f,self.kNv2f,self.kNv3f )
            # print('Rate',rC1f[-1],rNv1f[-1],rNv2f[-1],rNv3f[-1])
            ## SubsurfaceDissolution
            # rSD1f = self.kSD1f * Cacsub * CN
            # rSD1b = self.kSD1b * CNsub * Cac
            # rSD2f = self.kSD2f * Cacsub * CH
            # rSD2b = self.kSD2b * CHsub * Cac
            # print("\033[36m"+ f'rf = {rSD2f[-1]}'+'\033[39m')  
            # print("\033[36m"+ f'rb = {rSD2b[-1]}'+'\033[39m') 
            #########################
            dCNH3gdt_form = 0#rH3 + rER14 + 2*rER32 + 2*rLHW12 #+ rDISO4 + rDISO7 + rDISO8
            dCNH3gdt_des = rDes1f
            #coverage
            Con_c =  np.array([CN, CH, CNH, CNH2, Cac, CNH3, CN2, CH2, CNsub, CHsub, CNNH, CHNNH, CH2NNH2]).T
 
            return Con_c, dCNH3gdt_form, dCNH3gdt_des, rDes4f_NNH, rDes2f_HNNH, rDes3f_H2NNH2 #,CNNH, CHNNH, CH2NNH2



def makeplot(Conc, metal):
    """
    plot for surface coverages
    """
    #take log with the surface coverage
    #[CN, CH, CNH, CNH2, Cac, CNH3g]
    C_p = np.copy(Conc)
    
    filter4 = (C_p[0,:] == 0)
    C_p[0,:][filter4] = 1
    filter6 = (C_p[1:,:] < 0)
    C_p[1:,:][filter6] = 1E-90
    
    with errstate(divide='ignore'):
        
        Clog = np.log10(C_p)
   
    #take log with the data
    #CN, CH, CNH, CNH2, Cac, CNH3g, CN2, CH2, CNsub, CHsub, CNNH, CHNNH, CH2NNH2
    CNlog = Clog[:,0]
    CHlog = Clog[:,1] 
    CNHlog = Clog[:,2]
    CNH2log = Clog[:,3]
    Caclog = Clog[:,4]
    CNH3glog = Clog[:,5]
    CN2log = Clog[:,6]
    CH2log = Clog[:,7]
    CNsublog = Clog[:,8]
    CHsublog = Clog[:,9] 
    CNNHlog = Clog[:,10]
    CHNNHlog = Clog[:,11]
    CH2NNH2log = Clog[:,12]

    
    area = np.pi*0.1
    #plt.subplot(2,1,1)
    #plt.scatter(t, CNlog, s=area, alpha=0.1, c = 'g')
    fig, ax = plt.subplots()
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2.5) 
        ax.spines[axis].set_color('black') 
    ax.plot(t, CNlog, c = 'b')
    ax.plot(t, CHlog, c = 'limegreen')
    ax.plot(t, CNHlog, c = 'cornflowerblue')
    ax.plot(t, CNH2log, c = 'lightsteelblue')
    ax.plot(t, Caclog, '--',color='grey')
    ax.plot(t, CN2log,'+', color='darkblue')
    ax.plot(t, CH2log, 'g+')
    ax.plot(t, CNsublog,'x', color='gold')
    ax.plot(t, CHsublog, 'x', color = 'red')
    # ax.plot(t, CNNHlog, 'x', color='red')
    # ax.plot(t, CHNNHlog, '1', color='gold')
    # ax.plot(t, CH2NNH2log, '2', color='orange')
    ax.set_xlabel('Time'+' '+'$s^{-1}$')
    ax.set_ylabel('Coverage Ratior' +' '+'$log_{10}$')
    ax.set_ylim([-90, 0.5])
    #plt.ylim(-8, 5)
    ax.legend([r'$\Theta_{N}$',r'$\Theta_{H}$',r'$\Theta_{NH}$', r'$\Theta_{NH2}$',
                r'$\Theta_{*}$', r'$\Theta_{N_{2}}$',r'$\Theta_{H_{2}}$', 
                r'$\Theta_{N_{sub}}$',r'$\Theta_{H_{sub}}$'], loc = 0)
               # r'$\Theta_{NNH}$', r'$\Theta_{HNNH}$', r'$\Theta_{H2NNH2}$'], loc = 0)
    #plt.text(4, -5,'-10 = negative coverage')
    ax.annotate(metal, xy=(-8, -8), xycoords='axes points',
            size=10, ha='right', va='top',
            bbox=dict(boxstyle='round', fc='w'))
    
    plt.show()

 
def title_creater(label_list):
    # Head = 'Metal ΘN ΘH ΘNH ΘNH2 Θ* [NH3(g)] d[NH3(g)]/dt'
    # for c_head, label in enumerate(Head.split(' ')) : 
    for c_head, label in enumerate(label_list) : 
        writesheet.write(0, c_head+1, label) 
        
def excel_output(count, metal, TotalTOF):
    for counter, species in enumerate(TotalTOF.split(',')) : 
            writesheet.write(count+1, 0, str(metal))
            writesheet.write(count+1, counter+1, species) 
        
    
        
def title_creater_SP_ConcR(labname1,label_list1,labname2,label_list2):

    writesheet.write(0, 0, labname1)
    writesheet.write(1, 0, labname2)
    count = 0
    for c_head1, label1 in enumerate(label_list1) : 
        for c_head2, label2 in enumerate(label_list2) : 
            count += 1
            writesheet.write(0, count, label1) 
            writesheet.write(1, count, label2)
        
        
def excel_output_SP_ConcR(rowcount, colcount, metal, DesTOF):
    
    # for counter, species in enumerate(TotalTOF.split(',')) : 
    writesheet.write(rowcount+2,0, str(metal))
    writesheet.write(rowcount+2, colcount,  DesTOF)    
    
    
    
    
    
###################---Main---########################

odelist = ['CN', 'CH', 'CNH', 'CNH2', 'CNH3', 'CNsub', 'CHsub','CN2', 'CH2',\
            'CNNH', 'CNNH2','CHNNH', 'CNNH3' , 'CHNNH2' , 'CHNNH3' , 'CH2NNH2' ,'CH2NNH3'] 

analist = ['CN', 'CH','CNH', 'CNH2', 'Cac', 'CNH3', 'CN2', 'CH2', 'CNsub', 'CHsub', 'CNNH', 'CHNNH', 'CH2NNH2']
header = analist + ['rDesNH3', 'rDes4f_NNH', 'rDes2f_HNNH', 'rDes3f_H2NNH2']


C0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
t0, tf, tsteps = 0, 120, 100
T = 400

start = time.time()
ReadExlfile, ReadExlsheet='MKMInputSheet.xlsx', 'EnergyInput'
writebook = xlsxwriter.Workbook('MKMoutput.xlsx')
writesheet = writebook.add_worksheet()
DRC_label_list = ['H1f','H1b','H2f','H2b','H3f','H3b']
#[0,1e-25,1e-20,1e-15,1e-10,1e-5,1e-1,1]
#[1e-25,1e-24,1e-23,1e-22,1e-21]
#[1e-20,1e-19,1e-18,1e-17,1e-16]
#[1e-15,1e-14,1e-13,1e-12,1e-11]
#[1e-10,1e-9,1e-8,1e-7,1e-6]
#[1e-5,1e-4,1e-3,1e-2,1e-1, 1]
#################################
#[1,1e-1,1e-2,1e-3,1e-4,1e-5]
#[1e-6,1e-7,1e-8,1e-9,1e-10]
#[1e-11,1e-12,1e-13,1e-14,1e-15]
#[1e-16,1e-17,1e-18,1e-19,1e-20]
#[1e-21,1e-22,1e-23,1e-24,1e-25]
#[1e-26,1e-27,1e-28,1e-29,1e-30]
#[1e-31,1e-32,1e-33,1e-34,1e-35]
#[1e-36,1e-37,1e-38,1e-39,1e-40]
#[1e-41,1e-42,1e-43,1e-44,1e-45]
#[1e-46,1e-47,1e-48,1e-49,1e-50]


radicalconc =[1]
stickingprobability = [1]
Totalsurfaceratio = 10


All_metal = ['Co','Ni','Pd','Ga','Sn','Cu','Au',"Ag"]
SelectedMetals = ['Co']
# InputMetals = np.atleast_1d(All_metal[-2:])

title_creater_SP_ConcR('ConcR',radicalconc,'StickingP',stickingprobability)

columncount = 0
    
for CountR, concR in enumerate(radicalconc): 
    # M = 'Au'
    # SP = 1
        
        
    for CountSP, SP in enumerate(stickingprobability): #main cycle
        # M = 'Au'
        # conc = 1
        columncount += 1
    
        for CountM, M in enumerate(SelectedMetals): 
            # SP = 1E-30
            # conc = 1
            
            TOFs_AllMetal = []
            
            # for Count, label  in enumerate(DRC_label_list):
            # time steps
            t = np.linspace(t0, tf, tsteps)
            #solve ODE
            Ads = 'comp'
            ODE = Reaction_ODE(T,metal = M,Exlfile=ReadExlfile,Exlsheet=ReadExlsheet,\
                                InitialAds=Ads,N2vib='on',AAds='on',ER1='on', SubDso ='on',
                                    LHW ='on', ER2='on', Dissociation='on', Desorption = 'on',
                                    StickingProbability =SP, ConcRatio=concR,
                                    DRC_label = None,TSR = Totalsurfaceratio)
            #Integral the ODE
            # Ci = odeint(ODE.InitializeModel,t=t,y0=C0)
            #atol and rtol minimum  2.22e-14, 
            #Radau rtol=1.ER29012e-10, atol=1.49012e-14
            #Sn  BDF rtol=1.49012e-8, atol=1.49012e-8
            
            if M not in ['Fe','Ga','Sn', 'Pd'] :
                
                rtolGeneral=1e-10
                atolGeneral=1e-12
                methodGeneral = 'BDF'
                
                print("\033[1;35m"+f'Solver {methodGeneral}'+'\033[39m')
                print(f'rtol | {rtolGeneral}')
                print(f'atol | {atolGeneral}')
                C = solve_ivp(ODE.InitializeModel,t_span = [t0, tf], t_eval=t, y0 = C0, method=methodGeneral, dense_output=True,\
                              rtol=rtolGeneral, atol=atolGeneral, max_step=200000\
                              ,first_step=1e-30)
                Ci = C.sol(t).T
           
            elif M in ['Pd'] :
                
               rtolPd=1e-9
               atolPd=1e-10
               methodPd ='BDF'
               
               print("\033[1;35m"+f'Solver {methodPd}'+'\033[39m')
               print(f'rtol | {rtolPd}')
               print(f'atol | {atolPd}')
               C = solve_ivp(ODE.InitializeModel,t_span = [t0, tf], t_eval=t, y0 = C0, method=methodPd, dense_output=True,\
                             rtol=rtolPd, atol=atolPd, max_step=200000\
                             ,first_step=1e-20) #rtol=1e-8, atol=1e-12 with BDF
               Ci = C.sol(t).T
                
            elif M in ['Ga'] :
                rtolGa=1e-8
                atolGa=1e-9
                methodGa = 'BDF'
                print("\033[1;35m"+f'Solver {methodGa}'+'\033[39m')
                print(f'rtol | {rtolGa}')
                print(f'atol | {atolGa}')
                C = solve_ivp(ODE.InitializeModel,t_span = [t0, tf], t_eval=t, y0 = C0, method=methodGa, dense_output=True,\
                              rtol=rtolGa, atol=atolGa, max_step=300000\
                              ,first_step=1e-35) #rtol=1e-8, atol=1e-12 with BDF
                Ci = C.sol(t).T
                
            elif M in ['Sn'] :
                rtolSn=1e-8
                atolSn=1e-9
                methodSn = 'BDF'
                print("\033[1;35m"+f'Solver {methodSn}'+'\033[39m')
                print(f'rtol | {rtolSn}')
                print(f'atol | {atolSn}')
                C = solve_ivp(ODE.InitializeModel,t_span = [t0, tf], t_eval=t, y0 = C0, method=methodSn, dense_output=True,\
                              rtol=rtolSn, atol=atolSn, max_step=200000\
                              ,first_step=1e-20)
                Ci = C.sol(t).T
                
            elif M in ['Fe'] :
                print("\033[1;35m"+f'Solver BDF for Fe'+'\033[39m')
                C = solve_ivp(ODE.InitializeModel,t_span = [t0, tf], t_eval=t, y0 = C0, method='BDF', dense_output=True,\
                              rtol=1e-8, atol=1e-8, max_step=200000\
                              ,first_step=1e-20) #rtol=1e-4, atol=8e-8 with BDF
                Ci = C.sol(t).T
                
            Cf, rformNH3, rDesNH3, rDesNNH, rDesHNNH, rDes3fH2NNH2 = ODE.AnalyzeData(Ci)
            makeplot(Cf,M)
            
        
            
            # TOFs_AllMetal.append(rNH3g_TOF[-1])
            # DataOut =  f'{Cf[-1,:]}|{rDesNH3[-1]}|{rDesNNH[-1]}|{rDesHNNH[-1]}|{rDes3fH2NNH2[-1]}'.replace("'",'').replace("[",'').replace("]",'').replace(" ",',').replace("\n",'')
            # print(M,' ', Conc) 
            # DataOut = f'{dCNH3gdt_des[-1]}|'.replace("'",'').replace("[",'').replace("]",'')
            # excel_output(Count, M, DataOut.replace("|",','))
            excel_output_SP_ConcR(CountM, columncount, M, rDesNH3[-1])
            # print(M, '|' ,{dCNH3gdt_dir[-1]})
            print(M, 'NH3 Rform|' ,{rformNH3})
            print(M, 'NH3  Rdes|' ,{rDesNH3[-1]})
            print('----------------------------------')
            
writebook.close()

print('Running', time.time()-start, 'seconds.')

Ci = pd.DataFrame(Ci)
Ci.columns = odelist







