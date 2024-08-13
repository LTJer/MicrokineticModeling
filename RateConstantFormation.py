
from RateConstantEquation import RateConstant_Fomula
from EnergyTool import  EaTransformation



class RateConstant(RateConstant_Fomula) :  
    
    
    def __init__ (self,Exlfile,Exlsheet,metal,T,SN2,SH2,SN,SH,SNH3,SHNNH,SN2H4,
                  SNH,SNNH, DRC_n = None):
  
        
        '''
        Default entropy is at 1 atm in J/(K*mol)
        The reaction label should follow the excel sheet
        '''
        self.T = T
        self.SN2 = SN2
        self.SH2 = SH2
        self.SNH3 = SNH3
        self.SHNNH = SHNNH
        self.SN2H4 = SN2H4
        self.SN = SN
        self.SH = SH
        self.SNH = SNH
        self.SNNH = SNNH
        self.Exlfile = Exlfile
        self.Exlsheet = Exlsheet
        self.metal = metal
        
        
    def k_ThermalSimple (self):
        
        label_list = ['S1f','S1b','S2f','S2b']
        
        ES1f,ES1b,ES2f,ES2b = EaTransformation.Ea_FetchFromExcel(self.Exlfile,\
             self.Exlsheet, label_list , self.metal)[:,1].astype(float)
        
        kS1f =  self.k_Ey(-self.SN2, ES1f)  
        kS1b =  self.k_Ey((self.SN2*2), ES1b)  
        kS2f =  self.k_Ey(-self.SH2, ES2f) 
        kS2b =  self.k_Ey((self.SH2*2), ES2b)  
            
        return kS1f,kS1b,kS2f,kS2b
            
    def k_ThermalComplex (self):
        
        label_list = ['C1f','C1b','C2f','C2b','C3f','C3b','C4f','C4b']
        
        EC1f,EC1b,EC2f,EC2b,EC3f,EC3b,EC4f,EC4b = \
            EaTransformation.Ea_FetchFromExcel(self.Exlfile,\
             self.Exlsheet, label_list , self.metal)[:,1].astype(float)
        kC1f =  self.k_Ey(-self.SN2, EC1f)  
        kC1b =  self.k_Ey((self.SN2*2), EC1b)  
        kC2f =  self.k_Ey(-self.SH2, EC2f) 
        kC2b =  self.k_Ey((self.SH2*2), EC2b)
        kC3f =  self.k_sr(EC3f)
        kC3b =  self.k_sr(EC3b)
        kC4f =  self.k_sr(EC4f)
        kC4b =  self.k_sr(EC4b)
        
        return kC1f,kC1b,kC2f,kC2b,kC3f,kC3b,kC4f,kC4b
    
    def k_ThermalHydrogenation (self):
        
        label_list = ['H1f','H1b','H2f','H2b','H3f','H3b']
        
        EH1f,EH1b,EH2f,EH2b,EH3f,EH3b = \
            EaTransformation.Ea_FetchFromExcel(self.Exlfile,\
            self.Exlsheet, label_list , self.metal)[:,1].astype(float)
            
        kH1f =  self.k_sr(EH1f)
        kH1b =  self.k_sr(EH1b)
        kH2f =  self.k_sr(EH2f)
        kH2b =  self.k_sr(EH2b)
        kH3f =  self.k_sr(EH3f) 
        kH3b =  self.k_sr(EH3b)
        # kH3b =  self.k_sr(EH3b)
        
        # NH3* ⟶  NH3(g) + *
        # kH4b =  self.k_Ey((0), EH4b)
        # kH4f =  self.k_Ey(-self.SNH3, EH4f)  
          
        return kH1f,kH1b,kH2f,kH2b,kH3f,kH3b#,kH4b,kH4f
        
    def k_ActiveAdsorptions (self):#EAd3f,EAd3b,EAd4f,EAd4b):
        
        label_list = ['Ad1f','Ad1b','Ad2f','Ad2b', 'Ad3f','Ad3b']
        
        EAd1f,EAd1b,EAd2f,EAd2b,EAd3f,EAd3b = EaTransformation.Ea_FetchFromExcel(self.Exlfile,\
             self.Exlsheet, label_list , self.metal)[:,1].astype(float)
        
        kAd1f =  self.k_Ey(-self.SN, EAd1f)  
        kAd1b =  self.k_Ey((self.SN*2), EAd1b)  
        kAd2f =  self.k_Ey(-self.SH, EAd2f) 
        kAd2b =  self.k_Ey((self.SH*2), EAd2b)
        kAd3f =  self.k_Ey(-self.SNH, EAd3f)  
        kAd3b =  self.k_Ey((self.SNH*2), EAd3b)  
        # print(-self.SN, EAd1f)
        # print(kAd1f)
        # kAd4f =  self.k_Ey(-self.SH, EAd4f) 
        # kAd4b =  self.k_Ey((0), EAd4b)
        
        return kAd1f,kAd1b,kAd2f,kAd2b,kAd3f,kAd3b#,kAd4f,kAd4b
    
    def k_DissociativeAds_vib(self):
        
        label_list = ['S1f','S1b','S2f','S2b']
        
        ES1f,ES1b,ES2f,ES2b = EaTransformation.Ea_FetchFromExcel(self.Exlfile,\
             self.Exlsheet, label_list , self.metal)[:,1].astype(float)
        
        ENv1f,ENv1b,ENv2f,ENv2b,ENv3f,ENv3b,ENv4f,ENv4b,ENv5f,ENv5b,ENv6f,ENv6b,ENv7f,ENv7b,ENv8f,ENv8b,ENv9f,ENv9b,ENv10f,ENv10b,\
        EHv1f,EHv1b,EHv2f,EHv2b,EHv3f,EHv3b,EHv4f,EHv4b,EHv5f,EHv5b,EHv6f,EHv6b,EHv7f,EHv7b,EHv8f,EHv8b,EHv9f,EHv9b,EHv10f,EHv10b\
        =EaTransformation.Ea_DissociativeAds_vib(ES1f,ES1b,ES2f,ES2b)
        
        kNv1f =  self.k_Ey(-self.SN2, ENv1f)
        kNv1b =  self.k_Ey((self.SN2*2), ENv1b)
        kNv2f =  self.k_Ey(-self.SN2, ENv2f)
        kNv2b =  self.k_Ey((self.SN2*2), ENv2b) 
        kNv3f =  self.k_Ey(-self.SN2, ENv3f)
        kNv3b =  self.k_Ey((self.SN2*2), ENv3b)
        kNv4f =  self.k_Ey(-self.SN2, ENv4f)
        kNv4b =  self.k_Ey((self.SN2*2), ENv4b)
        kNv5f =  self.k_Ey(-self.SN2, ENv5f)
        kNv5b =  self.k_Ey((self.SN2*2), ENv5b)
        kNv6f =  self.k_Ey(-self.SN2, ENv6f)
        kNv6b =  self.k_Ey((self.SN2*2), ENv6b)
        kNv7f =  self.k_Ey(-self.SN2, ENv7f)
        kNv7b =  self.k_Ey((self.SN2*2), ENv7b) 
        kNv8f =  self.k_Ey(-self.SN2, ENv8f)
        kNv8b =  self.k_Ey((self.SN2*2), ENv8b)
        kNv9f =  self.k_Ey(-self.SN2, ENv9f)
        kNv9b =  self.k_Ey((self.SN2*2), ENv9b)
        kNv10f =  self.k_Ey(-self.SN2, ENv10f)
        kNv10b =  self.k_Ey((self.SN2*2), ENv10b)
        
        
        kHv1f =  self.k_Ey(-self.SH2, EHv1f)
        kHv1b =  self.k_Ey((self.SH2*2), EHv1b)
        kHv2f =  self.k_Ey(-self.SH2, EHv2f)
        kHv2b =  self.k_Ey((self.SH2*2), EHv2b) 
        kHv3f =  self.k_Ey(-self.SH2, EHv3f)
        kHv3b =  self.k_Ey((self.SH2*2), EHv3b)
        kHv4f =  self.k_Ey(-self.SH2, EHv4f)
        kHv4b =  self.k_Ey((self.SH2*2), EHv4b)
        kHv5f =  self.k_Ey(-self.SH2, EHv5f)
        kHv5b =  self.k_Ey((self.SH2*2), EHv5b)
        kHv6f =  self.k_Ey(-self.SH2, EHv6f)
        kHv6b =  self.k_Ey((self.SH2*2), EHv6b)
        kHv7f =  self.k_Ey(-self.SH2, EHv7f)
        kHv7b =  self.k_Ey((self.SH2*2), EHv7b) 
        kHv8f =  self.k_Ey(-self.SH2, EHv8f)
        kHv8b =  self.k_Ey((self.SH2*2), EHv8b)
        kHv9f =  self.k_Ey(-self.SH2, EHv9f)
        kHv9b =  self.k_Ey((self.SH2*2), EHv9b)
        kHv10f =  self.k_Ey(-self.SH2, EHv10f)
        kHv10b =  self.k_Ey((self.SH2*2), EHv10b)
    
        return kNv1f,kNv1b,kNv2f,kNv2b,kNv3f,kNv3b,kNv4f,kNv4b,kNv5f,kNv5b,kNv6f,kNv6b,kNv7f,kNv7b,kNv8f,kNv8b,kNv9f,kNv9b,kNv10f,kNv10b,\
            kHv1f,kHv1b,kHv2f,kHv2b,kHv3f,kHv3b,kHv4f,kHv4b,kHv5f,kHv5b,kHv6f,kHv6b,kHv7f,kHv7b,kHv8f,kHv8b,kHv9f,kHv9b,kHv10f,kHv10b
    
    
   
    
    def k_SubsurfaceDissolution(self):
        
        label_list = ['SD1f','SD1b','SD2f','SD2b']
        
        ESD1f,ESD1b,ESD2f,ESD2b\
            = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)
        
        kSD1f = self.k_Ey(0, ESD1f)
        kSD1b = self.k_Ey(0, ESD1b)
        kSD2f = self.k_Ey(0, ESD2f)
        kSD2b = self.k_Ey(0, ESD2b)
      
        return kSD1f,kSD1b,kSD2f,kSD2b
    
    
    def k_Desorption(self):
        
        label_list = ['Des1f','Des1b','Des2f','Des2b','Des3f','Des3b','Des4f','Des4b']
        
        EDes1f,EDes1b,EDes2f,EDes2b,EDes3f,EDes3b,EDes4f,EDes4b\
            = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)
            
        kDes1f = self.k_Ey(self.SNH3, EDes1f) # des
        kDes1b = self.k_Ey(-self.SNH3/2, EDes1b) # ads
        kDes2f = self.k_Ey(self.SHNNH, EDes2f)
        kDes2b = self.k_Ey(-self.SHNNH/2, EDes2b)
        kDes3f = self.k_Ey(self.SN2H4, EDes3f)
        kDes3b = self.k_Ey(-self.SN2H4/2, EDes3b)
        kDes4f = self.k_Ey(self.SNNH, EDes4f)
        kDes4b = self.k_Ey(-self.SNNH/2, EDes4b)
        # print(self.SNNH,EDes4f)
        return kDes1f,kDes1b,kDes2f,kDes2b,kDes3f,kDes3b,kDes4f,kDes4b
    
    
    def k_LangmuirHinshelwood(self):
        
        label_list = ['LHW1f','LHW1b','LHW2f','LHW2b','LHW3f','LHW3b','LHW4f','LHW4b','LHW5f','LHW5b','LHW6f'\
            ,'LHW6b','LHW7f','LHW7b','LHW8f','LHW8b','LHW9f','LHW9b','LHW10f','LHW10b','LHW11f','LHW11b'\
            ,'LHW12f','LHW12b']
        
        ELHW1f,ELHW1b,ELHW2f,ELHW2b,ELHW3f,ELHW3b,ELHW4f,ELHW4b,ELHW5f,ELHW5b,ELHW6f\
            ,ELHW6b,ELHW7f,ELHW7b,ELHW8f,ELHW8b,ELHW9f,ELHW9b,ELHW10f,ELHW10b,ELHW11f,ELHW11b\
            ,ELHW12f,ELHW12b\
            = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)
     
        kLHW1f = self.k_Ey(0, ELHW1f)
        kLHW1b = self.k_Ey(0, ELHW1b)
        kLHW2f = self.k_Ey(0, ELHW2f)
        kLHW2b = self.k_Ey(0, ELHW2b)
        kLHW3f = self.k_Ey(0, ELHW3f)
        kLHW3b = self.k_Ey(0, ELHW3b)
        kLHW4f = self.k_Ey(0, ELHW4f)
        kLHW4b = self.k_Ey(0, ELHW4b)
        kLHW5f = self.k_Ey(0, ELHW5f)
        kLHW5b = self.k_Ey(0, ELHW5b)
        kLHW6f = self.k_Ey(0, ELHW6f)
        kLHW6b = self.k_Ey(0, ELHW6b)
        kLHW7f = self.k_Ey(0, ELHW7f)
        kLHW7b = self.k_Ey(0, ELHW7b)
        kLHW8f = self.k_Ey(0, ELHW8f)
        kLHW8b = self.k_Ey(0, ELHW8b)
        kLHW9f = self.k_Ey(0, ELHW9f)
        kLHW9b = self.k_Ey(0, ELHW9b)
        kLHW10f = self.k_Ey(0, ELHW10f)
        kLHW10b = self.k_Ey(0, ELHW10b)
        kLHW11f = self.k_Ey(0, ELHW11f)
        kLHW11b = self.k_Ey(0, ELHW11b)
        kLHW12f = self.k_Ey(0, ELHW12f)
        kLHW12b = self.k_Ey(0, ELHW12b)
       
        return kLHW1f,kLHW1b,kLHW2f,kLHW2b,kLHW3f,kLHW3b,kLHW4f,kLHW4b,kLHW5f,kLHW5b,kLHW6f\
            ,kLHW6b,kLHW7f,kLHW7b,kLHW8f,kLHW8b,kLHW9f,kLHW9b,kLHW10f,kLHW10b,kLHW11f,kLHW11b\
            ,kLHW12f,kLHW12b
            
    def k_ER_Reactions1(self):
        
        label_list = ['ER11f','ER11b','ER12f','ER12b','ER13f','ER13b','ER14f','ER14b','ER15f','ER15b','ER16f','ER16b']
        
        EER11f,EER11b,EER12f,EER12b,EER13f,EER13b,EER14f,EER14b,EER15f,EER15b,EER16f,EER16b \
            = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)

        self.delSER16fnew = self.SH2 - self.SH
        
        kER11f = self.k_Ey(-self.SN, EER11f)   
        kER11b = self.k_Ey((self.SN*2), EER11b)  
        kER12f = self.k_Ey(-self.SH, EER12f)  
        kER12b = self.k_Ey((self.SH*2), EER12b)  
        kER13f = self.k_Ey(-self.SH, EER13f)  
        kER13b = self.k_Ey((self.SH*2), EER13b)  
        kER14f = self.k_Ey(-self.SH, EER14f)  
        kER14b = self.k_Ey((self.SH*2), EER14b)  
        kER15f = self.k_Ey(-self.SN, EER15f)  
        kER15b = self.k_Ey(self.SN*2, EER15b)  
        kER16f = self.k_Ey(-self.SH, EER16f)  
        kER16b = self.k_Ey(self.SH*2, EER16b)
        # kER16f = self.k_Ey(self.delSER16fnew, EER16f)  
        # kER16b = self.k_Ey(-self.delSER16fnew*2, EER16b)
        
        
        return kER11f,kER11b,kER12f,kER12b,kER13f,kER13b,kER14f,kER14b,kER15f,kER15b,kER16f,kER16b
    
            
    def k_ER_Reactions2(self):
        
            label_list = ['ER21f','ER21b','ER22f','ER22b','ER23f','ER23b','ER24f','ER24b','ER25f','ER25b',\
                        'ER26f','ER26b','ER27f','ER27b','ER28f','ER28b','ER29f','ER29b','ER30f','ER30b',\
                            'ER31f','ER31b','ER32f','ER32b','ER33f','ER33b','ER34f','ER34b','ER35f','ER35b']
       
            EER21f,EER21b,EER22f,EER22b,EER23f,EER23b,EER24f,EER24b,EER25f,EER25b,\
                        EER26f,EER26b,EER27f,EER27b,EER28f,EER28b,EER29f,EER29b,EER30f,EER30b,\
                            EER31f,EER31b,EER32f,EER32b,EER33f,EER33b,EER34f,EER34b,EER35f,EER35b\
                = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)
    
            
            kER21f = self.k_Ey(-self.SH, EER21f)   
            kER21b = self.k_Ey((self.SH*2), EER21b)  
            kER22f = self.k_Ey(-self.SH, EER22f)  
            kER22b = self.k_Ey((self.SH*2), EER22b)  
            kER23f = self.k_Ey(-self.SH, EER23f)  
            kER23b = self.k_Ey((self.SH*2), EER23b)  
            kER24f = self.k_Ey(-self.SH, EER24f)  
            kER24b = self.k_Ey((self.SH*2), EER24b)  
            kER25f = self.k_Ey(-self.SH, EER25f)  
            kER25b = self.k_Ey((self.SH*2), EER25b)  
            kER26f = self.k_Ey(-self.SH, EER26f)  
            kER26b = self.k_Ey((self.SH*2), EER26b)
            kER27f = self.k_Ey(-self.SH, EER27f)  
            kER27b = self.k_Ey((self.SH*2), EER27b)
            kER28f = self.k_Ey(-self.SH, EER28f)  
            kER28b = self.k_Ey((self.SH*2), EER28b)
            kER29f = self.k_Ey(-self.SH, EER29f)  
            kER29b = self.k_Ey((self.SH*2), EER29b)
            kER30f = self.k_Ey(-self.SH, EER30f)  
            kER30b = self.k_Ey((self.SH*2), EER30b)
            kER31f = self.k_Ey(-self.SH, EER31f)  
            kER31b = self.k_Ey((self.SH*2), EER31b)
            kER32f = self.k_Ey(-self.SH, EER32f)  
            kER32b = self.k_Ey((self.SH*2), EER32b)
            kER33f = self.k_Ey(-self.SN, EER33f)  
            kER33b = self.k_Ey((self.SN*2), EER33b)
            kER34f = self.k_Ey(-self.SN, EER34f)  
            kER34b = self.k_Ey((self.SN*2), EER34b)
            kER35f = self.k_Ey(-self.SN2, EER35f)  
            kER35b = self.k_Ey((self.SN2*2), EER35b)
            
            
            return kER21f,kER21b,kER22f,kER22b,kER23f,kER23b,kER24f,kER24b,kER25f,kER25b,\
                        kER26f,kER26b,kER27f,kER27b,kER28f,kER28b,kER29f,kER29b,kER30f,kER30b,\
                            kER31f,kER31b,kER32f,kER32b,kER33f,kER33b,kER34f,kER34b,kER35f,kER35b\
        
        
    def k_SurfaceDissociation(self):
            
            label_list = ['DISO1f','DISO1b','DISO2f','DISO2b','DISO3f','DISO3b'\
            ,'DISO4f','DISO4b','DISO5f','DISO5b','DISO6f','DISO6b'\
                ,'DISO7f','DISO7b','DISO8f','DISO8b']
       
            EDISO1f,EDISO1b,EDISO2f,EDISO2b,EDISO3f,EDISO3b,EDISO4f,EDISO4b,EDISO5f,EDISO5b,EDISO6f,EDISO6b\
                ,EDISO7f,EDISO7b,EDISO8f,EDISO8b\
                = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)
    
            
            kDISO1f = self.k_Ey(0, EDISO1f)   
            kDISO1b = self.k_Ey((0), EDISO1b)  
            kDISO2f = self.k_Ey(0, EDISO2f)  
            kDISO2b = self.k_Ey((0),EDISO2b)  
            kDISO3f = self.k_Ey(0, EDISO3f)  
            kDISO3b = self.k_Ey((0), EDISO3b)  
            kDISO4f = self.k_Ey(0, EDISO4f)  
            kDISO4b = self.k_Ey((0), EDISO4b)  
            kDISO5f = self.k_Ey(0, EDISO5f)  
            kDISO5b = self.k_Ey((0), EDISO5b)  
            kDISO6f = self.k_Ey(0, EDISO6f)  
            kDISO6b = self.k_Ey((0), EDISO6b)
            kDISO7f = self.k_Ey(0, EDISO7f)  
            kDISO7b = self.k_Ey((0), EDISO7b)
            kDISO8f = self.k_Ey(0, EDISO8f)  
            kDISO8b = self.k_Ey((0), EDISO8b)
            
            
            
            return kDISO1f,kDISO1b,kDISO2f,kDISO2b,kDISO3f,kDISO3b\
            ,kDISO4f,kDISO4b,kDISO5f,kDISO5b,kDISO6f,kDISO6b\
                ,kDISO7f,kDISO7b,kDISO8f,kDISO8b
        
        
    def k_NHradical(self):
                
                label_list = ['NHradical1f', 'NHradical1b','NHradical2f','NHradical2b']
           
                ENHradical1f,ENHradical1b,ENHradical2f,ENHradical2b,\
                    = EaTransformation.Ea_FetchFromExcel(self.Exlfile,self.Exlsheet,label_list,self.metal)[:,1].astype(float)
        
                
                kNHrad1f = self.k_Ey(-self.SNH, ENHradical1f)   
                kNHrad1b = self.k_Ey((self.SNH*2), ENHradical1b)  
                kNHrad2f = self.k_Ey(-self.SNH, ENHradical2f)  
                kNHrad2b = self.k_Ey((self.SNH*2),ENHradical2b)  
       
                         
                
                return kNHrad1f,kNHrad1b,kNHrad2f,kNHrad2b