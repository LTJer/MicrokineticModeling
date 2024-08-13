import numpy as np
from XlrdTool import WordFetchValue

class EaTransformation :
    
    def __init__(self) :
        
        pass
        
         
    def Ea_DissociativeAds_vib (E1f,E1b,E2f,E2b):
        
        #Barriers for N2 vibrational excited state
        ENvib = 0.288 * 96.48 * 1000 #Nitrogen vibrational excitated energy gap 
        EHvib = 0.561 * 96.48 * 1000 #Hydrogen vibrational excitated energy gap 
        #in j/mole
        #Nitrogen vibrational excited Ea
        ENv1f = E1f - ((1) * ENvib)
        if ENv1f < 0: ENv1f = 0
        ENv1b = E1b
        
        ENv2f = E1f - ((2) * ENvib)
        if ENv2f < 0: ENv2f = 0
        ENv2b = E1b
        
        ENv3f = E1f - ((3) * ENvib)
        if ENv3f < 0: ENv3f = 0
        ENv3b = E1b 
        
        ENv4f = E1f - ((4) * ENvib)
        if ENv4f < 0: ENv4f = 0
        ENv4b = E1b 
        
        ENv5f = E1f - ((5) * ENvib)
        if ENv5f < 0: ENv5f = 0
        ENv5b = E1b 
        
        ENv6f = E1f - ((6) * ENvib)
        if ENv6f < 0: ENv6f = 0
        ENv6b = E1b
        
        ENv7f = E1f - ((7) * ENvib)
        if ENv7f < 0: ENv7f = 0
        ENv7b = E1b
        
        ENv8f = E1f - ((8) * ENvib)
        if ENv8f < 0: ENv8f = 0
        ENv8b = E1b 
        
        ENv9f = E1f - ((9) * ENvib)
        if ENv9f < 0: ENv9f = 0
        ENv9b = E1b 
        
        ENv10f = E1f - ((10) * ENvib)
        if ENv10f < 0: ENv10f = 0
        ENv10b = E1b 
        
        #Hydrogen vibrational excited Ea
        EHv1f = E2f - ((1) * EHvib)
        if EHv1f < 0: EHv1f = 0
        EHv1b = E2b
        
        EHv2f = E2f - ((2) * EHvib)
        if EHv2f < 0: EHv2f = 0
        EHv2b = E2b
        
        EHv3f = E2f - ((3) * EHvib)
        if EHv3f < 0: EHv3f = 0
        EHv3b = E2b
        
        EHv4f = E2f - ((4) * EHvib)
        if EHv4f < 0: EHv4f = 0
        EHv4b = E2b
        
        EHv5f = E2f - ((5) * EHvib)
        if EHv5f < 0: EHv5f = 0
        EHv5b = E2b
        
        EHv6f = E2f - ((6) * EHvib)
        if EHv6f < 0: EHv6f = 0
        EHv6b = E2b
        
        EHv7f = E2f - ((7) * EHvib)
        if EHv7f < 0: EHv7f = 0
        EHv7b = E2b
        
        EHv8f = E2f - ((8) * EHvib)
        if EHv8f < 0: EHv8f = 0
        EHv8b = E2b
        
        EHv9f = E2f - ((9) * EHvib)
        if EHv9f < 0: EHv9f = 0
        EHv9b = E2b
        
        EHv10f = E2f - ((10) * EHvib)
        if EHv10f < 0: EHv10f = 0
        EHv10b = E2b
        
        return ENv1f,ENv1b,ENv2f,ENv2b,ENv3f,ENv3b,ENv4f,ENv4b,ENv5f,ENv5b,ENv6f,ENv6b,ENv7f,ENv7b,ENv8f,ENv8b,ENv9f,ENv9b,ENv10f,ENv10b,\
            EHv1f,EHv1b,EHv2f,EHv2b,EHv3f,EHv3b,EHv4f,EHv4b,EHv5f,EHv5b,EHv6f,EHv6b,EHv7f,EHv7b,EHv8f,EHv8b,EHv9f,EHv9b,EHv10f,EHv10b
    
    
    
    def Ea_FetchFromExcel (excelfile, sheet, labellist, metal):
        EaList = np.empty([0,2])
        
        for R_Label in labellist:
            Ea = WordFetchValue(excelfile, sheet, R_Label, metal)
            EaList = np.append(EaList,[[R_Label,Ea]], axis = 0)
            
        EaList[EaList == ''] = 0
        return EaList      