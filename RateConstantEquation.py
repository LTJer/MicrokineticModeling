import numpy as np


class RateConstant_Fomula :
    
    def __init__(self, T) :
        
        self.T = T

    
    def k_Ey(self, TS_S, TS_H): #Eyring equation
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
            
        TS_G = TS_H - self.T*TS_S
        
        if TS_G < 0:
            TS_G = 0
            
        k =  kb * self.T / h  * np.exp(-TS_G / (R * self.T)) 
        ksci  = np.float64("{:.1e}".format(k))
        # k2 = np.format_float_scientific(k, exp_digits=1)
        
        return ksci
           
    def k_sr(self,Eact):
        """
        Calculate reaction rate constant for a surface reaction in s-1
    
        T       - Temperature in K
        nu      - Pre-exponential factor in s^-1
        Eact    - Activation energy in J/mol
        """
        R = 8.3144598 # gas constant
        kb = 1.38064852e-23 # boltzmann constant
        h = 6.62607004e-34  # planck constant
        k =  kb * self.T / h * np.exp(-Eact / (R * self.T))
        ksci  = np.float64("{:.1e}".format(k))
        return ksci
    
    def k_ads(self, P, A, m):
        """
        Reaction rate constant for adsorption in s-1
        
        T           - Temperature in K
        P           - Pressure in Pa
        A           - Surface area in m^2
        m           - Mass of reactant in kg
        """
        kb = 1.38064852e-23 # boltzmann constant
        return P*A / np.sqrt(2 * np.pi * m * kb * self.T)
    
    def k_des_di(self, A, m, sigma, I, Edes):
        """
        Reaction rate constant for desorption with diatomic molecule in s-1
        
        T           - Temperature in K
        A           - Surface area in m^2
        m           - Mass of reactant in kg
        sigma       - Symmetry number
        theta_rot   - Rotational temperature in K
        Edes        - Desorption energy in J/mol
        I           - Momentum ineria 
        """
        kb = 1.38064852e-23 # boltzmann constant
        h = 6.62607004e-34  # planck constant
        R = 8.3144598       # gas constant
        theta_rot = h**2 / (8 * np.pi**2 * kb * I)
        
        return kb * self.T**3 / h**3 * A * (2 * np.pi * m * kb) / \
            (sigma * theta_rot) * np.exp(-Edes / (R*self.T))
            
    def k_des_poly(self, T, A, m, sigma, IA, IB, IC, Edes):
        """
        Reaction rate constant for desorption with polyatomic molecule in s-1
        
        T           - Temperature in K
        A           - Surface area in m^2
        m           - Mass of reactant in kg
        sigma       - Symmetry number
        IA,IB,IC    - Moment of inertia of the molecule in kg*m^2
        Edes        - Desorption energy in J/mol
        """
        kb = 1.38064852e-23     # boltzmann constant
        h = 6.62607004e-34      # planck constant
        R = 8.3144598           # gas constant
        c = 3e8                 #light speed
        B = h / (8 * np.pi**2 * c) #parameter to calculate Rotational constant
        RA = B / IA
        RB = B / IB
        RC = B / IC
        return kb**2.5 * self.T**3.5 / h**4.5 * A * (2 * np.pi * m * kb) / \
            (sigma * c**1.5) * np.pi**0.5 / (RA * RB * RC)**0.5 * np.exp(-Edes / (R*self.T))
            