import numpy as np

'''
calculte partition function and entropy 
assumption:
    1. ideal gas
    2. non-interacting particle
    3. High excited states are inaccessible
    4. N = 1
'''
#Entropy from partition function 

#deinfe constants 

R = 8.314 # J·K−1·mol−1
Sexp = R
#N2 gas entropy
#Constants input
mN2 = 28 * 1.66054e-27  #kg
mH2 = 2 * 1.66054e-27
mN = 14 * 1.66054e-27
mH = 1.66054e-27
mNH = mN + mH
mNNH = mN + mN + mH



Symmetry_Num_N2 = 2
Symmetry_Num_H2 = 2
Symmetry_Num_NH = 1
Symmetry_Num_NNH = 1

d_N2 = 1.09e-10 # m
d_H2 = 0.74e-10 # m

I_N2 = 2 * mN*(0.5*d_N2)**2  # kg m^2
I_H2 = 2 * mH*(0.5*d_H2)**2  # kg m^2
I_NH = 1.679e-47

# IA_NNH = ((mN*(0.633**2)) + (mN*(0.551**2)) + (mN*(1.481**2)))*(1e-10**2)
# IB_NNH = ((mN*(0.032**2)) + (mN*(0.032**2)) + (mN*(0.898**2)))*(1e-10**2)
# IB_NNH = ((mN*(0.633**2)) + (mN*(0.551**2)) + (mN*(1.439**2)))*(1e-10**2)

N2_VF = 2744  # cm-1
H2_VF = 4342 # cm-1
NH_VF = 3177 # cm-1


def NistEntropy(A,B,C,D,E,F,G,T):
   '''
     unit J * mol-1 * K-1 
   '''
   t = T/1000
   S = (A*np.log(t)) + (B*t) + (C*(t**2)/2) + (D*(t**3)/3) - (E/(2*(t**2))) + G
   
   return S

def Entropy2Ev(Entropy, T):
    '''
    Entropy input with unit [J/(mol*K)]
    '''
    JperMole2eV = 1.036e-5
    Entropy_eV = Entropy*T*JperMole2eV
    
    return Entropy_eV
def Translation(m, T, P):
    '''
    P is pressure in pa
    Three degree of freedoms
    '''
    kb = 1.38E-23 # J K-1
    h = 6.626e-34 # m2 kg / s  or J S
    R = 8.314 # J·K−1·mol−1
    
    Qt = (((2*np.pi*m*kb*T)/(h ** 2)) ** 1.5) * (kb*T/P)
    St = R*(np.log(Qt)+1.5)
    
    return Qt, St


def Rotation_Linear(I, sigma_r, T):
    '''
    Two degree of freedoms
    '''
    kb = 1.38E-23 # J K-1
    h = 6.626e-34 # m2 kg / s  or J S
    R = 8.314 # J·K−1·mol−1
    
    Theta_r = (h**2)/(8*(np.pi**2)*I*kb)
    Qr = (T/Theta_r)/sigma_r
    Sr = R*(1+np.log(Qr))
    
    return Qr, Sr

def Rotation_NonLinear(IA,IB,IC, sigma_r, T):
    '''
    Three degree of freedoms
    '''
    kb = 1.38E-23 # J K-1
    h = 6.626e-34 # m2 kg / s  or J S
    R = 8.314 # J·K−1·mol−1
    
    Theta_A = (h**2)/(8*(np.pi**2)*IA*kb)
    Theta_B = (h**2)/(8*(np.pi**2)*IB*kb)
    Theta_C = (h**2)/(8*(np.pi**2)*IC*kb)
    Qr = (((np.pi**0.5)*(T**1.5))/sigma_r)*(1/((Theta_A*Theta_B*Theta_C)**0.5))

    Sr = R*(1.5+np.log(Qr))
    
    return Qr, Sr
    
    
    

def Single_Vibration(v, T):
    '''
    ignore imaginary frequency 
    v is vibrational frequency with unit cm-1
    '''
    kb = 1.38E-23 # J K-1
    h = 6.626e-34 # m2 kg / s  or J S
    R = 8.314 # J·K−1·mol−1
    
    mu = v * 2.997E+10
    Theta_v = h*mu/kb
    T_ratio = Theta_v / T
    Qv_1mode = np.exp(-T_ratio/2)/(1-np.exp(-T_ratio))
    Sv_1mode = R*((T_ratio / (np.exp(T_ratio)-1))-np.log(1-np.exp(-T_ratio)))
    
    return Qv_1mode, Sv_1mode


def ElectronicMotion(w, T):
    '''
    w is the degeneracy of the energy level
    E is the energy of the level

    assumption:
    1. First electronic state is much greater than kb * T
    2. the energy of ground state is set to zero
    '''
    kb = 1.38E-23 # J K-1
    h = 6.626e-34 # m2 kg / s  or J S
    R = 8.314 # J·K−1·mol−1
    
    E = 0
    Qe = w * np.exp(-E/(kb*T))
    Se = R * np.log(Qe)
    
    return Qe, Se
    

def Entropy_Total(St, Sr, Sv, Se):
    '''
    Sum up entropies from 4 different partition function 
    in unit J * mol-1 * K-1 
    '''
    R = 8.314 # J·K−1·mol−1
    Sexp = R

    Stotal = Sexp + St + Sr + Sv + Se 
    
    return Stotal
     

def LatticeGasModel(VF_list, T):
    '''
    S [=] J * mol-1 * K-1
    Diatomic Gases 
    Treat all motions on the surface as stiff vibrational motions
    All partition functions on the surface are 1 or provide with vibrational
    frequency
    
    for diatomic molecule, six motions are vibrational motions
    '''
    Qe, Se = ElectronicMotion(1, T)
    Stot_Lattice = Sexp #+ Se
    S_list = []
    Q_list = []
    for VF in VF_list:
        if VF < 100: VF = 100 #VF=200 gives Qv = 1 at 298K
        Qv_1mode, Sv_1mode = Single_Vibration(VF, T)
        # print(Qv_1mode, Sv_1mode)
        Stot_Lattice += Sv_1mode
        S_list.append(Sv_1mode)
        Q_list.append(Qv_1mode)
    return Q_list, S_list, Stot_Lattice
    
def TwoDGasModel(VF_list_2D, sigma_r, m, I, P, T):
    
    '''
    Diatomic Gases (linear)
    The gas on the surface can x and y directions' translation,
    helicopter rotation, and vibration of z direction, car wheel, and normal mode'
    2 translation + 1 rotation + 3 vibration 
    '''
    
    Qt, St = Translation(m, T, P)
    Qr, Sr = Rotation_Linear(I, sigma_r, T)
    Qe, Se = ElectronicMotion(2, T)
    
    Stot_2D = Sexp + Se + ((2/3) * St) + ((0.5 * Sr) - 0.5*np.log(sigma_r)) 
    print(Stot_2D)
    Sv_list = []
    Qv_list = []
    for VF in VF_list_2D:
        if VF < 100: VF = 100 #VF=200 gives Qv = 1 at 298K
        Qv_1mode, Sv_1mode = Single_Vibration(VF, T)
        # print(Qv_1mode, Sv_1mode)
        Stot_2D += Sv_1mode
        Sv_list.append(Sv_1mode)
        Qv_list.append(Qv_1mode)
        
    
    
    return Qv_list, Sv_list, Stot_2D


def S_Natom_gas(T, P) :

    '''
    N atom
    T = temperature in K 
    P = pressure in Pa
    Ignore electronic contribution by w = 1
    '''
    mN = 14 * 1.66054e-27 #kg
    
    NQt, NSt = Translation(mN,T , P)
    
    NQe, NSe = ElectronicMotion(1, T)
    
    Stotal_N = NSt + Sexp + NSe
    
    print(f'Stotal_N = {Stotal_N}')
    
    return Stotal_N

def S_Hatom_gas(T, P) :
    '''
    H atom
    T = temperature in K 
    P = pressure in Pa
    Ignore electronic contribution by w = 1
    '''
    mH = 1.66054e-27 #kg
    
    HQt, HSt = Translation(mH, T, P) 
    
    HQe, HSe = ElectronicMotion(1, T)
    
    Stotal_H = HSt + Sexp + HSe
 
    print(f'Stotal_H = {Stotal_H}')
    
    return Stotal_H

def S_NH_gas_TransOnly(T, P) :
    '''
    NH(g) Translational Entropies
    '''
    mN = 14 * 1.66054e-27
    mH = 1.66054e-27
    mNH = mN + mH

    NHQt, NHSt = Translation(mNH,T, P)
    
    NHQe, NHSe = ElectronicMotion(1, T)
    
    Stotal_NH = NHSt + Sexp + NHSe
    
    print(f'Stotal_NH = {Stotal_NH}')
    print(NHSt, Sexp , NHSe)
    return Stotal_NH

def NNHgasEntropy(T,P):
    '''
    NNH
    P is pressure in Pa
    '''
    mN = 14 * 1.66054e-27
    mH = 1.66054e-27
    mNNH = mN + mN + mH
    sigma_r = 1
    vf_list = [2629.53, 1869.69, 1076.38] #cm-1
    IA_NNH = ((mN*(0.633**2)) + (mN*(0.551**2)) + (mN*(1.481**2)))*(1e-10**2)
    IB_NNH = ((mN*(0.032**2)) + (mN*(0.032**2)) + (mN*(0.898**2)))*(1e-10**2)
    IC_NNH = ((mN*(0.633**2)) + (mN*(0.551**2)) + (mN*(1.439**2)))*(1e-10**2)
    NNHQt, NNHSt = Translation(mNNH, T, P) 
    NNHQr, NNHSr = Rotation_NonLinear(IA_NNH,IB_NNH,IC_NNH, sigma_r, T)
    NNHQe, NNHSe = ElectronicMotion(1, T)
    
    NNHSv_sum = 0
    for vf in vf_list:   
        NNHQv, NNHSv = Single_Vibration(vf, T)
        NNHSv_sum += NNHSv

    Stotal_NNH = Entropy_Total(NNHSt, NNHSr, NNHSv_sum, NNHSe)
    
    return Stotal_NNH
    


def Entropy_N2H4(T):
    #1 bar
    Entropy_N2H4 = NistEntropy(35.1824, 96.0526,-40.5013,6.66807,-0.874233,77.9915,249.425,T)
    return Entropy_N2H4

def Entropy_HNNH_Z(T):
    #1 bar
    Entropy_HNNH_Z =  NistEntropy(8.899662,93.79691,-53.12132,12.11607,0.362722,207.8059,205.7009,T)
    return Entropy_HNNH_Z


def Entropy_N2(T):
    #1 bar
    if T <= 500:
        Entropy_N2 = NistEntropy(28.98641,1.853978,-9.647459,16.63537,0.000117,-8.671914,226.4168,T)
    if T > 500:
        Entropy_N2 = NistEntropy(19.50583,19.88705,-8.598535,1.369784,0.527601,-4.935202,212.39,T)
        
    return Entropy_N2

def Entropy_H2(T):
    #1 bar
    Entropy_H2 =  NistEntropy(33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,T)
    
    return Entropy_H2

def Entropy_NH3(T):
    #1 bar
    Entropy_NH3 =  NistEntropy(19.99563,49.77119,-15.37599,1.921168,0.189174,-53.30667,203.8591,T)
    
    return Entropy_NH3


def Entropy_Natom(T):
    #1 bar
    Entropy_Natom =  NistEntropy(21.13581,-0.388842,0.043545,0.024685,-0.025678,466.311,178.8263,T)
    
    return Entropy_Natom


def Entropy_Hatom(T):
    #1 bar
    Entropy_Hatom =  NistEntropy(20.78603,4.85E-10,-1.58292E-10,1.53E-11,3.20E-11,211.802,139.8711,T)
    
    return Entropy_Hatom


def Entropy_NH(T):
    #1 bar Imidogen
    if T <= 800:
        Entropy_NH = NistEntropy(31.31766,-8.931544,11.58629,-2.544099,-0.041861,367.382,221.0867,T)
    if T > 800:
        Entropy_NH = NistEntropy(28.0592,5.10084,-0.767827,0.056843,-1.05472,366.0623,211.8351,T)
        
    return Entropy_NH

def Entropy_N2_partition(T, PN2):
    '''
    N2
    T = temperature in K
    PN2 = total pressure in pa
    '''
    mN2 = 28 * 1.66054e-27  #kg
    mH2 = 2 * 1.66054e-27
    mN = 14 * 1.66054e-27
    mH = 1.66054e-27

    Symmetry_Num_N2 = 2
    Symmetry_Num_H2 = 2

    d_N2 = 1.09e-10 # m
    d_H2 = 0.74e-10 # m

    I_N2 = 2 * mN*(0.5*d_N2)**2  # kg m^2
    I_H2 = 2 * mH*(0.5*d_H2)**2  # kg m^2

    N2_VF = 2744  # cm-1
    H2_VF = 4342 # cm-1
    NH_VF = 3177 # cm-1
    
    N2Qt, N2St = Translation(mN2, T, PN2) 
    N2Qr, N2Sr = Rotation_Linear(I_N2, Symmetry_Num_N2, T)
    N2Qv, N2Sv = Single_Vibration(N2_VF, T)
    N2Qe, N2Se = ElectronicMotion(1, T)

    # print(f'N2 Tranlational  Entropy : {N2St}')   
    # print(f'N2 Rotational  Entropy   : {N2Sr}')
    # print(f'N2 Vibrational Entropy   : {N2Sv}')
    # print(f'N2 Electronic  Entropy   : {N2Se}')

    Stotal_N2 = Entropy_Total(N2St, N2Sr, N2Sv, N2Se)

    
    return Stotal_N2


def Entropy_H2_partition(T, PH2):
    '''
    N2
    T = temperature in K
    PH2 = total pressure in pa
    '''
    mN2 = 28 * 1.66054e-27  #kg
    mH2 = 2 * 1.66054e-27
    mN = 14 * 1.66054e-27
    mH = 1.66054e-27

    Symmetry_Num_N2 = 2
    Symmetry_Num_H2 = 2

    d_N2 = 1.09e-10 # m
    d_H2 = 0.74e-10 # m

    I_N2 = 2 * mN*(0.5*d_N2)**2  # kg m^2
    I_H2 = 2 * mH*(0.5*d_H2)**2  # kg m^2

    N2_VF = 2744  # cm-1
    H2_VF = 4342 # cm-1
    NH_VF = 3177 # cm-1
    
    H2Qt, H2St = Translation(mH2, T, PH2) 
    H2Qr, H2Sr = Rotation_Linear(I_H2, Symmetry_Num_H2, T)
    H2Qv, H2Sv = Single_Vibration(H2_VF, T)
    H2Qe, H2Se = ElectronicMotion(1, T)

    Stotal_H2 = Entropy_Total(H2St, H2Sr, H2Sv, H2Se)  
    return Stotal_H2

################################################################################
# T = 300
# comparasion entropies between partition function and NIST table 
# print(Entropy_N2_partition(T, 1e+5))
# print(Entropy_N2(T))

# print(Entropy_H2_partition(T, 1e+5))
# print(Entropy_H2(T))




## H_atom
# print(NistEntropy(20.78603, 4.850638e-10, -1.582916e-10, 1.525102e-11,3.196347e-11,211.802,139.8711, 400))
# # ## N_atom
# print(NistEntropy(21.13581, -0.388842, 0.043545, 0.024685,-0.025678,466.3110,178.8263, 400))
# ## H2
# print(NistEntropy(33.066178, -11.363417, 11.432816,-2.772874,-0.158558,-9.980797,172.707974, 400))
## N2
# print(NistEntropy(28.98641, 1.853978, -9.647459,16.63537,0.000117,-8.671914,226.4168, 400))




####################calculate

#########Gas phase###############

# '''
# H atom and N atom are in fuctions
# '''

# Pressure = [1E+5]#, 1E+4, 1E+3, 1E+2, 1E+1, 1E+0, 1E-1]
# T = 400

# Entropy_list = []
# for P in Pressure:
#     Entropy = S_Hatom_gas(T, P)
#     Entropy_list.append(Entropy)
    
# print(Entropy_list)
  

# '''
# NH atom and N atom are in fuctions
# '''

# Pressure = [1E+5, 1E+4, 1E+3, 1E+2, 1E+1, 1E+0, 1E-1]
# T =398

# Entropy_list = []
# for P in Pressure:
#     Entropy = S_NH_gas_TransOnly(T, P)
#     Entropy_list.append(Entropy)
    
# print(Entropy_list)
##############################
#non linear test
# T = 1000

# N2Qr, N2Sr = Rotation_Linear(I_N2, Symmetry_Num_N2, T)
# N2Qr, N2Sr_nonlinear = Rotation_NonLinear(I_N2,I_N2,I_N2, Symmetry_Num_N2, T)
################


# '''
# N2
# '''
# T = 400
# PN2 =  1e+5 #pa
# N2Qt, N2St = Translation(mN2, T, PN2) 
# N2Qr, N2Sr = Rotation_Linear(I_N2, Symmetry_Num_N2, T)
# N2Qv, N2Sv = Single_Vibration(N2_VF, T)
# N2Qe, N2Se = ElectronicMotion(1, T)

# print(f'N2 Tranlational  Entropy : {N2St}')   
# print(f'N2 Rotational  Entropy   : {N2Sr}')
# print(f'N2 Vibrational Entropy   : {N2Sv}')
# print(f'N2 Electronic  Entropy   : {N2Se}')

# Stotal_N2 = Entropy_Total(N2St, N2Sr, N2Sv, N2Se)
# print(f'Stotal_N2 = {Stotal_N2}')



# '''
# H2
# '''
# T =400
# PH2 = 1e+5 #pa
# H2Qt, H2St = Translation(mH2, T, PH2) 
# H2Qr, H2Sr = Rotation_Linear(I_H2, Symmetry_Num_H2, T)
# H2Qv, H2Sv = Single_Vibration(H2_VF, T)
# H2Qe, H2Se = ElectronicMotion(1, T)

# print(f'H2 Tranlational  Entropy : {H2St}')   
# print(f'H2 Rotational  Entropy   : {H2Sr}')
# print(f'H2 Vibrational Entropy   : {H2Sv}')
# print(f'H2 Electronic  Entropy   : {H2Se}')

# Stotal_H2 = Entropy_Total(H2St, H2Sr, H2Sv, H2Se)
# print(f'Stotal_H2 = {Stotal_H2}')

# '''
# NH
# '''
# T =298
# PNH = 1e+5 #pa
# NHQt, NHSt = Translation(mNH, T, PNH) 
# NHQr, NHSr = Rotation_Linear(I_NH, Symmetry_Num_NH, T)
# NHQv, NHSv = Single_Vibration(NH_VF, T)
# NHQe, NHSe = ElectronicMotion(3, T)

# Stotal_NH = Entropy_Total(NHSt, NHSr, NHSv, NHSe)
# print(f'Stotal_NH = {Stotal_NH}')


########################################################
#Surface entropy 

#Lattice gas
    
# '''
# N2 on Metal Surface 
# '''
# T = 673
# # VF_list = [2138.987654, 375.053821, 372.066372, 364.439390, 135.539465, 106.985859]
# # Q_list, S_list, Stot_Lattice = LatticeGasModel(VF_list, T)
# # print(f'Stot_Lattice = {Stot_Lattice}')
# # Qv_list, Sv_list, Stot_2D = TwoDGasModel(VF_list[0:3], Symmetry_Num_N2, mN2, I_N2, 6.932, T)
# # print(f'Stot_2D = {Stot_2D}')
    

# '''
# H2 on Pd surface 
# '''
# VF_list = [4251.228760, 245.242216, 221.663618, 185.691376]
# # Q_list, S_list, Stot_Lattice = LatticeGasModel(VF_list, T)
# # print(f'Stot_Lattice = {Stot_Lattice}')
# Qv_list, Sv_list, Stot_2D = TwoDGasModel(VF_list[0:3], Symmetry_Num_H2, mH2, I_H2, 27.728, T)
# print(f'Stot_2D = {Stot_2D}')



# '''
# N on surface 
# '''

# Temperature = [298, 673]
# Stot_Lattice_list = []
# for T in Temperature:
#     VF_list = [3350.281349, 3326.146610, 3173.247143, 2974.701236, 2785.325823, 1627.159238, 1540.277727, 1528.201107, 1460.097027, 1394.968429, 1196.454209, 1097.421178, 1056.367972, 900.699534, 296.809069, 258.395994, 229.106109, 160.523765, 122.760503, 81.592556, 54.844807
# ]
#     Q_list, S_list, Stot_Lattice = LatticeGasModel(VF_list, T)
#     Stot_Lattice_list.append(Stot_Lattice)
    
# print(f'Stot_Lattice = {Stot_Lattice_list}')



# '''
# H on surface 
# '''

# T_list = [398]
# for T in T_list :
#     VF_list = [1347.502459, 905.114137, 329.709668]
#     Q_list, S_list, Stot_Lattice = LatticeGasModel(VF_list, T)
#     print(f'Stot_Lattice = {Stot_Lattice}')
    
    
'''
lattice Surface species in ev
S [=] J * mol-1 * K-1
'''

# T_list = [398]
# for T in T_list :
#     VF_list = [3426.853417, 3358.303242, 3317.034201, 1561.097019, 1369.141389, 1187.555329, 1060.661934, 840.020217, 579.742553, 448.197486, 346.927920, 253.085196, 136.404371, 103.711472, 51.465279]
#     Q_list, S_list, Stot_Lattice = LatticeGasModel(VF_list, T)
#     Toev = 1.0364E-5
#     S_tot_ev = Stot_Lattice * Toev
   
#     print(f'Stot_Lattice = {Stot_Lattice}')
#     print(len(VF_list))

##########################
# print(NNHgasEntropy(1000,1e+5))
    