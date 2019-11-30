#%%
import numpy as np

total_vol = 140*140*100
no_atoms = 10**5

density_liquid = [0.74,0.68]
density_vapour = [0.012,0.03]
#%%
data_1 = np.genfromtxt('data_1.csv',delimiter=',',skip_header=1)
data_1 = data_1[:,1:]

data_2 = np.genfromtxt('data_2.csv',delimiter=',',skip_header=1)
data_2 = data_2[:,1:]
# %%

def ASL(theta, cos_theta, density_vapour,density_liquid):
    v_drop = (no_atoms-total_vol*density_vapour)\
            /(density_liquid-density_vapour)

    r_drop = pow(3*v_drop/(np.pi*(1+cos_theta)**2 * (2-cos_theta)),1/3)

    rsl = r_drop*np.sin(np.deg2rad(theta))

    return np.pi*rsl**2 

def ALV(theta,cos_theta, density_vapour,density_liquid):
    v_drop = (no_atoms-total_vol*density_vapour)\
            /(density_liquid-density_vapour)

    r_drop = pow(3*v_drop/(np.pi*(1+cos_theta)**2 * (2-cos_theta)),1/3)
    
    rsl = r_drop*np.sin(np.deg2rad(theta))

    hcent = np.sqrt(r_drop**2-rsl**2) + r_drop
    return np.pi*(hcent**2+rsl**2)

def GAMMA_LV(data):
    return -data[:,1]/data[:,2]

# %%
asl_1 = ASL(theta=data_1[:,-2], cos_theta=data_1[:,2], \
            density_vapour=density_vapour[0], \
            density_liquid=density_liquid[0])

asl_2 = ASL(theta=data_2[:,-2], cos_theta=data_2[:,2], \
            density_vapour=density_vapour[1], \
            density_liquid=density_liquid[1])

alv_1 = ALV(theta=data_1[:,-2], cos_theta=data_1[:,2], \
            density_vapour=density_vapour[0], \
            density_liquid=density_liquid[0])

alv_2 = ALV(theta=data_2[:,-2], cos_theta=data_2[:,2], \
            density_vapour=density_vapour[1], \
            density_liquid=density_liquid[1])

gamma_lv_1 = GAMMA_LV(data_1)
gamma_lv_2 = GAMMA_LV(data_2)

gibbs_1 = data_1[:,1]*asl_1 + gamma_lv_1*alv_1
gibbs_2 = data_2[:,1]*asl_2 + gamma_lv_2*alv_2

gibbs = gibbs_2-gibbs_1