#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
import matplotlib.pyplot as plt
import climlab
from climlab import constants as constants
import sys
import xarray as xr
import pdb
import copy as cp

# Pour fixer la taille de la police partout
import matplotlib as matplotlib
font = {'family' : 'monospace',
        'size'   : 15}
matplotlib.rc('font', **font)

outpath = './figures/'
os.makedirs(outpath, exist_ok=True)

units=r'W m$^{-2}$' # Unités puisance
alb=.25 # Albedo surface
levels=np.arange(200,330,20)
levels=[298]
Tlims=[180,310]
Nz=30 # nombre de niveaux verticaux

# Load the reference vertical temperature profile from the NCEP reanalysis
ncep_lev=np.load('npy/ncep_lev.npy')
ncep_T=np.load('npy/ncep_T.npy')+273.15

#  State variables (Air and surface temperature)
state = climlab.column_state(num_lev=30)

#  Fixed relative humidity
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)

#  Couple water vapor to radiation
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)

# Creation d'un modele couplé avec rayonemment et vapour d'eau
rcm = climlab.couple([rad,h2o], name='Radiative-Equilibrium Model')
rcm2 = climlab.process_like(rcm) # creation d'un clone du modele rcm

print('\n','\n','********************************************')
print('Control simulation ')
print('********************************************')
# Make the initial state isothermal
rcm.state.Tatm[:] = rcm.state.Ts
T=[]
q=[]
tr=[]
print(state)
# Plot temperature
for t in range(1000):
    T.append(cp.deepcopy(rcm.Tatm))
    q.append(cp.deepcopy(rcm.q))
    plt.plot(rcm.Tatm,rcm.lev)
    rcm.step_forward() #run the model forward one time step
    if abs(rcm.ASR - rcm.OLR)<1: # in W/m2
        tr.append(t)
plt.title('equilibrium reached at time t='+str(tr[0]))
plt.xlabel('temperature (K)')
plt.ylabel('pression (hPa)')
plt.gca().invert_yaxis()
fig_name=outpath+'fig1.png'
plt.savefig(fig_name,bbox_inches='tight')
plt.close()
print('output figure: ', fig_name)

#Plot humidity
for t in range(1000):
    plt.plot(q[t],rcm.lev)
plt.xlabel('specific humidity (kg/kg)')
plt.ylabel('pression (hPa)')
fig_name=outpath+'fig2.png'
plt.gca().invert_yaxis()
plt.savefig(fig_name,bbox_inches='tight')
plt.close()
print('output figure: ', fig_name)

# Quel est la sortie du modèle ?
#print('humidité:', rcm.humidity)
print('diagnostics: ',rcm.diagnostics,'\n')
print('tendencies',rcm.tendencies,'\n')
print('Tair: ',rcm.Tatm,'\n')
print('albedo',rcm.SW_flux_up[-1]/rcm.SW_flux_down[-1],'\n')
print('co2',rad.absorber_vmr['CO2'],'\n') #volumetric mixing ratio
print('ch4',rad.absorber_vmr['CH4'],'\n') #volumetric mixing ratio

print('\n','\n','********************************************')
print('Sensitivity to the concentration of gases in the atmosphere')
print('********************************************')

colors = ['k', 'r', 'g', 'orange', 'blue', 'purple']
linestyles = ['solid', 'dashed', 'dotted']

def field_to_array(f):
    try:
        a = np.asarray(f)
        if a.size != 0:
            return a
    except Exception:
        pass
    for attr in ('data', 'value', 'values', 'array'):
        if hasattr(f, attr):
            try:
                a = getattr(f, attr)
                return np.asarray(a)
            except Exception:
                pass
    try:
        a = f[:]
        return np.asarray(a)
    except Exception:
        pass
    try:
        return np.array([float(f)])
    except Exception:
        raise TypeError(f"Unable to convert object of type {type(f)} to numpy array.")

def mean_temperature(rcm_model):

    Tatm_arr = field_to_array(rcm_model.Tatm)
    Tatm_mean = float(np.mean(Tatm_arr))

    try:
        Ts_val = float(rcm_model.Ts)
    except Exception:
        Ts_arr = field_to_array(rcm_model.Ts)
        Ts_val = float(Ts_arr.flat[0])
    return Tatm_mean, Ts_val

#contrôle
Tatm_mean_ctrl, Ts_ctrl = mean_temperature(rcm)
label_ctrl = f'Contrôle – ⟨Tatm⟩={Tatm_mean_ctrl:.1f} K, Ts={Ts_ctrl:.1f} K'
plt.plot(field_to_array(rcm.Tatm)[::-1], rcm.lev[::-1], color=colors[0], label=label_ctrl, linestyle='solid')
plt.plot(Ts_ctrl, 1000, marker='s', color=colors[0])

#préindustrielle 280ppm
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)

rcm_preind = climlab.couple([rad, h2o], name='Radiative-Convective Model')
rcm_preind.absorber_vmr['CO2'] = 280e-6
rcm_preind.integrate_years(2)

Tatm_mean_preind, Ts_preind = mean_temperature(rcm_preind)
label_preind = f'CO₂ = 280 ppm – ⟨Tatm⟩={Tatm_mean_preind:.1f} K, Ts={Ts_preind:.1f} K'
plt.plot(field_to_array(rcm_preind.Tatm)[::-1], rcm_preind.lev[::-1], color='blue', label=label_preind)
plt.plot(Ts_preind, 1000, marker='s', color='blue')

# actuelle 420ppm
rcm_ajd = climlab.couple([rad, h2o], name='Radiative-Convective Model')
rcm_ajd.absorber_vmr['CO2'] = 420e-6
rcm_ajd.integrate_years(2)

Tatm_mean_ajd, Ts_ajd = mean_temperature(rcm_ajd)
label_ajd = f'CO₂ = 420 ppm – ⟨Tatm⟩={Tatm_mean_ajd:.1f} K, Ts={Ts_ajd:.1f} K'
plt.plot(field_to_array(rcm_ajd.Tatm)[::-1], rcm_ajd.lev[::-1], color='purple', label=label_ajd)
plt.plot(Ts_ajd, 1000, marker='s', color='purple')

# simulations double et quadruple pour chaque gaz
for gi, gg in enumerate(['O3', 'CO2', 'CH4']):
    for mult, style in zip([2, 4], ['--', ':']):
        state = climlab.column_state(num_lev=30)
        h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
        rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
        rcm_gas = climlab.couple([rad, h2o], name='Radiative-Convective Model')
        rcm_gas.absorber_vmr[gg] *= mult
        rcm_gas.integrate_years(2)

        Tatm_mean, Ts = mean_temperature(rcm_gas)
        label = f'{gg} ×{mult} – ⟨Tatm⟩={Tatm_mean:.1f} K, Ts={Ts:.1f} K'

        plt.plot(field_to_array(rcm_gas.Tatm)[::-1], rcm_gas.lev[::-1],
                 color=colors[gi+1], linestyle=style, label=label)
        plt.plot(Ts, 1000, marker='s', color=colors[gi+1])

# NCEP
plt.plot(ncep_T, ncep_lev, marker='x', color='k', label='NCEP reanalysis')

plt.gca().invert_yaxis()
plt.title('Sensibilité : gaz atmosphériques')
plt.ylabel('Pression (hPa)')
plt.xlabel('Température (K)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

fig_name = outpath + 'fig_sensitivity_gases_tempmean.png'
plt.savefig(fig_name, bbox_inches='tight', dpi=300)
plt.close()


print('\n\n********************************************')
print('albedo pour compenser co2')
print('********************************************')

#contrôle
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcm_control = climlab.couple([rad, h2o], name='Control')
rcm_control.integrate_years(2)
Ts_control = float(rcm_control.Ts)


#CO2 doublé
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcm_co2_2x = climlab.couple([rad, h2o], name='CO2_2x')
rcm_co2_2x.absorber_vmr['CO2'] = rcm_co2_2x.absorber_vmr['CO2'] * 2.0
rcm_co2_2x.integrate_years(2)
Ts_co2_2x = float(rcm_co2_2x.Ts)
dT_co2 = Ts_co2_2x - Ts_control


# Trouver l'albédo qui compense
# On cherche alpha_comp tel que Ts(CO2x2, alpha_comp) ≈ Ts_control
def Ts_for_albedo_with_co2_2x(alpha_value, integrate_years=2):
    st = climlab.column_state(num_lev=30)
    h = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=st)
    r = climlab.radiation.RRTMG(name='Radiation', state=st, specific_humidity=h.q, albedo=alpha_value)
    m = climlab.couple([r, h], name=f'co2_2x_alb_{alpha_value}')
    m.absorber_vmr['CO2'] = m.absorber_vmr['CO2'] * 2.0
    m.integrate_years(integrate_years)
    return float(m.Ts), m


a_low = 0.0
a_high = 0.9
tol_T = 0.05  # tolérance en K pour considérer Ts égal au contrôle
max_iter = 12

# Evaluations aux bornes
Ts_low, _ = Ts_for_albedo_with_co2_2x(a_low, integrate_years=2)
Ts_high, _ = Ts_for_albedo_with_co2_2x(a_high, integrate_years=2)
print(f"Ts at alb={a_low:.2f} -> {Ts_low:.3f} K; at alb={a_high:.2f} -> {Ts_high:.3f} K")


if not (Ts_high < Ts_control < Ts_low or Ts_low < Ts_control < Ts_high or True):
    pass

lo = a_low
hi = a_high
alpha_found = None
for it in range(max_iter):
    mid = 0.5 * (lo + hi)
    Ts_mid, _ = Ts_for_albedo_with_co2_2x(mid, integrate_years=2)
    resid = Ts_mid - Ts_control
    print(f"iter {it+1}: alpha={mid:.4f}, Ts_mid={Ts_mid:.3f}, resid={resid:.3f} K")
    if abs(resid) <= tol_T:
        alpha_found = mid
        break
    if Ts_mid > Ts_control:
        lo = mid
    else:
        hi = mid
else:
    alpha_found = mid
    Ts_mid, _ = Ts_for_albedo_with_co2_2x(alpha_found, integrate_years=2)
    resid = Ts_mid - Ts_control
    print("Bisection terminée sans tolérance stricte, retour dernier mid.")

alpha_comp = alpha_found
Ts_comp, rcm_comp = Ts_for_albedo_with_co2_2x(alpha_comp, integrate_years=2)
dT_comp = Ts_comp - Ts_control
delta_albedo = alpha_comp - alb

# Tracer les profils
plt.figure(figsize=(8,6))
plt.plot(rcm_control.Tatm[::-1], rcm_control.lev[::-1], color='k', label='Contrôle albedo = 0.2500')
plt.plot(rcm_control.Ts, 1000, marker='s', color='k')
plt.plot(rcm_co2_2x.Tatm[::-1], rcm_co2_2x.lev[::-1], color='r', linestyle='--', label='CO2 ×2')
plt.plot(rcm_co2_2x.Ts, 1000, marker='s', color='r')
plt.plot(rcm_comp.Tatm[::-1], rcm_comp.lev[::-1], color='b', linestyle='-.', label=f'CO2 ×2 + albedo={alpha_comp:.4f}')
plt.plot(rcm_comp.Ts, 1000, marker='s', color='b')
plt.gca().invert_yaxis()
plt.xlabel('Temperature (K)')
plt.ylabel('Pression (hPa)')
plt.title('Compensation du réchauffement par changement d\'albédo (CO2 doublé)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig_name = outpath + 'fig_albedo_compensation_co2_2x.png'
plt.savefig(fig_name, bbox_inches='tight')
plt.close()
print(f"Saved figure: {fig_name}")


print('\n','\n','********************************************')
print('Sensitivity to convection')
print('********************************************')
alb=.25
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcms={}
rcms['rcm0'] = climlab.couple([rad,h2o], name='Radiative-Convective Model')
conv = climlab.convection.ConvectiveAdjustment(name='Convection', state=state, adj_lapse_rate=3.6)
rcms['rcm1'] = climlab.couple([rad,conv,h2o], name='Radiative-Convective Model')
conv = climlab.convection.ConvectiveAdjustment(name='Convection', state=state, adj_lapse_rate=5)
rcms['rcm2'] = climlab.couple([rad,conv,h2o], name='Radiative-Convective Model')
conv = climlab.convection.ConvectiveAdjustment(name='Convection', state=state, adj_lapse_rate=6.5)
rcms['rcm3'] = climlab.couple([rad,conv,h2o], name='Radiative-Convective Model')
conv = climlab.convection.ConvectiveAdjustment(name='Convection', state=state, adj_lapse_rate=7.5) #lapse rate in degC per km
rcms['rcm4'] = climlab.couple([rad,conv,h2o], name='Radiative-Convective Model')
conv = climlab.convection.ConvectiveAdjustment(name='Convection', state=state, adj_lapse_rate=9.8)
rcms['rcm5'] = climlab.couple([rad,conv,h2o], name='Radiative-Convective Model')

mod_name=['control','conv-3.6 : très stable','conv-5 : stable', 'conv-6.5 : atmosphère standard, neutre', 'conv-7.5 : instable', 'conv-9.8 : gradient adiabatique sec, très instable']
for ai in range(6):
    rcms['rcm'+str(ai)].integrate_years(2)
    plt.plot(rcms['rcm'+str(ai)].Tatm[::-1], rcm.lev[::-1], marker='s', label=mod_name[ai],color=colors[ai])
    plt.plot(rcms['rcm'+str(ai)].Ts, 1000, marker='s',color=colors[ai])
plt.plot(ncep_T, ncep_lev, marker='x',color='k',label='NCEP reanalysis')
plt.gca().invert_yaxis()
plt.title('Sensitivity: convection')
plt.ylabel('Pression (hPa)')
plt.xlabel('Temperature (K)')
plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig_name=outpath+'lapse_rate.png'
print('output figure: ', fig_name)
plt.savefig(fig_name,bbox_inches='tight')
plt.close()


#bonus
#profil vertical de l'absorption solaire

def compute_SW_abs_profile(model):
    model.integrate_years(1)

    SW_down = np.array(model.SW_flux_down)
    SW_up = np.array(model.SW_flux_up)

    SW_net = SW_down - SW_up
    SW_abs_profile = -np.diff(SW_net)

    lev = np.array(model.lev)
    if len(SW_net) == len(lev) + 1:
        p_layer = lev
    elif len(SW_net) == len(lev):
        p_layer = 0.5 * (lev[:-1] + lev[1:])

    return SW_abs_profile, p_layer


#controle
alb = 0.2
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcm_std = climlab.couple([rad, h2o], name='Standard Atmosphere')

SW_std, p_layer = compute_SW_abs_profile(rcm_std)


#2x o3
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcm_2o3 = climlab.couple([rad, h2o], name='2xO3')
rcm_2o3.absorber_vmr['O3'] *= 2.0

SW_2o3, _ = compute_SW_abs_profile(rcm_2o3)


#5x o3
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcm_5o3 = climlab.couple([rad, h2o], name='5xO3')
rcm_5o3.absorber_vmr['O3'] *= 5.0

SW_5o3, _ = compute_SW_abs_profile(rcm_5o3)


#10x o3
state = climlab.column_state(num_lev=30)
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
rcm_10o3 = climlab.couple([rad, h2o], name='10xO3')
rcm_10o3.absorber_vmr['O3'] *= 10.0

SW_10o3, _ = compute_SW_abs_profile(rcm_10o3)

#graph
plt.figure(figsize=(6,7))
plt.plot(SW_std, p_layer, color='black', linewidth=2, label='Contrôle')
plt.plot(SW_2o3, p_layer, color='red', linestyle=':', linewidth=2, label='2×O₃')
plt.plot(SW_5o3, p_layer, color='blue', linestyle=':', linewidth=2, label='5×O₃')
plt.plot(SW_10o3, p_layer, color='green', linestyle=':', linewidth=2, label='10×O₃')

plt.gca().invert_yaxis()
plt.xlabel('Absorption solaire (W/m²)')
plt.ylabel('Pression (hPa)')
plt.title('Profil vertical de l’absorption solaire')
plt.legend(loc='lower right')
plt.grid(True, alpha=0.3)

fig_name = outpath + 'profil_absorption_o3.png'
plt.savefig(fig_name, bbox_inches='tight', dpi=300)
plt.close()






