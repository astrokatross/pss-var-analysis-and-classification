#!/usr/bin/python
# This script is to analyse the sources and find/define any variability. 
# By K.Ross 30/8/19

import numpy as np
import pandas as pd
import pprocess
import multiprocessing
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
from astropy.io import fits
import gpscssmodels
import read_variables_forfitting as reading
from fitparams_script import fitting_params, calc_curve_snr
import time
import color_color_plotting 
import sed_funct 
import scipy.stats as ss
from astropy.visualization import hist


freq = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
freq_mwa = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
freq_high = np.array([189,212,843,1400])
freq_low = np.array([107,115,122,130,143,151,158,166,174,181,197,204,220,227])
freq_xtra = [74,150,408,843,1400]
freq_gleam = np.array([76,84,92,99])

# load in the new csv as pandas to make things easy maybe change the thing above so that its all done here but youd have to include the fluxes etc 
save_dir = '/data/var_analysis/analysis/'
all_pop = pd.read_csv('/data/var_analysis/analysis/all_pop_yr1yr2_params.csv')

# Just a few sanity checks 
all_pop['alpha_low_yr1'].describe()
all_pop['alpha_low_yr2'].describe()
Name = all_pop['Name']

# This should already be done but if not just uncomment and run again 
# SN_cond_met_yr1 = np.array(np.zeros(np.shape(Name)[0]))
# for i in range(len(Name)):
#     src_sn_conds = len(np.where(np.array(all_pop.iloc[i,134:150])==True)[0])
#     if src_sn_conds>=8:
#         SN_cond_met_yr1[i] = 1
# all_pop['SN_cond_met_yr1']=SN_cond_met_yr1


# SN_cond_met_yr2 = np.array(np.zeros(np.shape(Name)[0]))
# for i in range(len(Name)):
#     src_sn_conds = len(np.where(np.array(all_pop.iloc[i,150:166])==True)[0])
#     if src_sn_conds>=8:
#         SN_cond_met_yr2[i] = 1
# all_pop['SN_cond_met_yr2']=SN_cond_met_yr2




all_pop['alpha_low_yr1'] = pd.to_numeric(all_pop['alpha_low_yr1'], errors='coerce')
all_pop['quad_turnover_yr1'] = pd.to_numeric(all_pop['quad_turnover_yr1'], errors='coerce')
all_pop['quad_turnover_error_yr1'] = pd.to_numeric(all_pop['quad_turnover_error_yr1'], errors='coerce')
all_pop['quad_curve_yr1'] = pd.to_numeric(all_pop['quad_curve_yr1'], errors='coerce')
all_pop['alpha_low_yr2'] = pd.to_numeric(all_pop['alpha_low_yr2'], errors='coerce')
all_pop['quad_turnover_yr2'] = pd.to_numeric(all_pop['quad_turnover_yr2'], errors='coerce')
all_pop['quad_turnover_error_yr2'] = pd.to_numeric(all_pop['quad_turnover_error_yr2'], errors='coerce')
all_pop['quad_curve_yr2'] = pd.to_numeric(all_pop['quad_curve_yr2'], errors='coerce')


# all_pop.to_csv(save_dir+'all_pop_yr1yr2_params.csv', index=None, header=True)

print('Alpha_low Year 1: '+str(np.nanmedian(all_pop['alpha_low_yr1']))+" +/- "+str(np.nanstd(all_pop['alpha_low_yr1'])))
print('Alpha_low Year 2: '+str(np.nanmedian(all_pop['alpha_low_yr2']))+" +/- "+str(np.nanstd(all_pop['alpha_low_yr2'])))


print('Alpha_high Year 1: '+str(np.nanmedian(all_pop['alpha_high_yr1']))+" +/- "+str(np.nanstd(all_pop['alpha_high_yr1'])))
print('Alpha_high Year 2: '+str(np.nanmedian(all_pop['alpha_high_yr2']))+" +/- "+str(np.nanstd(all_pop['alpha_high_yr2'])))



print('Quad_curve Year 1: '+str(np.nanmedian(all_pop['quad_curve_yr1']))+" +/- "+str(np.nanstd(all_pop['quad_curve_yr1'])))
print('Quad_curve Year 2: '+str(np.nanmedian(all_pop['quad_curve_yr2']))+" +/- "+str(np.nanstd(all_pop['quad_curve_yr2'])))
print('Quad_curve error Year 1: '+str(np.nanmedian(all_pop['quad_curve_error_yr1']))+" +/- "+str(np.nanstd(all_pop['quad_curve_error_yr1'])))
print('Quad_curve error Year 2: '+str(np.nanmedian(all_pop['quad_curve_error_yr2']))+" +/- "+str(np.nanstd(all_pop['quad_curve_error_yr2'])))

# STEP 1: Total
print(len(all_pop))
all_pop_unres = all_pop.query('(a_wide_yr1*b_wide_yr1)/(psf_a_wide_yr1*psf_b_wide_yr1) <= 1.1')
# print(len(all_pop_unres))

# STEP 2: already correct decs
# STEP 3: Unresolved
# Wide 200MHz >= 0.16Jy 
all_pop_unres_bright = all_pop_unres.query("S_150_yr1>= 0.16")# & int_flux_wide_yr2 >= 0.16")
# print(len(all_pop_unres_bright))

# STEP 4: Bringt (27181)--> Sources with 8 or more GLEAM points with S/N >=3 (27099)
# all_pop_unres_bright_sncut = all_pop_unres_bright.query('SN_cond_met_yr1==1')# &SN_cond_met_yr2==1')
# print(len(all_pop_unres_bright_sncut))


# STEP 5: SN Cut (27099) --> NVSS/SUMSS Counter Part (27091)
norm_pop_yr1 = all_pop[all_pop.S_sumss.notnull() | all_pop.S_nvss.notnull()].query("(a_wide_yr1*b_wide_yr2)/(psf_a_wide_yr1*psf_b_wide_yr1) <= 1.1& int_flux_wide_yr1 >= 0.16 & SN_cond_met_yr1==1")
norm_pop_yr2 = all_pop[all_pop.S_sumss.notnull() | all_pop.S_nvss.notnull()].query("(a_wide_yr2*b_wide_yr2)/(psf_a_wide_yr2*psf_b_wide_yr2) <= 1.1& int_flux_wide_yr2 >= 0.16 & SN_cond_met_yr2==1")
print(len(norm_pop_yr1))
print(len(norm_pop_yr2))
norm_pop_yr1.to_csv(save_dir+'norm_pop_yr1.csv')
norm_pop_yr2.to_csv(save_dir+'norm_pop_yr2.csv')



norm_pop_yr1 =pd.read_csv('/data/var_analysis/analysis/norm_pop_yr1.csv')
norm_pop_yr2 =pd.read_csv('/data/var_analysis/analysis/norm_pop_yr2.csv')

var_param_norm = norm_pop_yr2['var_param']
var_param_gps = norm_pop_yr2.query('gps_id==1')['var_param']

color_color_plotting.plt_var_hist('/data/var_analysis/Plots/var_param_hist.pdf', var_param_norm, var_param_gps)



# STEP 6a: High frequency soft sample YEAR 1
hs_conds_yr1 = "alpha_low_yr1 >= 0.1 &(alpha_low_yr1 - alpha_low_error_yr1 >= 0.1)& alpha_high_yr1 <=-0.5"
hs_pop_yr1 = norm_pop_yr1.query(hs_conds_yr1)
print("HS in Year1: "+ str(len(hs_pop_yr1)))
hs_csv = hs_pop_yr1.to_csv(save_dir+'hs_pop_yr1.csv', index=None, header=True)


# STEP 6b: High frequency hard sample 
hh_conds_yr1 =  "alpha_low_yr1 >= 0.1 &(alpha_low_yr1 - alpha_low_error_yr1 >= 0.1) & alpha_high_yr1 >-0.5 & alpha_high_yr1 <= 0 "
hh_pop_yr1 = norm_pop_yr1.query(hh_conds_yr1)
print("HH in Year1: "+ str(len(hh_pop_yr1)))
hh_csv = hh_pop_yr1.to_csv(save_dir+'hh_pop_yr1.csv', index=None, header=True)


# STEP 6c: GPS Sample
gps_conds_yr1 = "alpha_low_yr1 >= 0.1 & alpha_high_yr1 > 0 &(alpha_low_yr1 - alpha_low_error_yr1 >= 0.1)"
gps_pop_yr1 = norm_pop_yr1.query(gps_conds_yr1)
print("GPS in Year1: "+ str(len(gps_pop_yr1)))
gps_csv = gps_pop_yr1.to_csv(save_dir+'gps_pop_yr1.csv', index=None, header=True)


# STEP 6d: Low frequency soft sample 
ls_conds_yr1 = "quad_curve_yr1<=-0.2 & quad_curve_error_yr1<=0.2 & alpha_low_yr1<0.1 & alpha_high_yr1 <=-0.5 & quad_turnover_yr1>= 107. & quad_turnover_yr1 <= 231. & snr_curve_cut_yr1 == 1"
ls_pop_yr1 = norm_pop_yr1.query(ls_conds_yr1)
print("LS in Year1: "+ str(len(ls_pop_yr1)))
ls_csv = ls_pop_yr1.to_csv(save_dir+'ls_pop_yr1.csv', index=None, header=True)

# STEP 6e: Low frequency hard sample 
lh_conds_yr1 = "alpha_low_yr1<0.1 & alpha_high_yr1>-0.5  & alpha_high_yr1<= 0 & quad_turnover_yr1>= 107. & quad_turnover_yr1<= 231. & quad_curve_yr1<=-0.2 & quad_curve_error_yr1<=0.2 & snr_curve_cut_yr1 ==1"
lh_pop_yr1 = norm_pop_yr1.query(lh_conds_yr1)
print("LH in Year1: "+ str(len(lh_pop_yr1)))
lh_csv = lh_pop_yr1.to_csv(save_dir+'lh_pop_yr1.csv', index=None, header=True)


# STEP KAT: Classifying powlaw and flat 
flat_conds_yr1 = 'alpha_low_yr1<0.2 & alpha_low_yr1 >=-0.2 &(alpha_low_yr1 - alpha_low_error_yr1 <= 0.1) & quad_curve_yr1>=-0.2 '
flat_pop_yr1 = norm_pop_yr1.query(flat_conds_yr1)
print("Flat in Year1: "+ str(len(flat_pop_yr1)))


powlaw_conds_yr1 = 'quad_curve_yr1>-0.2 &  alpha_low_yr1<0.1 & alpha_high_yr1 <=-0.5 &(alpha_low_yr1 - alpha_low_error_yr1 <= 0.1) '
powlaw_pop_yr1 =  norm_pop_yr1[~norm_pop_yr1['Name'].isin(gps_pop_yr1['Name'])&~norm_pop_yr1['Name'].isin(hh_pop_yr1['Name'])&~norm_pop_yr1['Name'].isin(hs_pop_yr1['Name'])&~norm_pop_yr1['Name'].isin(ls_pop_yr1['Name'])&~norm_pop_yr1['Name'].isin(lh_pop_yr1['Name'])&~norm_pop_yr1['Name'].isin(flat_pop_yr1['Name'])]
print("Powlaw in Year1: " + str(len(powlaw_pop_yr1)))

color_color_plotting.plt_colour_color('/data/var_analysis/Plots/color-color-yr1.png', 'alpha_low_yr1', 'alpha_high_yr1', norm_pop_yr1, gps_pop_yr1, hh_pop_yr1, hs_pop_yr1, ls_pop_yr1, lh_pop_yr1)
color_color_plotting.plt_colour_curve('/data/var_analysis/Plots/color-curve-yr1.png', 'alpha_low_yr1', 'quad_curve_yr1', 'quad_curve_error_yr1<=0.3', norm_pop_yr1, gps_pop_yr1, hh_pop_yr1, hs_pop_yr1, ls_pop_yr1, lh_pop_yr1)





# STEP 6a: High frequency soft sample YEAR 1
hs_conds_yr2 = "alpha_low_yr2 >= 0.1 &(alpha_low_yr2 - alpha_low_error_yr2 >= 0.1)& alpha_high_yr2 <=-0.5"
hs_pop_yr2 = norm_pop_yr2.query(hs_conds_yr2)
print("HS in Year2: "+ str(len(hs_pop_yr2)))
hs_csv = hs_pop_yr2.to_csv(save_dir+'hs_pop_yr2.csv', index=None, header=True)


# STEP 6b: High frequency hard sample 
hh_conds_yr2 =  "alpha_low_yr2 >= 0.1 &(alpha_low_yr2 - alpha_low_error_yr2 >= 0.1) & alpha_high_yr2 >-0.5 & alpha_high_yr2 <= 0 "
hh_pop_yr2 = norm_pop_yr2.query(hh_conds_yr2)
print("HH in Year2: "+ str(len(hh_pop_yr2)))
hh_csv = hh_pop_yr2.to_csv(save_dir+'hh_pop_yr2.csv', index=None, header=True)


# STEP 6c: GPS Sample
gps_conds_yr2 = "alpha_low_yr2 >= 0.1 & alpha_high_yr2 > 0 &(alpha_low_yr2 - alpha_low_error_yr2 >= 0.1)"
gps_pop_yr2 = norm_pop_yr2.query(gps_conds_yr2)
print("GPS in Year2: "+ str(len(gps_pop_yr2)))
gps_csv = gps_pop_yr2.to_csv(save_dir+'gps_pop_yr2.csv', index=None, header=True)


# STEP 6d: Low frequency soft sample 
ls_conds_yr2 = "quad_curve_yr2<=-0.2 & quad_curve_error_yr2<=0.2 & alpha_low_yr2<0.1 & alpha_high_yr2 <=-0.5 & quad_turnover_yr2>= 107. & quad_turnover_yr2 <= 231. & snr_curve_cut_yr2 == 1"
ls_pop_yr2 = norm_pop_yr2.query(ls_conds_yr2)
print("LS in Year2: "+ str(len(ls_pop_yr2)))
ls_csv = ls_pop_yr2.to_csv(save_dir+'ls_pop_yr2.csv', index=None, header=True)

# STEP 6e: Low frequency hard sample 
lh_conds_yr2 = "alpha_low_yr2<0.1 & alpha_high_yr2>-0.5  & alpha_high_yr2<= 0 & quad_turnover_yr2>= 107. & quad_turnover_yr2<= 231. & quad_curve_yr2<=-0.2 & quad_curve_error_yr2<=0.2 & snr_curve_cut_yr2 ==1"
lh_pop_yr2 = norm_pop_yr2.query(lh_conds_yr2)
print("LH in Year2: "+ str(len(lh_pop_yr2)))
lh_csv = lh_pop_yr2.to_csv(save_dir+'lh_pop_yr2.csv', index=None, header=True)



color_color_plotting.plt_colour_color('/data/var_analysis/Plots/color-color-yr2.png', 'alpha_low_yr2', 'alpha_high_yr2', norm_pop_yr2, gps_pop_yr2, hh_pop_yr2, hs_pop_yr2, ls_pop_yr2, lh_pop_yr2)
color_color_plotting.plt_colour_curve('/data/var_analysis/Plots/color-curve-yr2.png', 'alpha_low_yr2', 'quad_curve_yr2', 'quad_curve_error_yr2<=0.3', norm_pop_yr2, gps_pop_yr2, hh_pop_yr2, hs_pop_yr2, ls_pop_yr2, lh_pop_yr2)
plt.clf()

data_dir = '/data/var_analysis/analysis'
save_dir = '/data/var_analysis/Plots/SEDs/'
name, ra_deg, dec_deg, flux_yr1, flux_err_yr1, flux_mwa_yr1, flux_mwa_err_yr1, flux_low_yr1, flux_err_low_yr1, flux_high_yr1, flux_err_high_yr1, S_white_yr1, local_rms_yr1, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err = reading.read_variables(data_dir+'/norm_pop_yr1.fits','_yr1', 1)
name, ra_deg, dec_deg, flux_yr2, flux_err_yr2, flux_mwa_yr2, flux_mwa_err_yr2, flux_low_yr2, flux_err_low_yr2, flux_high_yr2, flux_err_high_yr2, S_white_yr2, local_rms_yr2, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err = reading.read_variables(data_dir+'/norm_pop_yr1.fits','_yr2', 1)
models = [gpscssmodels.quad_plot]
directory = save_dir+'yr1_gps_non_vary/'
gleam=1
other_surveys=1
plt_df = norm_pop_yr1
var_param=np.array(plt_df['var_param'])
gps_id = np.array(plt_df['gps_id'])
norm_low_yr1_plt = np.array(plt_df['norm_low_yr1'].astype('float64'))
norm_low_yr2_plt = np.array(plt_df['norm_low_yr2'].astype('float64'))
alpha_low_yr1_plt = np.array(plt_df['alpha_low_yr1'].astype('float64'))*-1
alpha_low_yr2_plt = np.array(plt_df['alpha_low_yr2'].astype('float64'))*-1
norm_high_yr1_plt = np.array(plt_df['norm_high_yr1'].astype('float64'))
norm_high_yr2_plt = np.array(plt_df['norm_high_yr2'].astype('float64'))
alpha_high_yr1_plt = np.array(plt_df['alpha_high_yr1'].astype('float64'))*-1
alpha_high_yr2_plt = np.array(plt_df['alpha_high_yr2'].astype('float64'))*-1
norm_quad_yr1_plt = np.array(plt_df['norm_quad_yr1'].astype('float64'))
norm_quad_yr2_plt = np.array(plt_df['norm_quad_yr2'].astype('float64'))
alpha_quad_yr1_plt = np.array(plt_df['alpha_quad_yr1'].astype('float64'))
alpha_quad_yr2_plt = np.array(plt_df['alpha_quad_yr2'].astype('float64'))
quad_curve_yr1_plt = np.array(plt_df['quad_curve_yr1'].astype('float64'))
quad_curve_yr2_plt = np.array(plt_df['quad_curve_yr2'].astype('float64'))
redchisq_low_yr1 = np.array(plt_df['redchisq_low_yr1'].astype('float64'))
redchisq_low_yr2 = np.array(plt_df['redchisq_low_yr2'].astype('float64'))
redchisq_quad_yr1 = np.array(plt_df['redchisq_quad_yr1'].astype('float64'))
redchisq_quad_yr2 = np.array(plt_df['redchisq_quad_yr2'].astype('float64'))




powlaw_low_params_yr1 = np.stack((norm_low_yr1_plt, alpha_low_yr1_plt))
powlaw_low_params_yr2 = np.stack((norm_low_yr2_plt, alpha_low_yr2_plt))
powlaw_high_params_yr1 = np.stack((norm_high_yr1_plt, alpha_high_yr1_plt))
powlaw_high_params_yr2 = np.stack((norm_high_yr2_plt, alpha_high_yr2_plt))

quad_gleam_params_yr1 = np.stack((norm_quad_yr1_plt, alpha_quad_yr1_plt, quad_curve_yr1_plt))
quad_gleam_params_yr2 = np.stack((norm_quad_yr2_plt, alpha_quad_yr2_plt, quad_curve_yr2_plt))



for i in range(len(name)):
	if gps_id[i] == 1 and var_param[i] <= 20:
		if dec_deg[i] <= -20 and dec_deg[i] >=-35 and ra_deg[i]<=36 and ra_deg[i] >=18:
			print('I AM GROSS REGION PASSING: '+name[i])
		elif dec_deg[i] <= -22 and dec_deg[i] >= -35 and ra_deg[i] <=95 and ra_deg[i] >= 54:
			print('I AM GROSS REGION PASSING: '+name[i])
		else: 
			if redchisq_low_yr2[i] < redchisq_quad_yr2[i]:
				models_yr1 = [gpscssmodels.quad_plot]
				models_yr2 = [gpscssmodels.powlaw]
				print('Processing '+name[i])
				# print('Redchisq yr2 quad is: '+str(redchisq_quad_yr2[i]))
				# print('Redchisq yr2 plaw is: '+str(redchisq_low_yr2[i]))
				powlaw_low_yr1 = powlaw_low_params_yr1[:,i]
				powlaw_low_yr2 = powlaw_low_params_yr2[:,i]
				powlaw_high_yr1 = powlaw_high_params_yr1[:,i]
				powlaw_high_yr2 = powlaw_high_params_yr2[:,i]
				# powlaw_gleam_yr1 = powlaw_gleam_params_yr1[:,i]
				# powlaw_gleam_yr2 = powlaw_gleam_params_yr2[:,i]
				quad_gleam_yr1 = quad_gleam_params_yr1[:,i]
				quad_gleam_yr2 = quad_gleam_params_yr2[:,i]


				paras_yr1 = [quad_gleam_yr1]
				paras_yr2 = [powlaw_low_yr2]

				flux_plt_yr1 = flux_yr1[:,i]
				flux_plt_yr2 = flux_yr2[:,i]
				flux_err_plt_yr1 = flux_err_yr1[:,i]
				flux_err_plt_yr2 = flux_err_yr2[:,i]
				flux_plt_xtra = flux_xtra[i,:]
				flux_err_plt_xtra = flux_xtra_err[i,:]
				var_value = var_param[i]
				var_value = format(var_value,'.2f')
				name_plt = name[i]
				q_value = format(quad_curve_yr2_plt[i],'.2f')
				gps_val = np.str(gps_id[i])
				flux_gleam_plt = flux_gleam[:,i]
				flux_gleam_err_plt = flux_gleam_err[:,i]

				sed_funct.sed(directory, models_yr1, models_yr2,paras_yr1,paras_yr2,freq,flux_plt_yr1,flux_err_plt_yr1,flux_plt_yr2,flux_err_plt_yr2, name_plt, var_value, q_value, gps_val, freq_xtra, flux_plt_xtra, flux_err_plt_xtra, freq_gleam, flux_gleam_plt, flux_gleam_err_plt, freq_labels = True, log = True, resid = True, savefig=True, xtra_low=True, xtra_high=True)
			elif quad_curve_yr2_plt[i] >=0: 
				models_yr1 = [gpscssmodels.quad_plot]
				models_yr2 = [gpscssmodels.powlaw]
				print('Processing '+name[i]+' yr1 is powlaw')
				powlaw_low_yr1 = powlaw_low_params_yr1[:,i]
				powlaw_low_yr2 = powlaw_low_params_yr2[:,i]
				powlaw_high_yr1 = powlaw_high_params_yr1[:,i]
				powlaw_high_yr2 = powlaw_high_params_yr2[:,i]
				# powlaw_gleam_yr1 = powlaw_gleam_params_yr1[:,i]
				# powlaw_gleam_yr2 = powlaw_gleam_params_yr2[:,i]
				quad_gleam_yr1 = quad_gleam_params_yr1[:,i]
				quad_gleam_yr2 = quad_gleam_params_yr2[:,i]


				paras_yr1 = [quad_gleam_yr1]
				paras_yr2 = [powlaw_low_yr2]

				flux_plt_yr1 = flux_yr1[:,i]
				flux_plt_yr2 = flux_yr2[:,i]
				flux_err_plt_yr1 = flux_err_yr1[:,i]
				flux_err_plt_yr2 = flux_err_yr2[:,i]
				flux_plt_xtra = flux_xtra[i,:]
				flux_err_plt_xtra = flux_xtra_err[i,:]
				var_value = var_param[i]
				var_value = format(var_value,'.2f')
				name_plt = name[i]
				q_value = format(quad_curve_yr2_plt[i],'.2f')
				gps_val = np.str(gps_id[i])
				flux_gleam_plt = flux_gleam[:,i]
				flux_gleam_err_plt = flux_gleam_err[:,i]

				sed_funct.sed(directory, models_yr1, models_yr2,paras_yr1,paras_yr2,freq,flux_plt_yr1,flux_err_plt_yr1,flux_plt_yr2,flux_err_plt_yr2, name_plt, var_value, q_value, gps_val, freq_xtra, flux_plt_xtra, flux_err_plt_xtra, freq_gleam, flux_gleam_plt, flux_gleam_err_plt, freq_labels = True, log = True, resid = True, savefig=True, xtra_low=True, xtra_high=True)
			else: 
				models_yr1 = [gpscssmodels.quad_plot]
				models_yr2 = [gpscssmodels.quad_plot]
				print('Processing '+name[i]+' yr1 is quad')
				powlaw_low_yr1 = powlaw_low_params_yr1[:,i]
				powlaw_low_yr2 = powlaw_low_params_yr2[:,i]
				powlaw_high_yr1 = powlaw_high_params_yr1[:,i]
				powlaw_high_yr2 = powlaw_high_params_yr2[:,i]
				# powlaw_gleam_yr1 = powlaw_gleam_params_yr1[:,i]
				# powlaw_gleam_yr2 = powlaw_gleam_params_yr2[:,i]
				quad_gleam_yr1 = quad_gleam_params_yr1[:,i]
				quad_gleam_yr2 = quad_gleam_params_yr2[:,i]


				paras_yr1 = [quad_gleam_yr1]
				paras_yr2 = [quad_gleam_yr2]

				flux_plt_yr1 = flux_yr1[:,i]
				flux_plt_yr2 = flux_yr2[:,i]
				flux_err_plt_yr1 = flux_err_yr1[:,i]
				flux_err_plt_yr2 = flux_err_yr2[:,i]
				flux_plt_xtra = flux_xtra[i,:]
				flux_err_plt_xtra = flux_xtra_err[i,:]
				var_value = var_param[i]
				var_value = format(var_value,'.2f')
				name_plt = name[i]
				q_value = format(quad_curve_yr2_plt[i],'.2f')
				gps_val = np.str(gps_id[i])
				flux_gleam_plt = flux_gleam[:,i]
				flux_gleam_err_plt = flux_gleam_err[:,i]

				sed_funct.sed(directory, models_yr1, models_yr2,paras_yr1,paras_yr2,freq,flux_plt_yr1,flux_err_plt_yr1,flux_plt_yr2,flux_err_plt_yr2, name_plt, var_value, q_value, gps_val, freq_xtra, flux_plt_xtra, flux_err_plt_xtra, freq_gleam, flux_gleam_plt, flux_gleam_err_plt, freq_labels = True, log = True, resid = True, savefig=True, xtra_low=True, xtra_high=True)
