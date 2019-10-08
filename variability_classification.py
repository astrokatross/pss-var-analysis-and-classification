#!/usr/bin/python

# This is a script to classify the type of variability present by both using the BIC to determine change of fit and the fit of the residuals 
# K.Ross 8/10/19


import numpy as np
import matplotlib.pyplot as plt 
import gpscssmodels
import pandas as pd
import read_variables_forfitting as reading
import sed_funct 
plt.rcParams["font.family"] = "serif"


# Setting the models to compare (currently just wanna know if its curved or plaw) and setting the constant variables
data_dir = '/data/var_analysis/analysis/'
models_yr1 = [gpscssmodels.powlaw, gpscssmodels.quad_plot]
models_yr2 = [gpscssmodels.powlaw, gpscssmodels.quad_plot]
freq = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
freq_mwa = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
freq_high = np.array([189,212,843,1400])
freq_low = np.array([107,115,122,130,143,151,158,166,174,181,197,204,220,227])
freq_xtra = [74,150,408,843,1400]
freq_gleam = np.array([76,84,92,99])
sed_dir = '/data/var_analysis/Plots/SEDs/gps_variables/'


# Reading in the data 
norm_pop_yr2 =pd.read_csv(data_dir+'norm_pop_yr2.csv')
name, ra_deg, dec_deg, flux_yr1, flux_err_yr1, flux_mwa_yr1, flux_mwa_err_yr1, flux_low_yr1, flux_err_low_yr1, flux_high_yr1, flux_err_high_yr1, S_white_yr1, local_rms_yr1, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err = reading.read_variables(data_dir+'/norm_pop_yr2.fits','_yr1', 1)
name, ra_deg, dec_deg, flux_yr2, flux_err_yr2, flux_mwa_yr2, flux_mwa_err_yr2, flux_low_yr2, flux_err_low_yr2, flux_high_yr2, flux_err_high_yr2, S_white_yr2, local_rms_yr2, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err = reading.read_variables(data_dir+'/norm_pop_yr2.fits','_yr2', 1)
plt_df = norm_pop_yr2
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

# Setting the parameter values for the fits
powlaw_low_params_yr1 = np.stack((norm_low_yr1_plt, alpha_low_yr1_plt))
powlaw_low_params_yr2 = np.stack((norm_low_yr2_plt, alpha_low_yr2_plt))
powlaw_high_params_yr1 = np.stack((norm_high_yr1_plt, alpha_high_yr1_plt))
powlaw_high_params_yr2 = np.stack((norm_high_yr2_plt, alpha_high_yr2_plt))
quad_gleam_params_yr1 = np.stack((norm_quad_yr1_plt, alpha_quad_yr1_plt, quad_curve_yr1_plt))
quad_gleam_params_yr2 = np.stack((norm_quad_yr2_plt, alpha_quad_yr2_plt, quad_curve_yr2_plt))

chi_sq_plaw1, chi_sq_plaw2, chi_sq_quad1, chi_sq_quad2 = [np.zeros(len(name)) for dummy in range(4)]

for i in range(len(name)):
	powlaw_low_yr1 = powlaw_low_params_yr1[:,i]
	powlaw_low_yr2 = powlaw_low_params_yr2[:,i]
	quad_gleam_yr1 = quad_gleam_params_yr1[:,i]
	quad_gleam_yr2 = quad_gleam_params_yr2[:,i]
	paras_yr1 = [powlaw_low_yr1, quad_gleam_yr1]
	paras_yr2 = [powlaw_low_yr2, quad_gleam_yr2]
	model_points_plaw2 = models_yr2[0](freq,*paras_yr2[0])
	model_points_plaw1 = models_yr1[0](freq,*paras_yr1[0])
	model_points_quad2 = models_yr2[1](freq,*paras_yr2[1])
	model_points_quad1 = models_yr1[1](freq,*paras_yr1[1])

	# Calculating the chi squared for each model to be used in the BIC calculation
	chi_sq_plaw1[i] = np.sum(((flux_yr1[:,i]-model_points_plaw1)/flux_err_yr1[:,i])**2)
	chi_sq_plaw2[i] = np.sum(((flux_yr2[:,i]-model_points_plaw2)/flux_err_yr2[:,i])**2)
	chi_sq_quad1[i] = np.sum(((flux_yr1[:,i]-model_points_quad1)/flux_err_yr1[:,i])**2)
	chi_sq_quad2[i] = np.sum(((flux_yr2[:,i]-model_points_quad2)/flux_err_yr2[:,i])**2)



# Now to calculate BIC 
BIC_plaw_yr1 = chi_sq_plaw1+(2*np.log(16))
BIC_plaw_yr2 = chi_sq_plaw2+(2*np.log(16))
BIC_quad_yr1 = chi_sq_quad1+(3*np.log(16))
BIC_quad_yr2 = chi_sq_quad2+(3*np.log(16))

# Since quad - plaw means if delta_BIC is large then powerlaw is better fit by far if ~0 means powerlaw again 
delta_BIC_yr1 = BIC_quad_yr1 - BIC_plaw_yr1
delta_BIC_yr2 = BIC_quad_yr2 - BIC_plaw_yr2
print(np.where(delta_BIC_yr1<=-6))
print(np.where(delta_BIC_yr2<=-6))
curved_inds_yr1 = np.where(delta_BIC_yr1<=-6)
curved_inds_yr2 = np.where(delta_BIC_yr2<=-6)

curved_yr1_nms = name[curved_inds_yr1]
curved_yr2_nms = name[curved_inds_yr2]


for i in range(len(name)):
	if var_param[i] >=37 and gps_id[i]==1:
		if dec_deg[i] <= -20 and dec_deg[i] >=-35 and ra_deg[i]<=36 and ra_deg[i] >=18:
			print('I AM GROSS REGION PASSING: '+name[i])
		elif dec_deg[i] <= -22 and dec_deg[i] >= -35 and ra_deg[i] <=95 and ra_deg[i] >= 54:
			print('I AM GROSS REGION PASSING: '+name[i])
		else: 
			print('Processing '+name[i])
			powlaw_low_yr1 = powlaw_low_params_yr1[:,i]
			powlaw_low_yr2 = powlaw_low_params_yr2[:,i]
			powlaw_high_yr1 = powlaw_high_params_yr1[:,i]
			powlaw_high_yr2 = powlaw_high_params_yr2[:,i]
			quad_gleam_yr1 = quad_gleam_params_yr1[:,i]
			quad_gleam_yr2 = quad_gleam_params_yr2[:,i]
			if delta_BIC_yr1[i] <= -6 and delta_BIC_yr2[i] <=-6 and quad_curve_yr1_plt[i] <0 and quad_curve_yr2_plt[i] < 0:
				print('We\'re both curved!')
				models_yr1 = [gpscssmodels.quad_plot]
				models_yr2 = [gpscssmodels.quad_plot]
				paras_yr1 = [quad_gleam_yr1]
				paras_yr2 = [quad_gleam_yr2]

			elif delta_BIC_yr1[i] <= -6 and delta_BIC_yr2[i] > -6 and quad_curve_yr1_plt[i] <0:
				print('Only year1 curved') 
				models_yr1 = [gpscssmodels.quad_plot]
				models_yr2 = [gpscssmodels.powlaw]
				paras_yr1 = [quad_gleam_yr1]
				paras_yr2 = [powlaw_low_yr2]

			elif delta_BIC_yr2[i] <= -6 and delta_BIC_yr1[i] > 6 and quad_curve_yr2_plt[i] < 0: 
				print('Only year2 curved') 
				models_yr1 = [gpscssmodels.powlaw]
				models_yr2 = [gpscssmodels.quad_plot]
				paras_yr1 = [powlaw_low_yr1]
				paras_yr2 = [quad_gleam_yr2]

			else: 
				print('Neither of us is curved!') 
				models_yr1 = [gpscssmodels.powlaw]
				models_yr2 = [gpscssmodels.powlaw]
				paras_yr1 = [powlaw_low_yr1]
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

			sed_funct.sed(sed_dir, models_yr1, models_yr2,paras_yr1,paras_yr2,freq,flux_plt_yr1,flux_err_plt_yr1,flux_plt_yr2,flux_err_plt_yr2, name_plt, var_value, q_value, gps_val, freq_xtra, flux_plt_xtra, flux_err_plt_xtra, freq_gleam, flux_gleam_plt, flux_gleam_err_plt, freq_labels = True, log = True, resid = True, savefig=True, xtra_low=True, xtra_high=True)