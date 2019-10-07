# This script is to import both catalogues of GLEAM years and classify them. It should be run post colour-colour-fitparams.py for each year. 
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

# Creating the empty arrays to put in the fit values: NOTE YOU WILL NOT NEED SOME OF THESE AND SHOULD DISCARD THEM 
freq = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
freq_mwa = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
freq_high = np.array([189,212,843,1400])
freq_low = np.array([107,115,122,130,143,151,158,166,174,181,197,204,220,227])
freq_xtra = [74,150,408,843,1400]
freq_gleam = np.array([76,84,92,99])


# Setting up multi-porcessing stuff 
cores = multiprocessing.cpu_count()

name, ra_deg, dec_deg, flux_yr1, flux_err_yr1, flux_mwa_yr1, flux_mwa_err_yr1, flux_low_yr1, flux_err_low_yr1, flux_high_yr1, flux_err_high_yr1, S_white_yr1, local_rms_yr1, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err = reading.read_variables('/data/var_analysis/analysis/all_pop_yr1yr2_xmatch.fits','_yr1', 1)
name, ra_deg, dec_deg, flux_yr2, flux_err_yr2, flux_mwa_yr2, flux_mwa_err_yr2, flux_low_yr2, flux_err_low_yr2, flux_high_yr2, flux_err_high_yr2, S_white_yr2, local_rms_yr2, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err = reading.read_variables('/data/var_analysis/analysis/all_pop_yr1yr2_xmatch.fits','_yr2', 1)

print('fitting yr1 params')
params_yr1 = Parallel(n_jobs=cores)(delayed(fitting_params)(name[i], freq_low, freq_high, flux_high_yr1[:,i], flux_err_high_yr1[:,i], flux_low_yr1[:,i], flux_err_low_yr1[:,i], freq_mwa, flux_mwa_yr1[:,i], flux_mwa_err_yr1[:,i]) for i in range(len(name)))
print('fitting yr2 params')
params_yr2 = Parallel(n_jobs=cores)(delayed(fitting_params)(name[i], freq_low, freq_high, flux_high_yr2[:,i], flux_err_high_yr2[:,i], flux_low_yr2[:,i], flux_err_low_yr2[:,i], freq_mwa, flux_mwa_yr2[:,i], flux_mwa_err_yr2[:,i]) for i in range(len(name)))

alpha_low_yr1, alpha_low_error_yr1, alpha_high_yr1, alpha_high_error_yr1, norm_low_yr1, norm_high_yr1, quad_curve_yr1, quad_curve_error_yr1, alpha_quad_yr1, alpha_quad_error_yr1, quad_turnover_yr1, quad_turnover_error_yr1, norm_quad_yr1, redchisq_low_yr1, redchisq_high_yr1, redchisq_quad_yr1 = [np.zeros(len(name)) for dummy in range(16)]
alpha_low_yr2, alpha_low_error_yr2, alpha_high_yr2, alpha_high_error_yr2, norm_low_yr2, norm_high_yr2, quad_curve_yr2, quad_curve_error_yr2, alpha_quad_yr2, alpha_quad_error_yr2, quad_turnover_yr2, quad_turnover_error_yr2, norm_quad_yr2, redchisq_low_yr2, redchisq_high_yr2, redchisq_quad_yr2 = [np.zeros(len(name)) for dummy in range(16)]

for i in range(len(name)):
	alpha_low_yr1[i] = params_yr1[i][0]
	alpha_low_error_yr1[i] = params_yr1[i][1]
	alpha_high_yr1[i] = params_yr1[i][2]*-1
	alpha_high_error_yr1[i] = params_yr1[i][3]
	norm_low_yr1[i] = params_yr1[i][4]
	norm_high_yr1[i] = params_yr1[i][5]
	quad_curve_yr1[i] = params_yr1[i][6]
	quad_curve_error_yr1[i] = params_yr1[i][7]
	alpha_quad_yr1[i] = params_yr1[i][8]
	alpha_quad_error_yr1[i] = params_yr1[i][9]
	quad_turnover_yr1[i] = params_yr1[i][10]
	quad_turnover_error_yr1[i] = params_yr1[i][11]
	norm_quad_yr1[i] = params_yr1[i][12]
	alpha_low_yr2[i] = params_yr2[i][0]
	alpha_low_error_yr2[i] = params_yr2[i][1]
	alpha_high_yr2[i] = params_yr2[i][2]*-1
	alpha_high_error_yr2[i] = params_yr2[i][3]
	norm_low_yr2[i] = params_yr2[i][4]
	norm_high_yr2[i] = params_yr2[i][5]
	quad_curve_yr2[i] = params_yr2[i][6]
	quad_curve_error_yr2[i] = params_yr2[i][7]
	alpha_quad_yr2[i] = params_yr2[i][8]
	alpha_quad_error_yr2[i] = params_yr2[i][9]
	quad_turnover_yr2[i] = params_yr2[i][10]
	quad_turnover_error_yr2[i] = params_yr2[i][11]
	norm_quad_yr2[i] = params_yr2[i][12]
	redchisq_low_yr1[i] = params_yr1[i][13]
	redchisq_low_yr2[i] = params_yr2[i][13]
	redchisq_high_yr1[i] = params_yr1[i][14]
	redchisq_high_yr2[i] = params_yr2[i][14]
	redchisq_quad_yr1[i] = params_yr1[i][15]
	redchisq_quad_yr2[i] = params_yr2[i][15]

invalid_turnovers_yr1 = np.isnan(quad_turnover_yr1)
snr_curve_cut_yr1 = calc_curve_snr(freq_mwa, flux_mwa_yr1, local_rms_yr1, quad_turnover_yr1, invalid_turnovers_yr1, quad_curve_error_yr1)
invalid_turnovers_yr2 = np.isnan(quad_turnover_yr2)
snr_curve_cut_yr2 = calc_curve_snr(freq_mwa, flux_mwa_yr2, local_rms_yr2, quad_turnover_yr2, invalid_turnovers_yr2, quad_curve_error_yr2)

diff = flux_yr1-flux_yr2
diff_err = np.sqrt((flux_err_yr1)**2+(flux_err_yr2)**2) 

var_param = np.zeros(len(flux_yr1[0]))

for i in range(len(flux_yr1[0])):
    var_param[i]=(np.sum(((flux_yr1[:,i]-flux_yr2[:,i])**2)/(diff_err[:,i])**2))      # remove (1/16.) if you want normal chi2


all_pop = pd.read_csv('/data/var_analysis/analysis/all_pop_yr1yr2_params_jps.csv')
hl_pop = pd.read_csv('/data/var_analysis/analysis/hl_pop_yr1yr2.csv')
gps_pop = pd.read_csv('/data/var_analysis/analysis/gps_pop_yr1yr2.csv')

name_norm = np.array(all_pop['Name'])
name_gps = np.array(gps_pop['Name'])
name_hl = np.array(hl_pop['Name'])
gps_id = np.zeros(len(name_norm))


for i in range(len(name_norm)):
	if name_norm[i] in name_gps:
		print('I AM A GPS: '+name_norm[i])
		gps_id[i] = 1
	elif name_norm[i] in name_hl:
		print('I AM HL: '+name_norm[i])	
		gps_id[i] = 1
	else: 
		pass 

col1 = fits.Column(name='Name', format = '24A', array = name)
col2 = fits.Column(name='RA', format = 'E', array = ra_deg)
col3 = fits.Column(name='Dec', format = 'E', array = dec_deg)
col4 = fits.Column(name='alpha_low_yr1', format = 'E', array = alpha_low_yr1)
col5 = fits.Column(name='alpha_low_yr2', format = 'E', array = alpha_low_yr2)
col6 = fits.Column(name='alpha_high_yr1', format = 'E', array = alpha_high_yr1)
col7 = fits.Column(name='alpha_high_yr2', format = 'E', array = alpha_high_yr2)
col8 = fits.Column(name='norm_low_yr1', format = 'E', array = norm_low_yr1)
col9 = fits.Column(name='norm_low_yr2', format = 'E', array = norm_low_yr2)
col10 = fits.Column(name='norm_high_yr1', format = 'E', array = norm_high_yr1)
col11 = fits.Column(name='norm_high_yr2', format = 'E', array = norm_high_yr2)
col12 = fits.Column(name='quad_curve_yr1', format = 'E', array = quad_curve_yr1)
col13 = fits.Column(name='quad_curve_yr2', format = 'E', array = quad_curve_yr2)
col14 = fits.Column(name='alpha_quad_yr1', format = 'E', array = alpha_quad_yr1)
col15 = fits.Column(name='alpha_quad_yr2', format = 'E', array = alpha_quad_yr2)
col16 = fits.Column(name='quad_turnover_yr1', format = 'E', array = quad_turnover_yr1)
col17 = fits.Column(name='quad_turnover_yr2', format = 'E', array = quad_turnover_yr2)
col18 = fits.Column(name='norm_quad_yr1', format = 'E', array = norm_quad_yr1)
col19 = fits.Column(name='norm_quad_yr2', format = 'E', array = norm_quad_yr2)
col20 = fits.Column(name='alpha_low_error_yr1', format = 'E', array = alpha_low_error_yr1)
col21 = fits.Column(name='alpha_low_error_yr2', format = 'E', array = alpha_low_error_yr2)
col22 = fits.Column(name='alpha_high_error_yr1', format = 'E', array = alpha_high_error_yr1)
col23 = fits.Column(name='alpha_high_error_yr2', format = 'E', array = alpha_high_error_yr2)
col24 = fits.Column(name='quad_curve_error_yr1', format = 'E', array = quad_curve_error_yr1)
col25 = fits.Column(name='quad_curve_error_yr2', format = 'E', array = quad_curve_error_yr2)
col26 = fits.Column(name='alpha_quad_error_yr1', format = 'E', array = alpha_quad_error_yr1)
col27 = fits.Column(name='alpha_quad_error_yr2', format = 'E', array = alpha_quad_error_yr2)
col28 = fits.Column(name='quad_turnover_error_yr1', format = 'E', array = quad_turnover_error_yr1)
col29 = fits.Column(name='quad_turnover_error_yr2', format = 'E', array = quad_turnover_error_yr2)
col30 = fits.Column(name='snr_curve_cut_yr1', format = 'E', array = snr_curve_cut_yr1)
col31 = fits.Column(name='snr_curve_cut_yr2', format = 'E', array = snr_curve_cut_yr2)
col32 = fits.Column(name = 'redchisq_low_yr1', format='E', array = redchisq_low_yr1)
col33 = fits.Column(name = 'redchisq_low_yr2', format='E', array = redchisq_low_yr2)
col34 = fits.Column(name = 'redchisq_high_yr1', format='E', array = redchisq_high_yr1)
col35 = fits.Column(name = 'redchisq_high_yr2', format='E', array = redchisq_high_yr2)
col36 = fits.Column(name = 'redchisq_quad_yr1', format='E', array = redchisq_quad_yr1)
col37 = fits.Column(name = 'redchisq_quad_yr2', format='E', array = redchisq_quad_yr2)
col38 = fits.Column(name = 'var_param', format='E', array = var_param)
col39 = fits.Column(name = 'gps_id', format = 'E', array= gps_id)
cols = fits.ColDefs([col1, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  
tbhdu.writeto('/data/var_analysis/analysis/params_yr1_yr2.fits', overwrite = True)



