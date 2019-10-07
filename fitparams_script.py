# This code will calculate the plotting parameters for all sources in the given field for all the different models given by J.Callingham and then save them to a new catalogue 
print("#----------------------------------------------------------#")

import numpy as np
import scipy.optimize as opt
from astropy.io import fits
import gpscssmodels
import os
import time
import read_variables_forfitting as reading

# These are models you will use: FUTURE NOTE GET RID OF UNUSED ONES AND ALSO MOVE TO THEIR OWN SCRIPT CAUSE GROSS 

def redchisq(ydata,ymod,sd,deg):
    
    chisq=np.sum(((ydata-ymod)/sd)**2)

    nu=ydata.size-1-deg

    return [chisq, chisq/nu]

def quad(freq, a, b, c): # defining quadratic
    return c*np.power(np.log(freq),2) + b*np.log(freq) + a

def quad_plot(freq, a, b, c): # plotting quadratic in linear space
    return np.exp(c*np.power(np.log(freq),2) + b*np.log(freq) + a)

def curve(freq, freq_peak, peakfreq, alphathick, alphathin): # Model taken from Tschager et al. 2003. General fit not based on any physics.
    return freq_peak/(1 -np.exp(-1))*((freq/freq_peak)**alphathick)*(1 - np.exp(-(freq/freq_peak)**(alphathin-alphathick)))


def powlaw(freq,a,alpha): # defining powlaw as S = a*nu^-alpha. Important to have x value first in definition of function.
    return a*(freq**(-alpha))

# fit curvature models with rms cuts: use with literally everything flux and local_rms
def spectral_index_eval_curve(freq,flux,flux_err, local_rms):
    fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data
    if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
        return 'No flux?'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]
    # p0curve = [peakflux, peakfreq, -0.5, 0.7]
    p0curve = [peakflux, 0.5, 0.7, 300.]
    # p0quad = [0.6,-0.001,]
    p0gen = [peakflux,300.,0.8,-0.7]

    try:
        poptcurve, pcovcurve = opt.curve_fit(gpscssmodels.singinhomobremss, freqplot, fluxplot, p0 = p0curve, sigma = flux_errplot, maxfev = 10000)
        poptgen, pcovgen = opt.curve_fit(gpscssmodels.curve, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
        poptquad, pcovquad = opt.curve_fit(quad, freqplot, fluxplot, sigma = flux_errplot, maxfev = 10000)
        poptpowlaw_tot, pcovpowlaw_tot = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
        redchisq_curve = redchisq(fluxplot,gpscssmodels.singinhomobremss(freqplot,*poptcurve),flux_errplot,4)[1]
        redchisq_gen = redchisq(fluxplot,gpscssmodels.curve(freqplot,*poptgen),flux_errplot,4)[1]
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    return (poptcurve, pcovcurve, poptquad, pcovquad, poptpowlaw_tot, pcovpowlaw_tot, fluxplot, freqplot, flux_errplot, redchisq_curve, poptgen, pcovgen, redchisq_gen)


def spectral_index_eval_mwa(freq,flux,flux_err, local_rms_low):
    fluxposind = np.where((flux > 0))
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        return 'Curve_fit could not fit curve model.'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]
    p0quad = [peakflux,0.5,0.]

    try:
        poptquad, pcovquad = opt.curve_fit(quad, freqplot, np.log(fluxplot), p0 = p0quad, sigma = flux_errplot/fluxplot, maxfev = 10000)
        poptpowlaw_gleam, pcovpowlaw_gleam = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    redchisq_pow_gleam = redchisq(fluxplot,gpscssmodels.powlaw(freqplot,*poptpowlaw_gleam),flux_errplot,2)[1]
    redchisq_quad_gleam = redchisq(np.log(fluxplot),quad(freqplot,*poptquad),flux_errplot/fluxplot,3)[1]
    return (poptpowlaw_gleam, pcovpowlaw_gleam,poptquad, pcovquad, fluxplot, freqplot, flux_errplot,redchisq_pow_gleam,redchisq_quad_gleam, len(fluxposind[0]))


# Use this to fit the low 
def spectral_index_eval_rms(freq,flux,flux_err,local_rms_low):
    fluxposind = np.where((flux > 0) )
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        return 'Curve_fit could not fit powerlaw.'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]

    try:
        poptpowlaw, pcovpowlaw = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
        redchisq_pow = redchisq(fluxplot,gpscssmodels.powlaw(freqplot,*poptpowlaw),flux_errplot,2)[1]
    except (RuntimeError, TypeError): # If curve fit can not find a good fit, skip the source.    
        return 'Curve_fit could not fit powerlaw.'

    return (poptpowlaw, pcovpowlaw, fluxplot, freqplot, flux_errplot)


# Use this for the flux_high fits 
def spectral_index_eval_norms(freq,flux,flux_err):
    fluxposind = np.where((flux > 0))
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        return 'Curve_fit could not fit powerlaw.'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]

    try:
        poptpowlaw, pcovpowlaw = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
        redchisq_pow = redchisq(fluxplot,gpscssmodels.powlaw(freqplot,*poptpowlaw),flux_errplot,2)[1]
    except (RuntimeError, TypeError): # If curve fit can not find a good fit, skip the source.    
        return 'Curve_fit could not fit powerlaw.'

    return (poptpowlaw, pcovpowlaw, fluxplot, freqplot, flux_errplot)

# fit power-law without rms cuts
def fit_powlaw(freq,flux,flux_err):
    fluxposind = np.where((flux > 0))
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        return 'Curve_fit could not fit powerlaw.'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]

    try:
        poptpowlaw, pcovpowlaw = opt.curve_fit(powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
        redchisq_pow = redchisq(fluxplot,gpscssmodels.powlaw(freqplot,*poptpowlaw),flux_errplot,2)[1]
    except (RuntimeError, TypeError): # If curve fit can not find a good fit, skip the source.    
        return 'Curve_fit could not fit powerlaw.'

    return poptpowlaw, pcovpowlaw, fluxplot, freqplot, flux_errplot, redchisq_pow

# fit curvature models based on some physics NOTE ADD MORE HERE BUT USE WHEN MADE SOME CUTS NOT FOR ALL 
def fit_physical_curves(freq,flux,flux_err):
    fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data
    if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
        return 'No flux?'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]
    # p0curve = [peakflux, peakfreq, -0.5, 0.7]
    p0curve = [peakflux, 0.5, 0.7,]
    p0quad = [peakflux,-0.5,0.]
    p0gen = [peakflux,-0.8,0.]

    try:
    	#This is the actual curve fitting section!! Each one has a different curve they're trying to fit
        poptquad, pcovquad = opt.curve_fit(gpscssmodels.quad, freqplot, np.log(fluxplot), p0=p0quad, sigma = flux_errplot, maxfev = 10000)
        redchisq_quad = redchisq(np.log(fluxplot),quad(freqplot,*poptquad),flux_errplot/fluxplot,3)[1]
        # poptgen, pcovgen = opt.curve_fit(gpscssmodels.curve, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
        # poptcurve, pcovcurve = opt.curve_fit(gpscssmodels.curvepowlaw, freqplot, fluxplot, sigma = flux_errplot, maxfev = 10000, p0 = p0gen)
        # poptpowlaw_tot, pcovpowlaw_tot = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
        # redchisq_curve = redchisq(fluxplot,gpscssmodels.singinhomobremss(freqplot,*poptcurve),flux_errplot,4)[1]
        # redchisq_gen = redchisq(fluxplot,gpscssmodels.curve(freqplot,*poptgen),flux_errplot,4)[1]
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    return poptquad, pcovquad, redchisq_quad

# just fit to the MWA data, fits both a generic curve and powerlaw 
def fit_mwa(freq,flux,flux_err):
    fluxposind = np.where((flux > 0))
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        return 'Curve_fit could not fit curve model.'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]
    p0quad = [peakflux,0.5,0.]

    try:
        poptquad, pcovquad = opt.curve_fit(quad, freqplot, np.log(fluxplot), p0 = p0quad, sigma = flux_errplot/fluxplot, maxfev = 10000)
        poptpowlaw_gleam, pcovpowlaw_gleam = opt.curve_fit(powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    redchisq_pow_gleam = redchisq(fluxplot,gpscssmodels.powlaw(freqplot,*poptpowlaw_gleam),flux_errplot,2)[1]
    redchisq_quad_gleam = redchisq(np.log(fluxplot),quad(freqplot,*poptquad),flux_errplot/fluxplot,3)[1]
    return (poptpowlaw_gleam, pcovpowlaw_gleam,poptquad, pcovquad, fluxplot, freqplot, flux_errplot,redchisq_pow_gleam,redchisq_quad_gleam, len(fluxposind[0]))

def calc_curve_snr(freq_mwa, flux_mwa, local_rms, turnover_freq, invalid_turnovers, curve_error):
	freq_mwa = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
	snr_curve_limit = np.zeros(len(curve_error))
	snr_cut = np.zeros(len(curve_error))
	for i in range(len(snr_curve_limit)):
		if invalid_turnovers[i] == True:
			print("Im passing, invalid turnover") 
		else: 
		# Using i=7335 because thats the source in questio but change to specific pop and then general
			fluxpos_ind = np.where((flux_mwa[:,i]>0) & (flux_mwa[:,i]/local_rms[:,i]>=3))[0]
			if len(fluxpos_ind)==0:
				print('Passing because no flux???')
			else:
				# print(fluxpos_ind)
				# Only selecting the ones that are bright enough for testing 
				freq_src = freq_mwa[fluxpos_ind]
				flux_src = flux_mwa[:,i][fluxpos_ind]
				local_rms_src = local_rms[:,i][fluxpos_ind]
				turnover_freq_src = turnover_freq[i]
				# print(turnover_freq_src)
				# Calculating the peak and the points above/below 
				peak_point_ind = np.where(abs(turnover_freq_src-freq_src) == min(abs(turnover_freq_src - freq_src)))[0][0]
				lowest_freq_ind = np.where(freq_src == min(freq_src))[0][0]
				highest_freq_ind = np.where(freq_src == max(freq_src))[0][0]
				centre_freq_point = np.where(abs(freq_src - 151.) == min(abs(freq_src - 151.)))[0][0]
				# Here is the actual calculation
				if peak_point_ind >= centre_freq_point:
					if len(local_rms_src[peak_point_ind:highest_freq_ind]) == 0:
						local_rms_src_sum = local_rms_src[peak_point_ind]
						flux_src_sum = flux_src[peak_point_ind]
					else:
						local_rms_src_sum = local_rms_src[peak_point_ind:highest_freq_ind]
						flux_src_sum = flux_src[peak_point_ind:highest_freq_ind]
				if peak_point_ind < centre_freq_point:
					if len(local_rms_src[lowest_freq_ind:peak_point_ind]) == 0:
						local_rms_src_sum = local_rms_src[lowest_freq_ind]
						flux_src_sum = flux_src[lowest_freq_ind]
					else:
						local_rms_src_sum = local_rms_src[lowest_freq_ind:peak_point_ind]
						flux_src_sum = flux_src[lowest_freq_ind:peak_point_ind]
				snr_curve_limit[i] = -np.sqrt(np.sum(local_rms_src_sum**2)) / np.sum(flux_src_sum) + 0.2
				if curve_error[i]<=snr_curve_limit[i]:
					snr_cut[i] = 1
					# print("Sweet success")
			#         print(q_curve_error[i],snr_curve_limit[i])
				else: 
					# print("Parting is such sweet sorrow")
					continue
	return snr_cut 




def fitting_params(name, freq_low, freq_high, flux_high, flux_high_err, flux_low, flux_low_err, freq_mwa, flux_mwa, flux_mwa_err):
	print("Processing: "+name)
	# name, ra_deg, dec_deg, flux, flux_err, flux_mwa, flux_mwa_err, flux_low, flux_err_low, flux_high, flux_err_high, S_white, local_rms = reading.read_variables('/data/gleam/analysis/variability/code/all_pop.fits','_yr1')

	# invalid_turnovers_yr1 = np.isnan(quad_turnover)
	# snr_curve_cut_yr1 = calc_curve_snr(flux_mwa, local_rms, quad_turnover, invalid_turnovers, quad_curve_error)
	# FIRST FITTING THE POWER LAWS BOTH LOW AND HGIH 

	# high
	if fit_powlaw(freq_high,flux_high,flux_high_err) == 'Curve_fit could not fit powerlaw.' or fit_powlaw(freq_high,flux_high,flux_high_err) == 'No flux?':
		print('Could not fit high powerlaw to '+ name)
		alpha_high = np.nan
		norm_high = np.nan 
		alpha_high_error = np.nan
		redchisq_pow_high = np.nan
		pass
	else:    
		# print('Processing '+name)
		# popt = params of fit, pcov = covariance matrix of popt, others are just plotting stuff don't really need 
		poptpowlaw_high, pcovpowlaw_high, fluxplot_high, freqplot_high, flux_errplot_high, redchisq_pow_high = fit_powlaw(freq_high,flux_high,flux_high_err)
		# Note formula is: S=a*nu^-alpha
		alpha_high=poptpowlaw_high[1]
		norm_high=poptpowlaw_high[0]
		alpha_high_error = np.sqrt(pcovpowlaw_high[1][1])

	# low
	if fit_powlaw(freq_low,flux_low,flux_low_err) == 'Curve_fit could not fit powerlaw.' or fit_powlaw(freq_low,flux_low,flux_low_err) == 'No flux?':
		print('Could not fit high powerlaw to '+ name)
		alpha_low = np.nan
		norm_low = np.nan 
		alpha_low_error = np.nan
		redchisq_pow_low = np.nan
		pass
	else:
		# print('Processing '+name)
		# popt = params of fit, pcov = covariance matrix of popt, others are just plotting stuff don't really need but nice for later maybe 
		poptpowlaw_low, pcovpowlaw_low, fluxplot_low, freqplot_low, flux_errplot_low, redchisq_pow_low = fit_powlaw(freq_low,flux_low,flux_low_err)
		# Note formula is: S=a*nu^-alpha
		alpha_low=-poptpowlaw_low[1]
		norm_low=poptpowlaw_low[0]
		alpha_low_error = np.sqrt(pcovpowlaw_low[1][1])       


	# if spectral_index_eval_curve(freq,flux,flux_err,local_rms) == 'Curve_fit could not fit curve model.' or spectral_index_eval_curve(freq,flux,flux_err,local_rms) == 'No flux?':
	# 	print('Could not fit curve to '+ Name_bright[i])
	# 	poptcurve = np.array([np.nan,np.nan,np.nan,np.nan])
	# 	continue
	# else:    
	# 	poptcurve, pcovcurve, poptquad, pcovquad,poptpowlaw_tot, pcovpowlaw_tot, fluxplot, freqplot, flux_errplot,redchisq_curve, poptgen, pcovgen, redchisq_gen = spectral_index_eval_curve(freq,flux,flux_err,local_rms)
	# 	alpha_tot[i] = poptpowlaw_tot[1]
	# 	alpha_err_tot[i] = np.sqrt(pcovpowlaw_tot[1][1])
	# 	redchisq_curve_store[i] = redchisq_curve
	# 	redchisq_gen_store[i] = redchisq_gen

	#just fitting MWA band
	# if spectral_index_eval_mwa(freq_mwa,flux_mwa,flux_mwa_err,local_rms) == 'Curve_fit could not fit curve model.' or spectral_index_eval_mwa(freq_mwa,flux_mwa,flux_mwa_err,local_rms) == 'No flux?':
	# 	# print('Could not fit MWA data to '+ Name_bright[i])
	# 	# with open("quality_control/MWA_flux_not_fit_snrcut_less8pts.txt", "a") as myfile:
	# 	#     myfile.write(Name_bright[i]+'\n')
	# 	pass

	# 	# sedplots_gleam.sed([gpscssmodels.powlaw,gpscssmodels.powlaw],[poptpowlaw_high,poptpowlaw_low],np.concatenate((freqplot_low,freqplot_high)),np.concatenate((fluxplot_low,fluxplot_high)),np.concatenate((flux_errplot_low,flux_errplot_high)), freq_labels = True, savefig = True, title = Name_bright[i] )#+r' $\alpha$ = '+str(round(alpha[i],3))+r' Powlaw $ \chi_{\rm{red}}$  = '+str(round(redchisq_pow[i],2)))
	# else:
	# 	print('Processing '+name) 
	# 	poptpowlaw_gleam, pcovpowlaw_gleam, poptquad_gleam, pcovquad_gleam, fluxplot, freqplot, flux_errplot,redchisq_pow_gleam,redchisq_quad_gleam, number_flux_points = spectral_index_eval_mwa(freq_mwa,flux_mwa,flux_mwa_err,local_rms)
	# 	quad_curve_error = np.sqrt(pcovquad_gleam[2][2])
	# 	alpha_quad_error = np.sqrt(pcovquad_gleam[1][1])
	# 	alpha_quad = poptquad_gleam[1]
	# 	quad_curve = poptquad_gleam[2]
	# 	quad_turnover = np.exp(-alpha_quad/(2*quad_curve))
	# 	quad_turnover_error = quad_turnover * np.sqrt((alpha_quad_error / alpha_quad)**2 + (quad_curve_error / quad_curve)**2)
		# THIS HAS FIT MWA FREQUENCIES WITH GENERIC CURVE TO GET CURVATURE AS WELL AS A POWERLAW TO GET ALPHA KEY FOR NON GPS NORM POP 



	if fit_physical_curves(freq_mwa,flux_mwa,flux_mwa_err) == 'Curve_fit could not fit curve model.' or fit_physical_curves(freq_mwa,flux_mwa,flux_mwa_err) == 'No flux?':
		print('Could not fit MWA data to '+ name)
		alpha_low = np.nan
		alpha_low_error = np.nan
		alpha_high = np.nan
		alpha_high_error = np.nan
		norm_low = np.nan
		norm_high = np.nan 
		quad_curve = np.nan
		quad_curve_error = np.nan
		alpha_quad = np.nan
		alpha_quad_error = np.nan
		quad_turnover = np.nan
		quad_turnover_error = np.nan
		norm_quad = np.nan
		redchisq_quad = np.nan
		pass
	else:
		# print('Processing '+name) 
		poptquad, pcovquad, redchisq_quad = fit_physical_curves(freq_mwa,flux_mwa,flux_mwa_err)
		norm_quad = poptquad[0]
		quad_curve_error = np.sqrt(pcovquad[2][2])
		alpha_quad_error = np.sqrt(pcovquad[1][1])
		alpha_quad = poptquad[1]
		quad_curve = poptquad[2]
		quad_turnover = np.exp(-alpha_quad/(2*quad_curve))
		quad_turnover_error = quad_turnover * np.sqrt((alpha_quad_error / alpha_quad)**2 + (quad_curve_error / quad_curve)**2)

	# if alpha_low <= -1.5 and alpha_high >1: 
	# 	print('Plotting the SED')
	# 	name_plt = name 
	# 	paras=[poptpowlaw_low, poptpowlaw_high, poptquad]
	# 	seds_plot_func.sed(models,directory, paras,freq,flux,flux_err, name_plt, #alpha_low, alpha_high,
	# 	grid = False, freq_labels = False, log = True, bayes = False, resid = False, savefig=True)

	# print(alpha_high)
	# name, ra_deg, dec_deg, flux_yr2, flux_err_yr2, flux_mwa_yr2, flux_mwa_err_yr2, flux_low_yr2, flux_err_low_yr2, flux_high_yr2, flux_err_high_yr2, S_white_yr2, local_rms_yr2 = reading.read_variables('/data/gleam/analysis/variability/code/all_pop.fits','_yr2')



	# invalid_turnovers_yr2 = np.isnan(quad_turnover_yr2)
	# snr_curve_cut_yr2 = calc_curve_snr(flux_mwa_yr2, local_rms_yr2, quad_turnover_yr2, invalid_turnovers_yr2, quad_curve_error_yr2)
	return alpha_low, alpha_low_error, alpha_high, alpha_high_error, norm_low, norm_high, quad_curve, quad_curve_error, alpha_quad, alpha_quad_error, quad_turnover, quad_turnover_error, norm_quad, redchisq_pow_low, redchisq_pow_high, redchisq_quad #, alpha_curve, q_curve, norm_curve_err, alpha_curve_err, q_curve_err, turnover_quad, turnover_quad_error




