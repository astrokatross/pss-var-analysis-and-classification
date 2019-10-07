# Script function for plotting SEDs, GENERAL!
# K. Ross
# 21/06/19

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
from collections import OrderedDict
from matplotlib.ticker import StrMethodFormatter
import os
import re
from matplotlib import rc # To plot labels in serif rather than default.
from matplotlib.ticker import MaxNLocator
plt.rcParams["font.family"] = "serif"

def ticks_format_x(value, index):
	"""
	get the value and returns the value as:
		integer: [0,99]
		1 digit float: [0.1, 0.99]
		n*10^m: otherwise
	To have all the number of the same size they are all returned as latex strings
	"""
	if value in np.array([0.03,0.05,0.07,0.09,0.3,0.5,0.7,0.9,3,5,7,9,13,15,17,19,23,25,27,29]): # This will remove 90 and 900 MHz, replace number for anything you don't want to appear.
	# print value
	# print np.where(value == 0.5)
		return ''
	exp = np.floor(np.log10(value))
	base = value/10**exp
	if value == 900.0 or value == 90. or value == 700.: # This will remove 90 and 900 MHz, replace number for anything you don't want to appear.
		return ''
	if exp == 0 or exp == 1 or exp == 2 or exp ==3:   
		return '${0:d}$'.format(int(value))
	if exp == -1:
		return '${0:.1f}$'.format(value)
	else:
		return '${0:d}\\times10^{{{1:d}}}$'.format(int(base), int(exp))


# rc('text', usetex=True)
# rc('font',**{'family':'serif','serif':['serif'],'size':18})

cdict_models = {'singhomobremss':'red',
        'singhomobremsscurve':'maroon',
        'singhomobremssbreak': 'orangered',
        'singinhomobremss':'k',
        'singinhomobremsscurve':'#4682b4',
        'doubhomobremss':'saddlebrown',
        'doubhomobremsscurve':'dodgerblue',
        'doubhomobremssbreak':'olive',
        'doubhomobremssbreak':'DarkGoldenRod',
        'singSSA':'orchid',
        'singSSAcurve':'darkmagenta',
        'singSSAbreak':'indigo',
        'doubSSA':'navy',
        'doubSSAcurve':'sienna',
        'doubSSAbreak':'black',
        'powlaw': 'DarkOrange',
        'powlaw2': 'purple',
        'powlawbreak':'Chocolate',
        'internalbremss':'MediumSeaGreen',
        'curve':'k'
            } 

# Defining plotting routine

def sed(directory, models_yr1, models_yr2,paras_yr1,paras_yr2,freq, flux_yr1,flux_err_yr1,flux_yr2,flux_err_yr2, 
		name, var_param, q_value, gps_id, freq_xtra, flux_xtra,flux_err_xtra, freq_gleam, flux_gleam, flux_gleam_err, #alpha_low, alpha_high,
        grid = False, freq_labels = False, log = True, bayes = False, resid = True, savefig=False, xtra_high=False, xtra_low=False):


	# Ensuring that freq and flux are approability matched up.
	# ind = np.argsort(freq)
	# freq = freq[ind]
	# flux_yr1 = flux_yr1[ind]
	# flux_err_yr1 = flux_err_yr1[ind]
	# flux_yr2 = flux_yr2[ind]
	# flux_err_yr2 = flux_err_yr2[ind]
	max_freq = freq[-1]+100
	if resid == True:
		fig = plt.figure(1,figsize=(15, 10)) #(12,8)
		gs = plt.GridSpec(2,1, height_ratios = [3,1], hspace = 0)
		ax = plt.subplot(gs[0])
		ax1 = plt.subplot(gs[1])
		ax1.set_xlabel('Frequency (MHz)', fontsize = 20)
		ax1.set_ylabel('Residuals', fontsize = 30)
		ax1.xaxis.labelpad = 15
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(2)
			ax1.spines[axis].set_linewidth(2)
	else:
		fig = plt.figure(1,figsize=(15, 10))#(12, 8))
		gs = plt.GridSpec(1,1)
		ax = plt.subplot(gs[0])
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(2)
	freq_cont = np.array(list(range(60,1500)))
	ax.set_xlim(70.,1500)
	# if max_freq <= 1500.:
	# 	ax.set_xlim(70., 1500.)
	# 	print('xlims')
	# else:
	# 	ax.set_xlim(70., 21000.)


	# plt_flux_xtra = flux_xtra[flux_xtra !=0]
	# plt_flux_gleam =list(flux_gleam)

	# plt_flux_xtra = list(plt_flux_xtra)
	# plt_flux_yr1 = list(flux_yr1)
	# plt_flux_yr2 = list(flux_yr2)

	# plt_flux = plt_flux_yr1+plt_flux_yr2+plt_flux_xtra+plt_flux_gleam
	# plt_flux = np.array(plt_flux)
	# plt_flux = plt_flux[np.where(plt_flux>=0)]

	# ymin = min(plt_flux)-0.2*min(plt_flux)
	# ymax = max(plt_flux)+0.2*max(plt_flux)
	# print(ymin,ymax)
	# ax.set_ylim(ymin,ymax)


	if xtra_high == True and xtra_low == False:
		plt_flux_xtra = flux_xtra[flux_xtra !=0]
		plt_flux_xtra = list(plt_flux_xtra)
		plt_flux_yr1 = list(flux_yr1)
		plt_flux_yr2 = list(flux_yr2)
		plt_flux = plt_flux_yr1+plt_flux_yr2+plt_flux_xtra
		plt_flux = np.array(plt_flux)
		plt_flux = plt_flux[np.where(plt_flux>=0)]
		ymin = min(plt_flux)-0.2*min(plt_flux)
		ymax = max(plt_flux)+0.2*max(plt_flux)
		ax.set_ylim(ymin,ymax)
	elif xtra_high == True and xtra_low ==True: 
		plt_flux_xtra = flux_xtra[flux_xtra !=0]
		plt_flux_xtra = list(plt_flux_xtra)
		plt_flux_yr1 = list(flux_yr1)
		plt_flux_yr2 = list(flux_yr2)
		plt_flux_low = flux_gleam[flux_gleam != 0]
		plt_flux_low = list(plt_flux_low)
		plt_flux = plt_flux_yr1+plt_flux_yr2+plt_flux_xtra+plt_flux_low
		plt_flux = np.array(plt_flux)
		plt_flux = plt_flux[np.where(plt_flux>=0)]
		ymin = min(plt_flux)-0.2*min(plt_flux)
		ymax = max(plt_flux)+0.2*max(plt_flux)
		ax.set_ylim(ymin,ymax)
	elif xtra_high == False and xtra_low == True:
		plt_flux_xtra = flux_gleam[flux_gleam !=0]
		plt_flux_xtra = list(plt_flux_xtra)
		plt_flux_yr1 = list(flux_yr1)
		plt_flux_yr2 = list(flux_yr2)
		plt_flux = plt_flux_yr1+plt_flux_yr2+plt_flux_xtra
		plt_flux = np.array(plt_flux)
		plt_flux = plt_flux[np.where(plt_flux>=0)]
		ymin = min(plt_flux)-0.2*min(plt_flux)
		ymax = max(plt_flux)+0.2*max(plt_flux)
		ax.set_ylim(ymin,ymax)
	else:
		plt_flux_yr1 = list(flux_yr1)
		plt_flux_yr2 = list(flux_yr2)
		plt_flux = plt_flux_yr1+plt_flux_yr2
		ymin = min(plt_flux)-0.2*min(plt_flux)
		ymax = max(plt_flux)+0.2*max(plt_flux)
		ax.set_ylim(ymin,ymax)

	ax.set_xlabel(r'Frequency (MHz)', fontsize = 50,labelpad=10)
	ax.set_ylabel(r'Flux Density (Jy)', fontsize = 35)
	fig.suptitle(name, fontsize = 35, fontweight='bold')
	ax.set_title(r'Variability: '+ var_param + r', Curvature: '+q_value+ r', GPS_id: '+gps_id, fontsize=25)
	try:
		tt = len(models_yr1)
	except TypeError:
		tt = 1
		models_yr1 = [models_yr1]
		models_yr2 = [models_yr2]
		paras_yr1 = [paras_yr1]
		paras_yr2 = [paras_yr2]
	for i in range(tt):
		# Defining colours for models to make it easy to tell the difference
		try:
			color = cdict_models[models_yr1[i].__name__] # In case none of the models are listed here.        
		except KeyError:
			# print 'Model is not in colour dictionary. Defaulting to dark orange.'
			color = 'darkgreen'
		# print(models[i])
		yvals = models_yr1[i](freq_cont, *paras_yr1[i])
		ax.plot(freq_cont, yvals, color = 'C6', linestyle='-', lw=2)
		ax.plot(freq_cont, models_yr2[i](freq_cont, *paras_yr2[i]), color = 'purple', linestyle='-', lw=2)
		# ax.legend(loc='upper right')
		if resid == True:
			model_points_1 = models_yr1[0](freq,*paras_yr1[0])
			model_points_2 = models_yr2[0](freq,*paras_yr2[0])

			resid_yr1 = flux_yr1 - model_points_1
			resid_yr2 = flux_yr2 - model_points_2


			plt_resid_1 = resid_yr1/flux_err_yr1
			plt_resid_2 = resid_yr2/flux_err_yr2



			diff = flux_yr1-flux_yr2
			mean_resid = np.mean(diff)
			label_mean =format(mean_resid,'.6f')

			diff_err = np.sqrt((flux_err_yr1)**2+(flux_err_yr2)**2) 
			compx = np.linspace(min(freq_cont)-0.1*min(freq_cont),max(freq_cont)+0.1*max(freq_cont))
			compy = np.zeros(len(compx))
			

			ax1.plot(compx,compy,linestyle = '--',color = 'gray',linewidth = 2)
			ax1.set_xlim(min(freq_cont), max(freq_cont))
			for i in range(len(freq)):
					(_, caps, _) = ax1.errorbar(freq[i], diff[i], diff_err[i],marker = 'o',elinewidth=2,capsize=3,  markersize=5, color = 'k', linestyle='none')
					for cap in caps:
						cap.set_color('k') 
						cap.set_markeredgewidth(2) 

			# ax1.legend(loc='upper right')                     
		if log == True:
			ax.set_xscale('log')
			ax.set_yscale('log')

			subsx = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # ticks to show per decade
			ax.xaxis.set_minor_locator(ticker.LogLocator(subs=subsx)) #set the ticks position
			ax.xaxis.set_major_formatter(ticker.NullFormatter())   # remove the major ticks
			ax.xaxis.set_minor_formatter(ticker.FuncFormatter(ticks_format_x))

			# HERE IS THE Y-AXIS TICKS SETTINGS 
			# ax.set_ylim(ymin, ymax)
			locator=MaxNLocator(prune='both',nbins=8)
			ax.yaxis.set_minor_locator(locator)
			ax.yaxis.set_major_formatter(ticker.NullFormatter())   # remove the major ticks
			ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())

			if resid == True:
				ax.set_xticklabels('',minor = True)
				ax.set_xlabel('')
				ax1.set_xscale('log')
				ax1.xaxis.set_minor_locator(ticker.LogLocator(subs=subsx)) #set the ticks position
				ax1.xaxis.set_major_formatter(ticker.NullFormatter())   # remove the major ticks
				ax1.xaxis.set_minor_formatter(ticker.FuncFormatter(ticks_format_x))
				ax1.tick_params(axis='both',which='both',labelsize=22)
				ax1.tick_params(axis='both',which='major',length=8,width=1.5)
				ax1.tick_params(axis='both',which='minor',length=5,width=1.5)

	if freq_labels == True:
# COULD MAKE MORE GENERAL TO BE ABLE TO PLOT A GENERAL NUMBER OF EPOCHS i.e. ADD ANOTHER FOR I IN RANGE LENGTH OF EPOCHS AND THEN CHANGE FLUX TO BE ARRAAY OR ARRAYS 
		for i in range(len(freq)):
			if freq[i] in [76,84,92,99,115,123,130,143,151,158,166,174,181,189,197,204,212,219,227]:
				(_, caps, _) = ax.errorbar(freq[i], flux_yr1[i], flux_err_yr1[i],marker = 'o',capsize=3,elinewidth=1.5, markersize=5,  color = 'C6', linestyle='none',markeredgecolor='none')
				(_, caps, _) = ax.errorbar(freq[i], flux_yr2[i], flux_err_yr2[i],marker = 'o',capsize=2.8,elinewidth=1.5, markersize=5,  color = 'purple', linestyle='none',markeredgecolor='none')
				compx = np.linspace(min(freq)-0.1*min(freq),max(freq)+0.1*max(freq))

				ax.set_xlim(min(freq_cont), max(freq_cont))
				for cap in caps:
					# cap.set_color('pink') 
					cap.set_markeredgewidth(2)
			if freq[i] in [107]:
				(_, caps, _) = ax.errorbar(freq[i], flux_yr1[i], flux_err_yr1[i],marker = 'o',capsize=3,elinewidth=1.5, markersize=5,  color = 'C6', linestyle='none',markeredgecolor='none', label='Year 1')
				(_, caps, _) = ax.errorbar(freq[i], flux_yr2[i], flux_err_yr2[i],marker = 'o',capsize=2.8,elinewidth=1.5, markersize=5,  color = 'purple', linestyle='none',markeredgecolor='none', label='Year 2')
				compx = np.linspace(min(freq)-0.1*min(freq),max(freq)+0.1*max(freq))

				ax.set_xlim(min(freq_cont), max(freq_cont))
				for cap in caps:
					# cap.set_color('pink') 
					cap.set_markeredgewidth(2)
	if xtra_low ==True:
		for i in range(len(freq_gleam)):
			if freq_gleam[i] in [84,92,99]:
				(_, caps, _) = ax.errorbar(freq_gleam[i], flux_gleam[i], flux_gleam_err[i],marker = 'o',capsize=2.8,elinewidth=1.5, markersize=5,  color = 'k', linestyle='none',markeredgecolor='none')
				# compx = np.linspace(min(freq_mwa)-0.1*min(freq_mwa),max(freq_mwa)+0.1*max(freq_mwa))

				ax.set_xlim(min(freq_cont), max(freq_cont))
				for cap in caps:
					# cap.set_color('crimson') 
					cap.set_markeredgewidth(2)
			if freq_gleam[i] in [76]:
				(_, caps, _) = ax.errorbar(freq_gleam[i], flux_gleam[i], flux_gleam_err[i],marker = 'o',capsize=2.8,elinewidth=1.5, markersize=5,  color = 'k', linestyle='none',markeredgecolor='none', label='Gleam')
				# compx = np.linspace(min(freq_mwa)-0.1*min(freq_mwa),max(freq_mwa)+0.1*max(freq_mwa))

				ax.set_xlim(min(freq_cont), max(freq_cont))
				for cap in caps:
					# cap.set_color('crimson') 
					cap.set_markeredgewidth(2)
		ax.legend(loc='upper right')                                    
        # ax.errorbar(freq_xtra,flux_xtra,flux_err_xtra, marker = '>', color='k',label='Non-GLEAM')    
	else:   
		ax.errorbar(freq, flux_yr1, flux_err_yr1, marker = '.', color = 'C6', linestyle='none', markersize=7, label='GLEAM Year1')
		ax.errorbar(freq, flux_yr2, flux_err_yr2, marker = '.', color = 'purple', linestyle='none', markersize=7, label='GLEAM Year2')


		compx = np.linspace(min(freq)-0.1*min(freq),max(freq)+0.1*max(freq))
		ax.set_xlim(min(freq_cont), max(freq_cont))
		# ax.set_ylim(min)
	if xtra_high == True:
		for i in range(len(freq_xtra)):
			if flux_xtra[i] != 0.0:
				if freq_xtra[i] in [150]:  
					(_, caps, _) = ax.errorbar(freq_xtra[i], flux_xtra[i], flux_err_xtra[i], marker = 's',elinewidth=2,  markersize=10, color = 'k', linestyle='none', label = 'TGSS',markeredgecolor='none')
					for cap in caps:
						cap.set_color('dodgerblue') 
						cap.set_markeredgewidth(2)
				elif freq_xtra[i] in [408]:
					(_, caps, _) = ax.errorbar(freq_xtra[i], flux_xtra[i], flux_err_xtra[i], marker = '<',elinewidth=2, markersize=10,  color = 'forestgreen', linestyle='none', label = 'MRC',markeredgecolor='none')
					for cap in caps:
						cap.set_color('forestgreen') 
						cap.set_markeredgewidth(2)
				elif freq_xtra[i] in [1400]:
					(_, caps, _) = ax.errorbar(freq_xtra[i], flux_xtra[i], flux_err_xtra[i], marker = 'v',elinewidth=2, markersize=10,  color = 'navy', linestyle='none', label = 'NVSS',markeredgecolor='none')
					for cap in caps:
						cap.set_color('navy') 
						cap.set_markeredgewidth(2)            
				elif freq_xtra[i] in [74]:
					(_, caps, _) = ax.errorbar(freq_xtra[i], flux_xtra[i], flux_err_xtra[i], marker = 'v',elinewidth=2, markersize=10,  color = 'crimson', linestyle='none', label = 'VLSSR',markeredgecolor='none')
					for cap in caps:
						cap.set_color('navy') 
						cap.set_markeredgewidth(2)
				elif freq_xtra[i] in [843]:
					(_, caps, _) = ax.errorbar(freq_xtra[i], flux_xtra[i], flux_err_xtra[i], marker = '<',elinewidth=2, markersize=10,  color = 'orange', linestyle='none', label = 'SUMSS',markeredgecolor='none')
					for cap in caps:
						cap.set_color('orange') 
						cap.set_markeredgewidth(2)


	if grid == True:
		ax.grid(which='both')
		if resid == True:
			ax1.grid(axis = 'x',which = 'both')

	ax.legend(loc='upper right')

	if savefig == True:
		if not os.path.exists(directory):
			os.makedirs(directory)
			print('Creating directory '+ directory +' and saving figures in png format with title names.')
		# plt.legend(loc='upper right')
		plt.savefig(directory+str(name.replace(' ','_'))+'.png', bbox_inches='tight')
		plt.close()

# return resid_yr1, resid_yr2, plt_resid_1, plt_resi2_2
