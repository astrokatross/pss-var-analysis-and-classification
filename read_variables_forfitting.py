# This is a script entirely to just read the data so that it is external and not chunky in the main pipeline. 
# By K.Ross 30/8/19

from astropy.io import fits
import numpy as np


def read_variables(catalogue_dir,year_str, other_surveys):
	type_prefix = 'S_'
	# type_prefix = "e+'_err_'rr_int_flux_"
	# Reading in the raw catalogue of all sources in the field that have S_wide >0.16J
	hdulist = fits.open(catalogue_dir)
	tbdata=hdulist[1].data
	orig_cols=hdulist[1].columns
	hdulist.close()
	# models = [gpscssmodels.powlaw, gpscssmodels.powlaw, gpscssmodels.quad_plot]
	# Reading in important values 
	name = np.array(tbdata['Name'])
	ra_deg = np.array(tbdata['RA'])
	dec_deg = np.array(tbdata['Dec'])
	large_err_ind = np.where((dec_deg > 18.5) | (dec_deg < -72))			#Note probably dont need this cause moasaic yr2 not in this range
	small_err_ind = np.where((dec_deg <= 18.5) |  (dec_deg >= -72))
	ra_hms = np.array(tbdata['RA_hms'])
	dec_dms = np.array(tbdata['Dec_dms'])



	# This is just creating arrays of all the different frequencies you need for each type of fit, note: they are the same FREQ for each year but different FLUX
	freq = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
	freq_mwa = np.array([107,115,122,130,143,151,158,166,174,181,189,197,204,212,220,227])
	freq_high = np.array([189,212,843,1400])
	freq_low = np.array([107,115,122,130,143,151,158,166,174,181,197,204,220,227])
	freq_xtra = [74,150,408,843,1400]

	freq_list = ['S_vlssr','107', '115', '123', '130', '143','S_tgss','150', '158', '166', '174', '181', '189', '197', '204', '212', '220', '227','S_mrc','S_sumss','S_nvss']
	flux_mwa_list = ['107','115','123','130','143','150','158','166','174','181','189','197','204','212','220','227']
	freq_high_list = ['189', '212','S_sumss','S_nvss']
	freq_low_list = ['107','115','123','130','143','150','158','166','174','181','197','204','220','227']
	freq_xtra_surveys_list = ['S_mrc','S_sumss','S_nvss']


	# Creating the arrays to input all the flux measurements: need each individual year here 
	flux_err = []
	flux_temp = []
	flux_mwa = []
	flux_mwa_err = []
	flux_low = []
	flux_err_low = []
	flux_high = []
	flux_err_high = []
	local_rms_low = [] 
	local_rms = []


	# Reading in the flux values for each year but adding in the errors based on their positions etc
	# Making one big for loop that will set each of them, seems efficient right? Lol nope but ok. 
	# Ok I think Im going to have to go through them one by one because this is getting far too messy and I dont trust they'll stay in the correct order

	for flux in freq_list:
		if flux == 'S_vlssr':
			print("Adding: "+flux)
			# flux_mwa.append(np.array(tbdata[flux]))
			# flux_mwa_err.append(np.array(tbdata[flux])*0.05)
			# flux_temp.append(np.array(tbdata[flux]))
			# flux_err.append(np.array(tbdata[flux])*0.05)
		elif flux in ['107']:
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_low.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_low.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
			local_rms_low.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux in ['115', '123', '130','143', '150']:
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_low.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_low.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
			local_rms_low.append(np.array(tbdata['local_rms_'+flux+year_str]))		
		elif flux == "S_tgss":
			print("Adding: "+flux)
			tgss = np.array(tbdata[flux])
			tgss = np.where(np.isnan(tgss), 0., tgss)
			# flux_temp.append(tgss)
			tgss_err = np.array(tbdata['S_tgss_err'])
			tgss_err = np.where(np.isnan(tgss_err), 0., tgss_err)
			# flux_err.append(tgss_err)
		elif flux in ['158', '166']:#, '174', '181']:
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_low.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_low.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
			local_rms_low.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux in ['174', '181']:
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(2 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_low.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_low.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(2 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(2 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
			local_rms_low.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux == '189':
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_high.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_high.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux in ['197', '204']:
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_low.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_low.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(3 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
			local_rms_low.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux =='212':
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_high.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_high.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux in [ '220', '227']:
			print("Adding: "+flux)
			flux_temp.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_low.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_err_low.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			flux_mwa.append(np.array(tbdata[type_prefix+flux+year_str]))
			flux_mwa_err.append(np.sqrt(np.array(tbdata[type_prefix+flux+'_err'+year_str])**2 + np.array(4 / 100. * tbdata[type_prefix+flux+year_str])**2))
			local_rms.append(np.array(tbdata['local_rms_'+flux+year_str]))
			local_rms_low.append(np.array(tbdata['local_rms_'+flux+year_str]))
		elif flux == 'S_mrc':
			print("Adding: "+flux)
			mrc = np.array(tbdata[flux])
			mrc = np.where(np.isnan(mrc), 0., mrc)
			# flux_temp.append(mrc)
			# flux_err.append(mrc*0.05)
			# flux_high.append(mrc)
			# flux_err_high.append(mrc*0.05)
		elif flux in ['S_sumss']:
			print("Adding: "+flux)
			sumss = np.array(tbdata[flux])
			sumss = np.where(np.isnan(sumss), 0., sumss)
			# flux_temp.append(sumss)
			# flux_err.append(sumss*0.05)
			flux_high.append(sumss)
			flux_err_high.append(sumss*0.05)
		elif flux in ['S_nvss']:
			print("Adding: "+flux)
			nvss = np.array(tbdata[flux])
			nvss = np.where(np.isnan(nvss), 0., nvss)
			# flux_temp.append(nvss)
			# flux_err.append(nvss*0.05)
			flux_high.append(nvss)
			flux_err_high.append(nvss*0.05)

	flux=np.array(flux_temp)
	flux_err=np.array(flux_err)
	flux_mwa=np.array(flux_mwa)
	flux_mwa_err=np.array(flux_mwa_err)
	flux_low=np.array(flux_low)
	flux_err_low=np.array(flux_err_low)
	flux_high=np.array(flux_high)
	flux_err_high=np.array(flux_err_high)
	S_white = np.array(tbdata['int_flux_wide'+year_str]) #wide band image
	local_rms = np.array(local_rms)


	# Creating a continuous array for plotting/modelling 
	freq_cont = np.linspace(1,6000,20000)
	freq_cont = np.linspace(1,6000,20000)

	# Creating the empty arrays to put in the fit values: NOTE YOU WILL NOT NEED SOME OF THESE AND SHOULD DISCARD THEM 
	params = np.zeros(shape=(len(name),13))
	snr_curve_cut = np.zeros(len(name))
	gleam=1
	if gleam == 1:
		freq_gleam_list = ['076', '084', '092', '099']

		flux_gleam = []
		flux_gleam_err = []

		for freqs in freq_gleam_list:
			flux_gleam.append(np.array(tbdata['S_'+freqs]))
			flux_gleam_err.append(np.sqrt(np.array(tbdata['S_'+freqs+'_err'])**2 + np.array(3 / 100. * tbdata['S_'+freqs])**2))
		flux_gleam=np.array(flux_gleam)
		flux_gleam_err=np.array(flux_gleam_err)

	if other_surveys == 1:

		flux_list_extras = ['S_vlssr','S_tgss','S_mrc','S_sumss','S_nvss']
		freq_xtra = [74,150,408,843,1400]
		flux_xtra = np.zeros((np.shape(flux)[1],5))
		flux_xtra_err = np.zeros((np.shape(flux)[1],5))

		i=0
		for flux_temp in flux_list_extras: 
			S_flux = np.array(tbdata[flux_temp])
			S_flux = np.where(np.isnan(S_flux),0,S_flux)
			# S_flux = S_flux*0.001
			S_flux_err = S_flux*0.05
			# if flux == "S_sumss":
			# 	print('FIXING SUMMS!')
			# 	S_flux = S_flux*0.001
			# 	S_flux_err = S_flux_err*0.001
			flux_xtra[:,i] = S_flux
			flux_xtra_err[:,i] = S_flux_err
			i+=1



	return name, ra_deg, dec_deg, flux, flux_err, flux_mwa, flux_mwa_err, flux_low, flux_err_low, flux_high, flux_err_high, S_white, local_rms, flux_xtra, flux_xtra_err, flux_gleam, flux_gleam_err