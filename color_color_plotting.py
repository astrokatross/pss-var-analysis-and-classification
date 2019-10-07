import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as ss


def plt_colour_color(directory, alpha_low, alpha_high, all_pop, gps_pop, hh_pop, hs_pop, ls_pop, lh_pop):
	# Going to plot the alpha low vs alpha high plot to see if its a similar distribution 
	alpha_low_pts= np.array(all_pop[alpha_low])
	alpha_high_pts = np.array(all_pop[alpha_high])

	alpha_low_pts_gps = np.array(gps_pop[alpha_low])
	alpha_high_pts_gps = np.array(gps_pop[alpha_high])
	alpha_low_pts_hh = np.array(hh_pop[alpha_low])
	alpha_high_pts_hh = np.array(hh_pop[alpha_high])
	alpha_low_pts_hs = np.array(hs_pop[alpha_low])
	alpha_high_pts_hs = np.array(hs_pop[alpha_high])
	alpha_low_pts_ls = np.array(ls_pop[alpha_low])
	alpha_high_pts_ls = np.array(ls_pop[alpha_high])
	alpha_low_pts_lh = np.array(lh_pop[alpha_low])
	alpha_high_pts_lh = np.array(lh_pop[alpha_high])


	x_pts = np.linspace(-2,2,1000)
	y_pts = x_pts


	fig = plt.figure(1,figsize=(15, 10))
	gs = plt.GridSpec(1,1)
	ax = plt.subplot(gs[0])
	# Matching the cut off lines from Joes version 
	ax.axvline(0, color='r', linestyle='dashed', linewidth=2)
	ax.axhline(0, color='r', linestyle='dashed', linewidth=2)
	ax.axvline(0.1, color='b', linewidth=2)
	ax.axhline(-0.5, color='b', linewidth=2)
	ax.plot(x_pts,y_pts, color='r', linestyle='dashed', linewidth=2)


	ax.scatter(alpha_low_pts, alpha_high_pts,s=0.1, color='k', label='Norm pop')
	ax.scatter(alpha_low_pts_hh, alpha_high_pts_hh, s=10, color='r',label='My peaked sources')
	ax.scatter(alpha_low_pts_ls, alpha_high_pts_ls, s=10, color='r')
	ax.scatter(alpha_low_pts_hs, alpha_high_pts_hs, s=10, color='r')
	ax.scatter(alpha_low_pts_gps, alpha_high_pts_gps, s=10, color='r')
	ax.scatter(alpha_low_pts_lh, alpha_high_pts_lh, s=10, color='r')




	plt.xlabel('alpha_low')
	plt.ylabel('alpha_high')
	plt.title('Colour-Colour')
	plt.xlim([-2,2])
	plt.ylim([-3,2])
	plt.legend(loc='upper left')
	plt.savefig(directory,bbox_inches='tight')
	plt.clf()
	# plt.show()

def plt_colour_curve(directory, alpha_low, quad_curve, quad_curve_error_cond, all_pop, gps_pop, hh_pop, hs_pop, ls_pop, lh_pop):
	# Plotting alpha_low and curvature param 
	alpha_low_pts = np.array(all_pop.query(quad_curve_error_cond)[alpha_low])
	curvature_parameter = np.array(all_pop.query(quad_curve_error_cond)[quad_curve])

	print(str(len(alpha_low_pts))+" sources in the normal pop being plotted")

	alpha_low_pts_ls = np.array(ls_pop.query(quad_curve_error_cond)[alpha_low])
	curvature_param_ls = np.array(ls_pop.query(quad_curve_error_cond)[quad_curve])
	alpha_low_pts_gps = np.array(gps_pop.query(quad_curve_error_cond)[alpha_low])
	curvature_param_gps = np.array(gps_pop.query(quad_curve_error_cond)[quad_curve])
	alpha_low_pts_hh = np.array(hh_pop.query(quad_curve_error_cond)[alpha_low])
	curvature_param_hh = np.array(hh_pop.query(quad_curve_error_cond)[quad_curve])
	alpha_low_pts_hs = np.array(hs_pop.query(quad_curve_error_cond)[alpha_low])
	curvature_param_hs = np.array(hs_pop.query(quad_curve_error_cond)[quad_curve])
	alpha_low_pts_lh = np.array(lh_pop.query(quad_curve_error_cond)[alpha_low])
	curvature_param_lh = np.array(lh_pop.query(quad_curve_error_cond)[quad_curve])


	x_pts = np.linspace(-2,2,1000)
	y_pts = x_pts


	fig = plt.figure(1,figsize=(15, 10))
	gs = plt.GridSpec(1,1)
	ax = plt.subplot(gs[0])

	# Matching the cut off lines from Joes version 
	ax.axvline(np.nanmedian(alpha_low_pts), color='r', linestyle='dashed', linewidth=2)
	ax.axhline(-0.2, color='r', linestyle='dashed', linewidth=2)
	ax.axvline(0.1, color='b', linewidth=2)


	ax.scatter(alpha_low_pts, curvature_parameter,s=0.1, label='Normal pop')
	ax.scatter(alpha_low_pts_ls, curvature_param_ls, s=10,color='r')
	ax.scatter(alpha_low_pts_lh, curvature_param_lh, s=10, color='r')
	ax.scatter(alpha_low_pts_hh, curvature_param_hh, s=10, color='r')
	ax.scatter(alpha_low_pts_hs, curvature_param_hs, s=10, color='r')
	ax.scatter(alpha_low_pts_gps, curvature_param_gps, s=15, color='r',label='My peaked sources')
	ax.scatter(alpha_low_pts_hs, curvature_param_hs, s=10, color='r')

	plt.xlabel('alpha_low')
	plt.ylabel('Curvature Parameter')
	plt.title('Curvature vs Alpha_low')
	plt.xlim([-2,2])
	plt.ylim([-2,2])
	plt.legend(loc='upper left')
	plt.savefig(directory,bbox_inches='tight')
	plt.clf()
	# plt.show()

def plt_var_hist(directory, var_param_norm, var_param_gps):
	x_plt_log = np.linspace(-1,7,1000)
	minimum = 0
	maximum = 100
	no_bins_all = 100
	no_bins_gps = 75
	x_plt = np.linspace(0.0001, maximum, 1000)
	bins_norm = np.linspace(minimum,maximum,no_bins_all)
	bins_gps = np.linspace(minimum,maximum,no_bins_gps)

	n_norm, bins_norm = np.histogram(var_param_norm, bins_norm)
	n_gps, bins_gps = np.histogram(var_param_gps, bins_gps)
	plt.clf()
	width = np.diff(bins_norm)
	width_gps = np.diff(bins_gps)
	plt.bar(bins_norm[0:99],n_norm/max(n_norm),alpha=0.5, color='crimson', width=width, label='Normal Population')
	plt.bar(bins_gps[0:74],n_gps/max(n_gps),alpha=0.7, color='purple', width=width_gps, label='GPS Population')

	var_param_cut = var_param_norm[var_param_norm<=144]
	hist_fit = ss.lognorm.fit(var_param_cut)
	print(hist_fit)


	pdf_beta = ss.lognorm.pdf(x_plt,*hist_fit)
	pdf_plt = pdf_beta/max(pdf_beta)
	plt.plot(x_plt,pdf_plt, color='k', linestyle='--', linewidth=2, label='Lognormal PDF')
	plt.axvline(ss.lognorm.ppf(0.99,*hist_fit), color='k', linestyle=':', linewidth=2, label='99% Confidence Cut')


	print('The Var Param value to cut off significance is: ',ss.lognorm.ppf(0.995,*hist_fit))
	plt.title('Normalised Variability Parameter Histogram')
	plt.xlim(0,80)
	plt.legend(loc='upper right')
	plt.xlabel('Variability Parameter')
	plt.savefig(directory,bbox_inches='tight')
	plt.clf()