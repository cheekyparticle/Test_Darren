from __future__ import division
import matplotlib, os, sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from superplot import data_loader
from superplot import plot_options

from superplot.plotlib import plot_mod as pm 
import superplot.statslib.two_dim as two_dim
import superplot.statslib.one_dim as one_dim
import superplot.statslib.point as stats
from superplot.schemes import scheme_from_yaml
from superplot import plot_options as po

from superplot.plotlib import plots as plots

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

labels, data = data_loader.load(sys.argv[1], sys.argv[2])#"/mt/home/ereid/Downloads/AETERNA.txt""/mt/home/ereid/SB_MO_log_allpost.txt")
if len(sys.argv) == 5: 
	labels2, data2 = data_loader.load(sys.argv[1], sys.argv[4])
if len(sys.argv) == 6: 
	labels2, data2 = data_loader.load(sys.argv[1], sys.argv[4])
	labels3, data3 = data_loader.load(sys.argv[1], sys.argv[5])
if len(sys.argv) == 7: 
	labels2, data2 = data_loader.load(sys.argv[1], sys.argv[4])
	labels3, data3 = data_loader.load(sys.argv[1], sys.argv[5])
	labels4, data4 = data_loader.load(sys.argv[1], sys.argv[6])


#Plot options: didn't know exactly what all of these were for so I mostly guessed/ used the defaults

#BIN_LIMITS = [[-0.3, 2.6], [-7.0, -2.6]]
#ALPHA = [0.045500263896, 0.31731050786]
ALPHA = plot_options.default("alpha")



print len(data)
n_rows = len(data)-3 #Number of rows and columns in the corner plot (or number of variables other than mass being compared)

altcolors1 = plt.cm.seismic(np.linspace(0., 1, 256))# Alternative colormap
altcolors = altcolors1[128: 255]
altmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', altcolors)


colors1 = plt.cm.Reds(np.linspace(0., 1, 256))# custom colormap: basically Reds with the max value changed from off-white to pure white
#for i in range(10):
#	colors1[i,1] = 1.0-i*(1.0-colors1[10,1])
#	colors1[i,2] = 1.0-i*(1.0-colors1[10,2])
colors1[0,1] = 1.0
colors1[0,2] = 1.0
print colors1
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors1)

ax = dict() #Generate dictionary for subplots

linew = 1.0 #Set thickness for all contour lines

fig = plt.figure(figsize=(6.2,6))

plotpdfs = 1 # 1 or 0: if 1 will plot 1D pdfs

for i in range(n_rows):
	for j in range(n_rows-i):
		ax[i,j] = plt.subplot2grid((n_rows+plotpdfs, n_rows+plotpdfs), ((n_rows-1-j+plotpdfs),i))
		XINDEX = 2+i #First two entries are posterior data and chi-sq, so mass is data[2]
		YINDEX = (2+n_rows-j)
		pdf_data = two_dim.kde_posterior_pdf(data[XINDEX], data[YINDEX], data[0])
		pdf = pdf_data.pdf / pdf_data.pdf.sum()
		levels = [two_dim.critical_prof_like(aa) for aa in ALPHA]
		BIN_LIMITS =  [[min(data[XINDEX]), max(data[XINDEX])],[min(data[YINDEX]), max(data[YINDEX])]]
		#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
		plot_limits = [min(data[XINDEX]), max(data[XINDEX]),min(data[YINDEX]), max(data[YINDEX])]
		prof_data = two_dim.profile_like(data[XINDEX], data[YINDEX], data[1], 50)
		#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
		bin_LIMITS = np.array((BIN_LIMITS[0][0],BIN_LIMITS[0][1],BIN_LIMITS[1][0],BIN_LIMITS[1][1]))
		aspect = (plot_limits[1] - plot_limits[0] ) / (plot_limits[3] - plot_limits[2])
		print aspect, type(aspect)
		plt.im = plt.imshow( prof_data.prof_like.T,
							cmap = mymap,
							extent=bin_LIMITS,
							interpolation = 'bilinear',
							label = 'Profile Likelihood',
							origin = 'lower',
							aspect=aspect)
		
		cset = plt.contour(
				prof_data.prof_like.T, 
				levels,
				colors= ['Black', 'Black'],
				hold= 'on',
				extent=bin_LIMITS,
				interpolation = 'bilinear',
				origin=None,
				linestyles=['--', '-'],
				linewidths = linew+0.5
				)
				
		fmt = dict(zip(cset.levels, ['$2\sigma$ region', '$1\sigma$ region']))
		
		#plt.clabel(cset, inline=False, fmt=fmt, fontsize=12, hold='on')
		x_outside = 1E1 * abs(bin_LIMITS[1])
		y_outside = 1E1 * abs(bin_LIMITS[3])
		l2sg = plt.plot(x_outside, y_outside, '--', color='Black' , label='Xe+Ge+Ar ($2\sigma$)', alpha=1.0)
		l1sg = plt.plot(x_outside, y_outside, '-', color='Black' , label='Xe+Ge+Ar ($1\sigma$)', alpha=1.0) 
		#for name, style, color in zip(['Xe+Ge+Ar $2\sigma$ region', 'Xe+Ge+Ar $1\sigma$ region'], ['--', '-'], ['DarkOrange', 'Brown']):
		
			#plt.plot(x_outside, y_outside, style, color=color , label=name, alpha=1.0) 
							
		if len(sys.argv) == 5:
				BIN_LIMITS2 =  [[min(data2[XINDEX]), max(data2[XINDEX])],[min(data2[YINDEX]), max(data2[YINDEX])]]
				#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
				prof_data2 = two_dim.profile_like(data2[XINDEX], data2[YINDEX], data2[1], 50)
				#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
				bin_LIMITS2 = np.array((BIN_LIMITS2[0][0],BIN_LIMITS2[0][1],BIN_LIMITS2[1][0],BIN_LIMITS2[1][1]))
				cset = plt.contour(
									prof_data2.prof_like.T, 
									[levels[0]],
									colors= ['Blue'],
									hold= 'on',
									extent=bin_LIMITS2,
									interpolation = 'bilinear',
									origin=None,
									linestyles=['-'])
				fmt = dict(zip(cset.levels, [' Physics $1\sigma$ region']))
		
				plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')
				#x_outside = 1E1 * abs(bin_LIMITS2[1])
				#y_outside = 1E1 * abs(bin_LIMITS2[3])
				#for name, style in zip(['$Physics 2\sigma$ region'], [ '-']):
					#plt.plot(x_outside, y_outside, style, colors = 'Blues', label=name, alpha=1.0) 
	
		
		if len(sys.argv) == 6:
				BIN_LIMITS2 =  [[min(data2[XINDEX]), max(data2[XINDEX])],[min(data2[YINDEX]), max(data2[YINDEX])]]
				#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
				prof_data2 = two_dim.profile_like(data2[XINDEX], data2[YINDEX], data2[1], 50)
				#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
				bin_LIMITS2 = np.array((BIN_LIMITS2[0][0],BIN_LIMITS2[0][1],BIN_LIMITS2[1][0],BIN_LIMITS2[1][1]))
				cset = plt.contour(
									prof_data2.prof_like.T, 
									[levels[0]],
									colors= ['Blue'],
									hold= 'on',
									extent=bin_LIMITS2,
									interpolation = 'bilinear',
									origin=None,
									linestyles=['-'])
				fmt = dict(zip(cset.levels, [' Physics $1\sigma$ region']))
		
				#plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')
				x_outside = 1E1 * abs(bin_LIMITS2[1])
				y_outside = 1E1 * abs(bin_LIMITS2[3])
				for name, style in zip(['$Physics 2\sigma$ region', '$Physics 1\sigma$ region'], ['--', '-']):
					plt.plot(x_outside, y_outside, style, color = 'RoyalBlue', label=name, alpha=0.7) 
	
	
				BIN_LIMITS3 =  [[min(data3[XINDEX]), max(data3[XINDEX])],[min(data3[YINDEX]), max(data3[YINDEX])]]
				#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
				prof_data3 = two_dim.profile_like(data3[XINDEX], data3[YINDEX], data3[1], 50)
				#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
				bin_LIMITS3 = np.array((BIN_LIMITS3[0][0],BIN_LIMITS3[0][1],BIN_LIMITS3[1][0],BIN_LIMITS3[1][1]))
				cset = plt.contour(
									prof_data3.prof_like.T, 
									[levels[0]],
									colors= ['Green'],
									hold= 'on',
									extent=bin_LIMITS3,
									interpolation = 'bilinear',
									origin=None,
									linestyles=['-'])
				fmt = dict(zip(cset.levels, [' Physics $1\sigma$ region']))
		
				#plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')
				x_outside = 1E1 * abs(bin_LIMITS3[1])
				y_outside = 1E1 * abs(bin_LIMITS3[3])
				for name, style in zip(['$Physics 2\sigma$ region', '$Physics 1\sigma$ region'], ['--', '-']):
					plt.plot(x_outside, y_outside, style, color = 'RoyalBlue', label=name, alpha=0.7) 
	
		
		plt.xlim(BIN_LIMITS[0])
		plt.ylim(BIN_LIMITS[1])
		
		best_fit_x = stats.best_fit(data[1], data[XINDEX])
		best_fit_y = stats.best_fit(data[1], data[YINDEX])
		print best_fit_x, best_fit_y
		plt.plot([best_fit_x], [best_fit_y], marker = '*', zorder=2, color='white')
		
		
		
		if len(sys.argv) == 7:
				BIN_LIMITS2 =  [[min(data2[XINDEX]), max(data2[XINDEX])],[min(data2[YINDEX]), max(data2[YINDEX])]]
				#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
				prof_data2 = two_dim.profile_like(data2[XINDEX], data2[YINDEX], data2[1], 50)
				#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
				bin_LIMITS2 = np.array((BIN_LIMITS2[0][0],BIN_LIMITS2[0][1],BIN_LIMITS2[1][0],BIN_LIMITS2[1][1]))
				csetXe = plt.contour(
									prof_data2.prof_like.T, 
									[levels[0]],
									colors= ['Orange'],
									hold= 'on',
									extent=bin_LIMITS2,
									interpolation = 'bilinear',
									origin=None,
									linestyles=['--'],
									linewidths = linew)
				fmt = dict(zip(cset.levels, [' Physics $1\sigma$ region']))
		
				#plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')
				x_outside = 1E1 * abs(bin_LIMITS2[1])
				y_outside = 1E1 * abs(bin_LIMITS2[3])
				lXe = plt.plot(x_outside, y_outside, '--', color = 'Orange', label='Xe ($2\sigma$)', alpha=0.7) 
	
	
				BIN_LIMITS3 =  [[min(data3[XINDEX]), max(data3[XINDEX])],[min(data3[YINDEX]), max(data3[YINDEX])]]
				#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
				prof_data3 = two_dim.profile_like(data3[XINDEX], data3[YINDEX], data3[1], 50)
				#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
				bin_LIMITS3 = np.array((BIN_LIMITS3[0][0],BIN_LIMITS3[0][1],BIN_LIMITS3[1][0],BIN_LIMITS3[1][1]))
				csetGe = plt.contour(
									prof_data3.prof_like.T, 
									[levels[0]],
									colors= ['Blue'],
									hold= 'on',
									extent=bin_LIMITS3,
									interpolation = 'bilinear',
									origin=None,
									linestyles=['--'],
									linewidths = linew)
				fmt = dict(zip(cset.levels, [' Physics $1\sigma$ region']))
		
				#plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')
				x_outside = 1E1 * abs(bin_LIMITS3[1])
				y_outside = 1E1 * abs(bin_LIMITS3[3])
				lGe = plt.plot(x_outside, y_outside, '--', color = 'Blue', label='Ge ($2\sigma$)', alpha=0.7) 
				BIN_LIMITS4 =  [[min(data4[XINDEX]), max(data4[XINDEX])],[min(data4[YINDEX]), max(data4[YINDEX])]]
				#BIN_LIMITSXe =  [[min(dataXe[XINDEX]), max(dataXe[XINDEX])],[min(dataXe[YINDEX]), max(dataXe[YINDEX])]]
				prof_data4 = two_dim.profile_like(data4[XINDEX], data4[YINDEX], data4[1], 50)
				#prof_data_Xe = two_dim.profile_like(dataXe[XINDEX], dataXe[YINDEX], dataXe[1], 70)
				bin_LIMITS4 = np.array((BIN_LIMITS4[0][0],BIN_LIMITS4[0][1],BIN_LIMITS4[1][0],BIN_LIMITS4[1][1]))
				csetAr = plt.contour(
									prof_data4.prof_like.T, 
									[levels[0]],
									colors= ['Green'],
									hold= 'on',
									extent=bin_LIMITS4,
									interpolation = 'bilinear',
									origin=None,
									linestyles=['--'],
									linewidths = linew)
				fmt = dict(zip(cset.levels, [' Physics $1\sigma$ region']))
		
				#plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')
				x_outside = 1E1 * abs(bin_LIMITS4[1])
				y_outside = 1E1 * abs(bin_LIMITS4[3])
				lAr=plt.plot(x_outside, y_outside, '--', color = 'Green', label='Ar ($2\sigma$)', alpha=0.7) 
		plt.xlim(BIN_LIMITS[0])
		plt.ylim(BIN_LIMITS[1])
		
		best_fit_x = stats.best_fit(data[1], data[XINDEX])
		best_fit_y = stats.best_fit(data[1], data[YINDEX])
		print best_fit_x, best_fit_y
		plt.plot([best_fit_x], [best_fit_y], marker = '*', zorder=2,  color='white', mec='black')

		
		if j == 0:
			plt.xlabel(r"$%s$"%labels[XINDEX])#Only label axes on the edge
			#if i== n_rows-1:
				#plt.legend()
			
		else:
			plt.xticks([])#Remove axis ticks from central plots
		if i == 0:
			plt.ylabel(r"$%s$"%labels[YINDEX])
		 	ax[i,j].yaxis.set_label_coords(-0.4, 0.5) 
		else:
			plt.yticks([])

if plotpdfs == 1:
	axpdf = dict()
	for i in range(n_rows): #Generate 1D pdf plots for coefficients other than mass
		axpdf[n_rows-i] = plt.subplot2grid((n_rows+plotpdfs, n_rows+plotpdfs), ((n_rows-i),(n_rows-i)))
		YINDEX = (2+n_rows-i)
		prof_data1D = one_dim.prof_data(data[YINDEX], data[1], 50)#one_dim.posterior_pdf gives pdfs
		#print prof_data1D
		plt.plot(prof_data1D[2], prof_data1D[1], color='red')
		#plt.xlim([min(data[YINDEX]), max(data[YINDEX])])
		plt.yticks([])		
		if i == 0:
			plt.xlabel(r"$%s$"%labels[YINDEX])
		else:
			plt.xticks([])
		#plt.ylabel(r"$%s$"%labels[YINDEX])
	
	axpdf[0] = plt.subplot2grid((n_rows+plotpdfs, n_rows+plotpdfs), (0,0))
	prof_data1D = one_dim.prof_data(data[2], data[1], 50)
	plt.plot(prof_data1D[2], prof_data1D[1], color='red')
	#plt.xlim([min(data[2]), max(data[2])])
	plt.xticks([])
	plt.yticks([])	

#plt.figlegend((l1sg[0], l2sg[0], lXe[0], lGe[0], lAr[0]), ('Xe+Ge+Ar ($1\sigma$)','Xe+Ge+Ar ($2\sigma$)','Xe ($2\sigma$)','Ge ($2\sigma$)', 'Ar ($2\sigma$)' ), loc=[0.6, 0.6])
#plt.tight_layout()
plt.savefig(sys.argv[2].split('/')[0]+'newfig.pdf')
plt.show()
#plt.savefig('fig.png')
