
import mypy
import numpy
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages





data = mypy.readcat('../data/table_SFing_galaxies.dat')





##################################
###  Defining analysis parameters
##################################

zlo, zhi = 0.7, 1.3


dlmass = 0.5
lmassbins = numpy.arange(8.5-dlmass/2., 11.5+dlmass, dlmass)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

overdensbins          = numpy.array([-99., 0.5, 1.0, 99.])
overdensbins_plotting = numpy.array([-1., 0.5, 1.0, 2.5])
overdensbars = (overdensbins[1:] + overdensbins[:-1]) / 2.
overdens_labels = ['-$\infty$ < log(1+$\delta$) < %.1f' % overdensbins[1] + '   a.k.a. "field"',
                   '%.1f < log(1+$\delta$) < %.1f' % (overdensbins[1], overdensbins[2]) + '   a.k.a "groups"',
                   '%.1f < log(1+$\delta$) < $\infty$' % overdensbins[2] + '   a.k.a. "clusters"']
overdens_labels = ['"field"          -$\infty$ < log(1+$\delta$) < %.1f' % overdensbins[1],
                   '"groups"      %.1f < log(1+$\delta$) < %.1f' % (overdensbins[1], overdensbins[2]),
                   '"clusters"     %.1f < log(1+$\delta$) < $\infty$' % overdensbins[2]]
overdens_labels = ['"field"                    log(1+$\delta$) < %.1f' % overdensbins[1],
                   '"groups"      %.1f < log(1+$\delta$) < %.1f' % (overdensbins[1], overdensbins[2]),
                   '"clusters"    %.1f < log(1+$\delta$)         ' % overdensbins[2]]
overdens_labels = ['"field"', 
                   '"groups"', 
                   '"clusters"']

digi_lmass = numpy.digitize(data.lmass, lmassbins)
digi_overdens = numpy.digitize(data.log_overdens, overdensbins)

colors_overdens = ['b', 'g', 'r']

























#############################################
###  Plotting SFR histograms at fixed M*  ###
#############################################

if True:


	pdf = PdfPages('../figures/SFR_histograms.pdf')




	fig = pyplot.figure(figsize=(12.3, 10.6))
	sp1 = fig.add_subplot(311)
	sp2 = fig.add_subplot(312)
	sp3 = fig.add_subplot(313)
	sps = [sp3, sp2, sp1]

	fig.subplots_adjust(hspace=0)



	for i_lmassbin in range(1, len(lmassbins)):

		sp1.minorticks_on()
		sp2.minorticks_on()
		sp3.minorticks_on()

		#sp1.grid()
		#sp2.grid()
		#sp3.grid()

		sp3.set_xlabel('SFR$_{UV+IR}$ [M$_{\odot}$ / yr]')
		sp1.set_ylabel('Number')
		sp2.set_ylabel('Number')
		sp3.set_ylabel('Number')


		label_lmass = 'log(M) = %.2f - %.2f' % (lmassbins[i_lmassbin-1], lmassbins[i_lmassbin])

		###  list of the output from pyplot.hist for each overdens bin
		hist_data = []

		for i_overdensbin in range(1, len(overdensbins)):

			sp = sps[i_overdensbin-1]
			inds_lmass_overdens_bin = (digi_lmass == i_lmassbin) & (digi_overdens == i_overdensbin)

			subinds_negSFR = data.SFR_UVIR[inds_lmass_overdens_bin] <= 0
			subinds_posSFR = data.SFR_UVIR[inds_lmass_overdens_bin] > 0

			max_SFR = data.SFR_UVIR[inds_lmass_overdens_bin].max()
			min_SFR = data.SFR_UVIR[inds_lmass_overdens_bin].min()


			###  plotting histogram
			sfr_range = (-100, 100)
			nbins = int(numpy.count_nonzero(inds_lmass_overdens_bin)**0.5) * 6

			hist_data.append(sp.hist(data.SFR_UVIR[inds_lmass_overdens_bin],
				             histtype='stepfilled', color=colors_overdens[i_overdensbin-1], lw=3,
				             range=sfr_range, alpha=0.25,
				             bins=nbins))

			hist_outline = sp.hist(data.SFR_UVIR[inds_lmass_overdens_bin],
				             histtype='step', color=colors_overdens[i_overdensbin-1], lw=3,
				             range=sfr_range,
				             bins=nbins)


			sp.set_ylim(0, hist_outline[0].max()*1.15)
			sp.set_xlim(-45, 60)





			###  plotting median SFR
			sfr_medi = numpy.percentile(data.SFR_UVIR[inds_lmass_overdens_bin], 50)

			l = sp.axvline(sfr_medi, color='k', lw=5)
			l = sp.axvline(sfr_medi, color=colors_overdens[i_overdensbin-1], lw=1.)

			sp.axvline(0, color='gray', ls='--')






			text = overdens_labels[i_overdensbin-1]
			text += '\n%s' % label_lmass
			text += '\nN = %i' % numpy.count_nonzero(inds_lmass_overdens_bin)
			t = sp.text(0.03, 0.95, text, transform=sp.transAxes,
				        horizontalalignment='left', verticalalignment='top')
			#	           label=overdens_labels[i_overdensbin-1] + ' (N=%i)' % numpy.count_nonzero(inds_lmass_overdens_bin))


		pdf.savefig()
		sp1.clear()
		sp2.clear()
		sp3.clear()

	pdf.close()
	pyplot.close()













