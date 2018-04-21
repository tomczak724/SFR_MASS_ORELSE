
import mypy
import numpy
import pickle
from scipy import optimize
from astropy.io import fits
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages


def line(x, m, b):
	return m * x + b

def log_SFR_mass_formula(lmass, s0, m0, gamma):
	return s0 - numpy.log10(1 + (10.**lmass / m0)**-gamma)

def log_SFR_tomczak2016(z, lmass, s0_params=[0.448, 1.220, -0.174], m0_params=[9.458, 0.865, -0.132], gamma=1.091):
	'''
	Star-forming galaxies:
	  s0      = 0.448 + 1.220*z - 0.174*z**2
	  log(M0) = 9.458 + 0.865*z - 0.132*z**2
	  gamma   = 1.091

	All galaxies:
	  s0      = 0.195 + 1.157*z - 0.143*z**2
	  log(M0) = 9.244 + 0.753*z - 0.090*z**2
	  gamma   = 1.118

	'''
	s0 = s0_params[0] + s0_params[1] * z + s0_params[2] * z**2
	log_m0 = m0_params[0] + m0_params[1] * z + m0_params[2] * z**2
	return s0 - numpy.log10(1 + (10**lmass / 10**log_m0)**-gamma)



#######################################
###  plotting numbers per bin figure
#######################################

def plot_numbers_per_lmass_overdens_bin(n_field, n_groups, n_clusters, lmassbins, overdensbins):

	lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.
	overdensbars = (overdensbins[1:] + overdensbins[:-1]) / 2.

	fig = pyplot.figure(figsize=(12.7875, 6.175))
	fig.subplots_adjust(left=0.1)

	sp = fig.add_subplot(111)

	sp.minorticks_on()
	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp.set_ylabel('log( 1 + $\delta_{gal}$ )')

	sp.axis([lmassbins.min()-dlmass/3., lmassbins.max()+dlmass/3., -1.0, 2.5])
	sp.axis([lmassbins.min(), lmassbins.max(), overdensbins.min(), overdensbins.max()])


	###  plotting boundaries
	for lmass in lmassbins:
		sp.plot([lmass, lmass], [overdensbins.min(), overdensbins.max()], color='w', lw=5, zorder=1)
		sp.plot([lmass, lmass], [overdensbins.min(), overdensbins.max()], color='k', lw=2, zorder=2)

	for overdens in overdensbins:
		sp.plot([lmassbins.min(), lmassbins.max()], [overdens, overdens], color='w', lw=5, zorder=1)
		sp.plot([lmassbins.min(), lmassbins.max()], [overdens, overdens], color='k', lw=2, zorder=2)

	sp.fill_between([lmassbins.min(), lmassbins.max()], overdensbins[0], overdensbins[1],
		            color='b', alpha=0.15, zorder=1)
	sp.fill_between([lmassbins.min(), lmassbins.max()], overdensbins[1], overdensbins[2],
		            color='g', alpha=0.15, zorder=1)
	sp.fill_between([lmassbins.min(), lmassbins.max()], overdensbins[2], overdensbins[3],
		            color='r', alpha=0.15, zorder=1)



	###  plotting numbers
	for i in range(len(lmassbars)):
		x = lmassbars[i]
		y1, y2, y3 = overdensbars
		sp.text(x, y1, '%i' % n_field[i], horizontalalignment='center', verticalalignment='center')
		sp.text(x, y2, '%i' % n_groups[i], horizontalalignment='center', verticalalignment='center')
		sp.text(x, y3, '%i' % n_clusters[i], horizontalalignment='center', verticalalignment='center')


	figname = 'numbers_per_bin.pdf'
	pyplot.savefig('../figures/%s' % figname)
	print '\nwrote to:\n!open ../figures/%s' % figname
	pyplot.close()





'''
###   OLD/OBSOLETE VERSION 
class voronoi_catalog:
	def __init__(self, version, ngals, nzbins):
		self.version = version            # photometric catalog
		self.zlos = numpy.zeros(nzbins)
		self.zhis = numpy.zeros(nzbins)
		self.zbars = numpy.zeros(nzbins)
		self.ntot = numpy.zeros(nzbins)   # Median number of galaxies in bin from all iterat
		self.zspec_fraction = numpy.zeros(nzbins)   # Median fraction of galaxies in bin with spec-zs
		self.dens_matrix = numpy.zeros((ngals, nzbins)) - 99
		self.overdens_matrix = numpy.zeros((ngals, nzbins)) - 99
		self.bunit_dens = ''
		self.bunit_overdens = ''
		self.dens_at_z = numpy.zeros(ngals) - 99
		self.overdens_at_z = numpy.zeros(ngals) - 99
'''


class voronoi_catalog:
	def __init__(self, version, ngals, nzbins):
		self.version = version            # photometric catalog
		self.zlos = pylab.zeros(nzbins)
		self.zhis = pylab.zeros(nzbins)
		self.zbars = pylab.zeros(nzbins)
		self.ntot = pylab.zeros(nzbins)   # Median number of galaxies in bin from all iterat
		self.zspec_fraction = pylab.zeros(nzbins)   # Median fraction of galaxies in bin with spec-zs
		self.dens_matrix = pylab.zeros((ngals, nzbins))
		self.overdens_matrix = pylab.zeros((ngals, nzbins))
		self.bunit_dens = ''
		self.bunit_overdens = ''
		self.dens_at_z = pylab.zeros(ngals)       # Density at z_peak (z_spec when available)
		self.overdens_at_z = pylab.zeros(ngals)   # Overdensity at z_peak (z_spec when available)

		self.dens_Pz_weighted = pylab.zeros(ngals)       # Average density across all zslices weighted by photo P(z)
		                                                 # (z_spec when available)
		self.overdens_Pz_weighted = pylab.zeros(ngals)   # Average overdensity across all zslices weighted by photo P(z)
		                                                 # (z_spec when available)

		self.dens_matrix[:] = pylab.nan
		self.overdens_matrix[:] = pylab.nan
		self.dens_Pz_weighted[:] = pylab.nan
		self.overdens_Pz_weighted[:] = pylab.nan


github_dir = '/Users/atomczak/GitHub/ORELSE'
github_photo_cats = 'Catalogs/tomczak_catalogs'
github_radio_cats = 'Work_Lu/catalogs/AGN_z_0.55_1.3'
class field:
	def __init__(self, name, version, zclust, sigmaz, alpha=50, chi2red_thresh=10, 
		         uvj_slope=0.88, uvj_intercept=0.59, radio_AGN_filename=None, readcats=True):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz
		self.alpha = alpha

		if readcats:
			print '  reading: %s.cat' % version
			self.cat = mypy.readcat('%s/%s/%s/%s.cat.gz' % (github_dir, github_photo_cats, version, version))
			print '  reading: %s.zout' % version
			self.zout = mypy.readzout('%s/%s/%s/%s.zout.gz' % (github_dir, github_photo_cats, version, version))
			print '  reading: %s_ZFparams.fout' % version
			self.fout = mypy.readcat('%s/%s/%s/%s_ZFparams.fout.gz' % (github_dir, github_photo_cats, version, version))
			print '  reading: %s.fir' % version
			self.fir = mypy.readcat('%s/%s/%s/%s.fir.gz' % (github_dir, github_photo_cats, version, version))
			print '  reading: %s.restframe' % version
			self.restframe = mypy.readcat('%s/%s/%s/%s.restframe.gz' % (github_dir, github_photo_cats, version, version))
			#print '  reading: %s.restframe_colors' % version
			#self.rfcolors = mypy.readcat('%s/%s/%s/%s.restframe_colors.gz' % (github_dir, github_photo_cats, version, version))


			###  SETTING OBJECTS IDENTIFIED AS SECURE STARS FROM SPECTROSCOPY TO use=0
			self.crossmatch = mypy.readcat('%s/%s/%s/%s.crossmatch.gz' % (github_dir, github_photo_cats, version, version), dtype=str)
			self.star_inds = self.crossmatch.Q == '-1'
			for i_star in self.star_inds:
				id_phot_arr = self.crossmatch.id_phot[i_star].split(',')
				for id_phot in id_phot_arr:
					if id_phot != '-1':
						self.cat.use[int(id_phot)-1] *= 0


			print '  reading: %s_voronoi.pickle' % name
			self.voronoi = pickle.load(open('/Users/atomczak/Dropbox/ORELSE/DATA/VORONOI/%s_voronoi.pickle' % name, 'rb'))
			self.overdens_max = []
			for vi in range(len(self.voronoi.overdens_matrix)):
				self.overdens_max.append(self.voronoi.overdens_matrix[vi].max())
			self.overdens_max = numpy.array(self.overdens_max)


			###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
			chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
			cinds = (chi2red > chi2red_thresh) & (self.cat.z_spec < 0)
			self.cat.use[cinds] = 0.


			###  UVJ classification
			self.UV_color = -2.5 * numpy.log10(self.restframe.restflux_U / self.restframe.restflux_V)
			self.VJ_color = -2.5 * numpy.log10(self.restframe.restflux_V / self.restframe.restflux_J)
			self.uvj_class = mypy.uvj_select(self.UV_color, self.VJ_color, self.fout.z)
			#self.uvj_class = uvj_select(self.UV_color, self.VJ_color, self.fout.z, uvj_slope=uvj_slope, uvj_intercept=uvj_intercept)

			###  UVJ classification based on EAZY in 2-filter mode
			#self.UV_color_2filter = self.rfcolors.U_V_color
			#self.VJ_color_2filter = self.rfcolors.V_J_color
			#self.uvj_class_2filter = mypy.uvj_select(self.UV_color_2filter, self.VJ_color_2filter, self.fout.z)

			###  NUVrJ classification
			self.NUVr_color = -2.5 * numpy.log10(self.restframe.restflux_NUV / self.restframe.restflux_r)
			self.rJ_color = -2.5 * numpy.log10(self.restframe.restflux_r / self.restframe.restflux_J)
			self.nuvrj_class = mypy.nuvrj_select(self.NUVr_color, self.rJ_color, self.fout.z)


			###########################
			###  adding Radio AGN flag
			###########################

			self.radio_AGN_cat = None
			self.AGN_flag_radio = numpy.zeros(len(self.cat.id))
			if radio_AGN_filename != None:
				self.radio_AGN_cat = mypy.readcat('%s/%s/%s' % (github_dir, github_radio_cats, radio_AGN_filename), dtype=str)
				ids_radio_AGN = self.radio_AGN_cat.id_phot.astype(int)
				self.AGN_flag_radio[ids_radio_AGN - 1] = 1


			###########################
			###  adding Xray AGN flag
			###########################

			self.xray_AGN_cat = None
			self.AGN_flag_xray = numpy.zeros(len(self.cat.id))



			xyrd1 = self.cat.x[0], self.cat.y[0], self.cat.ra[0], self.cat.dec[0]
			xyrd2 = self.cat.x[1], self.cat.y[1], self.cat.ra[1], self.cat.dec[1]
			d_arcsec = mypy.radec_sep(xyrd1[2], xyrd1[3], xyrd2[2], xyrd2[3])
			d_pixel = ((xyrd1[0]-xyrd2[0])**2 + (xyrd1[1] - xyrd2[1])**2)**0.5 
			self.px_scale = d_arcsec / d_pixel

			### getting z band magnitude
			try: self.zmagnitude = 25 - 2.5 * numpy.log10(self.cat.fluxauto_z)
			except: pass
			if name != 'SC1324':
				try: self.zmagnitude = 25 - 2.5 * numpy.log10(self.cat.fluxauto_Zplus)
				except: pass

			###  setting spatial flag based on Convex Hull method
			#zspecinds = self.cat.z_spec > 0
			#self.spatial_flag = checkPoint(self.cat.ra[zspecinds], self.cat.dec[zspecinds], self.cat.ra, self.cat.dec)
			#self.inds_spatial = self.spatial_flag == 1

			###  setting spatial flag based on Concave Hull method (creates a "tighter" spatial selection than Cnvex Hull)
			zspecinds = self.cat.z_spec > 0
			points = numpy.array(zip(self.cat.x[zspecinds], self.cat.y[zspecinds]))
			self.concave_hull, self.edge_points = mypy.alpha_shape(points, alpha)
			self.inds_spatial = mypy.CheckPoints(self.concave_hull.buffer(10), self.cat.x, self.cat.y)
			self.area_pix2 = self.concave_hull.buffer(10).area
			self.area_arcmin2 = self.area_pix2 * (self.px_scale/60.)**2
			print ''

class baby_field:
	def __init__(self, name, version, zclust, sigmaz, alpha=50, chi2red_thresh=10, 
		         uvj_slope=0.88, uvj_intercept=0.59):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz
		self.alpha = alpha



#########################
###  Reading in catalogs
#########################
print ''

readcats = False

fields0 = []
fields0.append(field('RXJ1221', 'rxj1221+4918_v0.0.4', 0.700,  0.023, alpha=1./500,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('RCS0224', 'rcs0224-0002_v0.0.4', 0.772,  0.027, alpha=1./500,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('CL1350',  'cl1350+6007_v0.0.2',  0.804,  0.035, alpha=1./600,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('RXJ1716', 'rxj1716+6708_v0.0.9', 0.813,  0.017, alpha=1./500,  chi2red_thresh=8, readcats=readcats))
fields0.append(field('N5281',   'nep5281_v0.0.4',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10, radio_AGN_filename='FINAL_nep5281_no_avg_AGN_spec.cat', readcats=readcats))
fields0.append(field('SG0023',  'sg0023+0423_v0.2.0',  0.845,  0.025, alpha=1./500,  chi2red_thresh=14, radio_AGN_filename='FINAL_sc0023_no_avg_AGN_spec.cat', readcats=readcats))
fields0.append(field('SC1604',  'sc1604_v0.0.6',       0.910,  0.029, alpha=1./500,  chi2red_thresh=10, radio_AGN_filename='FINAL_sc1604_no_avg_AGN_spec.cat', readcats=readcats))
fields0.append(field('CL1429',  'cl1429+4241_v0.0.4',  0.920,  0.051, alpha=1./500,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('XLSS005', 'xlss005_v0.0.3',      1.000,  0.024, alpha=1./600,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('SC0910',  'cl0910+5422_v0.0.5',  1.110,  0.035, alpha=1./500,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('RXJ1053', 'rxj1053+5735_v0.0.3', 1.140,  0.031, alpha=1./500,  chi2red_thresh=10, readcats=readcats))
fields0.append(field('SC0849',  'sc0849+4452_v0.0.4',  1.261,  0.029, alpha=1./600,  chi2red_thresh=10, uvj_intercept=0.54, readcats=readcats))
print ''





fields = []
fields.append(fields0[0])    # RXJ1221
fields.append(fields0[1])    # RCS0224
fields.append(fields0[2])    # CL1350
fields.append(fields0[3])    # RXJ1716
fields.append(fields0[4])    # N5281
fields.append(fields0[5])    # SG0023
fields.append(fields0[6])    # SC1604
fields.append(fields0[7])    # CL1429
fields.append(fields0[8])    # XLSS005
fields.append(fields0[9])    # SC0910
fields.append(fields0[10])   # RXJ1053
fields.append(fields0[11])   # SC0849






##################################
###  Defining analysis parameters
##################################

zlo, zhi = 0.9, 1.3
zlo, zhi = 0.7, 0.9
zlo, zhi = 0.7, 1.3


dlmass = 0.5
lmassbins = numpy.arange(8.5-dlmass/2., 11.5+dlmass, dlmass)

#lmassbins = numpy.array([8.25, 9.5, 10.75, 12.0])

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
                   '"clusters"    %.1f < log(1+$\delta$)' % overdensbins[2]]








#######################
###  grabbing galaxies
#######################

if False:

	master_id_field = numpy.array([])
	master_id_gal   = numpy.array([])
	master_ra       = numpy.array([])
	master_dec      = numpy.array([])
	master_z        = numpy.array([])
	master_sfr      = numpy.array([])
	master_lmass    = numpy.array([])
	master_overdens = numpy.array([])



	for i in range(len(fields)):

		f = fields[i]


		###  binning galaxies by lmass, overdensity
		digi_lmass = numpy.digitize(f.fout.lmass, lmassbins)
		digi_overdens = numpy.digitize(f.voronoi.overdens_at_z, overdensbins)

		###  grabbing 24um SNRs for testing purposes
		snr24 = f.fir.ftot_mips24 / f.fir.etot_mips24
		snr24 = numpy.ones(len(f.fir.id)) * 99.

		inds = (f.cat.z_spec > zlo) & \
		       (f.cat.z_spec < zhi) & \
		       (f.cat.use == 1) & \
		       (f.AGN_flag_xray == 0) & \
		       (f.AGN_flag_radio == 0) & \
		       (f.fout.lmass > 8) & \
		       (f.uvj_class == 1) & \
		       (f.fir.weight_mips24 > 0.05) & \
		       (snr24 >= 3) & \
		       (f.voronoi.overdens_at_z > -90)


		master_id_field = numpy.append(master_id_field, [i] * numpy.count_nonzero(inds))
		master_id_gal = numpy.append(master_id_gal, f.cat.id[inds])
		master_ra = numpy.append(master_ra, f.cat.ra[inds])
		master_dec = numpy.append(master_dec, f.cat.dec[inds])
		master_z = numpy.append(master_z, f.cat.z_spec[inds])
		master_sfr = numpy.append(master_sfr, f.fir.SFR_UVIR[inds])
		master_lmass = numpy.append(master_lmass, f.fout.lmass[inds])
		master_overdens = numpy.append(master_overdens, f.voronoi.overdens_at_z[inds])

	digi_lmass = numpy.digitize(master_lmass, lmassbins)
	digi_overdens = numpy.digitize(master_overdens, overdensbins)








###########################################
###  grabbing galaxies from saved txt file
###########################################

if True:

	data = mypy.readcat('../data/table_SFing_galaxies.dat')

	master_id_field = data.ID_field
	master_id_gal   = data.ID_gal
	master_ra       = data.ra
	master_dec      = data.dec
	master_z        = data.z_spec
	master_sfr      = data.SFR_UVIR
	master_lmass    = data.lmass
	master_overdens = data.log_overdens

	digi_lmass = numpy.digitize(master_lmass, lmassbins)
	digi_overdens = numpy.digitize(master_overdens, overdensbins)









############################
###  plotting z histograms
############################

if True:

	fig = pyplot.figure(figsize=(16.5, 7))
	sp = fig.add_subplot(111)

	sp.minorticks_on()
	sp.set_xlabel('spectroscopic redshift')
	sp.set_ylabel('Number')
	sp.set_xlim(zlo-0.08, zhi+0.05)

	fig.subplots_adjust(left=0.07, bottom=0.1)

	hists = []

	colors = ['b', 'g', 'r']
	offsets = [-0.0015, 0.00, 0.0015]
	for i in range(len(overdensbars)):

		inds = (digi_overdens == i+1)
		print numpy.count_nonzero(inds)
		h_plot = sp.hist(master_z[inds], histtype='stepfilled', alpha=0.25, lw=3, color=colors[i],
			             bins=35, range=(zlo+offsets[i], zhi+offsets[i]))
		h_plot = sp.hist(master_z[inds], histtype='step', lw=3, color=colors[i],
			             bins=35, range=(zlo+offsets[i], zhi+offsets[i]))
		h = numpy.histogram(master_z[inds], bins=35, range=(zlo, zhi))
		hists.append(h)
		asdf = sp.axvline(-1, color=colors[i], lw=3, label=overdens_labels[i])

	sp.legend(loc=0, title='UVJ Star-forming galaxies with 24$\mu$m coverage')
	ymaxes = [h[0].max() for h in hists]
	sp.set_ylim(0, max(ymaxes)*1.2)


	###  adding label of fields used
	label = ''
	for f in fields: label += '%s\n' % f.name
	t = sp.text(0.01, 0.98, label, transform=sp.transAxes, 
	            horizontalalignment='left', verticalalignment='top')


	figname = 'redshift_histograms.pdf'
	pyplot.savefig('../figures/%s' % figname)
	print '\nwrote to:\n!open ../figures/%s' % figname
	pyplot.close()







#################################################
###  calculating weights for galaxies to map the
###  field/group z-dist. to the cluster z-dist.
#################################################

if True:

	zbins = hists[0][1]
	zbars = (zbins[1:] + zbins[:-1]) / 2.

	digi_z = numpy.digitize(master_z, zbins)


	weights_field2cluster = hists[2][0] * 1. / hists[0][0]
	weights_group2cluster = hists[2][0] * 1. / hists[1][0]


	###  fixing zbins where both cluster/field or group/field have zero galaxies
	weights_field2cluster[(hists[2][0]==0) & (hists[0][0]==0)] = 0.
	weights_group2cluster[(hists[2][0]==0) & (hists[1][0]==0)] = 0.


	###  fixing zbins where both field or group have zero galaxies but cluster has non-zero
	weights_field2cluster[(hists[2][0]!=0) & (hists[0][0]==0)] = 0.
	weights_group2cluster[(hists[2][0]!=0) & (hists[1][0]==0)] = 0.


	###  normalizing by lowest non-zero value
	weights_field2cluster /= weights_field2cluster[weights_field2cluster.nonzero()].min()
	weights_group2cluster /= weights_group2cluster[weights_group2cluster.nonzero()].min()







####################################
###  saving to intermediate catalog
####################################

if True:

	tablename = 'table_SFing_galaxies.dat'
	master_weights = numpy.array([])
	with open('../data/%s' % tablename, 'w') as outer:
		outer.write('#\n')
		outer.write('#  All galaxies in this table have been selected by these criteria:\n')
		outer.write('#\n')
		outer.write('#    1) Secure zspec between 0.75 < z < 1.3\n')
		outer.write('#    2) Rest-frame (U-V), (V-J) colors of star-forming galaxies\n')
		outer.write('#    3) MIPS 24um coverage (weight_mips24 > 0.05)\n')
		outer.write('#\n')
		outer.write('#  ORELSE ID_field:')
		outer.write('#\n')
		for i in range(len(fields)):
			outer.write('#    %i = %8s, %s\n' % (i, fields[i].name, fields[i].version))
		outer.write('#\n')
		outer.write('#  ID_field')
		outer.write('%15s' % 'ID_gal')
		outer.write('%15s' % 'ra')
		outer.write('%15s' % 'dec')
		outer.write('%15s' % 'z_spec')
		outer.write('%15s' % 'weight')
		outer.write('%15s' % 'lmass')
		outer.write('%15s' % 'SFR_UVIR')
		outer.write('%15s' % 'log_overdens')
		outer.write('\n')

		for i in range(len(master_z)):

			outer.write('%11i' % master_id_field[i])
			outer.write(' %14i' % master_id_gal[i])
			outer.write(' %14.7f' % master_ra[i])
			outer.write(' %14.7f' % master_dec[i])
			outer.write(' %14.7f' % master_z[i])

			if digi_overdens[i] == 1:
				master_weights = numpy.append(master_weights, weights_field2cluster[digi_z[i]-1])
				outer.write(' %14.4f' % weights_field2cluster[digi_z[i]-1])
			elif digi_overdens[i] == 2:
				master_weights = numpy.append(master_weights, weights_group2cluster[digi_z[i]-1])
				outer.write(' %14.4f' % weights_group2cluster[digi_z[i]-1])
			elif digi_overdens[i] == 3:
				master_weights = numpy.append(master_weights, 1.)
				outer.write(' %14.4f' % 1)

			outer.write(' %14.4f' % master_lmass[i])
			outer.write(' %14.4f' % master_sfr[i])
			outer.write(' %14.4f' % master_overdens[i])
			outer.write('\n')









##############################################################################
###  Creating "weighted" arrays of values for calculating z-matched medians.
###  I do this in a hand-wavey way where for each redshift bin in the z
###  histogram I multiply the numbers of field/groups galaxies (Nf, Ng) by 
###  the number of clusters galaxies (Nc) counted in that bin. If Nc = 0 then
###  I remove all field/groups galaxies within the z-slice
##############################################################################

if True:

	weighted_master_id_field = []
	weighted_master_id_gal   = []
	weighted_master_z        = []
	weighted_master_sfr      = []
	weighted_master_lmass    = []
	weighted_master_overdens = []

	zbins = hists[0][1]
	zbars = (zbins[1:] + zbins[:-1]) / 2.

	digi_z = numpy.digitize(master_z, zbins)
	digi_overdens = numpy.digitize(master_overdens, overdensbins)


	scale_factors_field = (weights_field2cluster + 0.5).astype(int)
	scale_factors_group = (weights_group2cluster + 0.5).astype(int)


	for i_zbin in range(1, len(zbins)):

		inds_field = (digi_z == i_zbin) & (digi_overdens == 1)
		inds_group = (digi_z == i_zbin) & (digi_overdens == 2)
		inds_cluster = (digi_z == i_zbin) & (digi_overdens == 3)


		###  adding field galaxies, in numbers matched to the z-distribution of the cluster sample
		scale_factor = scale_factors_field[i_zbin-1]
		weighted_master_id_field += scale_factor * master_id_field[inds_field].tolist()
		weighted_master_id_gal   += scale_factor * master_id_gal[inds_field].tolist()
		weighted_master_z        += scale_factor * master_z[inds_field].tolist()
		weighted_master_sfr      += scale_factor * master_sfr[inds_field].tolist()
		weighted_master_lmass    += scale_factor * master_lmass[inds_field].tolist()
		weighted_master_overdens += scale_factor * master_overdens[inds_field].tolist()

		###  adding group galaxies, in numbers matched to the z-distribution of the cluster sample
		scale_factor = scale_factors_group[i_zbin-1]
		weighted_master_id_field += scale_factor * master_id_field[inds_group].tolist()
		weighted_master_id_gal   += scale_factor * master_id_gal[inds_group].tolist()
		weighted_master_z        += scale_factor * master_z[inds_group].tolist()
		weighted_master_sfr      += scale_factor * master_sfr[inds_group].tolist()
		weighted_master_lmass    += scale_factor * master_lmass[inds_group].tolist()
		weighted_master_overdens += scale_factor * master_overdens[inds_group].tolist()

		###  adding cluster galaxies
		scale_factor = 1
		weighted_master_id_field += scale_factor * master_id_field[inds_cluster].tolist()
		weighted_master_id_gal   += scale_factor * master_id_gal[inds_cluster].tolist()
		weighted_master_z        += scale_factor * master_z[inds_cluster].tolist()
		weighted_master_sfr      += scale_factor * master_sfr[inds_cluster].tolist()
		weighted_master_lmass    += scale_factor * master_lmass[inds_cluster].tolist()
		weighted_master_overdens += scale_factor * master_overdens[inds_cluster].tolist()

	weighted_master_id_field = numpy.array(weighted_master_id_field)
	weighted_master_id_gal   = numpy.array(weighted_master_id_gal)
	weighted_master_z        = numpy.array(weighted_master_z)
	weighted_master_sfr      = numpy.array(weighted_master_sfr)
	weighted_master_lmass    = numpy.array(weighted_master_lmass)
	weighted_master_overdens = numpy.array(weighted_master_overdens)

	digi_lmass_weighted = numpy.digitize(weighted_master_lmass, lmassbins)
	digi_overdens_weighted = numpy.digitize(weighted_master_overdens, overdensbins)







########################################
###  fitting parameterization to SFR-M*
########################################


if True:
	pass











###############################
###  plotting SFR-M* relations
###############################

if True:

	fig          = pyplot.figure(figsize=(9.6875, 7.9749999999999996))
	fig_weighted = pyplot.figure(figsize=(9.6875, 7.9749999999999996))
	fig_resid          = pyplot.figure(figsize=(9.6875, 7.9749999999999996))
	fig_resid_weighted = pyplot.figure(figsize=(9.6875, 7.9749999999999996))

	sp           = fig.add_subplot(111)
	sp_weighted  = fig_weighted.add_subplot(111)
	sp_resid     = fig_resid.add_subplot(111)
	sp_resid_weighted  = fig_resid_weighted.add_subplot(111)

	'''
	xlo, xhi = 8., 12.
	ylo, yhi = 0., 2.
	h2d, xedges, yedges = numpy.histogram2d(master_lmass, master_lsfr, 
		                                    range=([xlo, xhi], [ylo, yhi]),
		                                    bins=(20, 20))
	cmap = pyplot.cm.gray_r
	vmin, vmax = h2d.min(), h2d.max() * 1.5
	cax = sp.imshow(h2d, interpolation='nearest', extent=(xlo, xhi, ylo, yhi), aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

	cbar = fig.colorbar(cax, label='Number')
	'''


	z_median = numpy.median(master_z[digi_overdens == 3])
	z_snaps = [zlo, z_median, zhi]

	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp.set_ylabel('log( SFR / [M$_{\odot}$ / yr] )')
	sp_weighted.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp_weighted.set_ylabel('log( SFR / [M$_{\odot}$ / yr] )')
	sp_resid.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp_resid.set_ylabel('log( SFR ) - log( SFR$_{z=%.2f}$ )' % z_median)
	sp_resid_weighted.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp_resid_weighted.set_ylabel('log( SFR ) - log( SFR$_{z=%.2f}$ )' % z_median)

	sp.axis([8.1, 12.4, -0.7, 2.05])
	sp.minorticks_on()
	sp_weighted.axis([8.1, 12.4, -0.7, 2.05])
	sp_weighted.minorticks_on()
	sp_resid.axis([8.1, 12.4, -0.4, 0.4])
	sp_resid.minorticks_on()
	sp_resid_weighted.axis([8.1, 12.4, -0.4, 0.4])
	sp_resid_weighted.minorticks_on()

	sp_resid.grid()
	sp_resid_weighted.grid()


	###  plotting ZFOURGE

	xmod = numpy.linspace(lmassbins.min(), lmassbins.max(), 100)
	ymod1 = log_SFR_tomczak2016(zlo, xmod)
	ymod2 = log_SFR_tomczak2016(z_median, xmod)
	ymod3 = log_SFR_tomczak2016(zhi, xmod)

	lsfr_zfourge = log_SFR_tomczak2016(z_median, lmassbars)

	ymods = [ymod1, ymod2, ymod3]
	alphas = [1., 1., 1.]
	lw = 1
	
	for i_z in range(3):
		sp.plot(xmod, ymods[i_z], color='k', ls=':',  lw=lw, alpha=alphas[i_z])
		sp.plot(xmod, ymods[i_z], color='k', ls='--', lw=lw, alpha=alphas[i_z])
		sp.text(xmod[-1]+dlmass/4., ymods[i_z][-1], '%.2f' % z_snaps[i_z],
			    horizontalalignment='left', verticalalignment='center')
		sp_weighted.plot(xmod, ymods[i_z], color='k', ls=':',  lw=lw, alpha=alphas[i_z])
		sp_weighted.plot(xmod, ymods[i_z], color='k', ls='--', lw=lw, alpha=alphas[i_z])
		sp_weighted.text(xmod[-1]+dlmass/4., ymods[i_z][-1], '%.2f' % z_snaps[i_z],
			             horizontalalignment='left', verticalalignment='center')
	#sp.axvline(0, color='k', lw=lw, ls='-', label='ZFOURGE')

	sp.fill_between(xmod, ymod1, ymod3, color='k', alpha=0.2, label='ZFOURGE')
	sp_weighted.fill_between(xmod, ymod1, ymod3, color='k', alpha=0.2, label='ZFOURGE')
	sp.axvline(0, color='k', alpha=0.2, lw=10, label='ZFOURGE')
	sp_weighted.axvline(0, color='k', alpha=0.2, lw=10, label='ZFOURGE')


	sp_resid.axhline(0, color='k', ls=':',  lw=lw*2, alpha=alphas[i_z])
	sp_resid.axhline(0, color='k', ls='--',  lw=lw*2, alpha=alphas[i_z])
	sp_resid_weighted.axhline(0, color='k', ls=':',  lw=lw*2, alpha=alphas[i_z])
	sp_resid_weighted.axhline(0, color='k', ls='--',  lw=lw*2, alpha=alphas[i_z])



	#######################################
	###  measuring median SFR-M* for field
	#######################################

	offset = -0.02
	n_field = []
	mean_sfr_field = []
	median_sfr_field = []
	median_sfr_weighted_field = []
	nmad_sfr_field = []
	median_lmasses_field = []
	for i_bin in range(1, len(lmassbins)):
		inds = (digi_lmass == i_bin) & (digi_overdens == 1)
		inds_weighted = (digi_lmass_weighted == i_bin) & (digi_overdens_weighted == 1)
		n_field.append(numpy.count_nonzero(inds))
		median_lmasses_field.append(numpy.median(master_lmass[inds]))
		try:
			mean = numpy.average(master_sfr[inds], weights=master_weights[inds])
			median = numpy.median(master_sfr[inds])
			median_weighted = numpy.median(weighted_master_sfr[inds_weighted])
			nmad = mypy.nmad(master_sfr[inds])

			mean_sfr_field.append(mean)
			median_sfr_field.append(median)
			median_sfr_weighted_field.append(median_weighted)
			nmad_sfr_field.append(nmad)
		except:
			mean_sfr_field.append(10**-4)
			median_sfr_field.append(10**-4)
			median_sfr_weighted_field.append(10**-4)
			nmad_sfr_field.append(10**-4)

	n_field = numpy.array(n_field)
	mean_sfr_field = numpy.array(mean_sfr_field)
	median_sfr_field = numpy.array(median_sfr_field)
	median_sfr_weighted_field = numpy.array(median_sfr_weighted_field)
	nmad_sfr_field = numpy.array(nmad_sfr_field)

	color = 'c'
	color = 'b'
	sp.errorbar(lmassbars+offset, numpy.log10(median_sfr_field), xerr=dlmass/2., yerr=numpy.log10((median_sfr_field + nmad_sfr_field/n_field**0.5) / median_sfr_field),
		        ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[0])
	#sp.errorbar(median_lmasses_field, numpy.log10(median_sfr_field), xerr=dlmass/2., yerr=numpy.log10((median_sfr_field + nmad_sfr_field/n_field**0.5) / median_sfr_field),
	#	        ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[0])
	sp_weighted.errorbar(lmassbars+offset, numpy.log10(median_sfr_weighted_field), xerr=dlmass/2., yerr=numpy.log10((median_sfr_field + nmad_sfr_field/n_field**0.5) / median_sfr_field),
		                 ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[0])

	###  residuals
	sp_resid.errorbar(lmassbars+offset, numpy.log10(median_sfr_field)-lsfr_zfourge, xerr=dlmass/2., yerr=numpy.log10((median_sfr_field + nmad_sfr_field/n_field**0.5) / median_sfr_field),
		              ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[0])
	sp_resid_weighted.errorbar(lmassbars+offset, numpy.log10(median_sfr_weighted_field)-lsfr_zfourge, xerr=dlmass/2., yerr=numpy.log10((median_sfr_field + nmad_sfr_field/n_field**0.5) / median_sfr_field),
		                       ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[0])




	########################################
	###  measuring median SFR-M* for groups
	########################################

	offset = 0.0
	n_groups = []
	mean_sfr_groups = []
	median_sfr_groups = []
	median_sfr_weighted_groups = []
	nmad_sfr_groups = []
	median_lmasses_groups = []
	for i_bin in range(1, len(lmassbins)):
		inds = (digi_lmass == i_bin) & (digi_overdens == 2)
		inds_weighted = (digi_lmass_weighted == i_bin) & (digi_overdens_weighted == 2)
		n_groups.append(numpy.count_nonzero(inds))
		median_lmasses_groups.append(numpy.median(master_lmass[inds]))
		try:
			mean = numpy.average(master_sfr[inds], weights=master_weights[inds])
			median = numpy.median(master_sfr[inds])
			median_weighted = numpy.median(weighted_master_sfr[inds_weighted])
			nmad = mypy.nmad(master_sfr[inds])

			mean_sfr_groups.append(mean)
			median_sfr_groups.append(median)
			median_sfr_weighted_groups.append(median_weighted)
			nmad_sfr_groups.append(nmad)
		except:
			mean_sfr_groups.append(10**-4)
			median_sfr_groups.append(10**-4)
			median_sfr_weighted_groups.append(10**-4)
			nmad_sfr_groups.append(10**-4)

	n_groups = numpy.array(n_groups)
	mean_sfr_groups = numpy.array(mean_sfr_groups)
	median_sfr_groups = numpy.array(median_sfr_groups)
	median_sfr_weighted_groups = numpy.array(median_sfr_weighted_groups)
	nmad_sfr_groups = numpy.array(nmad_sfr_groups)

	color = 'lime'
	color = 'g'
	sp.errorbar(lmassbars+offset, numpy.log10(median_sfr_groups), xerr=dlmass/2., yerr=numpy.log10((median_sfr_groups + nmad_sfr_groups/n_groups**0.5) / median_sfr_groups),
		        ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[1])
	#sp.errorbar(median_lmasses_groups, numpy.log10(median_sfr_groups), xerr=dlmass/2., yerr=numpy.log10((median_sfr_groups + nmad_sfr_groups/n_groups**0.5) / median_sfr_groups),
	#	        ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[1])
	sp_weighted.errorbar(lmassbars+offset, numpy.log10(median_sfr_weighted_groups), xerr=dlmass/2., yerr=numpy.log10((median_sfr_groups + nmad_sfr_groups/n_groups**0.5) / median_sfr_groups),
		                 ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[1])

	###  residuals
	sp_resid.errorbar(lmassbars+offset, numpy.log10(median_sfr_groups)-lsfr_zfourge, xerr=dlmass/2., yerr=numpy.log10((median_sfr_groups + nmad_sfr_groups/n_groups**0.5) / median_sfr_groups),
		              ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[1])
	sp_resid_weighted.errorbar(lmassbars+offset, numpy.log10(median_sfr_weighted_groups)-lsfr_zfourge, xerr=dlmass/2., yerr=numpy.log10((median_sfr_groups + nmad_sfr_groups/n_groups**0.5) / median_sfr_groups),
		                       ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[1])






	##########################################
	###  measuring median SFR-M* for clusters
	##########################################

	offset = 0.02
	n_clusters = []
	mean_sfr_clusters = []
	median_sfr_clusters = []
	median_sfr_weighted_clusters = []
	nmad_sfr_clusters = []
	median_lmasses_clusters = []
	for i_bin in range(1, len(lmassbins)):
		inds = (digi_lmass == i_bin) & (digi_overdens == 3)
		inds_weighted = (digi_lmass_weighted == i_bin) & (digi_overdens_weighted == 3)
		n_clusters.append(numpy.count_nonzero(inds))
		median_lmasses_clusters.append(numpy.median(master_lmass[inds]))
		try:
			mean = numpy.average(master_sfr[inds], weights=master_weights[inds])
			median = numpy.median(master_sfr[inds])
			median_weighted = numpy.median(weighted_master_sfr[inds_weighted])
			nmad = mypy.nmad(master_sfr[inds])

			mean_sfr_clusters.append(mean)
			median_sfr_clusters.append(median)
			median_sfr_weighted_clusters.append(median_weighted)
			nmad_sfr_clusters.append(nmad)
		except:
			mean_sfr_clusters.append(10**-4)
			median_sfr_clusters.append(10**-4)
			median_sfr_weighted_clusters.append(10**-4)
			nmad_sfr_clusters.append(10**-4)

	n_clusters = numpy.array(n_clusters)
	mean_sfr_clusters = numpy.array(mean_sfr_clusters)
	median_sfr_clusters = numpy.array(median_sfr_clusters)
	median_sfr_weighted_clusters = numpy.array(median_sfr_weighted_clusters)
	nmad_sfr_clusters = numpy.array(nmad_sfr_clusters)

	color = 'r'
	sp.errorbar(lmassbars+offset, numpy.log10(median_sfr_clusters), xerr=dlmass/2., yerr=numpy.log10((median_sfr_clusters + nmad_sfr_clusters/n_clusters**0.5) / median_sfr_clusters),
		        ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[2])
	#sp.errorbar(median_lmasses_clusters, numpy.log10(median_sfr_clusters), xerr=dlmass/2., yerr=numpy.log10((median_sfr_clusters + nmad_sfr_clusters/n_clusters**0.5) / median_sfr_clusters),
	#	        ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[2])
	sp_weighted.errorbar(lmassbars+offset, numpy.log10(median_sfr_weighted_clusters), xerr=dlmass/2., yerr=numpy.log10((median_sfr_clusters + nmad_sfr_clusters/n_clusters**0.5) / median_sfr_clusters),
		                 ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[2])

	###  residuals
	sp_resid.errorbar(lmassbars+offset, numpy.log10(median_sfr_clusters)-lsfr_zfourge, xerr=dlmass/2., yerr=numpy.log10((median_sfr_clusters + nmad_sfr_clusters/n_clusters**0.5) / median_sfr_clusters),
		              ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[2])
	sp_resid_weighted.errorbar(lmassbars+offset, numpy.log10(median_sfr_weighted_clusters)-lsfr_zfourge, xerr=dlmass/2., yerr=numpy.log10((median_sfr_clusters + nmad_sfr_clusters/n_clusters**0.5) / median_sfr_clusters),
		                       ls='-', marker='o', color=color, ms=9, mew=2, ecolor=color, elinewidth=1.5, label=overdens_labels[2])




	leg = sp.legend(loc=4, numpoints=1, fontsize=17, frameon=0)
	leg = sp_weighted.legend(loc=4, numpoints=1, fontsize=17, frameon=0,
		                     title='weighted to "clusters" z-dist.')
	leg = sp_resid.legend(loc=1, numpoints=1, fontsize=16, frameon=1, title='')
	leg = sp_resid_weighted.legend(loc=1, numpoints=1, fontsize=16, frameon=1,
		                           title='weighted to "clusters" z-dist.')



	###  plotting numbers of galaxies per lmass-overdens bin
	plot_numbers_per_lmass_overdens_bin(n_field, n_groups, n_clusters, lmassbins, overdensbins_plotting)


	###  adding label of fields used
	label = ''
	for f in fields: label += '%s\n' % f.name
	t = sp.text(0.03, 0.97, label, transform=sp.transAxes, 
	            horizontalalignment='left', verticalalignment='top')
	t = sp_weighted.text(0.03, 0.97, label, transform=sp_weighted.transAxes, 
	                     horizontalalignment='left', verticalalignment='top')
	#t = sp_resid.text(0.03, 0.03, label, transform=sp_resid.transAxes, 
	#                  horizontalalignment='left', verticalalignment='bottom')
	#t = sp_resid_weighted.text(0.03, 0.03, label, transform=sp_resid_weighted.transAxes, 
	#                           horizontalalignment='left', verticalalignment='bottom')



	figname          = 'sfr_mass0.pdf'
	figname_weighted = 'sfr_mass_weighted.pdf'
	figname_resid          = 'resid_sfr_mass0.pdf'
	figname_resid_weighted = 'resid_sfr_mass_weighted.pdf'

	print '\nwrote to:'
	fig.savefig('../figures/%s' % figname)
	print '!open ../figures/%s' % figname
	fig_weighted.savefig('../figures/%s' % figname_weighted)
	print '!open ../figures/%s' % figname_weighted
	fig_resid.savefig('../figures/%s' % figname_resid)
	print '!open ../figures/%s' % figname_resid
	fig_resid_weighted.savefig('../figures/%s' % figname_resid_weighted)
	print '!open ../figures/%s' % figname_resid_weighted
	print ''
	pyplot.close(fig)
	pyplot.close(fig_weighted)
	pyplot.close(fig_resid)
	pyplot.close(fig_resid_weighted)















