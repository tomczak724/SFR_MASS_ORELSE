
import mypy
import numpy
import pickle
from astropy import wcs
from astropy.io import fits
from matplotlib import pyplot


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
		self.dens_at_z = pylab.zeros(ngals) - 99
		self.overdens_at_z = pylab.zeros(ngals) - 99


github_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
class field:
	def __init__(self, name, version, zclust, sigmaz, alpha=50, chi2red_thresh=10, 
		         uvj_slope=0.88, uvj_intercept=0.59):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz
		self.alpha = alpha

		print '  reading: %s.cat' % version
		self.cat = mypy.readcat('%s/%s/%s.cat.gz' % (github_dir, version, version))
		print '  reading: %s.zout' % version
		self.zout = mypy.readzout('%s/%s/%s.zout.gz' % (github_dir, version, version))
		print '  reading: %s.fout' % version
		self.fout = mypy.readcat('%s/%s/%s.fout.gz' % (github_dir, version, version))
		print '  reading: %s.fir' % version
		self.fir = mypy.readcat('%s/%s/%s.fir.gz' % (github_dir, version, version))
		print '  reading: %s.restframe' % version
		self.restframe = mypy.readcat('%s/%s/%s.restframe.gz' % (github_dir, version, version))
		#print '  reading: %s.restframe_colors' % version
		#self.rfcolors = mypy.readcat('%s/%s/%s.restframe_colors.gz' % (github_dir, version, version))


		###  SETTING OBJECTS IDENTIFIED AS SECURE STARS FROM SPECTROSCOPY TO use=0
		self.crossmatch = mypy.readcat('%s/%s/%s.crossmatch.gz' % (github_dir, version, version), dtype=str)
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




###  NOTE: Not all fields have MIPS 24um data
print ''
fields = []
#fields.append(py_classes.field('N200',    'nep200_v0.0.6',       0.691,  0.027, alpha=1./600,  chi2red_thresh=7, uvj_intercept=0.49))
#fields.append(py_classes.field('RXJ1221', 'rxj1221+4918_v0.0.2', 0.700,  0.023, alpha=1./500,  chi2red_thresh=10))
#fields.append(py_classes.field('SC1324',  'sc1324_v0.0.3',       0.755,  0.033, alpha=1./600,  chi2red_thresh=10))
fields.append(py_classes.field('RCS0224', 'rcs0224-0002_v0.0.3', 0.772,  0.027, alpha=1./500,  chi2red_thresh=10))
#fields.append(py_classes.field('RXJ1716', 'rxj1716+6708_v0.0.8', 0.813,  0.021, alpha=1./500,  chi2red_thresh=8))
fields.append(py_classes.field('N5281',   'nep5281_v0.0.3',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10))
fields.append(py_classes.field('SG0023',  'sg0023+0423_v0.2.0',  0.845,  0.025, alpha=1./500,  chi2red_thresh=14))
fields.append(py_classes.field('SC1604',  'sc1604_v0.0.5',       0.910,  0.029, alpha=1./500,  chi2red_thresh=10))
fields.append(py_classes.field('CL1429',  'cl1429+4241_v0.0.3',  0.920,  0.051, alpha=1./500,  chi2red_thresh=10))
#fields.append(py_classes.field('CL1137',  'cl1137+3007_v0.0.1',  0.959,  0.032, alpha=1./500,  chi2red_thresh=10))
fields.append(py_classes.field('SC0910',  'cl0910+5422_v0.0.4',  1.110,  0.035, alpha=1./500,  chi2red_thresh=10))
fields.append(py_classes.field('RXJ1053', 'rxj1053+5735_v0.0.2', 1.140,  0.031, alpha=1./500,  chi2red_thresh=10))
fields.append(py_classes.field('SC0849',  'sc0849+4452_v0.0.3',  1.261,  0.029, alpha=1./600,  chi2red_thresh=10, uvj_intercept=0.54))
print ''






##################################
###  Defining analysis parameters
##################################

dlmass = 0.5
lmassbins = numpy.arange(8.5-dlmass/2., 11.5+dlmass, dlmass)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.
lmass_labels = ['%05.2f' % lmass_i for lmass_i in lmassbars]

overdensbins = numpy.array([-1.0, 0.5, 1.0, 2.5])
overdensbars = (overdensbins[1:] + overdensbins[:-1]) / 2.
overdens_labels = ['field', 'groups', 'clusters']

###  For a pixel scale of 0.2"/pixel
###  Flux [AB=25]   =   mips_zp * Flux [MJy/sr]
mips_zp = 10**6 * (1.19e-5)**2 * 10**((25 - 8.9) / 2.5)

table = mypy.readcat('../data/table_SFing_galaxies.dat')









######################################################
###  Generating postage stamps of objects with 
###  neighbors removed based on TPHOT's model-fitting
######################################################


###  reading in residual maps
resid_maps = []
wcs_maps = []
for f in fields:
	residname = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/tphot/images/%s_resid_pass2.fits' % f.name
	fopen = fits.open(residname)
	resid_maps.append(fopen)
	wcs_maps.append(wcs.WCS(residname))


radecs = []
modelnames = []
#for i in range(len(table.ID_gal)):
for i in range(50):

	mypy.progress_bar(i, len(table.ID_gal))


	###  grabbing various coordinates of galaxy
	f = fields[int(table.ID_field[i])]
	ra = f.cat.ra[int(table.ID_gal[i])-1]
	dec = f.cat.dec[int(table.ID_gal[i])-1]
	flux24 = f.fir.ftot_mips24[int(table.ID_gal[i])-1]
	weight = table.weight[i]
	radecs.append([ra, dec])


	###  reading in residual map
	resid = resid_maps[int(table.ID_field[i])][0].data
	wcsdat = wcs_maps[int(table.ID_field[i])]

	x,y = wcsdat.wcs_world2pix([ra], [dec], 1)
	x = int(x[0] + 0.5)
	y = int(y[0] + 0.5)



	###  adding flux-scaled model back to the residual map
	modelname = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/tphot/tphot_models/%s/mod-%i.fits' % (f.name, table.ID_gal[i])
	modelnames.append(modelname)
	try:
		model = fits.open(modelname)

		xlo = x - model[0].data.shape[1]/2
		xhi = xlo + model[0].data.shape[1]
		ylo = y - model[0].data.shape[0]/2
		yhi = ylo + model[0].data.shape[0]

		resid[ylo:yhi, xlo:xhi] += model[0].data * flux24 / mips_zp



		###  cutting out 25x25" postage stamp
		dx = 125
		x0_stamp = x - dx/2
		y0_stamp = y - dx/2
		stamp = resid[y0_stamp:y0_stamp+dx, x0_stamp:x0_stamp+dx]

		'''
		hdu = fits.Header()
		hdu['SIMPLE'] = "T" 
		hdu['BITPIX'] = -32
		hdu['NAXIS'] = 2
		hdu['NAXIS1'] = dx
		hdu['NAXIS2'] = dx
		hdu['CRPIX1'] = dx / 2
		hdu['CRPIX2'] = dx / 2
		hdu['CRVAL1'] = ra
		hdu['CRVAL2'] = dec
		hdu['CD1_1'] = model[0].header['CD1_1']
		hdu['CD1_2'] = model[0].header['CD2_1']
		hdu['CD2_1'] = model[0].header['CD1_2']
		hdu['CD2_2'] = model[0].header['CD2_2']
		hdu['WEIGHT'] = weight
		'''

		hdu = model[0].header.copy()
		hdu['CTYPE1'] = 'RA---TAN'
		hdu['CTYPE2'] = 'DEC--TAN'
		hdu['CRPIX1'] = dx / 2
		hdu['CRPIX2'] = dx / 2
		hdu['CRVAL1'] = ra
		hdu['CRVAL2'] = dec


		fits.writeto('../data/stamp_%s_%i.fits' % (f.name, table.ID_gal[i]), stamp, header=hdu, clobber=1)




		###  re-subtracting flux-scaled model back to the residual map
		resid[ylo:yhi, xlo:xhi] -= model[0].data * flux24 / mips_zp

	except: pass

















###  Stacking postage stamps in bins of 
###  lmass and overdensity

digi_lmass = numpy.digitize(table.lmass, lmassbins)
digi_overdens = numpy.digitize(table.log_overdens, overdensbins)


for i_overdens_bin in range(1, len(overdensbins)):
	overdens_label = overdens_labels[i_overdens_bin-1]
	for i_lmass_bin in range(1, len(lmassbins)):
		lmass_label = lmass_labels[i_lmass_bin-1]

		inds_in_bin = numpy.where((digi_overdens == i_overdens_bin) & (digi_lmass == i_lmass_bin))[0]

		stamps = []
		weights = []
		for ind in inds_in_bin:
			f = fields[int(table.ID_field[ind])]

			try:
				stamps.append(fits.getdata('../data/stamp_%s_%i.fits' % (f.name, table.ID_gal[ind])))
				weights.append(table.weight[ind])
			except:
				pass


		###  generating weighted stack
		stack = numpy.average(stamps, axis=0, weights=weights)

		outname = '../data/stack_%s_%s.fits' % (overdens_label, lmass_label)
		fits.writeto(outname, stack, clobber=1)
		print 'wrote to %s' % outname




















