
import mypy
import numpy
import pickle
from astropy.io import fits



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


