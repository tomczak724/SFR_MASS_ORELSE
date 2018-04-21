
###################################################
###################################################
###  This script reads in the ORELSE catalogs
###  and runs the selection algorithm to grab
###  the appropriate galaxies and saves a pickle
###  of all the relevant information for the
###  analysis. This is meant to save time by
###  not having to read in the catalogs every
###  time a change/update is made to the analysis
###################################################
###################################################

import time
import mypy
import numpy
import pickle
import pandas
import progressbar
from scipy import optimize
from astropy.io import fits, ascii
from matplotlib import pyplot
from astropy import coordinates, units, cosmology
from matplotlib.backends.backend_pdf import PdfPages

cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)

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
                 uvj_slope=0.88, uvj_intercept=0.59, radio_AGN_filename=None):
        self.name = name          # e.g. "NEP 200"
        self.version = version    # e.g. "nep200_v0.0.4"
        self.zclust = zclust      # cluster redshift
        self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
        self.zspec_lo = 0         # lower redshift bound for specz
        self.zspec_hi = 0         # upper redshift bound for specz
        self.alpha = alpha

        print '  reading: %s.cat' % version
        self.cat = mypy.readcat('%s/%s/%s/%s.cat.gz' % (github_dir, github_photo_cats, version, version))
        print '  reading: %s.zout' % version
        self.zout = mypy.readzout('%s/%s/%s/%s.zout.gz' % (github_dir, github_photo_cats, version, version))
        print '  reading: %s_ZFparams.fout' % version
        self.fout = mypy.readcat('%s/%s/%s/%s_ZFparams.fout.gz' % (github_dir, github_photo_cats, version, version))
        print '  reading: %s.fout' % version
        self.fout_old = mypy.readcat('%s/%s/%s/%s.fout.gz' % (github_dir, github_photo_cats, version, version))
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

fields0 = []
fields0.append(field('RXJ1221', 'rxj1221+4918_v0.0.4', 0.700,  0.023, alpha=1./500,  chi2red_thresh=10))
fields0.append(field('RCS0224', 'rcs0224-0002_v0.0.4', 0.772,  0.027, alpha=1./500,  chi2red_thresh=10))
fields0.append(field('CL1350',  'cl1350+6007_v0.0.2',  0.804,  0.035, alpha=1./600,  chi2red_thresh=10))
fields0.append(field('RXJ1716', 'rxj1716+6708_v0.0.9', 0.813,  0.017, alpha=1./500,  chi2red_thresh=8))
fields0.append(field('N5281',   'nep5281_v0.0.4',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10, radio_AGN_filename='FINAL_nep5281_no_avg_AGN_spec.cat'))
fields0.append(field('SG0023',  'sg0023+0423_v0.2.0',  0.845,  0.025, alpha=1./500,  chi2red_thresh=14, radio_AGN_filename='FINAL_sc0023_no_avg_AGN_spec.cat'))
fields0.append(field('SC1604',  'sc1604_v0.0.6',       0.910,  0.029, alpha=1./500,  chi2red_thresh=10, radio_AGN_filename='FINAL_sc1604_no_avg_AGN_spec.cat'))
fields0.append(field('CL1429',  'cl1429+4241_v0.0.4',  0.920,  0.051, alpha=1./500,  chi2red_thresh=10))
fields0.append(field('XLSS005', 'xlss005_v0.0.3',      1.000,  0.024, alpha=1./600,  chi2red_thresh=10))
fields0.append(field('SC0910',  'cl0910+5422_v0.0.5',  1.110,  0.035, alpha=1./500,  chi2red_thresh=10))
fields0.append(field('RXJ1053', 'rxj1053+5735_v0.0.3', 1.140,  0.031, alpha=1./500,  chi2red_thresh=10))
fields0.append(field('SC0849',  'sc0849+4452_v0.0.4',  1.261,  0.029, alpha=1./600,  chi2red_thresh=10, uvj_intercept=0.54))
print ''


###  X-ray catalog
xray_filename = '/Users/atomczak/GitHub/ORELSE/Catalogs/X-ray/X-ray_lum_cat.7.17.16.dat'
xray_catalog = mypy.readcat(xray_filename, dtype=str)

xray_lum_tot = xray_catalog.lum_full.astype(float)
xray_ra = xray_catalog.RA.astype(float)
xray_dec = xray_catalog.Dec.astype(float)








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
zlo, zhi = 0.6, 1.3


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







class combined_data_matrix(object):
    def __init__(self):
        pass

class combined_catalog(object):

    def __init__(self, header, n_galaxies, lmassbins, overdensbins, overdens_labels):

        self.time_stamp = time.ctime()
        self.header = header

        self.n_galaxies = n_galaxies

        self.lmassbins = lmassbins
        self.lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

        self.overdensbins = overdensbins
        self.overdensbars = (overdensbins[1:] + overdensbins[:-1]) / 2.

        self.overdens_labels = overdens_labels












#######################
###  grabbing galaxies
#######################

if True:

    field_names = numpy.array([f.name for f in fields])
    field_versions = numpy.array([f.version for f in fields])

    master_id_field    = numpy.array([])
    master_id_gal      = numpy.array([])
    master_ra          = numpy.array([])
    master_dec         = numpy.array([])
    master_zspec       = numpy.array([])
    master_ftot_mips24 = numpy.array([])
    master_etot_mips24 = numpy.array([])
    master_SFR_UVIR    = numpy.array([])
    master_SFR_FAST    = numpy.array([])
    master_lmass       = numpy.array([])
    master_ltau        = numpy.array([])
    master_lage        = numpy.array([])
    master_overdens    = numpy.array([])
    master_UV_color    = numpy.array([])
    master_VJ_color    = numpy.array([])
    master_UVJ_class   = numpy.array([])
    master_chi2        = numpy.array([])

    master_id_spec     = numpy.array([])
    master_slitmask    = numpy.array([])


    for i_field in range(len(fields)):

        f = fields[i_field]


        ###  binning galaxies by lmass, overdensity
        digi_lmass = numpy.digitize(f.fout.lmass, lmassbins)
        digi_overdens = numpy.digitize(f.voronoi.overdens_at_z, overdensbins)


        inds = numpy.where((f.cat.z_spec > zlo) & \
                           (f.cat.z_spec < zhi) & \
                           (f.cat.use == 1) & \
                           (f.AGN_flag_xray == 0) & \
                           (f.AGN_flag_radio == 0) & \
                           (f.fir.weight_mips24 > 0.05) & \
                           (f.voronoi.overdens_at_z > -90))[0]


        master_id_field = numpy.append(master_id_field, [i_field] * numpy.count_nonzero(inds))
        master_id_gal = numpy.append(master_id_gal, f.cat.id[inds])
        master_ra = numpy.append(master_ra, f.cat.ra[inds])
        master_dec = numpy.append(master_dec, f.cat.dec[inds])
        master_zspec = numpy.append(master_zspec, f.cat.z_spec[inds])
        master_ftot_mips24 = numpy.append(master_ftot_mips24, f.fir.ftot_mips24[inds])
        master_etot_mips24 = numpy.append(master_etot_mips24, f.fir.etot_mips24[inds])
        master_SFR_UVIR = numpy.append(master_SFR_UVIR, f.fir.SFR_UVIR[inds])
        master_SFR_FAST = numpy.append(master_SFR_FAST, 10**f.fout.lsfr[inds])
        master_lmass = numpy.append(master_lmass, f.fout.lmass[inds])
        master_ltau = numpy.append(master_ltau, f.fout.ltau[inds])
        master_lage = numpy.append(master_lage, f.fout.lage[inds])
        master_overdens = numpy.append(master_overdens, f.voronoi.overdens_at_z[inds])
        master_UV_color = numpy.append(master_UV_color, f.UV_color[inds])
        master_VJ_color = numpy.append(master_VJ_color, f.VJ_color[inds])
        master_UVJ_class = numpy.append(master_UVJ_class, f.uvj_class[inds])
        master_chi2 = numpy.append(master_chi2, f.zout.chi_p[inds])

        for i in inds:
            try:
                id_gal = f.cat.id[i]
                ind_crossmatch = numpy.where(f.crossmatch.id_phot == '%i' % id_gal)[0][0]

                master_id_spec = numpy.append(master_id_spec, f.crossmatch.id_spec[ind_crossmatch])
                master_slitmask = numpy.append(master_slitmask, f.crossmatch.slitmask[ind_crossmatch])
            except:
                pass

    digi_lmass = numpy.digitize(master_lmass, lmassbins)
    digi_overdens = numpy.digitize(master_overdens, overdensbins)




###  Creating Xray AGN flag in post
if True:
    master_AGN_flag_xray = []
    for i_gal in range(len(master_id_gal)):
        separations = mypy.radec_sep(xray_ra, xray_dec, master_ra[i_gal], master_dec[i_gal])
        ind_min_sep = numpy.where(separations == separations.min())[0][0]
        if separations.min() <= 2. and xray_lum_tot[ind_min_sep] > 10**42.5:
            master_AGN_flag_xray.append(True)
        else:
            master_AGN_flag_xray.append(False)
    master_AGN_flag_xray = numpy.array(master_AGN_flag_xray)


    master_id_field = master_id_field[-master_AGN_flag_xray]
    master_id_gal = master_id_gal[-master_AGN_flag_xray]
    master_ra = master_ra[-master_AGN_flag_xray]
    master_dec = master_dec[-master_AGN_flag_xray]
    master_zspec = master_zspec[-master_AGN_flag_xray]
    master_ftot_mips24 = master_ftot_mips24[-master_AGN_flag_xray]
    master_etot_mips24 = master_etot_mips24[-master_AGN_flag_xray]
    master_SFR_UVIR = master_SFR_UVIR[-master_AGN_flag_xray]
    master_SFR_FAST = master_SFR_FAST[-master_AGN_flag_xray]
    master_lmass = master_lmass[-master_AGN_flag_xray]
    master_ltau = master_ltau[-master_AGN_flag_xray]
    master_lage = master_lage[-master_AGN_flag_xray]
    master_overdens = master_overdens[-master_AGN_flag_xray]
    master_UV_color = master_UV_color[-master_AGN_flag_xray]
    master_VJ_color = master_VJ_color[-master_AGN_flag_xray]
    master_UVJ_class = master_UVJ_class[-master_AGN_flag_xray]
    master_chi2 = master_chi2[-master_AGN_flag_xray]

    digi_lmass = numpy.digitize(master_lmass, lmassbins)
    digi_overdens = numpy.digitize(master_overdens, overdensbins)





##############################################
###  Counting number of objects per slitmask
##############################################

non_ORELSE_slitmasks = ['GIOIA', 
                        'LRIS', 
                        'LRIS.scm1', 
                        'LRIS.scm2', 
                        'LRIS.scm4', 
                        'LRIS.scnm2', 
                        'Scl1', 
                        'Scl2', 
                        'Sgr1', 
                        'Sgr2', 
                        'Sgr3', 
                        'T08_FOCAS', 
                        'T08_LRIS', 
                        'VVDS', 
                        'oldLRIS', 
                        'oldLRISe']


if True:

    n_orelse_targets = 0
    n_non_orelse_targets = 0

    unique_slitmasks = numpy.unique(master_slitmask)

    for slitmask in unique_slitmasks:

        if slitmask in non_ORELSE_slitmasks:
            n_non_orelse_targets += master_slitmask.tolist().count(slitmask)
        else:
            n_orelse_targets += master_slitmask.tolist().count(slitmask)

        print('%11s  %i' % (slitmask, master_slitmask.tolist().count(slitmask)))











#################################################
###  calculating weights for galaxies to map the
###  field/group z-dist. to the cluster z-dist.
#################################################

if False:

    zbins = hists[0][1]
    zbars = (zbins[1:] + zbins[:-1]) / 2.

    digi_z = numpy.digitize(master_zspec, zbins)


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








#############################################
###  calculating my own global densities  ###
#############################################

if True:

    n_galaxies = len(master_id_field)
    log_eta = numpy.zeros(n_galaxies) + 5

    substructure_catalog = ascii.read('../data/ORELSE_substructures.dat')

    substructure_centroids = coordinates.SkyCoord(substructure_catalog['RA'], substructure_catalog['Dec'], unit=units.deg, frame='fk5')
    R200 = 3**0.5 * substructure_catalog['sigma_1Mpc']*(units.km/units.s) / \
           (10 * cosmo.H(substructure_catalog['z']))
    substructure_velocities = mypy.vel_recess(substructure_catalog['z'])


    ###  Initializing progress bar
    widgets = ['  Running: ', progressbar.Percentage(), 
               ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
               ' ', progressbar.ETA(), 
               ' ', progressbar.FileTransferSpeed()]
    pbar = progressbar.ProgressBar(widgets=widgets, maxval=n_galaxies+1).start()


    for i_gal in range(n_galaxies):

        pbar.update(i_gal)

        ###  calculate Rproj to every ORELSE substructure (that's right, every single one)
        galaxy_centroid = coordinates.SkyCoord(master_ra[i_gal], master_dec[i_gal], unit=units.degree, frame='fk5')

        separations = galaxy_centroid.separation(substructure_centroids)
        separations = separations.to(units.arcmin)

        Rproj = separations * cosmo.kpc_proper_per_arcmin(substructure_catalog['z'])
        Rproj = Rproj.to(units.Mpc)


        ###  calculating recessional velocities
        galaxy_velocity = mypy.vel_recess(master_zspec[i_gal])
        delta_velocities = abs(galaxy_velocity - substructure_velocities)


        ###  calculating all possible etas
        eta_all = (Rproj / R200).value * (delta_velocities / substructure_catalog['sigma_1Mpc'])


        ###  1) Find all substructures within Rproj < 5 Mpc and adopt smallest eta
        inds_5Mpc = numpy.where(Rproj < 5*units.Mpc)[0]

        if len(inds_5Mpc) > 0:
            log_eta[i_gal] = numpy.log10(eta_all[inds_5Mpc].min())









##########################
###  saving to pickle  ###
##########################

if True:

    tablename = 'table_all_galaxies_noXrayAGN.pickle'

    header_text = '\n'
    header_text += '#  All galaxies in this table have been selected by these criteria:\n'
    header_text += '#\n'
    header_text += '#    1) Secure zspec between %.3f < z < %.3f\n' % (zlo, zhi)
    header_text += '#    2) MIPS 24um coverage (weight_mips24 > 0.05)\n'
    header_text += '#    3) Not an identified radio AGN (only for N200, N5281, SG0023, SC1324, SC1604)\n'


    outer = combined_catalog(header_text, len(master_id_field), lmassbins, overdensbins, overdens_labels)

    outer.field_names = field_names
    outer.field_versions = field_versions

    outer.id_field  = master_id_field
    outer.id_gal    = master_id_gal
    outer.ra        = master_ra
    outer.dec       = master_dec
    outer.zspec     = master_zspec
    outer.ftot_mips24 = master_ftot_mips24
    outer.etot_mips24 = master_etot_mips24
    outer.SFR_UVIR  = master_SFR_UVIR
    outer.SFR_FAST  = master_SFR_FAST
    outer.lmass     = master_lmass
    outer.ltau      = master_ltau
    outer.lage      = master_lage
    outer.overdens  = master_overdens
    outer.UV_color  = master_UV_color
    outer.VJ_color  = master_VJ_color
    outer.UVJ_class = master_UVJ_class
    #outer.log_eta   = log_eta

    outer.digi_lmass = digi_lmass
    outer.digi_overdens = digi_overdens


    pickle.dump(outer, open('../data/%s' % tablename, 'wb'))




































#####################################
###  Adding "old" stellar masses  ###
#####################################

if False:

    class combined_catalog(object):
        def __init__(self, header, n_galaxies):

            self.time_stamp = time.ctime()
            self.header = header

            self.n_galaxies = n_galaxies



    galaxies_catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))



    lmass_old = numpy.zeros(galaxies_catalog.n_galaxies)
    lmass_new = numpy.zeros(galaxies_catalog.n_galaxies)

    for i_gal in range(galaxies_catalog.n_galaxies):

        i_field = int(galaxies_catalog.id_field[i_gal])
        f = fields0[i_field]

        id_gal = int(galaxies_catalog.id_gal[i_gal])
        lmass_old[i_gal] = f.fout_old.lmass[id_gal-1]
        lmass_new[i_gal] = f.fout.lmass[id_gal-1]



    galaxies_catalog.lmass_old = lmass_old

    pickle.dump(galaxies_catalog, open('../data/table_all_galaxies_withOldMasses.pickle', 'wb'))





















