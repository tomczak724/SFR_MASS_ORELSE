
import mypy
import numpy
import pickle
import progressbar
from scipy import optimize, stats
from astropy import coordinates, units, cosmology
from astropy.io import fits, ascii
from matplotlib import pyplot
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages

cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)


def log_SFR_tomczak2016(z, lmass, SF_galaxies=True):
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

    if SF_galaxies:
        s0_params = [0.448, 1.220, -0.174]
        m0_params = [9.458, 0.865, -0.132]
        gamma = 1.091
    else:
        s0_params = [0.195, 1.157, -0.143]
        m0_params = [9.244, 0.753, -0.090]
        gamma = 1.118

    s0 = s0_params[0] + s0_params[1] * z + s0_params[2] * z**2
    log_m0 = m0_params[0] + m0_params[1] * z + m0_params[2] * z**2
    return s0 - numpy.log10(1 + (10**lmass / 10**log_m0)**-gamma)


class combined_catalog(object):
    def __init__(self, header, n_galaxies):

        self.time_stamp = time.ctime()
        self.header = header

        self.n_galaxies = n_galaxies



galaxies_catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))


'''
###  inserting global densities, log(eta), into "catalog" object
galaxies_catalog.log_eta = numpy.zeros(galaxies_catalog.n_galaxies) + 5


###  grabbing global densities from Brian's catalog
if False:

    data_dir = '/Users/atomczak/Dropbox/ORELSE/DATA'
    data_file = 'allORELSE.globaldensity.eta.wIDRAdeczqflagqnimagSMge-100.0_MNUVrge-100.0_magcuts-100.0to50.0.cat'
    global_density_catalog = ascii.read('%s/%s' % (data_dir, data_file))

    n_matches = numpy.zeros(galaxies_catalog.n_galaxies)
    z_matches = numpy.zeros(galaxies_catalog.n_galaxies)

    outer = open('../data/missing_galaxies_from_eta_catalog.dat', 'w')
    outer.write('# field  photID  zspec\n')

    for i_gal in range(catalog.n_galaxies):

        f = catalog.field_names[catalog.id_field[i_gal]]
        id_gal = catalog.id_gal[i_gal]
        

        ###  Finding galaxies in the corresponding ORELSE field, and ...
        inds_field = numpy.where(global_density_catalog['field'] == f)[0]


        ###  ... finding the corresponding galaxy itself
        ind_gal = numpy.where(global_density_catalog['photID'][inds_field] == id_gal)[0]



        n_matches[i_gal] = len(ind_gal)
        try:
            z_matches[i_gal] = global_density_catalog['zspec'][inds_field][ind_gal[0]]
        except:
            pass


        if len(ind_gal) == 0:
            outer.write('%9s' % f)
            outer.write(' %7i' % id_gal)
            outer.write('   %.5f' % catalog.z[i_gal])
            outer.write('\n')

        outer.close()


###  calculating my own global densities
if True:

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
    pbar = progressbar.ProgressBar(widgets=widgets, maxval=galaxies_catalog.n_galaxies+1).start()


    for i_gal in range(galaxies_catalog.n_galaxies):

        pbar.update(i_gal)

        ###  calculate Rproj to every ORELSE substructure (that's right, every single one)
        galaxy_centroid = coordinates.SkyCoord(galaxies_catalog.ra[i_gal], galaxies_catalog.dec[i_gal], unit=units.degree, frame='fk5')

        separations = galaxy_centroid.separation(substructure_centroids)
        separations = separations.to(units.arcmin)

        Rproj = separations * cosmo.kpc_proper_per_arcmin(substructure_catalog['z'])
        Rproj = Rproj.to(units.Mpc)


        ###  calculating recessional velocities
        galaxy_velocity = mypy.vel_recess(galaxies_catalog.zspec[i_gal])
        delta_velocities = abs(galaxy_velocity - substructure_velocities)


        ###  calculating all possible etas
        eta_all = (Rproj / R200).value * (delta_velocities / substructure_catalog['sigma_1Mpc'])


        ###  1) Find all substructures within Rproj < 5 Mpc and adopt smallest eta
        inds_5Mpc = numpy.where(Rproj < 5*units.Mpc)[0]

        if len(inds_5Mpc) > 0:
            galaxies_catalog.log_eta[i_gal] = numpy.log10(eta_all[inds_5Mpc].min())
'''







###  binning data by log(eta)

eta_bins = numpy.array([0.1, 0.4, 2.0])
eta_labels = ['          $\eta$ < %.1f' % eta_bins[0],
              '%.1f < $\eta$ < %.1f' % (eta_bins[0], eta_bins[1]),
              '%.1f < $\eta$ < %.1f' % (eta_bins[1], eta_bins[2]),
              '%.1f < $\eta$' % eta_bins[2]]

digi_eta = numpy.digitize(10**galaxies_catalog.log_eta, eta_bins)







#############################################
###  Plotting SF and TOT panels together
###  Applying redshift-offsets:
###     I noticed that the redshift
###     distributions for the high-
###     density sample are sometimes
###     skewed to lower-z. Therefore
###     to account for this we will
###     be applying an offset based
###     on the ZFOURGE SFR-M* relations
#############################################

if True:



    pdf_zdists = PdfPages('../figures/redshift_distributions_weighted.pdf')



    fig_zdists = pyplot.figure(figsize=(16.5, 11.))

    sp2_zdists = fig_zdists.add_subplot(222)
    sp1_zdists = fig_zdists.add_subplot(221)
    sp4_zdists = fig_zdists.add_subplot(224)
    sp3_zdists = fig_zdists.add_subplot(223)


    fig_zdists.subplots_adjust(left=0.07, right=0.98, hspace=0, wspace=0)





    cmap = pyplot.cm.rainbow


    colors = ['r', 'orange', 'g', 'b']

    n_galaxies_tot = []
    z_galaxies_tot = []
    overdens_median_tot = []
    overdens_nmad_tot = []
    sfr_median_tot = []
    sfr_nmad_tot = []
    lmass_median_tot = []
    lmass_nmad_tot = []

    n_galaxies_sf = []
    z_galaxies_sf = []
    overdens_median_sf = []
    overdens_nmad_sf = []
    sfr_median_sf = []
    sfr_nmad_sf = []
    lmass_median_sf = []
    lmass_nmad_sf = []


    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        sp1_zdists.minorticks_on()
        sp2_zdists.minorticks_on()
        sp3_zdists.minorticks_on()
        sp4_zdists.minorticks_on()
        sp1_zdists.set_xlabel('spectroscopic redshift')
        sp2_zdists.set_xlabel('spectroscopic redshift')
        sp3_zdists.set_xlabel('spectroscopic redshift')
        sp4_zdists.set_xlabel('spectroscopic redshift')
        sp1_zdists.set_ylabel('Number')
        sp2_zdists.set_ylabel('Number')
        sp3_zdists.set_ylabel('CDF')
        sp4_zdists.set_ylabel('CDF')

        sp1_zdists.set_title('All Galaxies')
        sp2_zdists.set_title('Star-Forming Galaxies')

        sp3_zdists.axis([0.55, 1.35, 0, 1.1])
        sp4_zdists.axis([0.55, 1.35, 0, 1.1])


        n_galaxies_tot.append([])
        z_galaxies_tot.append([])
        overdens_median_tot.append([])
        overdens_nmad_tot.append([])
        sfr_median_tot.append([])
        sfr_nmad_tot.append([])
        lmass_median_tot.append([])
        lmass_nmad_tot.append([])

        n_galaxies_sf.append([])
        z_galaxies_sf.append([])
        overdens_median_sf.append([])
        overdens_nmad_sf.append([])
        sfr_median_sf.append([])
        sfr_nmad_sf.append([])
        lmass_median_sf.append([])
        lmass_nmad_sf.append([])

       #for i_overdens in range(len(eta_bins)+1)[::-1]:
        for i_overdens in range(len(eta_bins)+1):

            inds_lmass_overdens_tot = (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                      (digi_eta == i_overdens)

            inds_lmass_overdens_sf = (galaxies_catalog.UVJ_class == 1) & \
                                     (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                     (digi_eta == i_overdens)


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

            n_galaxies_tot[-1].append(ngals_tot)
            z_galaxies_tot[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist())
            sfr_median_tot[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            overdens_median_tot[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_tot]))
            overdens_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_tot]))
            lmass_median_tot[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf[-1].append(ngals_sf)
            z_galaxies_sf[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
            sfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            overdens_median_sf[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
            overdens_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
            lmass_median_sf[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))





            ###  plotting redshift distributions
            zlo, zhi = 0.7, 1.3
            offsets = [-0.002, -0.0007, 0.0007, 0.002]        
            color_overdens = colors[i_overdens]



            h_tot = sp1_zdists.hist(z_galaxies_tot[-1][-1], histtype='step', lw=3, color=color_overdens,
                             bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_tot[-1][-1])))

            h_sf = sp2_zdists.hist(z_galaxies_sf[-1][-1], histtype='step', lw=3, color=color_overdens,
                            bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_sf[-1][-1])))


            sp3_zdists.hist(z_galaxies_tot[-1][-1], histtype='step', lw=3, color=color_overdens,
                     bins=1000, range=(0.7, 1.3), cumulative=True, normed=True)

            sp4_zdists.hist(z_galaxies_sf[-1][-1], histtype='step', lw=3, color=color_overdens,
                     bins=1000, range=(0.7, 1.3), cumulative=True, normed=True)



        ###  organizing redshift distribution plots
        sp1_zdists.set_xlim(sp3_zdists.get_xlim())
        sp2_zdists.axis(sp1_zdists.axis())

        sp1_zdists.legend(loc=1)
        sp2_zdists.legend(loc=1)


        ###  KS tests
        ks_field_groups_tot = stats.ks_2samp(z_galaxies_tot[-1][0], z_galaxies_tot[-1][1])
        ks_field_clusters_tot = stats.ks_2samp(z_galaxies_tot[-1][0], z_galaxies_tot[-1][2])

        ks_field_groups_sf = stats.ks_2samp(z_galaxies_sf[-1][0], z_galaxies_sf[-1][1])
        ks_field_clusters_sf = stats.ks_2samp(z_galaxies_sf[-1][0], z_galaxies_sf[-1][2])



        ###  plotting text
        info = 'log(M$_*$) = %.1f' % galaxies_catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_tot[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_tot[1]

        sp3_zdists.text(0.03, 0.97, info, transform=sp3_zdists.transAxes, 
                 horizontalalignment='left', verticalalignment='top')

        info = 'log(M$_*$) = %.1f' % galaxies_catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_sf[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_sf[1]

        sp4_zdists.text(0.03, 0.97, info, transform=sp4_zdists.transAxes, 
                 horizontalalignment='left', verticalalignment='top')





        pdf_zdists.savefig()
        sp1_zdists.clear()
        sp2_zdists.clear()
        sp3_zdists.clear()
        sp4_zdists.clear()

    pdf_zdists.close()
    pyplot.close()

















    ###  estimating SFR offsets from median to fiducial redshifts
    z_fiducial = numpy.median(galaxies_catalog.zspec)

    lsfr_offsets_tot = []
    lsfr_offsets_sf = []

    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        lsfr_offsets_tot.append([])
        lsfr_offsets_sf.append([])

        for i_overdens in range(len(eta_bins)+1)[::-1]:

            z_medi_tot = numpy.median(z_galaxies_tot[i_lmass][i_overdens])
            z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])

            lsfr_fiducial_tot = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_z_medi_tot = log_SFR_tomczak2016(z_medi_tot, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_offsets_tot[-1].append(lsfr_fiducial_tot - lsfr_z_medi_tot)
            lsfr_offsets_sf[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)

    lsfr_offsets_tot = numpy.array(lsfr_offsets_tot)
    lsfr_offsets_sf = numpy.array(lsfr_offsets_sf)











    ###  plotting SFR vs. lmass
    lmass_median_tot = numpy.array(lmass_median_tot)
    lmass_nmad_tot = numpy.array(lmass_nmad_tot)
    sfr_median_tot = numpy.array(sfr_median_tot)
    sfr_nmad_tot = numpy.array(sfr_nmad_tot)
    overdens_median_tot = numpy.array(overdens_median_tot)
    overdens_nmad_tot = numpy.array(overdens_nmad_tot)
    n_galaxies_tot = numpy.array(n_galaxies_tot)

    lmass_median_sf = numpy.array(lmass_median_sf)
    lmass_nmad_sf = numpy.array(lmass_nmad_sf)
    sfr_median_sf = numpy.array(sfr_median_sf)
    sfr_nmad_sf = numpy.array(sfr_nmad_sf)
    overdens_median_sf = numpy.array(overdens_median_sf)
    overdens_nmad_sf = numpy.array(overdens_nmad_sf)
    n_galaxies_sf = numpy.array(n_galaxies_sf)



    ###  applying SFR offsets
    lsfr_median_tot = numpy.log10(sfr_median_tot) + lsfr_offsets_tot
    lsfr_median_sf = numpy.log10(sfr_median_sf) + lsfr_offsets_sf








    fig_sfrmass = pyplot.figure(figsize=(16.5, 8.1))

    sp2_sfrmass = fig_sfrmass.add_subplot(122)
    sp1_sfrmass = fig_sfrmass.add_subplot(121)

    sp1_sfrmass.minorticks_on()
    sp2_sfrmass.minorticks_on()
    sp1_sfrmass.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp2_sfrmass.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp1_sfrmass.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')
    sp2_sfrmass.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1_sfrmass.axis([8.6, 12.3, -0.4, 2.1])
    sp1_sfrmass.axis([8.1, 12.4, -0.85, 2.05])
    sp2_sfrmass.axis([8.1, 12.4, -0.85, 2.05])


    fig_sfrmass.subplots_adjust(wspace=0, left=0.07, right=0.98, )



    for i_overdens in range(len(eta_bins)+1)[::-1]:


        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_tot[:, i_overdens] > 8)
        inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf[:, i_overdens] > 8)


        elog_sfr_tot = numpy.log10((sfr_median_tot[:, i_overdens] + sfr_nmad_tot[:, i_overdens]/n_galaxies_tot[:, i_overdens]**0.5) / sfr_median_tot[:, i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf[:, i_overdens] + sfr_nmad_sf[:, i_overdens]/n_galaxies_sf[:, i_overdens]**0.5) / sfr_median_sf[:, i_overdens])

        color_overdens = colors[i_overdens]

        sp1_sfrmass.errorbar(lmass_median_tot[inds_plot_tot, i_overdens], lsfr_median_tot[inds_plot_tot, i_overdens], 
                     yerr=elog_sfr_tot[inds_plot_tot], xerr=lmass_nmad_tot[inds_plot_tot, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=eta_labels[i_overdens])

        sp2_sfrmass.errorbar(lmass_median_sf[inds_plot_sf, i_overdens], lsfr_median_sf[inds_plot_sf, i_overdens], 
                     yerr=elog_sfr_sf[inds_plot_sf], xerr=lmass_nmad_sf[inds_plot_sf, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=eta_labels[i_overdens])







    ###  plotting ZFOURGE

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 16), lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 84), lmass_model, SF_galaxies=False)

    sp1_sfrmass.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1_sfrmass.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec), zorder=1)

    print 'SF:  %.3f < z < %.3f' % (numpy.percentile(galaxies_catalog.zspec, 16), numpy.percentile(galaxies_catalog.zspec, 84))



    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 16), lmass_model, SF_galaxies=True)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 84), lmass_model, SF_galaxies=True)

    sp2_sfrmass.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2_sfrmass.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec), zorder=1)

    print 'ALL: %.3f < z < %.3f' % (numpy.percentile(galaxies_catalog.zspec, 16), numpy.percentile(galaxies_catalog.zspec, 84))





    t1 = sp1_sfrmass.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp1_sfrmass.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t2 = sp2_sfrmass.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp2_sfrmass.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_sfrmass.legend(loc=4, numpoints=1, fontsize=17, frameon=0)











#############################################
###  Plotting 4-panel of SFR vs. Mass:
###
###      All_galaxies         SF_galaxies
###      local_density        local_density
###
###      All_galaxies         SF_galaxies
###      global_density       global_density
###
#############################################

if True:


    figname = '../figures/redshift_distributions_weighted.pdf'
    pdf_zdists = PdfPages(figname)



    fig_zdists = pyplot.figure(figsize=(16.5, 11.))

    sp2_zdists = fig_zdists.add_subplot(222)
    sp1_zdists = fig_zdists.add_subplot(221)
    sp4_zdists = fig_zdists.add_subplot(224)
    sp3_zdists = fig_zdists.add_subplot(223)


    fig_zdists.subplots_adjust(left=0.07, right=0.98, hspace=0, wspace=0)





    cmap = pyplot.cm.rainbow


    colors = ['r', 'orange', 'g', 'b']

    n_galaxies_tot = []
    z_galaxies_tot = []
    overdens_median_tot = []
    overdens_nmad_tot = []
    sfr_median_tot = []
    sfr_nmad_tot = []
    lmass_median_tot = []
    lmass_nmad_tot = []

    n_galaxies_sf = []
    z_galaxies_sf = []
    overdens_median_sf = []
    overdens_nmad_sf = []
    sfr_median_sf = []
    sfr_nmad_sf = []
    lmass_median_sf = []
    lmass_nmad_sf = []



    ###  calculating median SFRs in GLOBAL overdensity bins
    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        sp1_zdists.minorticks_on()
        sp2_zdists.minorticks_on()
        sp3_zdists.minorticks_on()
        sp4_zdists.minorticks_on()
        sp1_zdists.set_xlabel('spectroscopic redshift')
        sp2_zdists.set_xlabel('spectroscopic redshift')
        sp3_zdists.set_xlabel('spectroscopic redshift')
        sp4_zdists.set_xlabel('spectroscopic redshift')
        sp1_zdists.set_ylabel('Number')
        sp2_zdists.set_ylabel('Number')
        sp3_zdists.set_ylabel('CDF')
        sp4_zdists.set_ylabel('CDF')

        sp1_zdists.set_title('All Galaxies')
        sp2_zdists.set_title('Star-Forming Galaxies')

        sp3_zdists.axis([0.55, 1.35, 0, 1.1])
        sp4_zdists.axis([0.55, 1.35, 0, 1.1])


        n_galaxies_tot.append([])
        z_galaxies_tot.append([])
        overdens_median_tot.append([])
        overdens_nmad_tot.append([])
        sfr_median_tot.append([])
        sfr_nmad_tot.append([])
        lmass_median_tot.append([])
        lmass_nmad_tot.append([])

        n_galaxies_sf.append([])
        z_galaxies_sf.append([])
        overdens_median_sf.append([])
        overdens_nmad_sf.append([])
        sfr_median_sf.append([])
        sfr_nmad_sf.append([])
        lmass_median_sf.append([])
        lmass_nmad_sf.append([])

        for i_overdens in range(len(eta_bins)+1):

            inds_lmass_overdens_tot = (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                      (digi_eta == i_overdens)

            inds_lmass_overdens_sf = (galaxies_catalog.UVJ_class == 1) & \
                                     (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                     (digi_eta == i_overdens)


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

            n_galaxies_tot[-1].append(ngals_tot)
            z_galaxies_tot[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist())
            sfr_median_tot[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            overdens_median_tot[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_tot]))
            overdens_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_tot]))
            lmass_median_tot[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf[-1].append(ngals_sf)
            z_galaxies_sf[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
            sfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            overdens_median_sf[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
            overdens_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
            lmass_median_sf[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))





            ###  plotting redshift distributions
            zlo, zhi = 0.7, 1.3
            offsets = [-0.002, -0.0007, 0.0007, 0.002]        
            color_overdens = colors[i_overdens]



            h_tot = sp1_zdists.hist(z_galaxies_tot[-1][-1], histtype='step', lw=3, color=color_overdens,
                             bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_tot[-1][-1])))

            h_sf = sp2_zdists.hist(z_galaxies_sf[-1][-1], histtype='step', lw=3, color=color_overdens,
                            bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_sf[-1][-1])))


            sp3_zdists.hist(z_galaxies_tot[-1][-1], histtype='step', lw=3, color=color_overdens,
                     bins=1000, range=(0.7, 1.3), cumulative=True, normed=True)

            sp4_zdists.hist(z_galaxies_sf[-1][-1], histtype='step', lw=3, color=color_overdens,
                     bins=1000, range=(0.7, 1.3), cumulative=True, normed=True)



        ###  organizing redshift distribution plots
        sp1_zdists.set_xlim(sp3_zdists.get_xlim())
        sp2_zdists.axis(sp1_zdists.axis())

        sp1_zdists.legend(loc=1)
        sp2_zdists.legend(loc=1)


        ###  KS tests
        ks_field_groups_tot = stats.ks_2samp(z_galaxies_tot[-1][0], z_galaxies_tot[-1][1])
        ks_field_clusters_tot = stats.ks_2samp(z_galaxies_tot[-1][0], z_galaxies_tot[-1][2])

        ks_field_groups_sf = stats.ks_2samp(z_galaxies_sf[-1][0], z_galaxies_sf[-1][1])
        ks_field_clusters_sf = stats.ks_2samp(z_galaxies_sf[-1][0], z_galaxies_sf[-1][2])



        ###  plotting text
        info = 'log(M$_*$) = %.1f' % galaxies_catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_tot[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_tot[1]

        sp3_zdists.text(0.03, 0.97, info, transform=sp3_zdists.transAxes, 
                 horizontalalignment='left', verticalalignment='top')

        info = 'log(M$_*$) = %.1f' % galaxies_catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_sf[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_sf[1]

        sp4_zdists.text(0.03, 0.97, info, transform=sp4_zdists.transAxes, 
                 horizontalalignment='left', verticalalignment='top')





        pdf_zdists.savefig()
        sp1_zdists.clear()
        sp2_zdists.clear()
        sp3_zdists.clear()
        sp4_zdists.clear()

    pdf_zdists.close()
    pyplot.close()
    print 'wrote to:\n%s' % figname
















    n_galaxies_tot_local = []
    z_galaxies_tot_local = []
    overdens_median_tot_local = []
    overdens_nmad_tot_local = []
    sfr_median_tot_local = []
    sfr_nmad_tot_local = []
    lmass_median_tot_local = []
    lmass_nmad_tot_local = []

    n_galaxies_sf_local = []
    z_galaxies_sf_local = []
    overdens_median_sf_local = []
    overdens_nmad_sf_local = []
    sfr_median_sf_local = []
    sfr_nmad_sf_local = []
    lmass_median_sf_local = []
    lmass_nmad_sf_local = []

    ###  calculating median SFRs in LOCAL overdensity bins
    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        n_galaxies_tot_local.append([])
        z_galaxies_tot_local.append([])
        overdens_median_tot_local.append([])
        overdens_nmad_tot_local.append([])
        sfr_median_tot_local.append([])
        sfr_nmad_tot_local.append([])
        lmass_median_tot_local.append([])
        lmass_nmad_tot_local.append([])

        n_galaxies_sf_local.append([])
        z_galaxies_sf_local.append([])
        overdens_median_sf_local.append([])
        overdens_nmad_sf_local.append([])
        sfr_median_sf_local.append([])
        sfr_nmad_sf_local.append([])
        lmass_median_sf_local.append([])
        lmass_nmad_sf_local.append([])

        for i_overdens in range(len(galaxies_catalog.overdensbars)):

            inds_lmass_overdens_tot = numpy.where((galaxies_catalog.digi_lmass == i_lmass+1) & \
                                                  (galaxies_catalog.digi_overdens == i_overdens+1))[0]

            inds_lmass_overdens_sf = numpy.where((galaxies_catalog.UVJ_class == 1) & \
                                                 (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                                 (galaxies_catalog.digi_overdens == i_overdens+1))[0]


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

            n_galaxies_tot_local[-1].append(ngals_tot)
            z_galaxies_tot_local[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist())
            sfr_median_tot_local[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot_local[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            overdens_median_tot_local[-1].append(numpy.median(galaxies_catalog.overdens[inds_lmass_overdens_tot]))
            overdens_nmad_tot_local[-1].append(mypy.nmad(galaxies_catalog.overdens[inds_lmass_overdens_tot]))
            lmass_median_tot_local[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot_local[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf_local[-1].append(ngals_sf)
            z_galaxies_sf_local[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
            sfr_median_sf_local[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf_local[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            overdens_median_sf_local[-1].append(numpy.median(galaxies_catalog.overdens[inds_lmass_overdens_sf]))
            overdens_nmad_sf_local[-1].append(mypy.nmad(galaxies_catalog.overdens[inds_lmass_overdens_sf]))
            lmass_median_sf_local[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf_local[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))

































    ###  estimating SFR offsets from median to fiducial redshifts
    z_fiducial = numpy.median(galaxies_catalog.zspec)

    lsfr_offsets_tot = []
    lsfr_offsets_sf = []

    lsfr_offsets_tot_local = []
    lsfr_offsets_sf_local = []

    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        lsfr_offsets_tot.append([])
        lsfr_offsets_sf.append([])

        lsfr_offsets_tot_local.append([])
        lsfr_offsets_sf_local.append([])

        for i_overdens in range(len(eta_bins)+1)[::-1]:

            ###  GLOBAL overdensity
            z_medi_tot = numpy.median(z_galaxies_tot[i_lmass][i_overdens])
            z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])

            lsfr_fiducial_tot = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_z_medi_tot = log_SFR_tomczak2016(z_medi_tot, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_offsets_tot[-1].append(lsfr_fiducial_tot - lsfr_z_medi_tot)
            lsfr_offsets_sf[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)


        for i_overdens in range(len(galaxies_catalog.overdensbars)):

            ###  LOCAL overdensity
            z_medi_tot = numpy.median(z_galaxies_tot_local[i_lmass][i_overdens])
            z_medi_sf = numpy.median(z_galaxies_sf_local[i_lmass][i_overdens])

            lsfr_fiducial_tot = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_z_medi_tot = log_SFR_tomczak2016(z_medi_tot, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_offsets_tot_local[-1].append(lsfr_fiducial_tot - lsfr_z_medi_tot)
            lsfr_offsets_sf_local[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)


    lsfr_offsets_tot = numpy.array(lsfr_offsets_tot)
    lsfr_offsets_sf = numpy.array(lsfr_offsets_sf)

    lsfr_offsets_tot_local = numpy.array(lsfr_offsets_tot_local)
    lsfr_offsets_sf_local = numpy.array(lsfr_offsets_sf_local)













    ###  plotting SFR vs. lmass
    lmass_median_tot = numpy.array(lmass_median_tot)
    lmass_nmad_tot = numpy.array(lmass_nmad_tot)
    sfr_median_tot = numpy.array(sfr_median_tot)
    sfr_nmad_tot = numpy.array(sfr_nmad_tot)
    overdens_median_tot = numpy.array(overdens_median_tot)
    overdens_nmad_tot = numpy.array(overdens_nmad_tot)
    n_galaxies_tot = numpy.array(n_galaxies_tot)

    lmass_median_sf = numpy.array(lmass_median_sf)
    lmass_nmad_sf = numpy.array(lmass_nmad_sf)
    sfr_median_sf = numpy.array(sfr_median_sf)
    sfr_nmad_sf = numpy.array(sfr_nmad_sf)
    overdens_median_sf = numpy.array(overdens_median_sf)
    overdens_nmad_sf = numpy.array(overdens_nmad_sf)
    n_galaxies_sf = numpy.array(n_galaxies_sf)

    lmass_median_tot_local = numpy.array(lmass_median_tot_local)
    lmass_nmad_tot_local = numpy.array(lmass_nmad_tot_local)
    sfr_median_tot_local = numpy.array(sfr_median_tot_local)
    sfr_nmad_tot_local = numpy.array(sfr_nmad_tot_local)
    overdens_median_tot_local = numpy.array(overdens_median_tot_local)
    overdens_nmad_tot_local = numpy.array(overdens_nmad_tot_local)
    n_galaxies_tot_local = numpy.array(n_galaxies_tot_local)

    lmass_median_sf_local = numpy.array(lmass_median_sf_local)
    lmass_nmad_sf_local = numpy.array(lmass_nmad_sf_local)
    sfr_median_sf_local = numpy.array(sfr_median_sf_local)
    sfr_nmad_sf_local = numpy.array(sfr_nmad_sf_local)
    overdens_median_sf_local = numpy.array(overdens_median_sf_local)
    overdens_nmad_sf_local = numpy.array(overdens_nmad_sf_local)
    n_galaxies_sf_local = numpy.array(n_galaxies_sf_local)



    ###  applying SFR offsets
    lsfr_median_tot = numpy.log10(sfr_median_tot) + lsfr_offsets_tot
    lsfr_median_sf = numpy.log10(sfr_median_sf) + lsfr_offsets_sf

    lsfr_median_tot_local = numpy.log10(sfr_median_tot_local) + lsfr_offsets_tot_local
    lsfr_median_sf_local = numpy.log10(sfr_median_sf_local) + lsfr_offsets_sf_local








    fig_sfrmass = pyplot.figure(figsize=(15.9, 15.3))

    sp2_sfrmass_global = fig_sfrmass.add_subplot(222)
    sp1_sfrmass_global = fig_sfrmass.add_subplot(221)
    sp2_sfrmass_local = fig_sfrmass.add_subplot(224)
    sp1_sfrmass_local = fig_sfrmass.add_subplot(223)


    sp1_sfrmass_global.grid()
    sp2_sfrmass_global.grid()
    sp1_sfrmass_local.grid()
    sp2_sfrmass_local.grid()

    sp1_sfrmass_global.minorticks_on()
    sp2_sfrmass_global.minorticks_on()
    sp1_sfrmass_local.minorticks_on()
    sp2_sfrmass_local.minorticks_on()

    sp1_sfrmass_local.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=22)
    sp2_sfrmass_local.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=22)
    sp1_sfrmass_local.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )', size=22)
    sp1_sfrmass_global.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )', size=22)

    sp1_sfrmass_global.axis([8.1, 12.4, -0.85, 2.05])
    sp2_sfrmass_global.axis([8.1, 12.4, -0.85, 2.05])
    sp1_sfrmass_local.axis([8.1, 12.4, -0.85, 2.05])
    sp2_sfrmass_local.axis([8.1, 12.4, -0.85, 2.05])


    fig_sfrmass.subplots_adjust(wspace=0, hspace=0, left=0.08, right=0.98, bottom=0.06)



    for i_overdens in range(len(eta_bins)+1)[::-1]:

        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_tot[:, i_overdens] > 8)
        inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf[:, i_overdens] > 8)


        elog_sfr_tot = numpy.log10((sfr_median_tot[:, i_overdens] + sfr_nmad_tot[:, i_overdens]/n_galaxies_tot[:, i_overdens]**0.5) / sfr_median_tot[:, i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf[:, i_overdens] + sfr_nmad_sf[:, i_overdens]/n_galaxies_sf[:, i_overdens]**0.5) / sfr_median_sf[:, i_overdens])

        color_overdens = colors[i_overdens]

        sp1_sfrmass_global.errorbar(lmass_median_tot[inds_plot_tot, i_overdens], lsfr_median_tot[inds_plot_tot, i_overdens], 
                     yerr=elog_sfr_tot[inds_plot_tot], xerr=lmass_nmad_tot[inds_plot_tot, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=eta_labels[i_overdens])

        sp2_sfrmass_global.errorbar(lmass_median_sf[inds_plot_sf, i_overdens], lsfr_median_sf[inds_plot_sf, i_overdens], 
                     yerr=elog_sfr_sf[inds_plot_sf], xerr=lmass_nmad_sf[inds_plot_sf, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=eta_labels[i_overdens])


    colors_local = ['b', '#7ead1a', 'r']
    colors_local = ['#001ddd', '#20c600', '#d11f1f']
    overdens_labels = ['          log(1+$\\delta$) < 0.5',
                       '0.5 < log(1+$\\delta$) < 1.0',
                       '1.0 < log(1+$\\delta$)']


    for i_overdens in range(len(galaxies_catalog.overdensbars)):

        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_tot_local[:, i_overdens] > 8)
        inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf_local[:, i_overdens] > 8)


        elog_sfr_tot = numpy.log10((sfr_median_tot_local[:, i_overdens] + sfr_nmad_tot_local[:, i_overdens]/n_galaxies_tot_local[:, i_overdens]**0.5) / sfr_median_tot_local[:, i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf_local[:, i_overdens] + sfr_nmad_sf_local[:, i_overdens]/n_galaxies_sf_local[:, i_overdens]**0.5) / sfr_median_sf_local[:, i_overdens])

        color_overdens = colors_local[i_overdens]

        sp1_sfrmass_local.errorbar(lmass_median_tot_local[inds_plot_tot, i_overdens], lsfr_median_tot_local[inds_plot_tot, i_overdens], 
                     yerr=elog_sfr_tot[inds_plot_tot], xerr=lmass_nmad_tot_local[inds_plot_tot, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=overdens_labels[i_overdens])

        sp2_sfrmass_local.errorbar(lmass_median_sf_local[inds_plot_sf, i_overdens], lsfr_median_sf_local[inds_plot_sf, i_overdens], 
                     yerr=elog_sfr_sf[inds_plot_sf], xerr=lmass_nmad_sf_local[inds_plot_sf, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=overdens_labels[i_overdens])




    ###  plotting ZFOURGE

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 16), lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 84), lmass_model, SF_galaxies=False)

    sp1_sfrmass_global.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1_sfrmass_global.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec), zorder=1)
    sp1_sfrmass_local.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1_sfrmass_local.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec), zorder=1)

    print 'SF:  %.3f < z < %.3f' % (numpy.percentile(galaxies_catalog.zspec, 16), numpy.percentile(galaxies_catalog.zspec, 84))



    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 16), lmass_model, SF_galaxies=True)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(galaxies_catalog.zspec, 84), lmass_model, SF_galaxies=True)

    sp2_sfrmass_global.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2_sfrmass_global.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec), zorder=1)
    sp2_sfrmass_local.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2_sfrmass_local.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec), zorder=1)

    print 'ALL: %.3f < z < %.3f' % (numpy.percentile(galaxies_catalog.zspec, 16), numpy.percentile(galaxies_catalog.zspec, 84))





    t1 = sp1_sfrmass_global.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp1_sfrmass_global.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t2 = sp2_sfrmass_global.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp2_sfrmass_global.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t3 = sp1_sfrmass_local.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp1_sfrmass_local.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t4 = sp2_sfrmass_local.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp2_sfrmass_local.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_sfrmass_global.legend(loc=4, title='Global Overdensity', numpoints=1, fontsize=17, frameon=0)
    leg = sp2_sfrmass_local.legend(loc=4, title='Local Overdensity', numpoints=1, fontsize=17, frameon=0)























    ##################################################
    ###  Plotting only for LOCAL overdensity bins  ###
    ##################################################

    fig_sfrmass = pyplot.figure(figsize=(14.5, 7.))

    sp2_sfrmass_local = fig_sfrmass.add_subplot(122)
    sp1_sfrmass_local = fig_sfrmass.add_subplot(121)


    sp1_sfrmass_local.grid()
    sp2_sfrmass_local.grid()

    sp1_sfrmass_local.minorticks_on()
    sp2_sfrmass_local.minorticks_on()

    sp1_sfrmass_local.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=22)
    sp2_sfrmass_local.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=22)
    sp1_sfrmass_local.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )', size=22)
    sp2_sfrmass_local.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )', size=22)

    sp1_sfrmass_local.axis([8.1, 12.25, -0.85, 2.05])
    sp2_sfrmass_local.axis([8.1, 12.25, -0.85, 2.05])


    fig_sfrmass.subplots_adjust(left=0.08, right=0.99, top=0.99, bottom=0.1, wspace=0)



    colors_local = ['#001ddd', '#20c600', '#d11f1f']
    colors_local = ['#3f4ea1', '#6db388', '#d92120']
    overdens_labels = ['          log(1+$\\delta$) < 0.5',
                       '0.5 < log(1+$\\delta$) < 1.0',
                       '1.0 < log(1+$\\delta$)']


    for i_overdens in range(len(galaxies_catalog.overdensbars)):

        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_tot_local[:, i_overdens] > 8)
        inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf_local[:, i_overdens] > 8)


        elog_sfr_tot = numpy.log10((sfr_median_tot_local[:, i_overdens] + sfr_nmad_tot_local[:, i_overdens]/n_galaxies_tot_local[:, i_overdens]**0.5) / sfr_median_tot_local[:, i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf_local[:, i_overdens] + sfr_nmad_sf_local[:, i_overdens]/n_galaxies_sf_local[:, i_overdens]**0.5) / sfr_median_sf_local[:, i_overdens])

        color_overdens = colors_local[i_overdens]

        sp1_sfrmass_local.errorbar(lmass_median_tot_local[inds_plot_tot, i_overdens], lsfr_median_tot_local[inds_plot_tot, i_overdens], 
                     yerr=elog_sfr_tot[inds_plot_tot], xerr=lmass_nmad_tot_local[inds_plot_tot, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=overdens_labels[i_overdens])

        sp2_sfrmass_local.errorbar(lmass_median_sf_local[inds_plot_sf, i_overdens], lsfr_median_sf_local[inds_plot_sf, i_overdens], 
                     yerr=elog_sfr_sf[inds_plot_sf], xerr=lmass_nmad_sf_local[inds_plot_sf, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=overdens_labels[i_overdens])




    ###  plotting ZFOURGE

    z16 = numpy.percentile(galaxies_catalog.zspec, 16)
    z84 = numpy.percentile(galaxies_catalog.zspec, 84)

    zfourge_label = 'ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec)
    zfourge_label = 'Tomczak et al. 2016\nz ~ %.2f' % numpy.median(galaxies_catalog.zspec)

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(z16, lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(z84, lmass_model, SF_galaxies=False)

    sp1_sfrmass_local.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')

    print 'SF:  %.3f < z < %.3f' % (numpy.percentile(galaxies_catalog.zspec, 16), numpy.percentile(galaxies_catalog.zspec, 84))



    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(z16, lmass_model, SF_galaxies=True)
    lsfr_upper = log_SFR_tomczak2016(z84, lmass_model, SF_galaxies=True)

    sp2_sfrmass_local.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2_sfrmass_local.axvline(0, color='#cccccc', lw=14, label=zfourge_label, zorder=1)

    print 'ALL: %.3f < z < %.3f' % (numpy.percentile(galaxies_catalog.zspec, 16), numpy.percentile(galaxies_catalog.zspec, 84))





    t3 = sp1_sfrmass_local.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp1_sfrmass_local.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t4 = sp2_sfrmass_local.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp2_sfrmass_local.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_sfrmass_local.legend(loc=4, numpoints=1, fontsize=17, frameon=0)










































