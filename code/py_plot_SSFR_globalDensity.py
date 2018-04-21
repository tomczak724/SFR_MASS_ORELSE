
import mypy
import numpy
import pandas
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
galaxies_catalog = pickle.load(open('../data/table_all_galaxies_noXrayAGN.pickle', 'rb'))









###  reading in Brian's environmental catalog
filedir = '/Users/atomczak/GitHub/ORELSE/Catalogs/Spec_z'
filename = 'allORELSE.localnglobaldensity.wqflag.noappmagcut.noSEDfitreq.w1stblendedserendip_NOHEADER.cat'
df_envData = pandas.read_csv('%s/%s' % (filedir, filename))

ldelta_tomczak_lemaux = []
inds_nomatch = []

for i_gal in range(len(galaxies_catalog.id_gal)):
    i_field = galaxies_catalog.id_field[i_gal]
    f_name = galaxies_catalog.field_names[i_field]
    photID = '%i' % galaxies_catalog.id_gal[i_gal]
    s = df_envData.loc[(df_envData['field'] == f_name) & (df_envData['photID'] == photID)]
    if s.shape[0] == 1:
        ldelta_tomczak_lemaux.append([galaxies_catalog.overdens[i_gal], float(s['log(1+delta_gal)'])])
    elif s.shape[0] > 1:
        ldelta_tomczak_lemaux.append([galaxies_catalog.overdens[i_gal], float(s['log(1+delta_gal)'].iloc[0])])
    else:
        ldelta_tomczak_lemaux.append([galaxies_catalog.overdens[i_gal], numpy.nan])
        inds_nomatch.append(i_gal)

ldelta_tomczak_lemaux = numpy.array(ldelta_tomczak_lemaux)
inds_nomatch = numpy.array(inds_nomatch)









###  binning data by log(eta)

eta_bins = numpy.array([0.1, 0.4, 2.0])
eta_labels = ['          $\eta$ < %.1f' % eta_bins[0],
              '%.1f < $\eta$ < %.1f' % (eta_bins[0], eta_bins[1]),
              '%.1f < $\eta$ < %.1f' % (eta_bins[1], eta_bins[2]),
              '%.1f < $\eta$' % eta_bins[2]]

digi_eta = numpy.digitize(10**galaxies_catalog.log_eta, eta_bins)










#########################################
###  Plotting log(SSFR) vs. log(eta)  ###
#########################################

###  binning in same M* bins as Figure 10 of Muzzin+2012
lmassbins2 = numpy.array([9.3, 10., 10.7])
lmassbins2 = numpy.array([9.4, 10.1, 10.8])
digi_lmass2 = numpy.digitize(galaxies_catalog.lmass, lmassbins2)

lmassbin2_labels = [' %4.1f < log(M$_*$) < %4.1f' % (lmassbins2[0], lmassbins2[1]),
                    '%4.1f < log(M$_*$) < %4.1f' % (lmassbins2[1], lmassbins2[2]),
                    '%4.1f < log(M$_*$)' % lmassbins2[2]]

if False:

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
    ssfr_median_tot = []
    ssfr_nmad_tot = []
    lmass_median_tot = []
    lmass_nmad_tot = []

    n_galaxies_sf = []
    z_galaxies_sf = []
    overdens_median_sf = []
    overdens_nmad_sf = []
    sfr_median_sf = []
    sfr_nmad_sf = []
    ssfr_median_sf = []
    ssfr_nmad_sf = []
    lmass_median_sf = []
    lmass_nmad_sf = []

    n_bootstraps = 1000
    sfr_median_tot_bootstraps = []
    sfr_median_sf_bootstraps = []


    for i_lmass in [1, 2, 3]:

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
        ssfr_median_tot.append([])
        ssfr_nmad_tot.append([])
        lmass_median_tot.append([])
        lmass_nmad_tot.append([])

        n_galaxies_sf.append([])
        z_galaxies_sf.append([])
        overdens_median_sf.append([])
        overdens_nmad_sf.append([])
        sfr_median_sf.append([])
        sfr_nmad_sf.append([])
        ssfr_median_sf.append([])
        ssfr_nmad_sf.append([])
        lmass_median_sf.append([])
        lmass_nmad_sf.append([])

        sfr_median_tot_bootstraps.append([])
        sfr_median_sf_bootstraps.append([])

        for i_overdens in range(len(eta_bins)+1)[::-1]:

            inds_lmass_overdens_tot = numpy.where((digi_lmass2 == i_lmass) & \
                                                  (digi_eta == i_overdens))[0]

            inds_lmass_overdens_sf = numpy.where((galaxies_catalog.UVJ_class == 1) & \
                                                 (digi_lmass2 == i_lmass) & \
                                                 (digi_eta == i_overdens))[0]


            ngals_tot = len(inds_lmass_overdens_tot)
            ngals_sf = len(inds_lmass_overdens_sf)

            n_galaxies_tot[-1].append(ngals_tot)
            z_galaxies_tot[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist())
            sfr_median_tot[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            ssfr_median_tot[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot] / 10**galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            ssfr_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot] / 10**galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            overdens_median_tot[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_tot]))
            overdens_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_tot]))
            lmass_median_tot[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf[-1].append(ngals_sf)
            z_galaxies_sf[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
            sfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            ssfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf] / 10**galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            ssfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf] / 10**galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            overdens_median_sf[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
            overdens_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
            lmass_median_sf[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))



            ###  running Bootstrap resamplings
            sfr_median_tot_bootstraps[-1].append([])
            sfr_median_sf_bootstraps[-1].append([])
        
            for i_bootstrap in range(n_bootstraps):

                inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_tot), len(inds_lmass_overdens_tot))
                inds_lmass_overdens_tot_bootstrap = inds_lmass_overdens_tot[inds_resamp]
                sfr_median_tot_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot_bootstrap]))

                inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_sf), len(inds_lmass_overdens_sf))
                inds_lmass_overdens_sf_bootstrap = inds_lmass_overdens_sf[inds_resamp]
                sfr_median_sf_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf_bootstrap]))





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

    for i_lmass in range(len(z_galaxies_tot)):

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
    ssfr_median_tot = numpy.array(ssfr_median_tot)
    ssfr_nmad_tot = numpy.array(ssfr_nmad_tot)
    overdens_median_tot = numpy.array(overdens_median_tot)
    overdens_nmad_tot = numpy.array(overdens_nmad_tot)
    n_galaxies_tot = numpy.array(n_galaxies_tot)

    lmass_median_sf = numpy.array(lmass_median_sf)
    lmass_nmad_sf = numpy.array(lmass_nmad_sf)
    sfr_median_sf = numpy.array(sfr_median_sf)
    sfr_nmad_sf = numpy.array(sfr_nmad_sf)
    ssfr_median_sf = numpy.array(ssfr_median_sf)
    ssfr_nmad_sf = numpy.array(ssfr_nmad_sf)
    overdens_median_sf = numpy.array(overdens_median_sf)
    overdens_nmad_sf = numpy.array(overdens_nmad_sf)
    n_galaxies_sf = numpy.array(n_galaxies_sf)


    ###  applying SFR offsets
    lsfr_median_tot = numpy.log10(sfr_median_tot) + lsfr_offsets_tot
    lsfr_median_sf = numpy.log10(sfr_median_sf) + lsfr_offsets_sf

    lssfr_median_tot = numpy.log10(ssfr_median_tot) + lsfr_offsets_tot
    lssfr_median_sf = numpy.log10(ssfr_median_sf) + lsfr_offsets_sf





    ###  Bootstrap scatters
    sfr_median_tot_bootstraps = numpy.array(sfr_median_tot_bootstraps)
    sfr_median_sf_bootstraps = numpy.array(sfr_median_sf_bootstraps)

    sfr_median_tot_bootstraps[sfr_median_tot_bootstraps <= 0] = 10**-5
    sfr_median_sf_bootstraps[sfr_median_sf_bootstraps <= 0] = 10**-5

    sfr_p16_tot_bootstrap = numpy.percentile(sfr_median_tot_bootstraps, 16, axis=2)
    sfr_p50_tot_bootstrap = numpy.percentile(sfr_median_tot_bootstraps, 50, axis=2)
    sfr_p84_tot_bootstrap = numpy.percentile(sfr_median_tot_bootstraps, 84, axis=2)

    sfr_p16_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 16, axis=2)
    sfr_p50_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 50, axis=2)
    sfr_p84_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 84, axis=2)

    lsfr_p16_tot_bootstrap = numpy.log10(sfr_p16_tot_bootstrap)
    lsfr_p50_tot_bootstrap = numpy.log10(sfr_p50_tot_bootstrap)
    lsfr_p84_tot_bootstrap = numpy.log10(sfr_p84_tot_bootstrap)

    lsfr_p16_sf_bootstrap = numpy.log10(sfr_p16_sf_bootstrap)
    lsfr_p50_sf_bootstrap = numpy.log10(sfr_p50_sf_bootstrap)
    lsfr_p84_sf_bootstrap = numpy.log10(sfr_p84_sf_bootstrap)






    fig_ssfrDens = pyplot.figure(figsize=(7.36, 7.17))

    sp2_ssfrDens = fig_ssfrDens.add_subplot(111)

    sp2_ssfrDens.grid()
    sp2_ssfrDens.minorticks_on()
    sp2_ssfrDens.set_xlabel('log( $\eta$ )')
    sp2_ssfrDens.set_ylabel('log( SSFR$_{UV+IR}$ / [yr$^{-1}$] )', size=22)

    sp2_ssfrDens.axis([0., 1.43, -10.12, -8.58])
    sp2_ssfrDens.axis([-1.9, 2.15, -10.22, -8.63])


    fig_ssfrDens.subplots_adjust(wspace=0, left=0.09, right=0.98)
    fig_ssfrDens.subplots_adjust(wspace=0, left=0.17, right=0.98)


    colors = ['b', 'yellow', 'orange']


    for i_lmass in range(len(z_galaxies_tot)):

        ###  error on the median
        elog_sfr_tot = numpy.log10((sfr_median_tot[i_lmass] + sfr_nmad_tot[i_lmass]/n_galaxies_tot[i_lmass]**0.5) / sfr_median_tot[i_lmass])
        elog_sfr_sf = numpy.log10((sfr_median_sf[i_lmass] + sfr_nmad_sf[i_lmass]/n_galaxies_sf[i_lmass]**0.5) / sfr_median_sf[i_lmass])

        ###  bootstrap error
        #elog_sfr_tot = lsfr_scatter_tot_bootstrap[i_lmass]
        #elog_sfr_sf = lsfr_scatter_sf_bootstrap[i_lmass]

        ###  bootstrap error
        elo_lsfr_tot = (lsfr_p50_tot_bootstrap - lsfr_p16_tot_bootstrap)[i_lmass]
        ehi_lsfr_tot = (lsfr_p84_tot_bootstrap - lsfr_p50_tot_bootstrap)[i_lmass]
        elog_sfr_tot = [elo_lsfr_tot, ehi_lsfr_tot]

        elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[i_lmass]
        ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[i_lmass]
        elog_sfr_sf = [elo_lsfr_sf, ehi_lsfr_sf]


        ###  making bootstrap errors symmetric by
        ###  taking average of lower, upper errors
        elog_sfr_tot = (elo_lsfr_tot + ehi_lsfr_tot) / 2.
        elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.



        color_overdens = colors[i_lmass]

        sp2_ssfrDens.plot(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                          ls='-', lw=3, color='k', zorder=1)



        sp2_ssfrDens.errorbar(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                     yerr=elog_sfr_sf, xerr=0,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor='k', 
                     label=lmassbin2_labels[i_lmass])


    t2 = sp2_ssfrDens.text(0.03, 0.97, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp2_ssfrDens.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_ssfrDens.legend(loc=4, numpoints=1, fontsize=18)





















############################################
###  Plotting log(SSFR) vs. log(1+dgal)  ###
############################################

zlo, zhi = 0.6, 1.3
#zlo, zhi = 0.6, 0.9  # low-z
#zlo, zhi = 0.9, 1.3  # high-z

###  binning in same M* bins as Figure 10 of Muzzin+2012
lmassbins2 = numpy.array([9.4, 10.1, 10.8])
digi_lmass2 = numpy.digitize(galaxies_catalog.lmass, lmassbins2)

lmassbin2_labels = [' %4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins2[0], lmassbins2[1]),
                    '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins2[1], lmassbins2[2]),
                    '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$)' % lmassbins2[2]]


overdensbins = galaxies_catalog.overdensbins
overdensbars = galaxies_catalog.overdensbars

###  Tomczak's log(1+dgal)
overdens = galaxies_catalog.overdens
digi_overdens = galaxies_catalog.digi_overdens

###  Lemaux's log(1+dgal)
#overdens = ldelta_tomczak_lemaux[:,1]
#digi_overdens = numpy.digitize(overdens, overdensbins)


if True:

    pdf_zdists = PdfPages('../figures/redshift_distributions_weighted.pdf')



    fig_zdists = pyplot.figure(figsize=(16.5, 11.))

    sp2_zdists = fig_zdists.add_subplot(222)
    sp1_zdists = fig_zdists.add_subplot(221)
    sp4_zdists = fig_zdists.add_subplot(224)
    sp3_zdists = fig_zdists.add_subplot(223)


    fig_zdists.subplots_adjust(left=0.07, right=0.98, hspace=0, wspace=0)





    cmap = pyplot.cm.rainbow
    colors = ['b', 'yellow', 'orange']



    n_galaxies_tot = []
    z_galaxies_tot = []
    overdens_median_tot = []
    overdens_nmad_tot = []
    sfr_median_tot = []
    sfr_nmad_tot = []
    ssfr_median_tot = []
    ssfr_nmad_tot = []
    lmass_median_tot = []
    lmass_nmad_tot = []

    n_galaxies_sf = []
    z_galaxies_sf = []
    overdens_median_sf = []
    overdens_nmad_sf = []
    sfr_median_sf = []
    sfr_nmad_sf = []
    ssfr_median_sf = []
    ssfr_nmad_sf = []
    lmass_median_sf = []
    lmass_nmad_sf = []

    n_bootstraps = 5000
    sfr_median_tot_bootstraps = []
    sfr_median_sf_bootstraps = []


    for i_lmass in range(1, len(lmassbins2)+1):

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
        ssfr_median_tot.append([])
        ssfr_nmad_tot.append([])
        lmass_median_tot.append([])
        lmass_nmad_tot.append([])

        n_galaxies_sf.append([])
        z_galaxies_sf.append([])
        overdens_median_sf.append([])
        overdens_nmad_sf.append([])
        sfr_median_sf.append([])
        sfr_nmad_sf.append([])
        ssfr_median_sf.append([])
        ssfr_nmad_sf.append([])
        lmass_median_sf.append([])
        lmass_nmad_sf.append([])

        sfr_median_tot_bootstraps.append([])
        sfr_median_sf_bootstraps.append([])

        for i_overdens in range(1, len(overdensbins)):

            inds_lmass_overdens_tot = numpy.where((galaxies_catalog.zspec >= zlo) & \
                                                  (galaxies_catalog.zspec <= zhi) & \
                                                  (digi_lmass2 == i_lmass) & \
                                                  (digi_overdens == i_overdens))[0]

            inds_lmass_overdens_sf = numpy.where((galaxies_catalog.zspec >= zlo) & \
                                                 (galaxies_catalog.zspec <= zhi) & \
                                                 (galaxies_catalog.UVJ_class == 1) & \
                                                 (digi_lmass2 == i_lmass) & \
                                                 (digi_overdens == i_overdens))[0]


            ngals_tot = len(inds_lmass_overdens_tot)
            ngals_sf = len(inds_lmass_overdens_sf)

            n_galaxies_tot[-1].append(ngals_tot)
            z_galaxies_tot[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist())
            sfr_median_tot[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            ssfr_median_tot[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot] / 10**galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            ssfr_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot] / 10**galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            overdens_median_tot[-1].append(numpy.median(overdens[inds_lmass_overdens_tot]))
            overdens_nmad_tot[-1].append(mypy.nmad(overdens[inds_lmass_overdens_tot]))
            lmass_median_tot[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf[-1].append(ngals_sf)
            z_galaxies_sf[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
            sfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            ssfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf] / 10**galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            ssfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf] / 10**galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            overdens_median_sf[-1].append(numpy.median(overdens[inds_lmass_overdens_sf]))
            overdens_nmad_sf[-1].append(mypy.nmad(overdens[inds_lmass_overdens_sf]))
            lmass_median_sf[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))



            ###  running Bootstrap resamplings
            sfr_median_tot_bootstraps[-1].append([])
            sfr_median_sf_bootstraps[-1].append([])
        
            for i_bootstrap in range(n_bootstraps):

                inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_tot), len(inds_lmass_overdens_tot))
                inds_lmass_overdens_tot_bootstrap = inds_lmass_overdens_tot[inds_resamp]
                sfr_median_tot_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot_bootstrap]))

                inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_sf), len(inds_lmass_overdens_sf))
                inds_lmass_overdens_sf_bootstrap = inds_lmass_overdens_sf[inds_resamp]
                sfr_median_sf_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf_bootstrap]))





            ###  plotting redshift distributions
            offsets = [-0.0013, 0., 0.0013]
            color_overdens = colors[i_overdens-1]



            h_tot = sp1_zdists.hist(z_galaxies_tot[-1][-1], histtype='step', lw=3, color=color_overdens,
                             bins=15, range=(0.7+offsets[i_overdens-1], 1.3+offsets[i_overdens-1]), 
                             label=str(len(z_galaxies_tot[-1][-1])))

            h_sf = sp2_zdists.hist(z_galaxies_sf[-1][-1], histtype='step', lw=3, color=color_overdens,
                            bins=15, range=(0.7+offsets[i_overdens-1], 1.3+offsets[i_overdens-1]), 
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

    z_1d_tot = []
    for zarr1 in z_galaxies_tot:
        for zarr2 in zarr1:
            z_1d_tot += zarr2

    z_1d_sf = []
    for zarr1 in z_galaxies_sf:
        for zarr2 in zarr1:
            z_1d_sf += zarr2

    z_fiducial = numpy.median(z_1d_tot)

    lsfr_offsets_tot = []
    lsfr_offsets_sf = []
    z_median_tot = []
    z_median_sf = []

    for i_lmass in range(len(z_galaxies_tot)):

        lsfr_offsets_tot.append([])
        lsfr_offsets_sf.append([])
        z_median_tot.append([])
        z_median_sf.append([])

        for i_overdens in range(len(overdensbars)):

            z_medi_tot = numpy.median(z_galaxies_tot[i_lmass][i_overdens])
            z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])

            lsfr_fiducial_tot = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_z_medi_tot = log_SFR_tomczak2016(z_medi_tot, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            z_median_tot[-1].append(z_medi_tot)
            z_median_sf[-1].append(z_medi_sf)

            lsfr_offsets_tot[-1].append(lsfr_fiducial_tot - lsfr_z_medi_tot)
            lsfr_offsets_sf[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)

    lsfr_offsets_tot = numpy.array(lsfr_offsets_tot)
    lsfr_offsets_sf = numpy.array(lsfr_offsets_sf)
    z_median_tot = numpy.array(z_median_tot)
    z_median_sf = numpy.array(z_median_sf)








    ###  plotting SFR vs. lmass
    lmass_median_tot = numpy.array(lmass_median_tot)
    lmass_nmad_tot = numpy.array(lmass_nmad_tot)
    sfr_median_tot = numpy.array(sfr_median_tot)
    sfr_nmad_tot = numpy.array(sfr_nmad_tot)
    ssfr_median_tot = numpy.array(ssfr_median_tot)
    ssfr_nmad_tot = numpy.array(ssfr_nmad_tot)
    overdens_median_tot = numpy.array(overdens_median_tot)
    overdens_nmad_tot = numpy.array(overdens_nmad_tot)
    n_galaxies_tot = numpy.array(n_galaxies_tot)

    lmass_median_sf = numpy.array(lmass_median_sf)
    lmass_nmad_sf = numpy.array(lmass_nmad_sf)
    sfr_median_sf = numpy.array(sfr_median_sf)
    sfr_nmad_sf = numpy.array(sfr_nmad_sf)
    ssfr_median_sf = numpy.array(ssfr_median_sf)
    ssfr_nmad_sf = numpy.array(ssfr_nmad_sf)
    overdens_median_sf = numpy.array(overdens_median_sf)
    overdens_nmad_sf = numpy.array(overdens_nmad_sf)
    n_galaxies_sf = numpy.array(n_galaxies_sf)


    ###  applying SFR offsets
    lsfr_median_tot = numpy.log10(sfr_median_tot) + lsfr_offsets_tot
    lsfr_median_sf = numpy.log10(sfr_median_sf) + lsfr_offsets_sf

    lssfr_median_tot = numpy.log10(ssfr_median_tot) + lsfr_offsets_tot
    lssfr_median_sf = numpy.log10(ssfr_median_sf) + lsfr_offsets_sf



    print '\nNgal for SF galaxies:'
    print n_galaxies_sf
    print ''



    ###  Bootstrap scatters
    sfr_median_tot_bootstraps = numpy.array(sfr_median_tot_bootstraps)
    sfr_median_sf_bootstraps = numpy.array(sfr_median_sf_bootstraps)

    sfr_median_tot_bootstraps[sfr_median_tot_bootstraps <= 0] = 10**-5
    sfr_median_sf_bootstraps[sfr_median_sf_bootstraps <= 0] = 10**-5

    sfr_p16_tot_bootstrap = numpy.percentile(sfr_median_tot_bootstraps, 16, axis=2)
    sfr_p50_tot_bootstrap = numpy.percentile(sfr_median_tot_bootstraps, 50, axis=2)
    sfr_p84_tot_bootstrap = numpy.percentile(sfr_median_tot_bootstraps, 84, axis=2)

    sfr_p16_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 16, axis=2)
    sfr_p50_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 50, axis=2)
    sfr_p84_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 84, axis=2)

    lsfr_p16_tot_bootstrap = numpy.log10(sfr_p16_tot_bootstrap)
    lsfr_p50_tot_bootstrap = numpy.log10(sfr_p50_tot_bootstrap)
    lsfr_p84_tot_bootstrap = numpy.log10(sfr_p84_tot_bootstrap)

    lsfr_p16_sf_bootstrap = numpy.log10(sfr_p16_sf_bootstrap)
    lsfr_p50_sf_bootstrap = numpy.log10(sfr_p50_sf_bootstrap)
    lsfr_p84_sf_bootstrap = numpy.log10(sfr_p84_sf_bootstrap)

    lsfr_median_tot_err = ((lsfr_p50_tot_bootstrap - lsfr_p16_tot_bootstrap) + \
                           (lsfr_p84_tot_bootstrap - lsfr_p50_tot_bootstrap)) / 2

    lsfr_median_sf_err = ((lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap) + \
                          (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)) / 2





    fig_ssfrDens = pyplot.figure(figsize=(7.36, 7.17))

    sp2_ssfrDens = fig_ssfrDens.add_subplot(111)

    sp2_ssfrDens.grid(color='#666666')
    sp2_ssfrDens.minorticks_on()
    sp2_ssfrDens.set_xlabel('log( 1 + $\delta_{gal}$ )', size=22)
    sp2_ssfrDens.set_ylabel('log( SSFR$_{UV+IR}$ / [yr$^{-1}$] )', size=22)

    sp2_ssfrDens.axis([0., 1.43, -10.12, -8.58])


    fig_ssfrDens.subplots_adjust(wspace=0, left=0.17, right=0.98)
    fig_ssfrDens.subplots_adjust(wspace=0, left=0.17, right=0.98, top=0.94)


    colors = ['b', 'yellow', 'orange']


    #for i_lmass in range(len(z_galaxies_tot))[::-1]:
    for i_lmass in range(len(z_galaxies_tot)):

        ###  error on the median
        elog_sfr_tot = numpy.log10((sfr_median_tot[i_lmass] + sfr_nmad_tot[i_lmass]/n_galaxies_tot[i_lmass]**0.5) / sfr_median_tot[i_lmass])
        elog_sfr_sf = numpy.log10((sfr_median_sf[i_lmass] + sfr_nmad_sf[i_lmass]/n_galaxies_sf[i_lmass]**0.5) / sfr_median_sf[i_lmass])


        ###  bootstrap error
        elo_lsfr_tot = (lsfr_p50_tot_bootstrap - lsfr_p16_tot_bootstrap)[i_lmass]
        ehi_lsfr_tot = (lsfr_p84_tot_bootstrap - lsfr_p50_tot_bootstrap)[i_lmass]
        elog_sfr_tot = [elo_lsfr_tot, ehi_lsfr_tot]

        elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[i_lmass]
        ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[i_lmass]
        elog_sfr_sf = [elo_lsfr_sf, ehi_lsfr_sf]



        ###  making bootstrap errors symmetric by
        ###  taking average of lower, upper errors
        elog_sfr_tot = (elo_lsfr_tot + ehi_lsfr_tot) / 2.
        elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.




        color_overdens = colors[i_lmass]

        sp2_ssfrDens.plot(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                          ls='-', lw=3.5, color='k', zorder=1)

        sp2_ssfrDens.errorbar(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                     yerr=elog_sfr_sf, xerr=overdens_nmad_sf[i_lmass] * 0,
                     ls='-', lw=1.5, marker='o', mew=1.8, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor='k', 
                     label=lmassbin2_labels[i_lmass])


    t2 = sp2_ssfrDens.text(0.03, 0.97, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp2_ssfrDens.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_ssfrDens.legend(loc=3, numpoints=1, fontsize=16)


    sp2_ssfrDens.set_title('Optically-selected galaxies at 0.6 < z$_{spec}$ < 1.3', verticalalignment='bottom')
    sp2_ssfrDens.set_title('UV/Optical-selected galaxies in all ORELSE fields', verticalalignment='bottom')



    ###  printing LaTeX table
    for i_lmass in range(len(z_galaxies_tot)):

        line = '%-40s' % lmassbin2_labels[i_lmass].strip().replace('<', '$<$')

        for i_overdens in range(len(lssfr_median_sf[i_lmass])):

            elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[i_lmass][i_overdens]
            ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[i_lmass][i_overdens]
            elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.

            line += ' &'
            line += ' $%5.2f' % lssfr_median_sf[i_lmass][i_overdens]
            line += '\\pm%4.2f$' % elog_sfr_sf

        line += ' \\\\'
        print line







    ###  saving as ascii tables
    #numpy.savetxt('../data/table_SSFR_density_Mstar_lowz___overdens_median_tot.dat', overdens_median_tot)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_lowz___lssfr_median_tot.dat', lsfr_median_tot)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_lowz___lsfr_median_tot_err.dat', lsfr_median_tot_err)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_lowz___overdens_median_sf.dat', overdens_median_sf)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_lowz___lssfr_median_sf.dat', lssfr_median_sf)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_lowz___lsfr_median_sf_err.dat', lsfr_median_sf_err)

    #numpy.savetxt('../data/table_SSFR_density_Mstar_highz___overdens_median_tot.dat', overdens_median_tot)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_highz___lssfr_median_tot.dat', lsfr_median_tot)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_highz___lsfr_median_tot_err.dat', lsfr_median_tot_err)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_highz___overdens_median_sf.dat', overdens_median_sf)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_highz___lssfr_median_sf.dat', lssfr_median_sf)
    #numpy.savetxt('../data/table_SSFR_density_Mstar_highz___lsfr_median_sf_err.dat', lsfr_median_sf_err)























#####################################################
###  printing LaTeX table for low- and high-z bins
#####################################################

lmassbins2 = numpy.array([9.4, 10.1, 10.8])
lmassbin3_labels = ['%.1f$-$%.1f' % (lmassbins2[0], lmassbins2[1]),
                    '%.1f$-$%.1f' % (lmassbins2[1], lmassbins2[2]),
                    '$>$%.1f' % lmassbins2[2]]

if True:

    lssfr_loz = numpy.loadtxt("../data/table_SSFR_density_Mstar_lowz___lssfr_median_sf.dat")
    lssfr_err_loz = numpy.loadtxt("../data/table_SSFR_density_Mstar_lowz___lsfr_median_sf_err.dat")

    lssfr_hiz = numpy.loadtxt("../data/table_SSFR_density_Mstar_highz___lssfr_median_sf.dat")
    lssfr_err_hiz = numpy.loadtxt("../data/table_SSFR_density_Mstar_highz___lsfr_median_sf_err.dat")


    ###  printing LaTeX table
    for i_lmass in range(len(lssfr_hiz)):

        line = '    %-11s' % lmassbin3_labels[i_lmass].strip().replace('<', '$<$')

        for i_overdens in range(len(lssfr_loz[i_lmass])):

            line += ' &'
            line += ' $%5.2f' % lssfr_loz[i_lmass][i_overdens]
            line += '\\pm%4.2f$' % lssfr_err_loz[i_lmass][i_overdens]

        line += ' \\\\'
        print line


    print('\n')


    ###  printing LaTeX table
    for i_lmass in range(len(lssfr_hiz)):

        line = '    %-11s' % lmassbin3_labels[i_lmass].strip().replace('<', '$<$')

        for i_overdens in range(len(lssfr_hiz[i_lmass])):

            line += ' &'
            line += ' $%5.2f' % lssfr_hiz[i_lmass][i_overdens]
            line += '\\pm%4.2f$' % lssfr_err_hiz[i_lmass][i_overdens]

        line += ' \\\\'
        print line


























###  plot delta_SFR of SF vs. ALL galaxies as function of density at fixed M*

if False:

    fig_deltaSFR = pyplot.figure(figsize=(9., 7.))

    sp_deltaSFR = fig_deltaSFR.add_subplot(111)

    sp_deltaSFR.minorticks_on()
    sp_deltaSFR.set_xlabel('log( $\eta$ )')
    sp_deltaSFR.set_ylabel('$\Delta$ log( SFR )     [SF - ALL]')

    sp_deltaSFR.axis([-1.7, 2.2, -0.1, 0.8])


    fig_deltaSFR.subplots_adjust(wspace=0, left=0.12, right=0.98)


    colors = ['b', 'yellow', 'orange']


    for i_lmass in range(len(z_galaxies_tot))[::-1]:

        elog_sfr_tot = numpy.log10((sfr_median_tot[i_lmass] + sfr_nmad_tot[i_lmass]/n_galaxies_tot[i_lmass]**0.5) / sfr_median_tot[i_lmass])
        elog_sfr_sf = numpy.log10((sfr_median_sf[i_lmass] + sfr_nmad_sf[i_lmass]/n_galaxies_sf[i_lmass]**0.5) / sfr_median_sf[i_lmass])

        color_overdens = colors[i_lmass]


        dSFR = lsfr_median_sf[i_lmass] - lsfr_median_tot[i_lmass]


        sp_deltaSFR.plot(overdens_median_sf[i_lmass], dSFR, 
                          ls='-', lw=3, color='k', zorder=1)


        sp_deltaSFR.errorbar(overdens_median_sf[i_lmass], dSFR, 
                     yerr=elog_sfr_sf, xerr=0,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor='k', 
                     label=lmassbin2_labels[i_lmass])

    leg = sp_deltaSFR.legend(loc=4, numpoints=1)















###  plot delta_SFR of SF vs. ALL galaxies as function of M* at fixed density

if False:

    fig_deltaSFR = pyplot.figure(figsize=(9., 7.))

    sp_deltaSFR = fig_deltaSFR.add_subplot(111)

    sp_deltaSFR.minorticks_on()
    sp_deltaSFR.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp_deltaSFR.set_ylabel('$\Delta$ log( SFR )     [SF - ALL]')

    sp_deltaSFR.axis([9.5, 11.1, -0.1, 0.7])


    fig_deltaSFR.subplots_adjust(wspace=0, left=0.12, right=0.98)


    colors = ['r', 'orange', 'g', 'b']


    for i_overdens in range(len(sfr_median_tot[0])):

        elog_sfr_tot = numpy.log10((sfr_median_tot[:, i_overdens] + sfr_nmad_tot[:, i_overdens]/n_galaxies_tot[:, i_overdens]**0.5) / sfr_median_tot[:, i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf[:, i_overdens] + sfr_nmad_sf[:, i_overdens]/n_galaxies_sf[:, i_overdens]**0.5) / sfr_median_sf[:, i_overdens])

        color_overdens = colors[i_overdens]


        dSFR = lsfr_median_sf[:, i_overdens] - lsfr_median_tot[:, i_overdens]


        sp_deltaSFR.plot(lmass_median_tot[:, i_overdens], dSFR, 
                          ls='-', lw=3, color='k', zorder=1)


        sp_deltaSFR.errorbar(lmass_median_tot[:, i_overdens], dSFR, 
                     yerr=elog_sfr_sf, xerr=0,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor='k', 
                     label=eta_labels[i_overdens])

    leg = sp_deltaSFR.legend(loc=4, numpoints=1)























