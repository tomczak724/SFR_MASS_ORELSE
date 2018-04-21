
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



###  Testing three UVJ selections:
###    1) The "normal" fiducial UVJ selection
###    2) More conserveative (i.e. smaller QU box)
###    3) Mre liberal (i.e. larger QU box)


def uvj_selection(uv, vj, delta_y=0.00):
    '''
    Star-forming -> 1
    Quiescent ----> 0

    Identifies galaxies as quiescent if both:

    (U-V) > slope1 * (V-J) + intercept1
    (U-V) > slope2 * (V-J) + intercept2

    The "delta_y" is to allow for testing of the boundary
    '''


    slope1, intercept1 = 0.00, 1.30
    slope2, intercept2 = 0.88, 0.59

    intercept1 += delta_y
    intercept2 += delta_y

    uvj_flag = numpy.ones(len(uv))

    inds = numpy.where((uv > slope1*vj + intercept1) & \
                       (uv > slope2*vj + intercept2))[0]

    uvj_flag[inds] = 0
    return uvj_flag



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



###  binning data by LOCAL overdensity
lmassbins2 = numpy.array([9.4, 10.1, 10.8])
digi_lmass2 = numpy.digitize(galaxies_catalog.lmass, lmassbins2)

lmassbin2_labels = [' %4.1f < log(M$_*$) < %4.1f' % (lmassbins2[0], lmassbins2[1]),
                    '%4.1f < log(M$_*$) < %4.1f' % (lmassbins2[1], lmassbins2[2]),
                    '%4.1f < log(M$_*$)' % lmassbins2[2]]


###  binning data by GLOBAL overdensity
eta_bins = numpy.array([0.1, 0.4, 2.0])
eta_labels = ['          $\eta$ < %.1f' % eta_bins[0],
              '%.1f < $\eta$ < %.1f' % (eta_bins[0], eta_bins[1]),
              '%.1f < $\eta$ < %.1f' % (eta_bins[1], eta_bins[2]),
              '%.1f < $\eta$' % eta_bins[2]]

digi_eta = numpy.digitize(10**galaxies_catalog.log_eta, eta_bins)






























dy = 0.15
vertices_conservative = [(-9., 1.3+dy), (0.807, 1.3+dy), (9., 8.51+dy)]
vertices_fiducial     = [(-9., 1.3), (0.807, 1.3), (9., 8.51)]
vertices_liberal      = [(-9., 1.3-dy), (0.807, 1.3-dy), (9., 8.51-dy)]
vertices_all = [vertices_conservative, vertices_fiducial, vertices_liberal]

delta_ys = [dy, 0., -dy]














###  setting up figure

fig = pyplot.figure(figsize=(21., 7.4))

sp1 = fig.add_subplot(131)
sp2 = fig.add_subplot(132)
sp3 = fig.add_subplot(133)



def fig_setup_LOCAL():

    sp1.grid()
    sp2.grid()
    sp3.grid()

    sp1.minorticks_on()
    sp2.minorticks_on()
    sp3.minorticks_on()

    sp1.set_xlabel('( V - J )$_{rest}$')
    sp1.set_ylabel('( U - V )$_{rest}$')
    sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp2.set_ylabel('log( SFR / [M$_{\odot}$ / yr] )')
    sp3.set_xlabel('log( 1 + $\delta_{gal}$ )')
    sp3.set_ylabel('log( SSFR / [yr$^{-1}$] )')

    fig.subplots_adjust(top=0.99, right=0.99, left=0.05)


    sp1.axis([-0.5, 2.2, -0.2, 2.5])
    sp2.axis([8.2, 12.3, -0.9, 2.0])
    sp3.axis([0., 1.45, -10.2, -8.6])



    for v in vertices_all:
        sp1.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='k', ls='--', lw=2.5, dashes=[4, 2], zorder=2)
        sp1.plot([v[1][0], v[2][0]], [v[1][1], v[2][1]], color='k', ls='--', lw=2.5, dashes=[4, 2], zorder=2)


    sp1.plot(galaxies_catalog.VJ_color, galaxies_catalog.UV_color, 
             ls='', marker='o', mfc='gray', mec='gray', ms=1, zorder=1)


def fig_setup_GLOBAL():

    sp1.grid()
    sp2.grid()
    sp3.grid()

    sp1.minorticks_on()
    sp2.minorticks_on()
    sp3.minorticks_on()

    sp1.set_xlabel('( V - J )$_{rest}$')
    sp1.set_ylabel('( U - V )$_{rest}$')
    sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp2.set_ylabel('log( SFR / [M$_{\odot}$ / yr] )')
    sp3.set_xlabel('log( $\eta$ )')
    sp3.set_ylabel('log( SSFR / [yr$^{-1}$] )')

    fig.subplots_adjust(top=0.99, right=0.99, left=0.05)


    sp1.axis([-0.5, 2.2, -0.2, 2.5])
    sp2.axis([8.2, 12.3, -0.9, 2.0])
    sp3.axis([-1.9, 2.1, -10.4, -8.7])



    for v in vertices_all:
        sp1.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='k', ls='--', lw=2.5, dashes=[4, 2], zorder=2)
        sp1.plot([v[1][0], v[2][0]], [v[1][1], v[2][1]], color='k', ls='--', lw=2.5, dashes=[4, 2], zorder=2)


    sp1.plot(galaxies_catalog.VJ_color, galaxies_catalog.UV_color, 
             ls='', marker='o', mfc='gray', mec='gray', ms=1, zorder=1)


















#####################################
###  For LOCAL evironment metric  ###
#####################################

pdf = PdfPages('../figures/uvj_boundary_variation_LOCAL.pdf')

for i_vertices in range(len(vertices_all)):

    fig_setup_LOCAL()

    v = vertices_all[i_vertices]

    sp1.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='#c40000', lw=4, zorder=3)
    sp1.plot([v[1][0], v[2][0]], [v[1][1], v[2][1]], color='#c40000', lw=4, zorder=3)



    ###  classifying SF vs. QU galaxies

    uvj_class = uvj_selection(galaxies_catalog.UV_color, 
                              galaxies_catalog.VJ_color,
                              delta_y=delta_ys[i_vertices])


    n_bootstraps = 1000












    ##############################
    ###  Plotting SFR-density  ###
    ##############################

    if True:

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

        sfr_median_sf_bootstraps = []


        for i_lmass in range(1, len(lmassbins2)+1):

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

            sfr_median_sf_bootstraps.append([])

            for i_overdens in range(1, len(galaxies_catalog.overdensbins)):

                inds_lmass_overdens_sf = numpy.where((uvj_class == 1) & \
                                                     (digi_lmass2 == i_lmass) & \
                                                     (galaxies_catalog.digi_overdens == i_overdens))[0]


                ngals_sf = len(inds_lmass_overdens_sf)

                n_galaxies_sf[-1].append(ngals_sf)
                z_galaxies_sf[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
                sfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
                sfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
                ssfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf] / 10**galaxies_catalog.lmass[inds_lmass_overdens_sf]))
                ssfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf] / 10**galaxies_catalog.lmass[inds_lmass_overdens_sf]))
                overdens_median_sf[-1].append(numpy.median(galaxies_catalog.overdens[inds_lmass_overdens_sf]))
                overdens_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.overdens[inds_lmass_overdens_sf]))
                lmass_median_sf[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
                lmass_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))



                ###  running Bootstrap resamplings
                sfr_median_sf_bootstraps[-1].append([])
            
                for i_bootstrap in range(n_bootstraps):

                    inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_sf), len(inds_lmass_overdens_sf))
                    inds_lmass_overdens_sf_bootstrap = inds_lmass_overdens_sf[inds_resamp]
                    sfr_median_sf_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf_bootstrap]))






        ###  estimating SFR offsets from median to fiducial redshifts
        z_fiducial = numpy.median(galaxies_catalog.zspec)

        lsfr_offsets_sf = []

        for i_lmass in range(len(z_galaxies_sf)):

            lsfr_offsets_sf.append([])

            for i_overdens in range(len(galaxies_catalog.overdensbars)):

                z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])
                lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_offsets_sf[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)

        lsfr_offsets_sf = numpy.array(lsfr_offsets_sf)



        ###  plotting SSFR vs. density

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
        lssfr_median_sf = numpy.log10(ssfr_median_sf) + lsfr_offsets_sf





        ###  Bootstrap scatters
        sfr_median_sf_bootstraps = numpy.array(sfr_median_sf_bootstraps)

        sfr_median_sf_bootstraps[sfr_median_sf_bootstraps <= 0] = 10**-5

        sfr_p16_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 16, axis=2)
        sfr_p50_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 50, axis=2)
        sfr_p84_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 84, axis=2)

        lsfr_p16_sf_bootstrap = numpy.log10(sfr_p16_sf_bootstrap)
        lsfr_p50_sf_bootstrap = numpy.log10(sfr_p50_sf_bootstrap)
        lsfr_p84_sf_bootstrap = numpy.log10(sfr_p84_sf_bootstrap)





        colors = ['b', 'yellow', 'orange']


        for i_lmass in range(len(z_galaxies_sf)):

            ###  error on the median
            elog_sfr_sf = numpy.log10((sfr_median_sf[i_lmass] + sfr_nmad_sf[i_lmass]/n_galaxies_sf[i_lmass]**0.5) / sfr_median_sf[i_lmass])


            ###  bootstrap error
            elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[i_lmass]
            ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[i_lmass]
            elog_sfr_sf = [elo_lsfr_sf, ehi_lsfr_sf]



            ###  making bootstrap errors symmetric by
            ###  taking average of lower, upper errors
            elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.



            color_overdens = colors[i_lmass]

            sp3.plot(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                              ls='-', lw=3.5, color='k', zorder=1)

            sp3.errorbar(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                         yerr=elog_sfr_sf, xerr=overdens_nmad_sf[i_lmass] * 0,
                         ls='-', lw=1.5, marker='o', mew=1.8, ms=9, elinewidth=2, zorder=2,
                         mfc=color_overdens, color=color_overdens, mec='k', ecolor='k', 
                         label=lmassbin2_labels[i_lmass])


        t2 = sp3.text(0.03, 0.97, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                      transform=sp3.transAxes, horizontalalignment='left', verticalalignment='top',
                      path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

        leg = sp3.legend(loc=3, numpoints=1, fontsize=18)
























    ####################################################
    ###  Plotting SFR-M* for LOCAL overdensity bins  ###
    ####################################################

    if True:



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

            n_galaxies_sf_local.append([])
            z_galaxies_sf_local.append([])
            overdens_median_sf_local.append([])
            overdens_nmad_sf_local.append([])
            sfr_median_sf_local.append([])
            sfr_nmad_sf_local.append([])
            lmass_median_sf_local.append([])
            lmass_nmad_sf_local.append([])

            for i_overdens in range(len(galaxies_catalog.overdensbars)):

                inds_lmass_overdens_sf = numpy.where((uvj_class == 1) & \
                                                     (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                                     (galaxies_catalog.digi_overdens == i_overdens+1))[0]


                ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

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

        lsfr_offsets_sf = []

        lsfr_offsets_sf_local = []

        for i_lmass in range(len(galaxies_catalog.lmassbars)):

            lsfr_offsets_sf_local.append([])

            for i_overdens in range(len(galaxies_catalog.overdensbars)):

                ###  LOCAL overdensity
                z_medi_sf = numpy.median(z_galaxies_sf_local[i_lmass][i_overdens])
                lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_offsets_sf_local[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)


        lsfr_offsets_sf_local = numpy.array(lsfr_offsets_sf_local)













        ###  plotting SFR vs. lmass

        lmass_median_sf_local = numpy.array(lmass_median_sf_local)
        lmass_nmad_sf_local = numpy.array(lmass_nmad_sf_local)
        sfr_median_sf_local = numpy.array(sfr_median_sf_local)
        sfr_nmad_sf_local = numpy.array(sfr_nmad_sf_local)
        overdens_median_sf_local = numpy.array(overdens_median_sf_local)
        overdens_nmad_sf_local = numpy.array(overdens_nmad_sf_local)
        n_galaxies_sf_local = numpy.array(n_galaxies_sf_local)



        ###  applying SFR offsets
        lsfr_median_sf_local = numpy.log10(sfr_median_sf_local) + lsfr_offsets_sf_local





        colors_local = ['#001ddd', '#20c600', '#d11f1f']
        overdens_labels = ['          log(1+$\\delta$) < 0.5',
                           '0.5 < log(1+$\\delta$) < 1.0',
                           '1.0 < log(1+$\\delta$)']


        for i_overdens in range(len(galaxies_catalog.overdensbars)):

            ###  don't plot lowest mass bin for intermediate- or high-density sample
            inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf_local[:, i_overdens] > 8)


            elog_sfr_sf = numpy.log10((sfr_median_sf_local[:, i_overdens] + sfr_nmad_sf_local[:, i_overdens]/n_galaxies_sf_local[:, i_overdens]**0.5) / sfr_median_sf_local[:, i_overdens])

            color_overdens = colors_local[i_overdens]

            sp2.errorbar(lmass_median_sf_local[inds_plot_sf, i_overdens], lsfr_median_sf_local[inds_plot_sf, i_overdens], 
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
        lsfr_lower = log_SFR_tomczak2016(z16, lmass_model, SF_galaxies=True)
        lsfr_upper = log_SFR_tomczak2016(z84, lmass_model, SF_galaxies=True)

        sp2.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
        sp2.axvline(0, color='#cccccc', lw=14, label=zfourge_label, zorder=1)





        t4 = sp2.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                      transform=sp2.transAxes, horizontalalignment='left', verticalalignment='top',
                      path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

        leg = sp2.legend(loc=4, numpoints=1, fontsize=17, frameon=0)



    pdf.savefig()
    sp1.clear()
    sp2.clear()
    sp3.clear()

pdf.close()




































######################################
###  For GLOBAL evironment metric  ###
######################################

pdf = PdfPages('../figures/uvj_boundary_variation_GLOBAL.pdf')

for i_vertices in range(len(vertices_all)):

    fig_setup_GLOBAL()

    v = vertices_all[i_vertices]

    sp1.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='#c40000', lw=4, zorder=3)
    sp1.plot([v[1][0], v[2][0]], [v[1][1], v[2][1]], color='#c40000', lw=4, zorder=3)



    ###  classifying SF vs. QU galaxies

    uvj_class = uvj_selection(galaxies_catalog.UV_color, 
                              galaxies_catalog.VJ_color,
                              delta_y=delta_ys[i_vertices])


    n_bootstraps = 1000











    
    ##############################
    ###  Plotting SFR-density  ###
    ##############################

    if True:

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

        sfr_median_sf_bootstraps = []


        for i_lmass in range(1, len(lmassbins2)+1):

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

            sfr_median_sf_bootstraps.append([])

            for i_overdens in range(len(eta_bins)+1):

                inds_lmass_overdens_sf = numpy.where((uvj_class == 1) & \
                                                     (digi_lmass2 == i_lmass) & \
                                                     (digi_eta == i_overdens))[0]


                ngals_sf = len(inds_lmass_overdens_sf)

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
                sfr_median_sf_bootstraps[-1].append([])
            
                for i_bootstrap in range(n_bootstraps):

                    inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_sf), len(inds_lmass_overdens_sf))
                    inds_lmass_overdens_sf_bootstrap = inds_lmass_overdens_sf[inds_resamp]
                    sfr_median_sf_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf_bootstrap]))






        ###  estimating SFR offsets from median to fiducial redshifts
        z_fiducial = numpy.median(galaxies_catalog.zspec)

        lsfr_offsets_sf = []

        for i_lmass in range(len(z_galaxies_sf)):

            lsfr_offsets_sf.append([])

            for i_overdens in range(len(eta_bins)+1):

                z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])
                lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_offsets_sf[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)

        lsfr_offsets_sf = numpy.array(lsfr_offsets_sf)



        ###  plotting SSFR vs. density

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
        lssfr_median_sf = numpy.log10(ssfr_median_sf) + lsfr_offsets_sf





        ###  Bootstrap scatters
        sfr_median_sf_bootstraps = numpy.array(sfr_median_sf_bootstraps)

        sfr_median_sf_bootstraps[sfr_median_sf_bootstraps <= 0] = 10**-5

        sfr_p16_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 16, axis=2)
        sfr_p50_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 50, axis=2)
        sfr_p84_sf_bootstrap = numpy.percentile(sfr_median_sf_bootstraps, 84, axis=2)

        lsfr_p16_sf_bootstrap = numpy.log10(sfr_p16_sf_bootstrap)
        lsfr_p50_sf_bootstrap = numpy.log10(sfr_p50_sf_bootstrap)
        lsfr_p84_sf_bootstrap = numpy.log10(sfr_p84_sf_bootstrap)





        colors = ['b', 'yellow', 'orange']


        for i_lmass in range(len(z_galaxies_sf)):

            ###  error on the median
            elog_sfr_sf = numpy.log10((sfr_median_sf[i_lmass] + sfr_nmad_sf[i_lmass]/n_galaxies_sf[i_lmass]**0.5) / sfr_median_sf[i_lmass])


            ###  bootstrap error
            elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[i_lmass]
            ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[i_lmass]
            elog_sfr_sf = [elo_lsfr_sf, ehi_lsfr_sf]



            ###  making bootstrap errors symmetric by
            ###  taking average of lower, upper errors
            elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.



            color_overdens = colors[i_lmass]

            sp3.plot(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                              ls='-', lw=3.5, color='k', zorder=1)

            sp3.errorbar(overdens_median_sf[i_lmass], lssfr_median_sf[i_lmass], 
                         yerr=elog_sfr_sf, xerr=overdens_nmad_sf[i_lmass] * 0,
                         ls='-', lw=1.5, marker='o', mew=1.8, ms=9, elinewidth=2, zorder=2,
                         mfc=color_overdens, color=color_overdens, mec='k', ecolor='k', 
                         label=lmassbin2_labels[i_lmass])


        t2 = sp3.text(0.03, 0.97, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                      transform=sp3.transAxes, horizontalalignment='left', verticalalignment='top',
                      path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

        leg = sp3.legend(loc=4, numpoints=1, fontsize=18)
    























    #####################################################
    ###  Plotting SFR-M* for GLOBAL overdensity bins  ###
    #####################################################

    if True:

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

        sfr_median_sf_bootstraps = []


        for i_lmass in range(len(galaxies_catalog.lmassbars)):

            n_galaxies_sf.append([])
            z_galaxies_sf.append([])
            overdens_median_sf.append([])
            overdens_nmad_sf.append([])
            sfr_median_sf.append([])
            sfr_nmad_sf.append([])
            lmass_median_sf.append([])
            lmass_nmad_sf.append([])

            for i_overdens in range(len(eta_bins)+1):

                inds_lmass_overdens_sf = (uvj_class == 1) & \
                                         (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                         (digi_eta == i_overdens)


                ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

                n_galaxies_sf[-1].append(ngals_sf)
                z_galaxies_sf[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
                sfr_median_sf[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
                sfr_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
                overdens_median_sf[-1].append(numpy.median(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
                overdens_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.log_eta[inds_lmass_overdens_sf]))
                lmass_median_sf[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
                lmass_nmad_sf[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))












        ###  estimating SFR offsets from median to fiducial redshifts
        z_fiducial = numpy.median(galaxies_catalog.zspec)

        lsfr_offsets_sf = []

        lsfr_offsets_sf = []

        for i_lmass in range(len(galaxies_catalog.lmassbars)):

            lsfr_offsets_sf.append([])

            for i_overdens in range(len(eta_bins)+1)[::-1]:

                ###  GLOBAL overdensity
                z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])
                lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)
                lsfr_offsets_sf[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)


        lsfr_offsets_sf = numpy.array(lsfr_offsets_sf)












        ###  plotting SFR vs. lmass

        lmass_median_sf = numpy.array(lmass_median_sf)
        lmass_nmad_sf = numpy.array(lmass_nmad_sf)
        sfr_median_sf = numpy.array(sfr_median_sf)
        sfr_nmad_sf = numpy.array(sfr_nmad_sf)
        overdens_median_sf = numpy.array(overdens_median_sf)
        overdens_nmad_sf = numpy.array(overdens_nmad_sf)
        n_galaxies_sf = numpy.array(n_galaxies_sf)



        ###  applying SFR offsets
        lsfr_median_sf = numpy.log10(sfr_median_sf) + lsfr_offsets_sf



        colors_global = ['r', 'orange', 'g', 'b']


        for i_overdens in range(len(eta_bins)+1)[::-1]:

            ###  don't plot lowest mass bin for intermediate- or high-density sample
            inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf[:, i_overdens] > 8)


            elog_sfr_sf = numpy.log10((sfr_median_sf[:, i_overdens] + sfr_nmad_sf[:, i_overdens]/n_galaxies_sf[:, i_overdens]**0.5) / sfr_median_sf[:, i_overdens])

            color_overdens = colors_global[i_overdens]

            sp2.errorbar(lmass_median_sf[inds_plot_sf, i_overdens], lsfr_median_sf[inds_plot_sf, i_overdens], 
                         yerr=elog_sfr_sf[inds_plot_sf], xerr=lmass_nmad_sf[inds_plot_sf, i_overdens],
                         ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                         mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                         label=eta_labels[i_overdens])




        ###  plotting ZFOURGE

        z16 = numpy.percentile(galaxies_catalog.zspec, 16)
        z84 = numpy.percentile(galaxies_catalog.zspec, 84)

        zfourge_label = 'ZFOURGE:  z ~ %.2f' % numpy.median(galaxies_catalog.zspec)
        zfourge_label = 'Tomczak et al. 2016\nz ~ %.2f' % numpy.median(galaxies_catalog.zspec)

        lmass_model = numpy.linspace(8.3, 11.7, 100)
        lsfr_lower = log_SFR_tomczak2016(z16, lmass_model, SF_galaxies=True)
        lsfr_upper = log_SFR_tomczak2016(z84, lmass_model, SF_galaxies=True)

        sp2.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
        sp2.axvline(0, color='#cccccc', lw=14, label=zfourge_label, zorder=1)





        t4 = sp2.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                      transform=sp2.transAxes, horizontalalignment='left', verticalalignment='top',
                      path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

        leg = sp2.legend(loc=4, numpoints=1, fontsize=17, frameon=0)































    pdf.savefig()
    sp1.clear()
    sp2.clear()
    sp3.clear()



pdf.close()




