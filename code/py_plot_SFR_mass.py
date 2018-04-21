
import mypy
import numpy
import pandas
import pickle
from scipy import optimize, stats
from astropy.io import fits
from matplotlib import pyplot
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages




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



catalog =          pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))
galaxies_catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))
galaxies_catalog = pickle.load(open('../data/table_all_galaxies_noXrayAGN.pickle', 'rb'))


inds = (catalog.zspec >= 0.92) & \
       (catalog.zspec <= 1.50)

inds = (catalog.zspec >= 0.70) & \
       (catalog.zspec <= 0.92)

catalog.SFR_FAST = catalog.SFR_FAST[inds]
catalog.SFR_UVIR = catalog.SFR_UVIR[inds]
catalog.UVJ_class = catalog.UVJ_class[inds]
catalog.UV_color = catalog.UV_color[inds]
catalog.VJ_color = catalog.VJ_color[inds]
catalog.dec = catalog.dec[inds]
catalog.digi_lmass = catalog.digi_lmass[inds]
catalog.digi_overdens = catalog.digi_overdens[inds]
catalog.id_field = catalog.id_field[inds]
catalog.id_gal = catalog.id_gal[inds]
catalog.lmass = catalog.lmass[inds]
catalog.overdens = catalog.overdens[inds]
catalog.ra = catalog.ra[inds]
catalog.zspec = catalog.zspec[inds]












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


















#################################################################
###  Plotting for Star-Forming subsample in overdensity bins  ###
#################################################################

if False:

    fig = pyplot.figure(figsize=(9.6875, 7.9749999999999996))

    sp1 = fig.add_subplot(111)

    sp1.minorticks_on()
    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp1.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1.axis([8.6, 12.3, -0.4, 2.1])
    sp1.axis([8.1, 12.4, -0.85, 2.05])


    z_all = []
    cmap = pyplot.cm.rainbow


    colors = ['b', 'g', 'r']

    for i_overdens in range(len(catalog.overdensbars)):


        n_galaxies = []
        z_galaxies = []
        overdens_median = []
        overdens_nmad = []
        sfr_median = []
        sfr_nmad = []
        lmass_median = []
        lmass_nmad = []

        for i_lmass in range(len(catalog.lmassbars)):

            inds_lmass_overdens = (catalog.UVJ_class == 1) & \
                                  (catalog.digi_lmass == i_lmass+1) & \
                                  (catalog.digi_overdens == i_overdens+1)

            ngals = numpy.count_nonzero(inds_lmass_overdens)

            ###  skipping log(M*)~8.5
            if catalog.lmassbars[i_lmass] < 8.51 and ngals < 8:
                continue

            n_galaxies.append(ngals)
            z_galaxies += catalog.zspec[inds_lmass_overdens].tolist()
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))
            overdens_nmad.append(mypy.nmad(catalog.overdens[inds_lmass_overdens]))
            lmass_median.append(numpy.median(catalog.lmass[inds_lmass_overdens]))
            lmass_nmad.append(mypy.nmad(catalog.lmass[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        z_all += z_galaxies
        z_galaxies = numpy.array(z_galaxies)
        overdens_median = numpy.array(overdens_median)
        overdens_nmad = numpy.array(overdens_nmad)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        lmass_median = numpy.array(lmass_median)
        lmass_nmad = numpy.array(lmass_nmad)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_overdens = colors[i_overdens]


        ###  plotting SFR vs. lmass
        sp1.errorbar(lmass_median, numpy.log10(sfr_median), yerr=elog_sfr, xerr=lmass_nmad,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=catalog.overdens_labels[i_overdens])




    ###  plotting ZFOURGE

    z_all = numpy.array(z_all)

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(z_all, 16), lmass_model)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(z_all, 84), lmass_model)

    sp1.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(z_all), zorder=1)









    t = sp1.text(0.03, 0.97, 'Star-Forming', fontweight='normal', fontsize=24, color='b', 
                 transform=sp1.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp1.legend(loc=4, numpoints=1, fontsize=17, frameon=0)






    ###  plotting full sample

    if False:

        n_galaxies = []
        overdens_median = []
        sfr_median = []
        sfr_nmad = []
        lmass_median = []

        for i_lmass in range(len(catalog.lmassbars)):


            ###  skipping log(M*)~8.5
            if catalog.lmassbars[i_lmass] < 8.51:
                continue

            inds_lmass_overdens = (catalog.UVJ_class == 1) & \
                                  (catalog.digi_lmass == i_lmass+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))
            lmass_median.append(numpy.median(catalog.lmass[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        lmass_median = numpy.array(lmass_median)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_overdens = 'gray'


        ###  plotting SFR vs. lmass
        sp1.errorbar(lmass_median, numpy.log10(sfr_median), yerr=elog_sfr, ls='-', lw=1.5, marker='o', 
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     mew=1.5, ms=9, elinewidth=2, zorder=2)
























#######################################################
###  Plotting for Total sample in overdensity bins  ###
#######################################################

if False:

    fig = pyplot.figure(figsize=(9.6875, 7.9749999999999996))

    sp1 = fig.add_subplot(111)

    sp1.minorticks_on()
    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp1.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1.axis([8.6, 12.3, -0.4, 2.1])
    sp1.axis([8.1, 12.4, -0.85, 2.05])


    z_all = []
    cmap = pyplot.cm.rainbow


    colors = ['b', 'g', 'r']

    for i_overdens in range(len(catalog.overdensbars)):


        n_galaxies = []
        z_galaxies = []
        overdens_median = []
        overdens_nmad = []
        sfr_mean = []
        sfr_median = []
        sfr_nmad = []
        lmass_median = []
        lmass_nmad = []

        for i_lmass in range(len(catalog.lmassbars)):

            inds_lmass_overdens = (catalog.digi_lmass == i_lmass+1) & \
                                  (catalog.digi_overdens == i_overdens+1)

            ngals = numpy.count_nonzero(inds_lmass_overdens)

            ###  skipping log(M*)~8.5
            if catalog.lmassbars[i_lmass] < 8.51 and ngals < 8:
                continue

            n_galaxies.append(ngals)
            z_galaxies += catalog.zspec[inds_lmass_overdens].tolist()
            sfr_mean.append(numpy.mean(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))
            overdens_nmad.append(mypy.nmad(catalog.overdens[inds_lmass_overdens]))
            lmass_median.append(numpy.median(catalog.lmass[inds_lmass_overdens]))
            lmass_nmad.append(mypy.nmad(catalog.lmass[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        z_all += z_galaxies
        z_galaxies = numpy.array(z_galaxies)
        overdens_median = numpy.array(overdens_median)
        overdens_nmad = numpy.array(overdens_nmad)
        sfr_mean = numpy.array(sfr_mean)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        lmass_median = numpy.array(lmass_median)
        lmass_nmad = numpy.array(lmass_nmad)
        n = len(sfr_median)

        print sfr_mean

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_overdens = colors[i_overdens]


        ###  plotting SFR vs. lmass
        sp1.errorbar(lmass_median, numpy.log10(sfr_median), yerr=elog_sfr, xerr=lmass_nmad,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=catalog.overdens_labels[i_overdens])




    ###  plotting ZFOURGE

    z_all = numpy.array(z_all)

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(z_all, 16), lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(z_all, 84), lmass_model, SF_galaxies=False)

    sp1.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(z_all), zorder=1)









    t = sp1.text(0.03, 0.97, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                 transform=sp1.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp1.legend(loc=4, numpoints=1, fontsize=17, frameon=0)






    ###  plotting full sample

    if False:

        n_galaxies = []
        overdens_median = []
        sfr_median = []
        sfr_nmad = []
        lmass_median = []

        for i_lmass in range(len(catalog.lmassbars)):


            ###  skipping log(M*)~8.5
            if catalog.lmassbars[i_lmass] < 8.51:
                continue

            inds_lmass_overdens = (catalog.UVJ_class == 1) & \
                                  (catalog.digi_lmass == i_lmass+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))
            lmass_median.append(numpy.median(catalog.lmass[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        lmass_median = numpy.array(lmass_median)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_overdens = 'gray'


        ###  plotting SFR vs. lmass
        sp1.errorbar(lmass_median, numpy.log10(sfr_median), yerr=elog_sfr, ls='-', lw=1.5, marker='o', 
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     mew=1.5, ms=9, elinewidth=2, zorder=2)
























#############################################
###  Plotting SF and TOT panels together  ###
#############################################

if False:

    fig = pyplot.figure(figsize=(16.5, 8.1))

    sp1 = fig.add_subplot(122)
    sp2 = fig.add_subplot(121)

    sp1.minorticks_on()
    sp2.minorticks_on()
    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp1.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')
    sp2.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1.axis([8.6, 12.3, -0.4, 2.1])
    sp1.axis([8.1, 12.4, -0.85, 2.05])
    sp2.axis([8.1, 12.4, -0.85, 2.05])


    fig.subplots_adjust(wspace=0, left=0.07, right=0.98, )


    z_all_tot = []
    z_all_sf = []
    cmap = pyplot.cm.rainbow


    colors = ['b', 'g', 'r']

    for i_overdens in range(len(catalog.overdensbars)):


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

        for i_lmass in range(len(catalog.lmassbars)):

            inds_lmass_overdens_tot = (catalog.digi_lmass == i_lmass+1) & \
                                      (catalog.digi_overdens == i_overdens+1)

            inds_lmass_overdens_sf = (catalog.UVJ_class == 1) & \
                                     (catalog.digi_lmass == i_lmass+1) & \
                                     (catalog.digi_overdens == i_overdens+1)


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

            ###  skipping log(M*)~8.5
            if catalog.lmassbars[i_lmass] < 8.51 and ngals_tot < 8:
                continue

            n_galaxies_tot.append(ngals_tot)
            z_galaxies_tot += catalog.zspec[inds_lmass_overdens_tot].tolist()
            sfr_median_tot.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            overdens_median_tot.append(numpy.median(catalog.overdens[inds_lmass_overdens_tot]))
            overdens_nmad_tot.append(mypy.nmad(catalog.overdens[inds_lmass_overdens_tot]))
            lmass_median_tot.append(numpy.median(catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot.append(mypy.nmad(catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf.append(ngals_sf)
            z_galaxies_sf += catalog.zspec[inds_lmass_overdens_sf].tolist()
            sfr_median_sf.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            overdens_median_sf.append(numpy.median(catalog.overdens[inds_lmass_overdens_sf]))
            overdens_nmad_sf.append(mypy.nmad(catalog.overdens[inds_lmass_overdens_sf]))
            lmass_median_sf.append(numpy.median(catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf.append(mypy.nmad(catalog.lmass[inds_lmass_overdens_sf]))

        n_galaxies_tot = numpy.array(n_galaxies_tot)
        z_all_tot += z_galaxies_tot
        z_galaxies_tot = numpy.array(z_galaxies_tot)
        overdens_median_tot = numpy.array(overdens_median_tot)
        overdens_nmad_tot = numpy.array(overdens_nmad_tot)
        sfr_median_tot = numpy.array(sfr_median_tot)
        sfr_nmad_tot = numpy.array(sfr_nmad_tot)
        lmass_median_tot = numpy.array(lmass_median_tot)
        lmass_nmad_tot = numpy.array(lmass_nmad_tot)
        ntot = len(sfr_median_tot)
        elog_sfr_tot = numpy.log10((sfr_median_tot + sfr_nmad_tot/n_galaxies_tot**0.5) / sfr_median_tot)

        n_galaxies_sf = numpy.array(n_galaxies_sf)
        z_all_sf += z_galaxies_sf
        z_galaxies_sf = numpy.array(z_galaxies_sf)
        overdens_median_sf = numpy.array(overdens_median_sf)
        overdens_nmad_sf = numpy.array(overdens_nmad_sf)
        sfr_median_sf = numpy.array(sfr_median_sf)
        sfr_nmad_sf = numpy.array(sfr_nmad_sf)
        lmass_median_sf = numpy.array(lmass_median_sf)
        lmass_nmad_sf = numpy.array(lmass_nmad_sf)
        nsf = len(sfr_median_sf)
        elog_sfr_sf = numpy.log10((sfr_median_sf + sfr_nmad_sf/n_galaxies_sf**0.5) / sfr_median_sf)






        ###  plotting SFR vs. lmass
        color_overdens = colors[i_overdens]

        sp1.errorbar(lmass_median_sf, numpy.log10(sfr_median_sf), yerr=elog_sfr_sf, xerr=lmass_nmad_sf,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=catalog.overdens_labels[i_overdens])

        sp2.errorbar(lmass_median_tot, numpy.log10(sfr_median_tot), yerr=elog_sfr_tot, xerr=lmass_nmad_tot,
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=catalog.overdens_labels[i_overdens])




    ###  plotting ZFOURGE

    z_all_sf = numpy.array(z_all_sf)

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(z_all_sf, 16), lmass_model)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(z_all_sf, 84), lmass_model)

    sp1.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(z_all_sf), zorder=1)

    print 'SF:  %.3f < z < %.3f' % (numpy.percentile(z_all_sf, 16), numpy.percentile(z_all_sf, 84))



    z_all_tot = numpy.array(z_all_tot)

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(z_all_tot, 16), lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(z_all_tot, 84), lmass_model, SF_galaxies=False)

    sp2.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(z_all_tot), zorder=1)

    print 'ALL: %.3f < z < %.3f' % (numpy.percentile(z_all_tot, 16), numpy.percentile(z_all_tot, 84))





    t1 = sp1.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp1.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t2 = sp2.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp2.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp1.legend(loc=4, numpoints=1, fontsize=17, frameon=0)
























#############################################
###  Checking redshift distributions for  ### 
###  each mass-overdensity bin            ###
#############################################

if False:



    pdf = PdfPages('../figures/redshift_distributions.pdf')



    fig = pyplot.figure(figsize=(16.5, 11.))

    sp2 = fig.add_subplot(222)
    sp1 = fig.add_subplot(221)
    sp4 = fig.add_subplot(224)
    sp3 = fig.add_subplot(223)


    fig.subplots_adjust(left=0.07, right=0.98, hspace=0, wspace=0)


    z_all_tot = []
    z_all_sf = []
    hists_tot = []
    hists_sf = []
    cmap = pyplot.cm.rainbow


    colors = ['b', 'g', 'r']

    for i_lmass in range(len(catalog.lmassbars)):


        sp1.minorticks_on()
        sp2.minorticks_on()
        sp3.minorticks_on()
        sp4.minorticks_on()
        sp1.set_xlabel('spectroscopic redshift')
        sp2.set_xlabel('spectroscopic redshift')
        sp3.set_xlabel('spectroscopic redshift')
        sp4.set_xlabel('spectroscopic redshift')
        sp1.set_ylabel('Number')
        sp2.set_ylabel('Number')
        sp3.set_ylabel('CDF')
        sp4.set_ylabel('CDF')

        sp1.set_title('All Galaxies')
        sp2.set_title('Star-Forming Galaxies')

        sp3.axis([0.55, 1.35, 0, 1.1])
        sp4.axis([0.55, 1.35, 0, 1.1])


        z_galaxies_tot = []
        z_galaxies_sf = []

        for i_overdens in range(len(catalog.overdensbars)):
    
            inds_lmass_overdens_tot = (catalog.digi_lmass == i_lmass+1) & \
                                      (catalog.digi_overdens == i_overdens+1)

            inds_lmass_overdens_sf = (catalog.UVJ_class == 1) & \
                                     (catalog.digi_lmass == i_lmass+1) & \
                                     (catalog.digi_overdens == i_overdens+1)


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)


            z_galaxies_tot.append(catalog.zspec[inds_lmass_overdens_tot].tolist())
            z_galaxies_sf.append(catalog.zspec[inds_lmass_overdens_sf].tolist())


            ###  plotting
            zlo, zhi = 0.7, 1.3
            offsets = [-0.0015, 0.00, 0.0015]        
            color_overdens = colors[i_overdens]



            h_tot = sp1.hist(z_galaxies_tot[-1], histtype='step', lw=3, color=color_overdens,
                             bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_tot[-1])))

            h_sf = sp2.hist(z_galaxies_sf[-1], histtype='step', lw=3, color=color_overdens,
                            bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_sf[-1])))

            hists_tot.append(h_tot[0])
            hists_sf.append(h_sf[0])


            sp3.hist(z_galaxies_tot[-1], histtype='step', lw=3, color=color_overdens,
                     bins=1000, range=(0.7, 1.3), cumulative=True, normed=True)

            sp4.hist(z_galaxies_sf[-1], histtype='step', lw=3, color=color_overdens,
                     bins=1000, range=(0.7, 1.3), cumulative=True, normed=True)


        sp1.set_xlim(sp3.get_xlim())
        sp2.axis(sp1.axis())

        sp1.legend(loc=1)
        sp2.legend(loc=1)


        ###  KS tests
        ks_field_groups_tot = stats.ks_2samp(z_galaxies_tot[0], z_galaxies_tot[1])
        ks_field_clusters_tot = stats.ks_2samp(z_galaxies_tot[0], z_galaxies_tot[2])

        ks_field_groups_sf = stats.ks_2samp(z_galaxies_sf[0], z_galaxies_sf[1])
        ks_field_clusters_sf = stats.ks_2samp(z_galaxies_sf[0], z_galaxies_sf[2])



        ###  plotting text
        info = 'log(M$_*$) = %.1f' % catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_tot[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_tot[1]

        sp3.text(0.03, 0.97, info, transform=sp3.transAxes, 
                 horizontalalignment='left', verticalalignment='top')

        info = 'log(M$_*$) = %.1f' % catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_sf[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_sf[1]

        sp4.text(0.03, 0.97, info, transform=sp4.transAxes, 
                 horizontalalignment='left', verticalalignment='top')









        pdf.savefig()
        sp1.clear()
        sp2.clear()
        sp3.clear()
        sp4.clear()

    pdf.close()
    pyplot.close()
























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

if False:



    pdf_zdists = PdfPages('../figures/redshift_distributions_weighted.pdf')



    fig_zdists = pyplot.figure(figsize=(16.5, 11.))

    sp2_zdists = fig_zdists.add_subplot(222)
    sp1_zdists = fig_zdists.add_subplot(221)
    sp4_zdists = fig_zdists.add_subplot(224)
    sp3_zdists = fig_zdists.add_subplot(223)


    fig_zdists.subplots_adjust(left=0.07, right=0.98, hspace=0, wspace=0)




    z_all_tot = []
    z_all_sf = []
    hists_tot = []
    hists_sf = []
    cmap = pyplot.cm.rainbow


    colors = ['b', 'g', 'r']

    overdens_labels = ['          log(1+$\\delta_{gal}$) < 0.5',
                       '0.5 < log(1+$\\delta_{gal}$) < 1.0',
                       '1.0 < log(1+$\\delta_{gal}$)']


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

    n_bootstraps = 1000
    sfr_median_tot_bootstraps = []
    sfr_median_sf_bootstraps = []


    for i_lmass in range(len(catalog.lmassbars)):

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

        sfr_median_tot_bootstraps.append([])
        sfr_median_sf_bootstraps.append([])

        for i_overdens in range(len(catalog.overdensbars)):

            inds_lmass_overdens_tot = numpy.where((catalog.digi_lmass == i_lmass+1) & \
                                                  (catalog.digi_overdens == i_overdens+1))[0]

            inds_lmass_overdens_sf = numpy.where((catalog.UVJ_class == 1) & \
                                                 (catalog.digi_lmass == i_lmass+1) & \
                                                 (catalog.digi_overdens == i_overdens+1))[0]


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

            n_galaxies_tot[-1].append(ngals_tot)
            z_galaxies_tot[-1].append(catalog.zspec[inds_lmass_overdens_tot].tolist())
            sfr_median_tot[-1].append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot[-1].append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            overdens_median_tot[-1].append(numpy.median(catalog.overdens[inds_lmass_overdens_tot]))
            overdens_nmad_tot[-1].append(mypy.nmad(catalog.overdens[inds_lmass_overdens_tot]))
            lmass_median_tot[-1].append(numpy.median(catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot[-1].append(mypy.nmad(catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf[-1].append(ngals_sf)
            z_galaxies_sf[-1].append(catalog.zspec[inds_lmass_overdens_sf].tolist())
            sfr_median_sf[-1].append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf[-1].append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            overdens_median_sf[-1].append(numpy.median(catalog.overdens[inds_lmass_overdens_sf]))
            overdens_nmad_sf[-1].append(mypy.nmad(catalog.overdens[inds_lmass_overdens_sf]))
            lmass_median_sf[-1].append(numpy.median(catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf[-1].append(mypy.nmad(catalog.lmass[inds_lmass_overdens_sf]))


            ###  running Bootstrap resamplings
            sfr_median_tot_bootstraps[-1].append([])
            sfr_median_sf_bootstraps[-1].append([])
        
            for i_bootstrap in range(n_bootstraps):

                inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_tot), len(inds_lmass_overdens_tot))
                inds_lmass_overdens_tot_bootstrap = inds_lmass_overdens_tot[inds_resamp]
                sfr_median_tot_bootstraps[-1][-1].append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens_tot_bootstrap]))

                inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_sf), len(inds_lmass_overdens_sf))
                inds_lmass_overdens_sf_bootstrap = inds_lmass_overdens_sf[inds_resamp]
                sfr_median_sf_bootstraps[-1][-1].append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens_sf_bootstrap]))





            ###  plotting redshift distributions
            zlo, zhi = 0.7, 1.3
            offsets = [-0.0015, 0.00, 0.0015]        
            color_overdens = colors[i_overdens]



            h_tot = sp1_zdists.hist(z_galaxies_tot[-1][-1], histtype='step', lw=3, color=color_overdens,
                             bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_tot[-1][-1])))

            h_sf = sp2_zdists.hist(z_galaxies_sf[-1][-1], histtype='step', lw=3, color=color_overdens,
                            bins=15, range=(0.7+offsets[i_overdens], 1.3+offsets[i_overdens]), 
                             label=str(len(z_galaxies_sf[-1][-1])))

            hists_tot.append(h_tot[0])
            hists_sf.append(h_sf[0])


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
        info = 'log(M$_*$) = %.1f' % catalog.lmassbars[i_lmass]
        info += '\nKS: p$_{f-g}$ = %.3f' % ks_field_groups_tot[1]
        info += '\nKS: p$_{f-c}$ = %.3f' % ks_field_clusters_tot[1]

        sp3_zdists.text(0.03, 0.97, info, transform=sp3_zdists.transAxes, 
                 horizontalalignment='left', verticalalignment='top')

        info = 'log(M$_*$) = %.1f' % catalog.lmassbars[i_lmass]
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
    z_fiducial = numpy.median(catalog.zspec)

    lsfr_offsets_tot = []
    lsfr_offsets_sf = []

    for i_lmass in range(len(catalog.lmassbars)):

        lsfr_offsets_tot.append([])
        lsfr_offsets_sf.append([])

        for i_overdens in range(len(catalog.overdensbars)):

            z_medi_tot = numpy.median(z_galaxies_tot[i_lmass][i_overdens])
            z_medi_sf = numpy.median(z_galaxies_sf[i_lmass][i_overdens])

            lsfr_fiducial_tot = log_SFR_tomczak2016(z_fiducial, catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_z_medi_tot = log_SFR_tomczak2016(z_medi_tot, catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, catalog.lmassbars[i_lmass], SF_galaxies=True)

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



    for i_overdens in range(len(catalog.overdensbars)):


        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (catalog.lmassbars > 8.51) | (n_galaxies_tot[:,i_overdens] > 8)
        inds_plot_sf = (catalog.lmassbars > 8.51) | (n_galaxies_sf[:,i_overdens] > 8)

        ###  error on the median
        elog_sfr_tot = numpy.log10((sfr_median_tot[inds_plot_tot,i_overdens] + sfr_nmad_tot[inds_plot_tot,i_overdens]/n_galaxies_tot[inds_plot_tot,i_overdens]**0.5) / sfr_median_tot[inds_plot_tot,i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf[inds_plot_sf,i_overdens] + sfr_nmad_sf[inds_plot_sf,i_overdens]/n_galaxies_sf[inds_plot_sf,i_overdens]**0.5) / sfr_median_sf[inds_plot_sf,i_overdens])

        ###  bootstrap error
        elo_lsfr_tot = (lsfr_p50_tot_bootstrap - lsfr_p16_tot_bootstrap)[:,i_overdens]
        ehi_lsfr_tot = (lsfr_p84_tot_bootstrap - lsfr_p50_tot_bootstrap)[:,i_overdens]
        elog_sfr_tot = [elo_lsfr_tot[inds_plot_tot], ehi_lsfr_tot[inds_plot_tot]]

        elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[:,i_overdens]
        ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[:,i_overdens]
        elog_sfr_sf = [elo_lsfr_sf[inds_plot_sf], ehi_lsfr_sf[inds_plot_sf]]


        color_overdens = colors[i_overdens]

        sp1_sfrmass.errorbar(lmass_median_tot[inds_plot_tot, i_overdens], lsfr_median_tot[inds_plot_tot, i_overdens], 
                     yerr=elog_sfr_tot, xerr=lmass_nmad_tot[inds_plot_tot, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=overdens_labels[i_overdens])

        sp2_sfrmass.errorbar(lmass_median_sf[inds_plot_sf, i_overdens], lsfr_median_sf[inds_plot_sf, i_overdens], 
                     yerr=elog_sfr_sf, xerr=lmass_nmad_sf[inds_plot_sf, i_overdens],
                     ls='-', lw=1.5, marker='o', mew=1.5, ms=9, elinewidth=2, zorder=2,
                     mfc=color_overdens, color=color_overdens, mec='k', ecolor=color_overdens, 
                     label=overdens_labels[i_overdens])




    ###  plotting ZFOURGE

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(catalog.zspec, 16), lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(catalog.zspec, 84), lmass_model, SF_galaxies=False)

    sp1_sfrmass.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp1_sfrmass.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(catalog.zspec), zorder=1)

    print 'SF:  %.3f < z < %.3f' % (numpy.percentile(catalog.zspec, 16), numpy.percentile(catalog.zspec, 84))



    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(numpy.percentile(catalog.zspec, 16), lmass_model, SF_galaxies=True)
    lsfr_upper = log_SFR_tomczak2016(numpy.percentile(catalog.zspec, 84), lmass_model, SF_galaxies=True)

    sp2_sfrmass.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2_sfrmass.axvline(0, color='#cccccc', lw=10, label='ZFOURGE:  z ~ %.2f' % numpy.median(catalog.zspec), zorder=1)

    print 'ALL: %.3f < z < %.3f' % (numpy.percentile(catalog.zspec, 16), numpy.percentile(catalog.zspec, 84))





    t1 = sp1_sfrmass.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp1_sfrmass.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t2 = sp2_sfrmass.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='b', 
                  transform=sp2_sfrmass.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_sfrmass.legend(loc=4, numpoints=1, fontsize=17, frameon=0)



























#############################################
###  FIGURE FOR PAPER v01
###
###  Plotting SFR-M* in LOCAL overdensity
###  bins for ALL and SF galaxies
#############################################

zlo, zhi = 0.6, 1.3
#zlo, zhi = 0.6, 0.9  # low-z
#zlo, zhi = 0.9, 1.3  # high-z


overdensbins = galaxies_catalog.overdensbins
overdensbars = galaxies_catalog.overdensbars

###  Tomczak's log(1+dgal)
overdens = galaxies_catalog.overdens
digi_overdens = galaxies_catalog.digi_overdens

###  Lemaux's log(1+dgal)
#overdens = ldelta_tomczak_lemaux[:,1]
#digi_overdens = numpy.digitize(overdens, overdensbins)


if True:



    n_galaxies_tot_local = []
    z_galaxies_tot_local = []
    z_galaxies_tot_local_flat = []
    overdens_median_tot_local = []
    overdens_nmad_tot_local = []
    sfr_mean_tot_local = []
    sfr_median_tot_local = []
    sfr_nmad_tot_local = []
    lmass_median_tot_local = []
    lmass_nmad_tot_local = []

    n_galaxies_sf_local = []
    z_galaxies_sf_local = []
    z_galaxies_sf_local_flat = []
    overdens_median_sf_local = []
    overdens_nmad_sf_local = []
    sfr_mean_sf_local = []
    sfr_median_sf_local = []
    sfr_nmad_sf_local = []
    lmass_median_sf_local = []
    lmass_nmad_sf_local = []

    n_bootstraps = 5000
    sfr_median_tot_bootstraps = []
    sfr_median_sf_bootstraps = []


    ###  calculating median SFRs in LOCAL overdensity bins
    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        n_galaxies_tot_local.append([])
        z_galaxies_tot_local.append([])
        overdens_median_tot_local.append([])
        overdens_nmad_tot_local.append([])
        sfr_mean_tot_local.append([])
        sfr_median_tot_local.append([])
        sfr_nmad_tot_local.append([])
        lmass_median_tot_local.append([])
        lmass_nmad_tot_local.append([])

        n_galaxies_sf_local.append([])
        z_galaxies_sf_local.append([])
        overdens_median_sf_local.append([])
        overdens_nmad_sf_local.append([])
        sfr_mean_sf_local.append([])
        sfr_median_sf_local.append([])
        sfr_nmad_sf_local.append([])
        lmass_median_sf_local.append([])
        lmass_nmad_sf_local.append([])

        sfr_median_tot_bootstraps.append([])
        sfr_median_sf_bootstraps.append([])

        for i_overdens in range(len(overdensbars)):

            inds_lmass_overdens_tot = numpy.where((galaxies_catalog.zspec >= zlo) & \
                                                  (galaxies_catalog.zspec <= zhi) & \
                                                  (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                                  (digi_overdens == i_overdens+1))[0]

            inds_lmass_overdens_sf = numpy.where((galaxies_catalog.zspec >= zlo) & \
                                                 (galaxies_catalog.zspec <= zhi) & \
                                                 (galaxies_catalog.UVJ_class == 1) & \
                                                 (galaxies_catalog.digi_lmass == i_lmass+1) & \
                                                 (digi_overdens == i_overdens+1))[0]


            ngals_tot = numpy.count_nonzero(inds_lmass_overdens_tot)
            ngals_sf = numpy.count_nonzero(inds_lmass_overdens_sf)

            n_galaxies_tot_local[-1].append(ngals_tot)
            z_galaxies_tot_local[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist())
            z_galaxies_tot_local_flat += galaxies_catalog.zspec[inds_lmass_overdens_tot].tolist()
            sfr_median_tot_local[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_mean_tot_local[-1].append(numpy.mean(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            sfr_nmad_tot_local[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot]))
            overdens_median_tot_local[-1].append(numpy.median(overdens[inds_lmass_overdens_tot]))
            overdens_nmad_tot_local[-1].append(mypy.nmad(overdens[inds_lmass_overdens_tot]))
            lmass_median_tot_local[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_tot]))
            lmass_nmad_tot_local[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_tot]))

            n_galaxies_sf_local[-1].append(ngals_sf)
            z_galaxies_sf_local[-1].append(galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist())
            z_galaxies_sf_local_flat += galaxies_catalog.zspec[inds_lmass_overdens_sf].tolist()
            sfr_mean_sf_local[-1].append(numpy.mean(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_median_sf_local[-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            sfr_nmad_sf_local[-1].append(mypy.nmad(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf]))
            overdens_median_sf_local[-1].append(numpy.median(overdens[inds_lmass_overdens_sf]))
            overdens_nmad_sf_local[-1].append(mypy.nmad(overdens[inds_lmass_overdens_sf]))
            lmass_median_sf_local[-1].append(numpy.median(galaxies_catalog.lmass[inds_lmass_overdens_sf]))
            lmass_nmad_sf_local[-1].append(mypy.nmad(galaxies_catalog.lmass[inds_lmass_overdens_sf]))




            ###  running Bootstrap resamplings
            sfr_median_tot_bootstraps[-1].append([])
            sfr_median_sf_bootstraps[-1].append([])
        
            for i_bootstrap in range(n_bootstraps):

                if ngals_tot > 0:
                    inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_tot), len(inds_lmass_overdens_tot))
                    inds_lmass_overdens_tot_bootstrap = inds_lmass_overdens_tot[inds_resamp]
                    sfr_median_tot_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_tot_bootstrap]))
                else:
                    sfr_median_tot_bootstraps[-1][-1].append(numpy.nan)


                if ngals_sf > 0:
                    inds_resamp = numpy.random.randint(0, len(inds_lmass_overdens_sf), len(inds_lmass_overdens_sf))
                    inds_lmass_overdens_sf_bootstrap = inds_lmass_overdens_sf[inds_resamp]
                    sfr_median_sf_bootstraps[-1][-1].append(numpy.median(galaxies_catalog.SFR_UVIR[inds_lmass_overdens_sf_bootstrap]))
                else:
                    sfr_median_sf_bootstraps[-1][-1].append(numpy.nan)


































    ###  estimating SFR offsets from median to fiducial redshifts
    z_galaxies_tot_local_flat = numpy.array(z_galaxies_tot_local_flat)
    z_galaxies_sf_local_flat = numpy.array(z_galaxies_sf_local_flat)

    z_fiducial = (zlo + zhi) / 2.
    z_fiducial = numpy.median(galaxies_catalog.zspec)
    z_fiducial = numpy.median(z_galaxies_tot_local_flat)

    lsfr_offsets_tot_local = []
    lsfr_offsets_sf_local = []

    for i_lmass in range(len(galaxies_catalog.lmassbars)):

        lsfr_offsets_tot_local.append([])
        lsfr_offsets_sf_local.append([])

        for i_overdens in range(len(overdensbars)):

            ###  LOCAL overdensity
            z_medi_tot = numpy.median(z_galaxies_tot_local[i_lmass][i_overdens])
            z_medi_sf = numpy.median(z_galaxies_sf_local[i_lmass][i_overdens])

            lsfr_fiducial_tot = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_fiducial_sf = log_SFR_tomczak2016(z_fiducial, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_z_medi_tot = log_SFR_tomczak2016(z_medi_tot, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=False)
            lsfr_z_medi_sf = log_SFR_tomczak2016(z_medi_sf, galaxies_catalog.lmassbars[i_lmass], SF_galaxies=True)

            lsfr_offsets_tot_local[-1].append(lsfr_fiducial_tot - lsfr_z_medi_tot)
            lsfr_offsets_sf_local[-1].append(lsfr_fiducial_sf - lsfr_z_medi_sf)

    lsfr_offsets_tot_local = numpy.array(lsfr_offsets_tot_local)
    lsfr_offsets_sf_local = numpy.array(lsfr_offsets_sf_local)













    ###  plotting SFR vs. lmass
    lmass_median_tot_local = numpy.array(lmass_median_tot_local)
    lmass_nmad_tot_local = numpy.array(lmass_nmad_tot_local)
    sfr_mean_tot_local = numpy.array(sfr_mean_tot_local)
    sfr_median_tot_local = numpy.array(sfr_median_tot_local)
    sfr_nmad_tot_local = numpy.array(sfr_nmad_tot_local)
    overdens_median_tot_local = numpy.array(overdens_median_tot_local)
    overdens_nmad_tot_local = numpy.array(overdens_nmad_tot_local)
    n_galaxies_tot_local = numpy.array(n_galaxies_tot_local)

    lmass_median_sf_local = numpy.array(lmass_median_sf_local)
    lmass_nmad_sf_local = numpy.array(lmass_nmad_sf_local)
    sfr_mean_sf_local = numpy.array(sfr_mean_sf_local)
    sfr_median_sf_local = numpy.array(sfr_median_sf_local)
    sfr_nmad_sf_local = numpy.array(sfr_nmad_sf_local)
    overdens_median_sf_local = numpy.array(overdens_median_sf_local)
    overdens_nmad_sf_local = numpy.array(overdens_nmad_sf_local)
    n_galaxies_sf_local = numpy.array(n_galaxies_sf_local)



    ###  applying SFR offsets
    lsfr_median_tot_local = numpy.log10(sfr_median_tot_local) + lsfr_offsets_tot_local
    lsfr_median_sf_local  = numpy.log10(sfr_median_sf_local) + lsfr_offsets_sf_local

    lsfr_mean_tot_local = numpy.log10(sfr_mean_tot_local) + lsfr_offsets_tot_local
    lsfr_mean_sf_local  = numpy.log10(sfr_mean_sf_local) + lsfr_offsets_sf_local




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













    fig_sfrmass = pyplot.figure(figsize=(14.5, 7.))

    sp2_sfrmass_local = fig_sfrmass.add_subplot(122)
    sp1_sfrmass_local = fig_sfrmass.add_subplot(121)


    #sp1_sfrmass_local.grid()
    #sp2_sfrmass_local.grid()

    sp1_sfrmass_local.minorticks_on()
    sp2_sfrmass_local.minorticks_on()

    sp1_sfrmass_local.set_xlabel('log( $\mathit{M}_*$ / $\mathit{M}_{\odot}$ )', size=22)
    sp2_sfrmass_local.set_xlabel('log( $\mathit{M}_*$ / $\mathit{M}_{\odot}$ )', size=22)
    sp1_sfrmass_local.set_ylabel('log( SFR$_{UV+IR}$ / [$\mathit{M}_{\odot}$ / yr] )', size=22)
    sp2_sfrmass_local.set_ylabel('log( SFR$_{UV+IR}$ / [$\mathit{M}_{\odot}$ / yr] )', size=22)

    sp1_sfrmass_local.axis([8.1, 12.25, -0.85, 2.05])
    sp2_sfrmass_local.axis([8.1, 12.25, -0.85, 2.05])


    fig_sfrmass.subplots_adjust(left=0.08, right=0.99, top=0.99, bottom=0.1, wspace=0)



    colors_local = ['#001ddd', '#20c600', '#d11f1f']
    colors_local = ['#3f4ea1', '#6db388', '#d92120']
    overdens_labels = ['          log(1+$\\delta_{gal}$) < 0.5',
                       '0.5 < log(1+$\\delta_{gal}$) < 1.0',
                       '1.0 < log(1+$\\delta_{gal}$)']


    for i_overdens in range(len(overdensbars)):

        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_tot_local[:, i_overdens] > 8)
        inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf_local[:, i_overdens] > 8)

        ###  throwing out lowest-M* bin of intermediate-density bin due to noisy data
        if i_overdens == 1:
            inds_plot_tot[0] = False
            inds_plot_sf[0] = False

        ###  throwing out lowest-M* two bins of high-density bin due to noisy data
        if i_overdens == 2:
            inds_plot_tot[0:2] = False
            inds_plot_sf[0:2] = False


        ###  errors on median
        #elog_sfr_tot = numpy.log10((sfr_median_tot_local[:, i_overdens] + sfr_nmad_tot_local[:, i_overdens]/n_galaxies_tot_local[:, i_overdens]**0.5) / sfr_median_tot_local[:, i_overdens])
        #elog_sfr_sf = numpy.log10((sfr_median_sf_local[:, i_overdens] + sfr_nmad_sf_local[:, i_overdens]/n_galaxies_sf_local[:, i_overdens]**0.5) / sfr_median_sf_local[:, i_overdens])



        ###  bootstrap error
        elo_lsfr_tot = (lsfr_p50_tot_bootstrap - lsfr_p16_tot_bootstrap)[:, i_overdens]
        ehi_lsfr_tot = (lsfr_p84_tot_bootstrap - lsfr_p50_tot_bootstrap)[:, i_overdens]
        elog_sfr_tot = [elo_lsfr_tot, ehi_lsfr_tot]

        elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[:, i_overdens]
        ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[:, i_overdens]
        elog_sfr_sf = [elo_lsfr_sf, ehi_lsfr_sf]


        ###  making bootstrap errors symmetric by
        ###  taking average of lower, upper errors
        elog_sfr_tot = (elo_lsfr_tot + ehi_lsfr_tot) / 2.
        elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.






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

    z16 = numpy.percentile(z_galaxies_tot_local_flat, 16)
    z84 = numpy.percentile(z_galaxies_tot_local_flat, 84)

    zfourge_label = 'ZFOURGE:  z ~ %.2f' % numpy.median(z_galaxies_tot_local_flat)
    zfourge_label = 'Tomczak et al. 2016\nz ~ %.2f' % numpy.median(z_galaxies_tot_local_flat)

    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(z16, lmass_model, SF_galaxies=False)
    lsfr_upper = log_SFR_tomczak2016(z84, lmass_model, SF_galaxies=False)

    sp1_sfrmass_local.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')

    print 'SF:  %.3f < z < %.3f' % (numpy.percentile(z_galaxies_tot_local_flat, 16), numpy.percentile(z_galaxies_tot_local_flat, 84))



    lmass_model = numpy.linspace(8.3, 11.7, 100)
    lsfr_lower = log_SFR_tomczak2016(z16, lmass_model, SF_galaxies=True)
    lsfr_upper = log_SFR_tomczak2016(z84, lmass_model, SF_galaxies=True)

    sp2_sfrmass_local.fill_between(lmass_model, lsfr_lower, lsfr_upper, color='#cccccc')
    sp2_sfrmass_local.axvline(0, color='#cccccc', lw=14, label=zfourge_label, zorder=1)

    print 'ALL: %.3f < z < %.3f' % (numpy.percentile(z_galaxies_tot_local_flat, 16), numpy.percentile(z_galaxies_tot_local_flat, 84))





    t3 = sp1_sfrmass_local.text(0.04, 0.95, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp1_sfrmass_local.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t4 = sp2_sfrmass_local.text(0.04, 0.95, 'Star-Forming\nGalaxies', fontweight='normal', fontsize=24, color='k', 
                  transform=sp2_sfrmass_local.transAxes, horizontalalignment='left', verticalalignment='top',
                  path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    leg = sp2_sfrmass_local.legend(loc=4, numpoints=1, fontsize=17, frameon=0)







    ###  printing LaTeX-reday table

    print '\\begin{table}'
    print '    \\begin{center}'
    print '    \\caption{SFR-$M_*$ Relations}'
    print '    \\label{tab:sfrmass}'
    print '    \\begin{tabular}{c|rrrc|rrc}'
    print '    \\hline \\\\[-5mm]'
    print '    \\hline \\\\[-3mm]'
    #print '    Overdensity  &   log($M_*$)  &  log$(\widetilde{M_*})$  &  log$(\widetilde{\\rm SFR})$  &  $N$ \\\\'
    #print '        Bin      & [$M_{\odot}$] &      [$M_{\odot}$]       &         [$M_{\odot}$/yr] & \\\\'
    print '    Overdensity  &   \\multicolumn{1}{c}{log($M_*$)}  &  \\multicolumn{1}{c}{log$(\\widetilde{M_*})_{\\rm ALL}$}  &  \\multicolumn{1}{c}{log$(\\widetilde{\\rm SFR})_{\\rm ALL}$}  &  $N_{\\rm ALL}$  &  \\multicolumn{1}{c}{log$(\\widetilde{M_*})_{\\rm SF}$}  &  \\multicolumn{1}{c}{log$(\\widetilde{\\rm SFR})_{\\rm SF}$}  &  $N_{\\rm SF}$ \\\\'
    print '        Bin      & \\multicolumn{1}{c}{[$M_{\\odot}$]} & \\multicolumn{1}{c}{[$M_{\\odot}$]} & \\multicolumn{1}{c}{[$M_{\\odot}$/yr]} & & \\multicolumn{1}{c}{[$M_{\\odot}$]} & \\multicolumn{1}{c}{[$M_{\\odot}$/yr]} & \\\\'
    print '    \hline \\\\[-3mm]'

    for i_overdens in range(len(overdensbars)):

        ###  don't plot lowest mass bin for intermediate- or high-density sample
        inds_plot_tot = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_tot_local[:, i_overdens] > 15)
        inds_plot_sf = (galaxies_catalog.lmassbars > 8.51) | (n_galaxies_sf_local[:, i_overdens] > 15)

        ###  errors on the median
        elog_sfr_tot = numpy.log10((sfr_median_tot_local[:, i_overdens] + sfr_nmad_tot_local[:, i_overdens]/n_galaxies_tot_local[:, i_overdens]**0.5) / sfr_median_tot_local[:, i_overdens])
        elog_sfr_sf = numpy.log10((sfr_median_sf_local[:, i_overdens] + sfr_nmad_sf_local[:, i_overdens]/n_galaxies_sf_local[:, i_overdens]**0.5) / sfr_median_sf_local[:, i_overdens])


        ###  bootstrap errors
        elo_lsfr_tot = (lsfr_p50_tot_bootstrap - lsfr_p16_tot_bootstrap)[:, i_overdens]
        ehi_lsfr_tot = (lsfr_p84_tot_bootstrap - lsfr_p50_tot_bootstrap)[:, i_overdens]
        elog_sfr_tot = (elo_lsfr_tot + ehi_lsfr_tot) / 2.

        elo_lsfr_sf = (lsfr_p50_sf_bootstrap - lsfr_p16_sf_bootstrap)[:, i_overdens]
        ehi_lsfr_sf = (lsfr_p84_sf_bootstrap - lsfr_p50_sf_bootstrap)[:, i_overdens]
        elog_sfr_sf = (elo_lsfr_sf + ehi_lsfr_sf) / 2.



        s = ''
        for i_lmass in range(len(inds_plot_sf)):

            if inds_plot_sf[i_lmass]:

                s += '                 & '

                # lower/upper bounds of massbin
                s += '$%5.2f' % galaxies_catalog.lmassbins[i_lmass]
                s += ' - '
                s += '%5.2f$' % galaxies_catalog.lmassbins[i_lmass+1]
                s += '  &  '



                # median/nmad mass of massbin ALL
                s += '$%5.2f' % lmass_median_tot_local[i_lmass, i_overdens]
                s += ' \pm '
                s += '%4.2f$' % lmass_nmad_tot_local[i_lmass, i_overdens]
                s += '  &  '

                # median/nmad SFR of massbin ALL
                s += '$%4.2f' % lsfr_median_tot_local[i_lmass, i_overdens]
                s += ' \pm '
                s += '%4.2f$' % elog_sfr_tot[i_lmass]
                s += '  &  '

                # number galaxies of bin ALL
                s += '%3i' % n_galaxies_tot_local[i_lmass, i_overdens]
                s += ' &\n'
                s += ' ' * 39


                # median/nmad mass of massbin SF
                s += '$%5.2f' % lmass_median_sf_local[i_lmass, i_overdens]
                s += ' \pm '
                s += '%4.2f$' % lmass_nmad_sf_local[i_lmass, i_overdens]
                s += '  &  '

                # median/nmad SFR of massbin SF
                s += '$%4.2f' % lsfr_median_sf_local[i_lmass, i_overdens]
                s += ' \pm '
                s += '%4.2f$' % elog_sfr_sf[i_lmass]
                s += '  &  '

                # number galaxies of bin SF
                s += '%3i' % n_galaxies_sf_local[i_lmass, i_overdens]
                s += ' \\\\'
                s += '\n'


        print s, '    \\hline \\\\[-3mm]'









