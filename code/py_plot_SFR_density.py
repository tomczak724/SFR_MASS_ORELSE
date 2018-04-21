
import mypy
import numpy
import pickle
from scipy import optimize
from astropy.io import fits
from matplotlib import pyplot
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages


def line(x, m, b):
    return m * x + b


class combined_catalog(object):

    def __init__(self, header, n_galaxies):

        self.time_stamp = time.ctime()
        self.header = header

        self.n_galaxies = n_galaxies



catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))

inds = (catalog.z >= 0.70) & \
       (catalog.z <= 0.92)

inds = (catalog.z >= 0.92) & \
       (catalog.z <= 1.50)

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
catalog.z = catalog.z[inds]













#############################################
###  Plotting for Star-Forming subsample  ###
#############################################

if True:

    fig = pyplot.figure(figsize=(8.5, 11.))

    sp1 = fig.add_subplot(111)

    sp1.minorticks_on()
    sp1.set_xlabel('log( 1 + $\delta_{gal}$ )')
    sp1.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1.axis([0, 2, -0.1, 1.9])



    cmap = pyplot.cm.rainbow



    for i_lmass in range(len(catalog.lmassbars)):


        ###  skipping log(M*)~8.5
        if catalog.lmassbars[i_lmass] < 8.51:
            continue


        n_galaxies = []
        overdens_median = []
        sfr_median = []
        sfr_nmad = []

        for i_overdens in range(len(catalog.overdensbars)):

            inds_lmass_overdens = (catalog.UVJ_class > -1) & \
                                  (catalog.digi_lmass == i_lmass+1) & \
                                  (catalog.digi_overdens == i_overdens+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_lmass = cmap(i_lmass / (len(catalog.lmassbars)-1.))


        ###  plotting SFR vs. overdensity
        sp1.plot(overdens_median, numpy.log10(sfr_median), ls='-', marker='o', 
                 lw=4, mew=4, ms=9, color='#cccccc', mec='#cccccc', zorder=1)
        sp1.errorbar(overdens_median, numpy.log10(sfr_median), yerr=elog_sfr, ls='-', lw=1.5, marker='o', 
                     mfc=color_lmass, color=color_lmass, mec=color_lmass, ecolor=color_lmass, 
                     mew=1, ms=8, elinewidth=2, zorder=2)


        ###  adding text indicating M* bin
        t = '10$^{%s}$ M$_{\odot}$' % catalog.lmassbars[i_lmass]
        y_text = sp1.transLimits.transform((0, numpy.log10(sfr_median[-1])))[1]
        sp1.text(0.97, y_text, t, color=color_lmass, transform=sp1.transAxes,
                 horizontalalignment='right', verticalalignment='center', fontweight='bold')



























#####################################################
###  Plotting for Star-Forming and Total samples  ###
#####################################################

if True:

    #fig = pyplot.figure(figsize=(13.6, 9.6))
    fig = pyplot.figure(figsize=(12.24, 8.89))

    sp2 = fig.add_subplot(122)
    sp1 = fig.add_subplot(121)

    sp1.minorticks_on()
    sp2.minorticks_on()
    sp1.set_xlabel('log( 1 + $\delta_{gal}$ )')
    sp2.set_xlabel('log( 1 + $\delta_{gal}$ )')
    sp1.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1.axis([-0.05, 1.7, -0.15, 1.9])
    sp2.axis([-0.05, 1.7, -0.15, 1.9])

    fig.subplots_adjust(wspace=0, left=0.08, bottom=0.09, top=0.98, right=0.98)




    font = FontProperties()
    font.set_family('sans-serif')


    t = sp1.text(0.03, 0.97, 'Star-Forming', fontweight='normal', fontsize=24, color='b', 
                 transform=sp1.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t = sp2.text(0.03, 0.97, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                 transform=sp2.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])






    cmap = pyplot.cm.rainbow

    for i_lmass in range(len(catalog.lmassbars))[::-1]:


        ###  skipping log(M*)~8.5
        if catalog.lmassbars[i_lmass] < 8.51:
            continue




        ###  plotting SF galaxies

        n_galaxies = []
        overdens_median = []
        sfr_median = []
        sfr_nmad = []

        for i_overdens in range(len(catalog.overdensbars)):

            inds_lmass_overdens = (catalog.UVJ_class == 1) & \
                                  (catalog.digi_lmass == i_lmass+1) & \
                                  (catalog.digi_overdens == i_overdens+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_lmass = cmap(i_lmass / (len(catalog.lmassbars)-1.))


        ###  plotting SFR vs. overdensity
        label = '10$^{%s}$ M$_{\odot}$' % catalog.lmassbars[i_lmass]
        sp1.errorbar(overdens_median, numpy.log10(sfr_median), yerr=elog_sfr, ls='-', lw=1.5, marker='o', 
                     mfc=color_lmass, color=color_lmass, mec='k', ecolor=color_lmass, 
                     mew=1., ms=9, elinewidth=2, zorder=2, label=label)








        ###  plotting Total galaxies

        n_galaxies = []
        overdens_median = []
        sfr_median = []
        sfr_nmad = []

        for i_overdens in range(len(catalog.overdensbars)):

            inds_lmass_overdens = (catalog.digi_lmass == i_lmass+1) & \
                                  (catalog.digi_overdens == i_overdens+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_lmass = cmap(i_lmass / (len(catalog.lmassbars)-1.))


        ###  plotting SFR vs. overdensity
        label = '10$^{%s}$ M$_{\odot}$' % catalog.lmassbars[i_lmass]
        label = '%s' % catalog.lmassbars[i_lmass]
        sp2.errorbar(overdens_median, numpy.log10(sfr_median), yerr=elog_sfr, ls='-', lw=1.5, marker='o', 
                     mfc=color_lmass, color=color_lmass, mec='k', ecolor=color_lmass, 
                     mew=1., ms=9, elinewidth=2, zorder=2, label=label)


    sp2.legend(loc=1, fontsize=18, title='log( M$_*$ / M$_{\odot}$ )', numpoints=1)



























#####################################################
###  Plotting for Star-Forming and Total samples  ###
###          with 4 overdensity bins              ###
#####################################################

if True:


    overdensbins2 = [-99, 0.3, 0.75, 1.2, 99]
    digi_overdens2 = numpy.digitize(catalog.overdens, overdensbins2)





    #fig = pyplot.figure(figsize=(13.6, 9.6))
    fig = pyplot.figure(figsize=(12.24, 8.89))

    sp2 = fig.add_subplot(122)
    sp1 = fig.add_subplot(121)

    sp1.minorticks_on()
    sp2.minorticks_on()
    sp1.set_xlabel('log( 1 + $\delta_{gal}$ )')
    sp2.set_xlabel('log( 1 + $\delta_{gal}$ )')
    sp1.set_ylabel('log( SFR$_{UV+IR}$ / [M$_{\odot}$ / yr] )')

    sp1.axis([-0.05, 1.7, -0.15, 1.9])
    sp2.axis([-0.05, 1.7, -0.15, 1.9])

    fig.subplots_adjust(wspace=0, left=0.08, bottom=0.09, top=0.98, right=0.98)




    font = FontProperties()
    font.set_family('sans-serif')


    t = sp1.text(0.03, 0.97, 'Star-Forming', fontweight='normal', fontsize=24, color='b', 
                 transform=sp1.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

    t = sp2.text(0.03, 0.97, 'All Galaxies', fontweight='normal', fontsize=24, color='k', 
                 transform=sp2.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])






    cmap = pyplot.cm.rainbow

    for i_lmass in range(len(catalog.lmassbars))[::-1]:


        ###  skipping log(M*)~8.5
        if catalog.lmassbars[i_lmass] < 8.51:
            continue




        ###  plotting SF galaxies

        n_galaxies = []
        overdens_median = []
        overdens_nmad = []
        sfr_median = []
        sfr_nmad = []

        for i_overdens in range(len(overdensbins2)-1):

            inds_lmass_overdens = (catalog.UVJ_class == 1) & \
                                  (catalog.digi_lmass == i_lmass+1) & \
                                  (digi_overdens2 == i_overdens+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))
            overdens_nmad.append(mypy.nmad(catalog.overdens[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        overdens_nmad = numpy.array(overdens_nmad)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_lmass = cmap(i_lmass / (len(catalog.lmassbars)-1.))


        ###  plotting SFR vs. overdensity
        label = '10$^{%s}$ M$_{\odot}$' % catalog.lmassbars[i_lmass]
        sp1.errorbar(overdens_median, numpy.log10(sfr_median), yerr=elog_sfr, xerr=overdens_nmad,
                     ls='-', lw=1.5, marker='o', 
                     mfc=color_lmass, color=color_lmass, mec='k', ecolor=color_lmass, 
                     mew=1., ms=9, elinewidth=2, zorder=2, label=label)








        ###  plotting Total galaxies

        n_galaxies = []
        overdens_median = []
        overdens_nmad = []
        sfr_median = []
        sfr_nmad = []

        for i_overdens in range(len(overdensbins2)-1):

            inds_lmass_overdens = (catalog.digi_lmass == i_lmass+1) & \
                                  (digi_overdens2 == i_overdens+1)

            n_galaxies.append(numpy.count_nonzero(inds_lmass_overdens))
            sfr_median.append(numpy.median(catalog.SFR_UVIR[inds_lmass_overdens]))
            sfr_nmad.append(mypy.nmad(catalog.SFR_UVIR[inds_lmass_overdens]))
            overdens_median.append(numpy.median(catalog.overdens[inds_lmass_overdens]))
            overdens_nmad.append(mypy.nmad(catalog.overdens[inds_lmass_overdens]))

        n_galaxies = numpy.array(n_galaxies)
        overdens_median = numpy.array(overdens_median)
        overdens_nmad = numpy.array(overdens_nmad)
        sfr_median = numpy.array(sfr_median)
        sfr_nmad = numpy.array(sfr_nmad)
        n = len(sfr_median)

        elog_sfr = numpy.log10((sfr_median + sfr_nmad/n_galaxies**0.5) / sfr_median)



        color_lmass = cmap(i_lmass / (len(catalog.lmassbars)-1.))


        ###  plotting SFR vs. overdensity
        label = '10$^{%s}$ M$_{\odot}$' % catalog.lmassbars[i_lmass]
        label = '%s' % catalog.lmassbars[i_lmass]
        sp2.errorbar(overdens_median, numpy.log10(sfr_median), yerr=elog_sfr, xerr=overdens_nmad,
                     ls='-', lw=1.5, marker='o', 
                     mfc=color_lmass, color=color_lmass, mec='k', ecolor=color_lmass, 
                     mew=1., ms=9, elinewidth=2, zorder=2, label=label)






        ###  fitting line to data points
        fit, cov = optimize.curve_fit(line, overdens_median, numpy.log10(sfr_median), p0=[-0.1, 0.9], sigma=elog_sfr)
        print '%4.1f :   (%s, %s)' % (catalog.lmassbars[i_lmass], fit[0], cov[0][0]**0.5)









    sp2.legend(loc=1, fontsize=18, title='log( M$_*$ / M$_{\odot}$ )', numpoints=1)


































