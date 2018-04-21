
import mypy
import numpy
import pickle
from scipy import optimize, stats
from astropy.io import fits
from matplotlib import pyplot
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages




class combined_catalog(object):

    def __init__(self, header, n_galaxies):

        self.time_stamp = time.ctime()
        self.header = header

        self.n_galaxies = n_galaxies



zlo, zhi = 0.6, 1.3

galaxies_catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))













###############################
###  plotting zspec histograms
###############################

if True:

    fig = pyplot.figure(figsize=(15.3, 6.))
    sp = fig.add_subplot(111)

    sp.minorticks_on()
    sp.set_xlabel('spectroscopic redshift', size=20)
    sp.set_ylabel('Number', size=20)
    sp.set_xlim(zlo-0.04, zhi+0.04)

    fig.subplots_adjust(left=0.06, bottom=0.1)

    hists = []

    colors = ['#3f4ea1', '#6db388', '#d92120']
    overdens_labels = ['          log(1+$\\delta_{gal}$) < 0.5',
                       '0.5 < log(1+$\\delta_{gal}$) < 1.0',
                       '1.0 < log(1+$\\delta_{gal}$)']
    offsets = [-0.0015, 0.00, 0.0015]
    for i in range(len(galaxies_catalog.overdensbars)):

        inds = (galaxies_catalog.digi_overdens == i+1)
        print numpy.count_nonzero(inds)
        h_plot = sp.hist(galaxies_catalog.zspec[inds], histtype='stepfilled', alpha=0.35, lw=3, color=colors[i],
                         bins=35, range=(zlo+offsets[i], zhi+offsets[i]), label=overdens_labels[i])
        h_plot = sp.hist(galaxies_catalog.zspec[inds], histtype='step', lw=3, color=colors[i],
                         bins=35, range=(zlo+offsets[i], zhi+offsets[i]))
        h = numpy.histogram(galaxies_catalog.zspec[inds], bins=35, range=(zlo, zhi))
        hists.append(h)
        #asdf = sp.axvline(-1, color=colors[i], lw=3, label=overdens_labels[i])

        print 'median zspec: %.2f' % numpy.median(galaxies_catalog.zspec[inds])

    sp.legend(loc=0, fontsize=20)
    ymaxes = [h[0].max() for h in hists]
    sp.set_ylim(0, max(ymaxes)*1.2)



    fig.savefig('../figures/redshift_histograms.pdf')
    pyplot.close()

















###############################
###  plotting zspec histograms
###############################

overdens_labels = ['          log(1+$\\delta_{gal}$) < 0.5',
                   '0.5 < log(1+$\\delta_{gal}$) < 1.0',
                   '1.0 < log(1+$\\delta_{gal}$)          ']

if True:

    fig = pyplot.figure(figsize=(15.3, 6.))
    sp = fig.add_subplot(111)

    sp.minorticks_on()
    sp.set_xlabel('spectroscopic redshift', size=20)
    sp.set_ylabel('Number', size=20)
    sp.set_xlim(zlo-0.04, zhi+0.04)

    fig.subplots_adjust(left=0.06, bottom=0.105)

    hists = []

    colors       = ['#3f4ea1', '#6db388', '#d92120']
    colors_dark  = ['#323e81', '#4e976a', '#b11b1b']
    colors_faded = ['#919bd4', '#abd4ba', '#ee9090']
    offsets = [-0.0015, 0.00, 0.0015]

    n_hist_bins = 40

    for i in range(len(galaxies_catalog.overdensbars)):

        inds = (galaxies_catalog.digi_overdens == i+1)
        inds_sf = (galaxies_catalog.digi_overdens == i+1) & (galaxies_catalog.UVJ_class == 1)
        print numpy.count_nonzero(inds), numpy.count_nonzero(inds_sf)
        h_plot = sp.hist(galaxies_catalog.zspec[inds], histtype='stepfilled', alpha=0.35, lw=3, color=colors[i],
                         bins=n_hist_bins, range=(zlo+offsets[i], zhi+offsets[i]))
        h_plot = sp.hist(galaxies_catalog.zspec[inds], histtype='step', lw=3, color=colors_dark[i],
                         bins=n_hist_bins, range=(zlo+offsets[i], zhi+offsets[i]))
        h = numpy.histogram(galaxies_catalog.zspec[inds], bins=35, range=(zlo, zhi))
        hists.append(h)
        #asdf = sp.axvline(-1, color=colors[i], lw=3, label=overdens_labels[i])

        print 'median zspec: %.2f' % numpy.median(galaxies_catalog.zspec[inds])

        # add a rectangle
        label = overdens_labels[i]
        #label += '(%i/%i)' % (numpy.count_nonzero(inds), numpy.count_nonzero(inds_sf))
        rect = mpatches.Rectangle((0, 0), 0.05, 0.1, 
                                  color=colors_faded[i], 
                                  ec=colors_dark[i], lw=3,
                                  label=label)

        sp.add_patch(rect)



    leg = sp.legend(loc=0, fontsize=19, frameon=False)
    ymaxes = [h[0].max() for h in hists]
    sp.set_ylim(0, max(ymaxes)*1.2)



    fig.savefig('../figures/redshift_histograms.pdf')
    pyplot.close()








